import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from ete3 import NCBITaxa
from gatherer.bioinfo import *
from gatherer.util import *

def make_transcriptome_file(annotation_path, genome_path, outdir):
    print("Loading annotation")
    annotation = pd.read_csv(annotation_path, sep="\t", header=None,
                names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    print("Loading genome: " + genome_path)
    genome_dict = seqListToDict(readSeqsFromFasta(genome_path), header_to_name = header_to_id)
    transcriptome = []
    lncRNA_seqs = []
    print("Creating transcriptome file")
    for index, row in annotation.iterrows():
        #print(fasta_header)
        if str(row["seqname"]) in genome_dict:
            s = genome_dict[str(row["seqname"])] #cant find key PGUA01000001.1 #TODO
            new_header = get_gff_attributes(row["attribute"])["ID"]
            from_seq = int(row["start"])
            to_seq = int(row["end"])
            begin = min(from_seq,to_seq)-1
            up_to = max(from_seq,to_seq)
            new_seq = s[begin:up_to]
            if "lncRNA" in row["attribute"]:
                lncRNA_seqs.append((new_header, new_seq))
            transcriptome.append((new_header, new_seq))
    print("Writing transcriptome")
    writeFastaSeqs(transcriptome, outdir + "/transcriptome.fasta")
    writeFastaSeqs(lncRNA_seqs, outdir + "/lncRNA.fasta")

def evol_sim(taxid1, taxid2):
    l1 = set(ncbi.get_lineage(taxid1))
    l2 = set(ncbi.get_lineage(taxid2))
    common_taxid = l1.intersection(l2)
    return len(common_taxid)

def read_annotation(gff_path):
    details = {}
    with open(gff_path,'r') as stream:
        for line in stream:
            cells = line.rstrip("\n").split("\t")
            source = cells[1]
            attrs = get_gff_attributes(cells[-1])
            attrs['source'] = source
            details[attrs['ID']] = attrs
    return details

def get_best_homolog(df):
    sorted = df.sort_values(["identity","qcovs"],
            ascending=[False,False])
    return sorted.iloc[0]

def align(query, db, outdir, threads, aligner):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outfile = outdir+"/alignments.paf"
    if not os.path.exists(outfile):
        cmd = " ".join([aligner,
            "-x splice:hq -uf -t", str(threads),
            db, query,
            ">", outfile])
        code = runCommand(cmd)
        if code != 0:
            print("Error: Alignment unsucessful.")
            os.remove(outfile)
            return None
    return outfile

def get_tax_name(id, ncbi):
    ids = ncbi.translate_to_names([id])
    return str(ids[0])

def get_rna_type_str(rna_name, gene_details):
    tp = gene_details[rna_name]['type'].split(";")
    if len(tp) == 1 or tp[0] == 'Cis-reg':
        return tp[0]
    else:
        return tp[1]

def count_type(tp_str, gene_details):
    count = 0
    for ID, details in gene_details.items():
        if get_rna_type_str(ID) == tp_str:
            count += 1
    return count

def classify(value, th):
    return "VALID" if value >= th else "INVALID"

def contaminant_removal(args, confs, tmpDir, stepDir):
    gff_annotation = stepDir["remove_redundancies"] + "/annotation.gff"
    if "contaminant_db" in confs:
        make_transcriptome_file(stepDir["remove_redundancies"] + "/annotation.gff", 
                            args['genome_link'], tmpDir)
        contaminant_db = confs["contaminant_db"]
        query = tmpDir+"/transcriptome.fasta"
        db = contaminant_db
        outdir = tmpDir
        ncbi = NCBITaxa()
        paf_file = align(query, db, outdir, int(confs['threads']), confs['minimap2'])
        if paf_file == None:
            return False
        
        print("Reading gff annotation for ncRNA details")
        gene_details = read_annotation(gff_annotation)

        print('Reading paf file')
        minimap_df = pd.read_csv(paf_file, sep='\t', header=None, index_col=False,
                    names=["qseqid","qseq_len","qstart","qend","strand",
                        "sseqid","sseq_len","sstart","send","matchs",
                        "block_len","quality","13th","14th","15th","16th","17th","18th"])
        minimap_df = minimap_df.astype({"qstart": 'int32', "qend": 'int32', "qseq_len": "int32",
                    "sstart": 'int32', "send": 'int32', "sseq_len": "int32",
                    "quality": 'int32', "block_len": "int32", "matchs": "int32"})

        print("Making new columns")

        minimap_df["type"] = minimap_df.apply(
            lambda row:  get_rna_type_str(row['qseqid'], gene_details), axis=1)
        minimap_df["id"] = np.arange(len(minimap_df))
        
        print("Filtering...")
        minimap_df["qcovs"] = minimap_df.apply(
            lambda row: (row["qend"]-row["qstart"]) / row["qseq_len"], axis=1)
        minimap_df["identity"] = minimap_df.apply(
            lambda row: row["matchs"] / row["block_len"], axis=1)
        print(str(minimap_df.head()))
        print(str(len(minimap_df)) + " alignments")
        high_th = 0.90
        minimap_df["result"] = minimap_df.apply(
            lambda row: classify(min(row['identity'], row['qcovs']), high_th), axis=1)
        minimap_filtered = minimap_df[minimap_df['result'] != "INVALID"]

        print("Finding best hits")
        best_hits = set()
        for name, hits in minimap_filtered.groupby(["qseqid"]):
            hit = get_best_homolog(hits)
            best_hits.add(hit['id'])
        minimap_filtered["best_hit"] = minimap_filtered.apply(
            lambda row: row['id'] in best_hits, axis=1)
        contaminated_df = minimap_filtered[minimap_filtered['best_hit'] == True]
        
        type_counts = {rna_type: len(type_rows) 
                        for rna_type, type_rows in contaminated_df.groupby(['type'])}
        contaminant_list = contaminated_df['qseqid'].tolist()
        report = "Contaminant types:\n"
        for tp, count in type_counts.items():
            report += "\t"+tp+"\t"+str(count)+"\n"
        report += "Contaminant list (removed from annotation):\n"
        for cont in contaminant_list:
            report += "\t"+cont+"\n"
        
        open(tmpDir+"/report.txt",'w').write(report)
        print(report)
        cleaned = 0
        with open(gff_annotation, 'r') as in_stream:
            with open(tmpDir+"/annotation.gff",'w') as out_stream:
                for line in in_stream:
                    cells = line.rstrip("\n").split("\t")
                    attrs = get_gff_attributes(cells[-1])
                    if "ID" in attrs:
                        if attrs['ID'] in contaminant_list:
                            cleaned += 1
                        else:
                            out_stream.write(line)
        print("Cleaned ", cleaned, " lines from annotation")
    else:
        runCommand("cp " + gff_annotation + " " + tmpDir+"/annotation.gff")
    return True

def read_ids2go(filepath):
    gos_dict = {}
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split("\t")
            rfam_id = cols[0]
            gos_str = cols[1]
            go_ids = gos_str.split(";")
            gos_dict[rfam_id] = set(go_ids)
    return gos_dict

def read_rfam2go(filepath):
    gos_dict = {}
    with open(filepath, 'r') as stream:
        for raw_line in stream.readlines():
            cols = raw_line.rstrip("\n").split()
            rfam_id = cols[0].split(":")[-1]
            if not rfam_id in gos_dict:
                gos_dict[rfam_id] = set()
            go_str = cols[-1]
            gos_dict[rfam_id].add(go_str)
    return gos_dict

def write_id2go(filepath, gos_dict):
    with open(filepath, 'w') as stream:
        for key, gos in gos_dict.items():
            for go in gos:
                stream.write(key + "\t" + go + "\n")

def has_rfam_alignment(df):
    count = 0
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "genbank" in attrs:
            if attrs["genbank"] != "None":
                count += 1
    return count

def has_rfam_id(df):
    count = 0
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "rfam" in attrs:
            if attrs["rfam"] != "None":
                count += 1
    return count

def number_of_genbank_ids(df):
    count = set()
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "genbank" in attrs:
            if attrs["genbank"] != "None":
                count.add(attrs["genbank"])
    return len(count)

def get_rfam_ids(df):
    count = set()
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "rfam" in attrs:
            if attrs["rfam"] != "None":
                count.add(attrs["rfam"])
    return count

def get_ncrna_type(attrs_str):
    attrs = get_gff_attributes(attrs_str)
    if "type" in attrs:
        return attrs["type"]
    else:
        return "other"

def get_subtype(tp):
    parts = tp.split(";")
    if len(parts) == 1:
        return "No subtype"
    else:
        return parts[1]

def group_rows(input_rows):
    higher_level = 100
    for row in input_rows:
        level = len(row[0].split(";"))
        if level < higher_level:
            higher_level = level

    head_rows = []
    for row in input_rows:
        level = len(row[0].split(";"))
        if level == higher_level:
            head_rows.append(row)
    if len(head_rows) == len(input_rows):
        input_rows.sort(key=lambda row: row[1], reverse=True)
        return [[row, []] for row in input_rows]
    else:
        row_groups = []
        for head_row in head_rows:
            head_row_id = ";".join(head_row[0].split(";")[0:higher_level])
            group = []
            ids = []
            for row in input_rows:
                if head_row_id in row[0] and head_row_id != row[0]:
                    group.append(row)
                    ids.append(row[0])
            #print("Group of " + head_row_id + ": " + str(ids))
            grouped_group = group_rows(group)
            row_groups.append([head_row, grouped_group])
        row_groups.sort(key=lambda row: row[0][1], reverse=True)
        
        return row_groups

def expand_groups(groups):
    rows = []
    for head_row, sub_rows in groups:
        rows.append(head_row)
        for sub_row, sub_group in sub_rows:
            rows.append(sub_row)
            if len(sub_group) > 0:
                expanded = expand_groups(sub_group)
                rows += expanded
    return rows

def sort_by_genes(input_rows):
    input_rows.sort(key=lambda row: row[0], reverse=False)
    grouped = group_rows(input_rows)
    new_rows = expand_groups(grouped)
    return new_rows

def review_df(df, sources):
    by_source = {src: len(src_group) for src, src_group in df.groupby(["source"])}
    for src in sources:
        if not src in by_source:
            by_source[src] = 0
    n_rnas = len(df)
    families = get_rfam_ids(df)
    return n_rnas, families, by_source

def review_annotations(args, confs, tmpDir, stepDir):
    annotation = pd.read_csv(stepDir["contaminant_removal"] + "/annotation.gff", sep="\t", header=None,
        names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    print("Reviewing annotation sources")
    source_list = annotation["source"].unique().tolist()
    
    type_lists = []
    print("Retrieving rna_type")
    annotation["rna_type"]    = annotation.apply(lambda row: get_ncrna_type(row["attribute"]).split(";")[0],
                                                axis = 1)
    annotation["rna_subtype"] = annotation.apply(lambda row: get_subtype(get_ncrna_type(row["attribute"])),
                                                axis = 1)
    sources = annotation["source"].tolist()
    sources.sort()
    for rna_type, type_df in annotation.groupby(['rna_type']):
        type_rnas, families, by_source = review_df(type_df,sources)
        type_line = [rna_type, type_rnas, families] + [by_source[src] for src in sources]

        subtype_lines = []
        for subtype, subtype_df in type_df.groupby(['rna_subtype']):
            subtype_rnas, subtype_families, subtype_by_source = review_df(subtype_df, sources)
            subtype_line = [subtype, subtype_rnas, subtype_families] + [subtype_by_source[src] for src in sources]
            subtype_lines.append(subtype_line)
        subtype_lines.sort(key=lambda x: x[1], reverse=True)
        if len(subtype_lines) >= 2:
            type_lists.append([type_line, subtype_lines])
        else:
            type_lists.append([type_line, []])
    type_lists.sort(key=lambda x: x[0][1], reverse=True)

    all_rnas, all_families, all_by_source = review_df(annotation, sources)
    type_lists.append(["All", all_rnas, all_families] + [all_by_source[src] for src in sources])
    review_rows = []
    for tp, subtps in type_lists:
        tp[0] += "\t"
        review_rows.append(tp)
        for subtp_line in subtps:
            subtp_line[0] = "\t"+subtp_line[0]
            review_rows.append(subtp_line)
    for i in range(len(sources)):
        sources[i] = sources[i].replace("reference_mapping","Reference Mapping")
        sources[i] = sources[i].replace("rnasamba","RNA Samba")
        sources[i] = sources[i].replace("db_alignment","Database Alignment")
    with open(tmpDir+"/type_review.tsv",'w') as stream:
        stream.write("\t".join(['ncRNA Types\t', 'Number of ncRNAs', 'RFAM Families'] + sources)+"\n")
        for cells in review_rows:
            stream.write("\t".join([str(x) for x in cells])+"\n")

    '''print("Counting RNA Types...")
    def count_rna_types(df):
        type_counts = {rna_type: len(type_annotation) 
                    for rna_type, type_annotation in df.groupby(["rna_type"])}
        expanded_type_counts = {}
        for rna_type in type_counts.keys():
            total = 0
            for another_type, counts in type_counts.items():
                if rna_type in another_type:
                    total += counts
            expanded_type_counts[rna_type] = total
        expanded_type_counts["All"] = len(df)
        return expanded_type_counts
        
    counts_by_sources = {"All": count_rna_types(annotation)}
    for source_name, source_annotation in tqdm(annotation.groupby(["source"])):
        counts_by_sources[source_name] = count_rna_types(source_annotation)
    
    print("Counting RFAM families")
    def count_rfam_families(df):
        rfam_ids = {rna_type: get_rfam_ids(type_annotation) 
                    for rna_type, type_annotation in df.groupby(["rna_type"])}
        expanded_rfam_counts = {}
        for rna_type in tqdm(rfam_ids.keys()):
            rfam_total = set()
            for another_type, ids in rfam_ids.items():
                if rna_type in another_type:
                    rfam_total.update(ids)
            expanded_rfam_counts[rna_type] = len(rfam_total)
        expanded_rfam_counts["All"] = len(get_rfam_ids(df))
        return expanded_rfam_counts
    rfam_counts = count_rfam_families(annotation)

    def make_row(type_name, type_annotation):
        new_row = [type_name]
        new_row.append(counts_by_sources["All"][type_name])
        new_row.append(rfam_counts[type_name])
        for source in source_list:
            n = 0
            source_counts = counts_by_sources[source]
            if type_name in source_counts:
                n = source_counts[type_name]
            new_row.append(n)
        return new_row

    print("Grouping by rna_type")
    for rna_type, type_annotation in tqdm(annotation.groupby(["rna_type"])):
        review_rows.append(make_row(rna_type, type_annotation))
    review_rows = sort_by_genes(review_rows)
    review_rows = [make_row("All", annotation)] + review_rows
    print("\n".join([str(row) for row in review_rows]))
    type_review_df = pd.DataFrame(review_rows, 
                            columns=["ncRNA Type", "Total", "Families"] + source_list)
    type_review_df.to_csv(tmpDir+"/type_review.tsv", sep="\t", index = False)'''
    return True

def write_transcriptome(args, confs, tmpDir, stepDir):
    make_transcriptome_file(stepDir["contaminant_removal"] + "/annotation.gff", 
                            args['genome_link'], tmpDir)
    return True

def make_id2go(args, confs, tmpDir, stepDir):
    id2go_path = confs["rfam2go"]
    if os.path.exists(id2go_path):
        print("Loading ids2go associations")
        global_ids2go = read_rfam2go(id2go_path)
        print("Loading annotation")
        annotation = pd.read_csv(stepDir["contaminant_removal"] + "/annotation.gff", sep="\t", header=None,
            names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])

        print("Associating IDs to GO terms")
        local_ids2go = {}
        ids = []
        for index, row in annotation.iterrows():
            attrs = get_gff_attributes(row["attribute"])
            ID = attrs["ID"]
            ids.append(ID)
            if "rfam" in attrs:
                rfam_id = attrs["rfam"]
                if rfam_id in global_ids2go:
                    go_list = global_ids2go[rfam_id]
                    local_ids2go[ID] = set(go_list)
        api_annotation = stepDir["get_functional_info"] + "/retrieved_functions.id2go"
        if os.path.exists(api_annotation):
            with open(api_annotation, 'r') as stream:
                for line in stream:
                    cells = line.rstrip("\n").split("\t")
                    if not cells[0] in local_ids2go:
                        local_ids2go[cells[0]] = set()
                    local_ids2go[cells[0]].add(cells[1])
        
        write_id2go(tmpDir + "/id2go.tsv", local_ids2go)
        print("Writing population: " + str(len(ids)) + " ids")
        with open(tmpDir + "/ids.txt", 'w') as stream:
            for ID in ids:
                stream.write(ID + "\n")
        return True
    else:
        print(id2go_path + " does not exist.")
        return False
