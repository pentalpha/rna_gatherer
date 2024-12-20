import json
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
    tp = gene_details[rna_name]['type'].split("|")
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
    no_id = 0
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        if "rfam" in attrs:
            if attrs["rfam"] != "None":
                count.add(attrs["rfam"])
            else:
                no_id += 1
        else:
            no_id += 1
    return count, no_id

def get_ncrna_type(attrs_str):
    attrs = get_gff_attributes(attrs_str)
    if attrs['ID'] == 'URS0001BBEC3B.0':
        attrs['type'] = 'Gene|miRNA'
    if "type" in attrs:
        tp = attrs["type"].replace(';', '|')
        if not tp.startswith('Gene'):
            tp = 'other|'+tp
        return tp
    else:
        return "other"
    
def get_ncrna_ids(df):
    ids = []
    for index, row in df.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        ids.append(attrs['ID'])
    return ids

def get_subtype(tp):
    parts = tp.replace(';', '|').split("|")
    if len(parts) == 1:
        return "No subtype"
    else:
        return parts[1].replace(';No subtype', '').rstrip(';')

def review_df(df, sources):
    by_source = {src: len(src_group) for src, src_group in df.groupby(["source"])}
    for src in sources:
        if not src in by_source:
            by_source[src] = 0
    families, no_family = get_rfam_ids(df)
    ncrna_ids = get_ncrna_ids(df)
    return families, by_source, ncrna_ids, len(df)-no_family

def review_annotations(args, confs, tmpDir, stepDir):
    annotation = pd.read_csv(stepDir["contaminant_removal"] + "/annotation.gff", sep="\t", header=None,
        names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    print("Reviewing annotation sources")
    source_list = annotation["source"].unique().tolist()
    
    type_lists = []
    print("Retrieving rna_type")
    annotation["rna_type"]    = annotation.apply(lambda row: get_ncrna_type(row["attribute"]).split("|")[0],
                                                axis = 1)
    annotation["rna_subtype"] = annotation.apply(lambda row: get_subtype(get_ncrna_type(row["attribute"])),
                                                axis = 1)
    sources = source_list
    sources.sort()
    df_columns = ['ncRNA Type', 'Number of ncRNAs', 'RFAM Families', 'With RFAM Family'] + sources
    print(df_columns)
    for rna_type, type_df in annotation.groupby(['rna_type']):
        print('Main type:', rna_type)
        if rna_type == 'Gene;sRNA':
            attrs = type_df['attribute'].tolist()
            for t in attrs:
                print(t)
        print('Secondary types:', type_df['rna_subtype'].unique().tolist())
        families, by_source, ncrna_ids, with_family = review_df(type_df,sources)
        by_source['ncRNA Type'] = rna_type
        by_source['RFAM Families'] = len(families)
        by_source['Number of ncRNAs'] = len(ncrna_ids)
        by_source['With RFAM Family'] = with_family / len(ncrna_ids)
        for c in sources:
            by_source[c] = by_source[c] / len(ncrna_ids)
        by_source['ncRNA IDs'] = ncrna_ids
        by_source['RFAM IDs'] = list(families)
        print([by_source[c] for c in df_columns])
        type_line = by_source

        subtype_lines = []
        for subtype, subtype_df in type_df.groupby(['rna_subtype']):
            subtype_families, subtype_by_source, subtype_ncrna_ids, with_family = review_df(subtype_df, sources)
            if subtype != 'No subtype':
                subtype_by_source['ncRNA Type'] = rna_type+';'+subtype
                subtype_by_source['RFAM Families'] = len(subtype_families)
                subtype_by_source['Number of ncRNAs'] = len(subtype_ncrna_ids)
                subtype_by_source['With RFAM Family'] = (with_family / len(subtype_ncrna_ids))*100
                for c in sources:
                    subtype_by_source[c] = (subtype_by_source[c] / len(subtype_ncrna_ids))*100
                subtype_by_source['ncRNA IDs'] = subtype_ncrna_ids
                subtype_by_source['RFAM IDs'] = list(subtype_families)
                print([subtype_by_source[c] for c in df_columns])
                subtype_lines.append(subtype_by_source)
        subtype_lines.sort(key=lambda x: x['Number of ncRNAs'], reverse=True)
        if len(subtype_lines) >= 2:
            type_lists.append([type_line, subtype_lines])
        else:
            type_lists.append([type_line, []])
    type_lists.sort(key=lambda x: x[0]['Number of ncRNAs'], reverse=True)

    all_families, all_by_source, all_ncrna_ids, with_family = review_df(annotation, sources)
    all_by_source['ncRNA Type'] = 'All'
    all_by_source['RFAM Families'] = len(all_families)
    all_by_source['Number of ncRNAs'] = len(all_ncrna_ids)
    all_by_source['With RFAM Family'] = (with_family / len(all_ncrna_ids))*100
    for c in sources:
        all_by_source[c] = (all_by_source[c] / len(all_ncrna_ids))*100
    all_by_source['ncRNA IDs'] = all_ncrna_ids
    all_by_source['RFAM IDs'] = list(all_families)
    type_lists.append([all_by_source, []])
    review_rows = []
    for tp, subtps in type_lists:
        review_rows.append(tp)
        for subtp_line in subtps:
            review_rows.append(subtp_line)

    df_header = '\t'.join(df_columns)
    df_header = df_header.replace("reference_mapping","Reference Mapping")
    df_header = df_header.replace("rnasamba","RNA Samba")
    df_header = df_header.replace("db_alignment","Database Alignment")

    with open(tmpDir+"/type_review.tsv",'w') as stream:
        stream.write(df_header+"\n")
        for cells in review_rows:
            newline = "\t".join([str(cells[col]) for col in df_columns])
            stream.write(newline+'\n')
            print(newline)
    
    types_json = json.dumps(review_rows, indent=4)
    open(tmpDir+"/type_review.json",'w').write(types_json)
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
