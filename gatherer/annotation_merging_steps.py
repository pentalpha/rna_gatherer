import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from gatherer.bioinfo import *
from gatherer.util import *
from gatherer.netutils import *
from gatherer.rna_type import *

def filter_non_transcripts(gff_path, gff_output_path):
    ref_annotation = pd.read_csv(gff_path, sep="\t", header=None, names=["seqname", 
            "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    ref_annotation = ref_annotation[ref_annotation["feature"] == "transcript"]
    ref_annotation_file = gff_output_path
    ref_annotation.to_csv(ref_annotation_file, sep="\t", index=False, header=False)

def get_info(args, confs, tmpDir, stepDir):
    map_gff_path = stepDir["map_to_genome"] + "/reference_mapped.gff"
    ref_gff_path = stepDir["prepare_ref_annotation"] + "/reference.gff"
    aln_gff_path = stepDir["ncrna_alignment_parsing"] + "/alignment_annotation.gff"

    all_gff_path = tmpDir + "/all_references.gff"
    wrote = join_gffs_in_one([map_gff_path, ref_gff_path, aln_gff_path], 
                                all_gff_path)

    if wrote:
        seqs = {}
        if "reference_fasta" in args:
            reference_fasta_path = args["reference_fasta"]
            seqs = readSeqsFromFasta(args["reference_fasta"])
            seqs = {header_to_id(header): seq for header,seq in seqs}
        retrieval_stats = update_with_info(all_gff_path, 
                                        tmpDir + "/annotation_with_meta.gff",
                                        confs, seqs_dict = seqs)
        print("Information retrieval results: " + "\n\t".join([stat_name+": "+str(value)
                                    for stat_name, value in retrieval_stats.items()]))
        
        '''retrieve_func_annotation_rnacentral(all_gff_path, tmpDir + "/retrieved_functions.id2go",
                                confs, args['taxon_id'])'''
    return True

def get_functional_info(args, confs, tmpDir, stepDir):
    all_gff_path = stepDir["get_info"] + "/annotation_with_meta.gff"
    retrieve_func_annotation_rnacentral(all_gff_path, tmpDir + "/retrieved_functions.id2go",
                                confs, args['taxon_id'])
    return True

def run_gffcompare(args, confs, tmpDir, stepDir):
    rfam_mappings = stepDir["parse_infernal"] + "/rfam_annotation_genome.gff"
    tRNA_mappings = stepDir["parse_trna"] + "/tRNAs.gff"
    novel_mappings = tmpDir + "/novel_mappings.gff"
    filter_non_transcripts(rfam_mappings, novel_mappings)
    gffs = [novel_mappings,tRNA_mappings]

    lnc_mappings = stepDir["lnc_alignment_parsing"] + "/lncRNA_annotation.gff"
    if os.path.exists(lnc_mappings):
        gffs = gffs + [lnc_mappings]

    ref = stepDir["get_info"] + "/annotation_with_meta.gff"
    if os.path.exists(ref):
        gffs = [os.path.abspath(ref)] + gffs

    gffs = remove_redundant_names(gffs, tmpDir)

    all_lines = []
    for gff in gffs:
        with open(gff,'r') as stream:
            all_lines = all_lines + [line.rstrip("\n") 
                                    for line in stream.readlines()]
    all_gffs = "\n".join(all_lines)+"\n"
    all_mappings_path = tmpDir+"/all_mappings.gff"
    with open(all_mappings_path, 'w') as stream:
        stream.write(all_gffs)
    cmd = " ".join(["cd ", tmpDir, " && ", confs["gffcompare"], "-o gffcmp"] + gffs)
    ret = runCommand(cmd)
    return ret == 0

def best_id_in_source(ids, hits, source):
    best = ids[0]
    for i in ids:
        best_hit = hits[best]
        hit = hits[i]
        best_len = int(best_hit['end']) - int(best_hit['start'])
        current_len = int(hit['end']) - int(hit['start'])
        if current_len > best_len:
            best = i
    return best

def best_id(ids, hits):
    id_by_source = {}
    for id_str in ids:
        hit = hits[id_str]
        source = hit["source"]
        if not source in id_by_source:
            id_by_source[source] = list()
        id_by_source[source].append(id_str)
    best_by_source = {}
    for source in id_by_source.keys():
        best_by_source[source] = best_id_in_source(id_by_source[source], hits, source)
    if "reference" in best_by_source:
        return best_by_source["reference"]
    elif "reference_mapping" in best_by_source:
        return best_by_source["reference_mapping"]
    elif "db_alignment" in best_by_source:
        return best_by_source["db_alignment"]
    elif "tRNAscan-SE" in best_by_source:
        return best_by_source["tRNAscan-SE"]
    elif "cmscan" in best_by_source:
        return best_by_source["cmscan"]
    elif "rnasamba" in best_by_source:
        return best_by_source["rnasamba"]
    else:
        print("Error: no known source in " + str(id_by_source))
        return None

unkown_rna = ["other", "misc_rna", "misc_RNA"]

def update_attrs(attr_str):
    attrs = get_gff_attributes(attr_str)
    if "family" in attrs:
        attrs["rfam"] = attrs["family"]
    
    if "rfam" in attrs:
        new_type = get_rna_type(attrs["rfam"])
        attrs["type"] = new_type
    else:
        if not "type" in attrs:
            attrs["type"] = "other"
        #Replace 'misc_rna' with 'other'
        if attrs["type"] in unkown_rna:
            attrs["type"] = "other"
        attrs["type"] = get_full_type(attrs["type"])

    return get_gff_attributes_str(attrs)

def remove_redundancies(args, confs, tmpDir, stepDir):
    print("Reading raw annotation")
    raw_annotation = pd.read_csv(stepDir["run_gffcompare"] + "/all_mappings.gff", sep="\t", header=None)
    raw_annotation.columns = ["seqname", "source", "feature", "start", "end", "score", "strand",
                                "frame", "attribute"]
    hits = {}
    for index, row in raw_annotation.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        hits[attrs["ID"]] = row
    raw_annotation["attribute"] = raw_annotation.apply(
        lambda row: update_attrs(row["attribute"]),axis=1
    )
    redundant_ids = []
    loci_file = stepDir["run_gffcompare"] + "/gffcmp.loci"
    print("Reading loci file: " + loci_file)
    loci = pd.read_csv(loci_file, sep="\t")
    new_gff_lines = []
    novel_ids = 0
    not_found_desc = 0
    with open(loci_file, "r") as loci_stream:
        for raw_line in loci_stream:
            line = raw_line.rstrip("\n")
            cells = line.split("\t")
            ids_cells = cells[2:len(cells)]
            ids = []
            for cell in ids_cells:
                for subcell in cell.split(","):
                    if subcell != "-":
                        ids.append(subcell)
            redundant_ids.append(ids)

    print("Calculating best ids")
    non_redundant_ids = [best_id(id_list, hits) for id_list in redundant_ids]
    def is_filtered(row, valid_ids):
        attrs = get_gff_attributes(row["attribute"])
        ID = attrs["ID"]
        filtered = not ID in valid_ids
        return filtered
    raw_annotation["filtered"] = raw_annotation.apply(
        lambda row: is_filtered(row,non_redundant_ids),
        axis=1
    )

    print("All annotations: " + str(len(raw_annotation)))
    raw_annotation = raw_annotation[raw_annotation["filtered"] == False]
    print("Without redundancies: " + str(len(raw_annotation)))
    #annotation.is_copy = False
    raw_annotation = raw_annotation.drop('filtered', axis=1)
    print("Writing final annotation")
    raw_annotation.to_csv(tmpDir + "/annotation.gff", sep="\t", index=False, header=False)

    return True
