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
    print(gff_path, 'has', len(ref_annotation), 'lines')
    ref_annotation = ref_annotation[ref_annotation["feature"] == "transcript"]
    ref_annotation.to_csv(gff_output_path, sep="\t", index=False, header=False)
    print('Saved', len(ref_annotation), 'transcripts to', gff_output_path)

def get_info(args, confs, tmpDir, stepDir):
    map_gff_path = stepDir["map_to_genome"] + "/reference_mapped.gff"
    ref_gff_path = stepDir["prepare_ref_annotation"] + "/reference.gff"
    aln_gff_path = stepDir["ncrna_alignment_parsing"] + "/alignment_annotation.gff"
    infernal_gff_path = stepDir["parse_infernal"] + "/rfam_annotation_genome.gff"

    all_gff_path = tmpDir + "/all_references.gff"
    wrote = join_gffs_in_one([map_gff_path, ref_gff_path, aln_gff_path, infernal_gff_path], 
                                all_gff_path)

    if wrote:
        retrieval_stats = update_with_info(all_gff_path, 
                                        tmpDir + "/annotation_with_meta.gff")
        print("Information retrieval results: " + "\n\t".join([stat_name+": "+str(value)
                                    for stat_name, value in retrieval_stats.items()]))
        json.dump(retrieval_stats, open(tmpDir + '/rnacentral_consulting_stats.json', 'w'), indent=4)
    return True

def get_functional_info(args, confs, tmpDir, stepDir):
    all_gff_path = stepDir["get_info"] + "/annotation_with_meta.gff"
    retrieval_stats = retrieve_func_annotation([all_gff_path], tmpDir + "/retrieved_functions.id2go", confs)
    json.dump(retrieval_stats, open(tmpDir + '/retrieval_stats.json', 'w'), indent=4)
    return True

def run_gffcompare(args, confs, tmpDir, stepDir):
    #rfam_mappings = stepDir["parse_infernal"] + "/rfam_annotation_genome.gff"
    tRNA_mappings = stepDir["parse_trna"] + "/tRNAs.gff"
    #novel_mappings = tmpDir + "/novel_mappings.gff"
    #filter_non_transcripts(rfam_mappings, novel_mappings)
    gffs = [tRNA_mappings]

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
    
def merge_attrs(best_id, other_ids, all_hits):
    best_attrs = get_gff_attributes(all_hits[best_id]['attribute'])
    other_attrs = [get_gff_attributes(all_hits[other]['attribute']) for other in other_ids if other != best_id]

    for other in other_attrs:
        for attr, val in other.items():
            if not attr in best_attrs:
                best_attrs[attr] = val
    
    rna_types = list(set([d['type'] if 'type' in d else '' for d in [best_attrs] + other_attrs]))
    rna_types.sort(key = lambda tp: (tp.count('|'), not 'other' in tp.lower(), len(tp)))
    best_attrs['type'] = rna_types[-1]
    if len(rna_types) > 1:
        print("from", rna_types, 'choose', rna_types[-1])

    return get_gff_attributes_str(best_attrs)

unkown_rna = ["other", "misc_rna", "misc_RNA"]

def update_attrs(attr_str, type_tree):
    attrs = get_gff_attributes(attr_str)
    if "family" in attrs:
        attrs["rfam"] = attrs["family"]
    
    if not "type" in attrs:
        attrs["type"] = "other"
    #Replace 'misc_rna' with 'other'
    if attrs["type"] in unkown_rna:
        attrs["type"] = "other"
    attrs["type"] = standardize_rna_type(attrs["type"], type_tree)
    if attrs["type"] == None:
        print(attr_str)
        print('Invalid type in', attr_str)
        quit(1)
    return get_gff_attributes_str(attrs)

def remove_redundancies(args, confs, tmpDir, stepDir):
    print("Reading raw annotation")
    all_mappints_path = stepDir["run_gffcompare"] + "/all_mappings.gff"
    raw_annotation = pd.read_csv(all_mappints_path, sep="\t", header=None)
    raw_annotation.columns = ["seqname", "source", "feature", "start", "end", "score", "strand",
                                "frame", "attribute"]
    
    print('Updating attributes')
    type_tree = load_rna_types2()
    raw_annotation["attribute"] = raw_annotation.apply(
        lambda row: update_attrs(row["attribute"], type_tree),axis=1
    )
    hits = {}
    print('Reading rows from', all_mappints_path)
    for index, row in raw_annotation.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        hits[attrs["ID"]] = row
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
    '''print("Merging attrs")
    merged_attrs = {non_redundant_ids[i]: merge_attrs(non_redundant_ids[i], redundant_ids[i], hits) 
        for i in range(len(non_redundant_ids))}'''
    
    '''def is_filtered(row, valid_ids):
        ID = row["attribute"].split('ID=')[-1].split(';')[0]
        return not ID in valid_ids'''
    print("Calculating filter")
    raw_annotation["filtered"] = raw_annotation.apply(
        lambda row: not row["attribute"].split('ID=')[-1].split(';')[0] in non_redundant_ids,
        axis=1
    )

    print("All annotations: " + str(len(raw_annotation)))
    raw_annotation = raw_annotation[raw_annotation["filtered"] == False]
    print("Without redundancies: " + str(len(raw_annotation)))
    #annotation.is_copy = False
    raw_annotation = raw_annotation.drop('filtered', axis=1)
    '''def get_merged(row, merged_attrs):
        attrs = get_gff_attributes(row["attribute"])
        return merged_attrs[attrs["ID"]]'''
    print("Applying merged attributes")
    '''raw_annotation["attribute"] = raw_annotation.apply(
        lambda row: merged_attrs[get_gff_attributes(row["attribute"])['ID']],
        axis=1
    )'''
    print("Writing final annotation")
    raw_annotation.to_csv(tmpDir + "/annotation.gff", sep="\t", index=False, header=False)

    return True
