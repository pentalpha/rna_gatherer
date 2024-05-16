import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from gatherer.bioinfo import readSeqsFromFasta, filterSeqs, writeFastaSeqs, getFastaHeaders, seqListToDict
from gatherer.bioinfo import cluster_all_ranges, shortFastaHeader
from gatherer.bioinfo import read_plast_extended
from gatherer.bioinfo import get_gff_attributes, get_gff_attributes_str
from gatherer.util import runCommand, write_file, getFilesWith

def infernal(fasta, cmscan, rfam, threads):
    output_name = fasta.rstrip(".fasta") + ".tsv"
    output_name2 = fasta.rstrip(".fasta") + ".out"
    new_fasta = fasta + "_running"
    runCommand("mv " + fasta + " " + new_fasta)
    cmd = (cmscan + " -o " + output_name2 + " --tblout " + output_name + " -E 0.01 --acc --cpu "
        + str(threads) + " " + rfam + " " + new_fasta)
    runCommand("rm " + output_name)
    runCommand("rm " + fasta.rstrip(".fasta") + ".out")
    code = runCommand(cmd)
    if(code != 0):
        runCommand("mv " + new_fasta + " " + fasta)
    else:
        runCommand("mv " + new_fasta + " " + fasta + "_done")

def check_progress(output_dir):
    total = 0
    done = 0
    for f in getFilesWith(output_dir, ".fasta"):
        if "_done" in f:
            done += 1
        total += 1
    complete = (done == total)
    return complete, total, done

def print_progress(output_dir):
    complete, total, done = check_progress(output_dir)
    print("Infernal processing: " + str(done) + "/" + str(total))
    return complete

def start(output_dir, cmscan, rfam, threads):
    while(True):
        fasta = ""
        for f in getFilesWith(output_dir, ".fasta"):
            if not(("_running" in f) or ("_done" in f)):
                fasta = f
                break
        if fasta == "":
            break
        infernal(fasta, cmscan, rfam, threads)
        print_progress(output_dir)
        
def erase_comments(path):
    tmp = path+".tmp"
    runCommand("grep -v '^#' " + path + " > " + tmp)
    runCommand("mv " + tmp + " " + path)

def get_infernal_output(output_dir, output_file):
    paths = getFilesWith(output_dir, ".tsv")
    paths2 = getFilesWith(output_dir, ".out")
    if len(paths) > 0:
        if len(paths2) > 0:
            results2 = " ".join(paths2)
            runCommand("cat " + results2 + " > " + output_file.rstrip(".tsv")+".out")
            erase_comments(output_file.rstrip(".tsv")+".out")
        results = " ".join(paths)
        runCommand("cat " + results + " > " + output_file + ".tsv")
        erase_comments(output_file)
        return True
    else:
        print("No results ready yet")
        return False

def read_infernal_output(infernal_path):
    infernal_tsv = infernal_path
    full_names = {}
    lines = []
    print("Reading lines from infernal tabular output")
    with open(infernal_tsv) as f:
        line = f.readline()
        to_print = 0
        while line != "":
            if line[0] != "#":
                elements = [x for x in line.split() if x != ""]
                cols = {"rna_name": elements[0],"rfam":elements[1], "seqname": elements[2], 
                        "qstart": elements[7], "qend": elements[8], "strand": elements[9],
                        "evalue": elements[15]}
                if to_print > 0:
                    print("Line:\n\t"+line)
                    print("\t"+str(cols))
                to_print -= 1
                lines.append(cols)
            line = f.readline()
    #print(str(lines))

    print("Creating gff file about seqs identified")
    rows = []
    name_used = {}
    for hit in lines:
        base_name = hit['seqname']+"_"+hit['rna_name']
        if not base_name in name_used:
            name_used[base_name] = 0
        unique_name = base_name+"_"+str(name_used[base_name])
        name_used[base_name] += 1
        start = min(int(hit["qstart"]),int(hit["qend"]))
        end = max(int(hit["qstart"]),int(hit["qend"]))
        row = {"seqname": hit['seqname'], "source": "cmscan",
        "feature": "transcript", "start": start,
        "end":end, "score": hit["evalue"],
        "strand": hit["strand"],
        "frame": ".", "attribute":"ID="+unique_name+";name="+hit['rna_name']
            +";evalue="+hit['evalue']+";rfam="+hit['rfam']}
        rows.append(row)
    gff = pd.DataFrame(rows, columns = ["seqname", "source",
        "feature", "start", "end", "score", "strand", 
        "frame", "attribute"])
    return gff

def split_genome(args, confs, tmpDir, stepDir):
    output_dir = args["data_dir"] + "/genome_parts"
    if not "genome" in args:
        print("Cannot run annotation without a path to a genome.")
        return False
    fasta_path = os.path.abspath(args["genome"])
    n = 100
    #creating shortcut to genome fasta
    runCommand("ln -s " + fasta_path + " " + args["genome_link"])
    fasta_path = args["genome_link"]
    print("Reading input fasta")
    seqs = readSeqsFromFasta(fasta_path)
    total_length = sum([len(entry[1]) for entry in seqs])
    print("Total length:" + str(total_length))
    max_length = int(total_length / n)
    
    current_length = 0
    part = []
    parts = []
    print("Spliting parts of fasta")
    cumulative = 0
    parts_done = 0
    for seq in seqs:
        if n > parts_done:
            max_length = int((total_length-cumulative) / (n-parts_done))
        if ((current_length >= max_length)):
            parts.append(part)
            parts_done += 1
            part = []
            current_length = 0
        part.append(shortFastaHeader(seq))
        cumulative += len(seq[1])
        current_length += len(seq[1])
    if len(part) > 0:
        parts.append(part)
    
    file_names = [output_dir + "/" + str(i) + ".fasta" for i in range(len(parts))]
    runCommand("mkdir " + output_dir)
    print("Writing fasta files")
    for i in range(len(parts)):
        writeFastaSeqs(parts[i], file_names[i])
    return True

def run_infernal(args, confs, tmpDir, stepDir):
    output_dir = args["data_dir"] + "/genome_parts"
    if print_progress(output_dir):
        return True
    else:
        start(output_dir, confs["cmscan"], confs["rfam_cm"], confs["threads"])
        return print_progress(output_dir)

def merge_infernal_outs(args, confs, tmpDir, stepDir):
    output_dir = args["data_dir"] + "/genome_parts"
    output_file = tmpDir + "/infernal"
    success = get_infernal_output(output_dir,output_file)
    return success

def parse_infernal(args, confs, tmpDir, stepDir):
    infernal_tsv = stepDir["merge_infernal_outs"] + "/infernal.tsv"
    gff = read_infernal_output(infernal_tsv)
    print("Writing .gff file")
    gff.to_csv(tmpDir + "/rfam_annotation_genome.gff", sep="\t", index=False, header=False)
    return True

def analyze(args, confs, tmpDir, stepDir):
    print("Analyzing cmscan annotation")
    annotation = pd.read_csv(stepDir["parse_infernal"] + "/rfam_annotation_genome.gff", sep="\t", header=None)
    annotation.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    rfams = dict()
    n_families = 0
    ids = set()
    exons = 0
    for index, row in annotation.iterrows():
        attrs = get_gff_attributes(row["attribute"])
        rfam = attrs["rfam"]
        if row["feature"] == "noncoding_exon":
            exons += 1
        if rfam in rfams:
            rfams[rfam] += 1
        else:
            rfams[rfam] = 1
            n_families += 1
        ids.add(attrs["transcript_id"])
    with open(tmpDir + "/stats.txt", 'w') as stream: 
        stream.write(str(len(ids)) + " transcript_id\n")
        stream.write(str(exons) + " exons\n")
        stream.write(str(n_families) + " families\n")

    pairs = []
    with open(tmpDir + "/family_sizes.tsv", 'w') as stream:
        for key in rfams.keys():
            pairs.append((key, rfams[key]))
            stream.write(key + "\t" + str(rfams[key]) + "\n")
    pairs.sort(key=lambda x: x[1])
    sorted_keys = [x[0] for x in pairs]
    sorted_vals = [x[1] for x in pairs]
    import matplotlib as mpl
    mpl.rcParams['ytick.labelsize'] = 6
    import matplotlib.pyplot as plt
    plt.figure(figsize=(5, len(sorted_keys)/10))
    #plt.xticks(fontsize=8)
    plt.barh([4*x for x in range(len(sorted_keys))], sorted_vals, height=3.0,align='center', tick_label=sorted_keys)
    for i, v in enumerate(sorted_vals):
        plt.text(v, 4*i, " " + str(v), fontsize=6, color='blue', va='center', fontweight='bold')
    #lt.xticks(range(len(sorted_keys)), sorted_keys)
    plt.savefig(args, confs, tmpDir + "/family_sizes.png", dpi = 300, format='png')
    
    return True