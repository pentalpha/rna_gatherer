import os
from gatherer.util import runCommand

def run_trnascan(args, confs, tmpDir, stepDir):
    genome_path = args["genome_link"]
    tRNAscan = confs["tRNAscan-SE"]
    output_file = tmpDir + "/trna_raw.txt"
    stats_file = tmpDir + "/trna_stats.txt"
    cmd = " ".join([tRNAscan, "-o", output_file, "-m", stats_file, genome_path])
    code = runCommand(cmd)
    return code == 0

def to_gff(trna_line, names):
    seqname = trna_line[0]
    begin = trna_line[2]
    end = trna_line[3]
    posA = min(int(begin),int(end))
    posB = max(int(begin),int(end))
    trna_type = trna_line[4]
    anticodon = trna_line[5]
    score = trna_line[8]
    strand = "+"
    if int(begin) > int(end):
        strand = "-"
    rna_name = "tRNA_"+trna_type
    name = seqname + "_" + rna_name
    if not name in names:
        names[name] = -1
    names[name] += 1
    unique_name = name + "_" + str(names[name])
    attrs = ["ID="+unique_name,"type=tRNA","tRNA_type="+rna_name,"anticodon="+anticodon]
    return "\t".join(
        [seqname,"tRNAscan-SE","transcript",str(posA),str(posB),score,strand,".",
            ";".join(attrs)]
    )

def valid_trna_type(gff_line):
    return not("tRNA_type=tRNA_Undet" in gff_line or "tRNA_type=tRNA_Pseudo" in gff_line)

def parse_trna(args, confs, tmpDir, stepDir):
    raw_trnas = open(stepDir["run_trnascan"] + "/trna_raw.txt", 'r')
    raw_lines = raw_trnas.readlines()
    raw_trnas.close()
    raw_lines = raw_lines[3:]
    lines = [raw_line.split() for raw_line in raw_lines]
    names = {}
    gff_lines = [to_gff(line,names) for line in lines]
    valid_gff = []
    for line in gff_lines:
        if valid_trna_type(line):
            valid_gff.append(line)
    gff_content = "\n".join(valid_gff)
    with open(tmpDir+"/tRNAs.gff", 'w') as stream:
        stream.write(gff_content)
    print("tRNAs of type 'Undet' or 'Pseudo' filtered: " 
            + str(len(gff_lines) - len(valid_gff)))
    return True
