from gatherer.bioinfo import *
from gatherer.util import *
from gatherer.netutils import *
import pandas as pd
import os
from tqdm import tqdm

def ncrna_alignment_minimap(args, confs, tmpDir, stepDir):
    print("Aligning to ncRNA databases")
    rna_dbs = None
    if "rna_dbs" in confs:
        if len(confs["rna_dbs"]) > 0:
            rna_dbs = confs["rna_dbs"]
    if rna_dbs == None:
        print("No ncRNA databases configured, not aligning.")
        return True
    for db_name, db_path in tqdm(rna_dbs.items()):
        print("\tAligning to " + db_name + " (" + db_path + ")")
        paf_output = tmpDir + "/"+db_name+"-mapping.paf"
        cmd = " ".join([confs["minimap2"],
            "-x splice:hq -uf -t", str(confs["threads"]),
            args["genome_index"], db_path, 
            ">", paf_output])
        code = runCommand(cmd)
        if code != 0:
            print("Error: Alignment unsucessful.")
            os.remove(paf_output)
    return True

def ncrna_alignment_parsing(args, confs, tmpDir, stepDir):
    print("Parsing database alignments")
    paf_files = getFilesWith(stepDir["ncrna_alignment_minimap"], "-mapping.paf")
    gffs = []
    fastas = []
    for paf_file in tqdm(paf_files):
        db_name = paf_file.split("/")[-1].split("-")[0]
        print("\tParsing " + db_name)
        genome_alignment = paf_file
        
        output_gff = tmpDir + "/"+db_name+"-annotation.gff"
        annotated_fasta = tmpDir + "/"+db_name+"-in_genome.fasta"
        mapped_fasta = minimap_annotation(genome_alignment, output_gff,
                        annotated_fasta, source="db_alignment", mol_type=None, 
                        db_name=db_name)
        gffs.append(output_gff)
        fastas.append(annotated_fasta)
    
    if len(gffs) > 0:
        all_gff = tmpDir + "/alignment_annotation.gff"
        all_fasta = tmpDir + "/alignment_sequences.fasta"
        join_files_in_one(gffs, all_gff)
        join_files_in_one(fastas, all_fasta)
    return True
