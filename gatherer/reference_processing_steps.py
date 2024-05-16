import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from gatherer.bioinfo import *
from gatherer.netutils import *
from gatherer.util import *
import hashlib
import requests
import time
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

def prepare_ref_annotation(args, confs, tmpDir, stepDir):
    annotation_lines = 0
    no_parent_lines = 0
    if "reference_gff" in args:
        ref = args["reference_gff"]
        prepared_gff = tmpDir + "/reference.gff"
        with open(prepared_gff, 'w') as out_stream:
            with open(ref, 'r') as in_stream:
                for line in in_stream:
                    if not (line.startswith("!") or line.startswith("#")):
                        annotation_lines += 1
                        cells = line.rstrip("\n").split("\t")
                        if (not "Parent=" in cells[8]) and (not "parent=" in cells[8]):
                            no_parent_lines += 1
                            cells[8] = cells[8].rstrip(";") + ";ref_db=" + cells[1]
                            cells[1] = "reference"
                            out_stream.write("\t".join(cells)+"\n")
    print(annotation_lines, " annotation lines in reference gff.")
    print(no_parent_lines, " lines in processed reference.")
    return True

def map_to_genome(args, confs, tmpDir, stepDir):
    if "reference_fasta" in args:
        to_map_path = args["reference_fasta"]

        genome_alignment = tmpDir + "/reference_mapped.paf"
        cmd = " ".join([confs["minimap2"],
            "-x splice:hq -uf -t", str(confs["threads"]),
            args["genome_index"], to_map_path, 
            ">", genome_alignment])
        code = runCommand(cmd)
        if code != 0:

            return False
        output_gff = tmpDir + "/reference_mapped.gff"
        annotated_fasta = tmpDir + "/reference_mapped.fasta"
        mapped_fasta = minimap_annotation(genome_alignment, output_gff,
                        annotated_fasta, source="reference_mapping", 
                        query_file=to_map_path, mol_type=None)
    return True