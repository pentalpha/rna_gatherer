"""

Author: Pitágoras Alves (github.com/pentalpha)
"""

import gzip
import random
import sys
import obonet
import numpy as np
import argparse
import multiprocessing
import pandas as pd

from os import path, mkdir

from gatherer.netutils import QuickGoRetriever

'''
Loading configurations and command line arguments
'''
from config import configs, require_files
from gatherer.util import chunks, runCommand
mandatory_files = ["go_obo"]
require_files(mandatory_files)

#set configuration values
confs = {}
for conf in configs:
    confs[conf] = configs[conf]

def getArgs():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-txid", "--taxon-id", required=True,
        help=("Base taxon to download protein sequences and annotations"))
    ap.add_argument("-g", "--genome", required=True,
        help=("Path to genome fasta"))
    ap.add_argument("-o", "--output-dir", required=True, help=("Output directory."))
    
    default_threads = max(2, multiprocessing.cpu_count()-1)
    ap.add_argument("-p", "--processes", required=False,
        default=default_threads, help=("CPUs to use. Default: " + str(default_threads)+"."))
    
    return vars(ap.parse_args())

if __name__ == '__main__':
    cmdArgs = getArgs()
    genome_path = cmdArgs["genome"]
    parent_taxon_id = cmdArgs["taxon_id"]
    output_dir = cmdArgs["output_dir"]
    go_path = confs["go_obo"]

    threads = int(cmdArgs["processes"])
    if not path.exists(output_dir):
        mkdir(output_dir)

    output_proteins = path.join(output_dir, 'proteins.fasta.gz')
    if not path.exists(output_proteins):
        cmd = "wget -O OUTPUT https://rest.uniprot.org/uniprotkb/stream\?compressed\=true\&format\=fasta\&includeIsoform\=true\&query\=%28%28taxonomy_id%3A27723%29%29"
        cmd = cmd.replace('27723', parent_taxon_id).replace('OUTPUT', output_proteins)
        runCommand(cmd)

    protein_list = []
    for rawline in gzip.open(output_proteins, 'rt'):
        if rawline.startswith('>'):
            protein_list.append(rawline.lstrip('>').split('|')[1])
    
    proteins_downloaded_path = path.join(output_dir, 'proteins.txt')
    open(proteins_downloaded_path, 'w').write("\n".join(protein_list))

    mpi_db_path = genome_path.replace('.fasta', '.mpi')
    if not path.exists(mpi_db_path):
        miniprot_cmd = [
            confs['miniprot'], '-t8', '-d', mpi_db_path, genome_path
        ]
        code = runCommand(' '.join(miniprot_cmd))

    annotations_path = path.join(output_dir, 'annotations.tsv')
    #protein_list = protein_list[:700]
    annotations_paths = []
    if not path.exists(annotations_path):
        attempts = 3
        for i in range(1, attempts+1):
            waittime = i
            chunk_size = int(240 / i)
            random.shuffle(protein_list)
            id_chunks = list(chunks(protein_list, chunk_size))
            
            print(len(protein_list), 'proteins divided into', len(id_chunks), 'chunks')
            quickgo_retriever = QuickGoRetriever('https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch', 
                parent_taxon_id, waittime)
            annotation_lines = []
            annotation_attempt_path = path.join(output_dir, 'annotations_'+str(i)+'.tsv')
            with multiprocessing.Pool(3) as pool:
                annotation_chunks = pool.map(quickgo_retriever.retrieve_quickgo_prot_annotations, id_chunks)
                
                annotations_file = open(annotation_attempt_path, 'w')
                for chunk in annotation_chunks:
                    for line in chunk:
                        rawline = '\t'.join(line)+'\n'
                        annotations_file.write(rawline)
                annotations_file.close()
                annotations_paths.append(annotation_attempt_path)

        annotations = set()
        for p in annotations_paths:
            for rawline in open(p, 'r'):
                cells = rawline.rstrip('\n').split('\t')
                if len(cells) == 4:
                    annotations.add((cells[0], cells[1], cells[2], cells[3]))
            print(len(annotations), 'protein annotations loaded')
        annotations = sorted(annotations)
        single_ann_output = open(annotations_path, 'w')
        #single_ann_output.write('uniprot\ttaxon\tgo\taspect\n')
        for cells in annotations:
            single_ann_output.write('\t'.join(cells)+'\n')
        single_ann_output.close()
    
    annotations = pd.read_csv(annotations_path, sep='\t', header=None)
    annotations.rename(columns={0 :'uniprot', 1: 'taxon', 2: 'go', 3: 'aspect'}, 
        inplace=True)
    print(len(annotations), 'protein annotations loaded')
    print(annotations.head())
    annotated_proteins = set(annotations['uniprot'].tolist())
    ann_proteins_path = path.join(output_dir, 'annotated_proteins.txt')
    open(ann_proteins_path, 'w').write("\n".join(annotated_proteins))

    annotated_proteins_fasta = path.join(output_dir, 'annotated_proteins.fasta')
    if not path.exists(annotated_proteins_fasta):
        print_protein = False
        annotated_proteins_fasta_file = open(annotated_proteins_fasta, 'w')
        for rawline in gzip.open(output_proteins, 'rt'):
            if rawline.startswith('>'):
                protid = rawline.lstrip('>').rstrip('\n').split('|')[1]
                print_protein = protid in annotated_proteins
            if print_protein:
                annotated_proteins_fasta_file.write(rawline)
        annotated_proteins_fasta_file.close()

    protein_alignment_path = path.join(output_dir, 'annotated_proteins.gff')
    #if not path.exists(protein_alignment_path):
    miniprot_cmd2 = [
        confs['miniprot'], '-t', str(int(threads/2)), 
        '--gff', genome_path,
        annotated_proteins_fasta, '>', protein_alignment_path
    ]
    code = runCommand(' '.join(miniprot_cmd2))