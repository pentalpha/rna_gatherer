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
import gffutils
import pyfaidx
from tqdm import tqdm
from os import path, mkdir

from gatherer.bioinfo import cov_id_from_minimap_line
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

def create_mrna_annotation(genome_path, parent_taxon_id, output_dir, 
        threads, min_coverage, min_identity):
    
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
    if not path.exists(protein_alignment_path):
        miniprot_cmd2 = [
            confs['miniprot'], '-t', str(int(threads/2)), 
            '--gff', genome_path,
            annotated_proteins_fasta, '>', protein_alignment_path
        ]
        code = runCommand(' '.join(miniprot_cmd2))

    protein_alignment_filtered_path = path.join(output_dir, 'annotated_proteins.filtered.gff')
    paf_id = None
    paf_cov = None

    protein_alignment_filtered_file = open(protein_alignment_filtered_path, 'w')
    mRNA_count = 0
    mRNA_kept = 0
    mRNA_small_cov = 0
    mRNA_small_id = 0
    last_mrna_kept = False
    for rawline in open(protein_alignment_path, 'r'):
        
        if rawline.startswith('##PAF\t'):
            rawline = rawline.replace('##PAF\t', '')
            cells = rawline.rstrip('\n').split('\t')
            paf_cov, paf_id = cov_id_from_minimap_line(cells)
        else:
            cells = rawline.rstrip('\n').split('\t')
            if len(cells) >= 9:
                if cells[2] == 'mRNA':
                    mRNA_count += 1
                    if paf_cov and paf_id:
                        cells[8] = ('PAFIdentity='+str(round(paf_id, 4))
                                    +';PAFCoverage='+str(round(paf_cov, 4))
                                    +';'+cells[8])
                    info_fields = {part.split('=')[0]: part.split('=')[1] for part in cells[8].split(';')}
                    if (float(info_fields['Identity']) >= min_identity
                        and float(info_fields['PAFCoverage']) >= min_coverage):
                        protein_alignment_filtered_file.write('\t'.join(cells)+'\n')
                        mRNA_kept += 1
                        last_mrna_kept = True

                    elif float(info_fields['Identity']) < min_identity:
                        mRNA_small_id += 1
                        last_mrna_kept = False
                    elif float(info_fields['PAFCoverage']) < min_coverage:
                        mRNA_small_cov += 1
                        last_mrna_kept = False
                    else:
                        last_mrna_kept = False
                else:
                    if last_mrna_kept:
                        protein_alignment_filtered_file.write(rawline)
            else:
                protein_alignment_filtered_file.write(rawline)

            paf_cov = None
            paf_id = None
    protein_alignment_filtered_file.close()
    
    print('mRNA_count', mRNA_count)
    print('mRNA_kept', mRNA_kept)
    print('mRNA_kept/mRNA_count', mRNA_kept/mRNA_count)
    print('mRNA_small_id', mRNA_small_id)
    print('mRNA_small_cov', mRNA_small_cov)

    protein_alignment_filtered_fasta_path = path.join(output_dir, 
        'annotated_proteins.filtered.mRNA.fasta')
    if not path.exists(protein_alignment_filtered_fasta_path):
        print('Loading', protein_alignment_filtered_path)
        db = gffutils.create_db(protein_alignment_filtered_path, 
            protein_alignment_filtered_path+'.db', 
            force=True, )
        print('Loading', genome_path)
        fasta = pyfaidx.Fasta(genome_path)
        print('Loaded inputs, creating', protein_alignment_filtered_fasta_path)

        print('Collecting mRNA names')
        cds_by_parent = {}
        parent_to_protein_name = {}
        targets = set()
        for mRNA in tqdm(db.features_of_type('mRNA', order_by='start')):
            target = mRNA.attributes['Target'][0]
            targets.add(target.split(' ')[0])
            mrna_id = mRNA.attributes['ID'][0]
            #print(mrna_id, target)
            parent_to_protein_name[mrna_id] = target
            cds_by_parent[mrna_id] = ''
        
        print('Collecting cds sequences')
        for cds in tqdm(db.features_of_type('CDS', order_by='start')):
            cds_parent = cds.attributes['Parent'][0]
            cds_seq = cds.sequence(fasta)
            if cds_parent in cds_by_parent:
                cds_by_parent[cds_parent] += cds_seq
        
        print('Writing mRNA sequences')
        fasta_output = open(protein_alignment_filtered_fasta_path, 'w')
        for mrna_id, protein_name in parent_to_protein_name.items():
            coding_sequence = cds_by_parent[mrna_id]
            fasta_output.write('>'+protein_name+' '+mrna_id+'\n')
            fasta_output.write(coding_sequence+'\n')
        
        fasta_output.close()

        print(len(parent_to_protein_name.keys()), 'mRNAs writen')
        print('Of', len(targets), 'proteins')
    
    return (protein_alignment_filtered_fasta_path, 
        protein_alignment_filtered_path,
        annotations_path)

def create_all_rnas_fasta(output_dir, ncrna_fasta_path, mrna_fasta_path):
    output_path = output_dir + '/coding_and_noncoding.fasta'
    output = open(output_path, 'w')
    for f in [ncrna_fasta_path, mrna_fasta_path]:
        input_f = open(f, 'r')
        for line in input_f:
            if line.startswith(">"):
                output.write(line.rstrip("\n")+'\n')
            else:
                output.write(line.rstrip("\n")+'\n')
        input_f.close()
    output.close()

    index_path = output_path+'.salmon_idx'
    if not path.exists(index_path):
        index_cmd = ['salmon', 'index', '-t', output_path, 
            '-i', index_path, '-k 31']
        print(index_cmd)
        runCommand(' '.join(index_cmd))
    return output_path, index_path

def run_salmon(output_dir, index_path, fastqs_table, procs):
    samples = {}
    fastqs_dir = path.dirname(fastqs_table)
    for rawline in open(fastqs_table, 'r'):
        cells = rawline.rstrip('\n').split('\t')
        if len(cells) == 2:
            gname = cells[0]
            r1, r2 = cells[1].split(',')
            samples[gname] = (fastqs_dir+'/'+r1, fastqs_dir+'/'+r2)
            print(gname, samples[gname])
        else:
            print('not two columns:')
            print(rawline)
            quit(1)
    
    salmon_quant_dir = output_dir + '/salmon_quant'
    if not path.exists(salmon_quant_dir):
        mkdir(salmon_quant_dir)
    salmon_results = {}
    print('Running salmon quant for samples')
    for sample_name, fastqs in tqdm(samples.items()):
        print(sample_name)
        #salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq -o transcripts_quant
        sample_quant_dir = salmon_quant_dir+'/'+sample_name
        sample_quant_file = sample_quant_dir+'/quant.sf'
        stdout_path = sample_quant_dir+'.stdout'
        stderr_path = sample_quant_dir+'.stderr'
        if not path.exists(sample_quant_file):
            quant_cmd = ['salmon', 'quant', '-i', index_path, '-l IU', 
                        '-1', fastqs[0], '-2', fastqs[1], '-o', sample_quant_dir,
                        '-p', str(procs), '2>'+stderr_path, '1>'+stdout_path]
            print(quant_cmd)
            runCommand(' '.join(quant_cmd))
        else:
            print(sample_quant_file, 'already created')
        salmon_results[sample_name] = sample_quant_file

def getArgs():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-txid", "--taxon-id", required=True,
        help=("Base taxon to download protein sequences and annotations"))
    ap.add_argument("-g", "--genome", required=True,
        help=("Path to genome fasta"))
    ap.add_argument("-nc", "--nc-rna-fasta", required=True,
        help=("Path to fasta file with ncRNA sequences"))
    ap.add_argument("-fq", "--fastq-table", required=True,
        help=("Path to .TSV file where first column is the tissue name and"
              +" second column is a ',' separated list of .fastq.gz files."))
    ap.add_argument("-o", "--output-dir", required=True, help=("Output directory."))
    ap.add_argument("-cov", "--coverage", required=True, help=("Minimum protein coverage"))
    ap.add_argument("-id", "--identity", required=True, help=("Minimum protein identity"))
    
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
    min_coverage = float(cmdArgs["coverage"])
    min_identity = float(cmdArgs["identity"])
    fastq_table_path = cmdArgs['fastq_table']
    ncrna_fasta_path = cmdArgs['nc_rna_fasta']
    if not path.exists(output_dir):
        mkdir(output_dir)

    protein_annot = create_mrna_annotation(genome_path, 
        parent_taxon_id, output_dir, threads, min_coverage, 
        min_identity)
    
    mrna_fasta_path = protein_annot[0]
    protein_alignment_gff_path = protein_annot[1]
    protein_go_annotations_path = protein_annot[2]

    all_rnas_fasta, all_rnas_index = create_all_rnas_fasta(output_dir, ncrna_fasta_path, 
        mrna_fasta_path)

    count_reads_table = run_salmon(output_dir, all_rnas_index, 
        fastq_table_path, threads)

    