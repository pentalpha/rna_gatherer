import os
import subprocess
import pandas as pd
import numpy as np
from tqdm import tqdm

def runCommand(cmd, print_cmd=True):
    if print_cmd:
        print("\t> " + cmd)
    process = subprocess.call(cmd, shell=True)
    return process

def align_with_minimap2(genomes, analysis_dir, transcripts_fasta):
    align_results = []
    for family, name, taxid, fasta_db in genomes:
        result_paf = analysis_dir + '/align_to_'+name.replace(' ', '_')+'.paf'
        print('Aligning to', fasta_db)
        created = os.path.exists(result_paf)
        if not created:
            '''index_1 = fasta_db + '.mmi'
            index_2 = os.path.dirname(fasta_db) + '/index.mmi'
            correct_indexes = [i for i in [index_1, index_2] if os.path.exists(i)]
            if len(correct_indexes) == 0:
                print('No index found for', fasta_db)
                index_cmd = 'minimap2 -d ' + index_2 + ' ' + fasta_db
                runCommand(index_cmd)
                correct_indexes.append(index_2)
            
            index_path = correct_indexes[0]'''
            cmd = " ".join(["minimap2",
                "-x sr -k19 -w9 -c -t", str(8),
                fasta_db, transcripts_fasta,
                ">", result_paf])
            code = runCommand(cmd)

            cmd2 = " ".join(["minimap2",
                "-x splice -c -t", str(8),
                fasta_db, transcripts_fasta,
                ">>", result_paf])
            code = runCommand(cmd2)
            
        align_results.append([family, name, taxid, result_paf])
    return align_results

def load_paf(paf_path):
    filters = {'is_similar':  {'identity': 0.8, 'qcovs': 0.7},
               'is_homology': {'identity': 0.9, 'qcovs': 0.9}}
    
    lines = []
    header = ["qseqid","qseq_len","qstart","qend","strand",
            "sseqid","sseq_len","sstart","send","matchs",
            "block_len","quality"]
    for rawline in open(paf_path):
        cells = rawline.rstrip('\n').split('\t')
        row = {header[i]: cells[i] for i in range(len(header))}
        as_value = None
        for x in cells:
            if x.startswith('AS:i:'):
                as_value = int(x.replace('AS:i:', ''))
                break
        row['AS'] = as_value
        lines.append(row)
    
    df = pd.DataFrame(lines)
    df = df.astype({"qstart": 'int32', "qend": 'int32', "qseq_len": "int32",
            "sstart": 'int32', "send": 'int32', "sseq_len": "int32",
            "quality": 'int32', "block_len": "int32", "matchs": "int32"})
    df["id"] = np.arange(len(df))
    df["qcovs"] = df.apply(
            lambda row: (row["qend"]-row["qstart"]) / row["qseq_len"], axis=1)
    df["identity"] = df.apply(
            lambda row: row["matchs"] / row["block_len"], axis=1)
    for filter_name, th_dict in filters.items():
        ths = [(filter_name+'_min_'+key, key, val) for key, val in th_dict.items()]
        filter_names = []
        for filter_col, attr, min_val in ths:
            print(filter_col, attr, min_val)
            df[filter_col] = df[attr] >= min_val
            filter_names.append(filter_col)
        df[filter_name] = df[filter_names[0]] & df[filter_names[1]]
        del df[filter_names[0]]
        del df[filter_names[1]]
    
    lines2 = []
    for qseqid, alignment_lines in tqdm(df.groupby('qseqid')):
        #alignment_lines = alignment_lines[alignment_lines['qcovs'] >= min_cov]
        alignment_lines.sort_values(['is_homology', 'is_similar', 'identity'])
        top_row = dict(alignment_lines.iloc[-1])
        lines2.append(top_row)
    df = pd.DataFrame(lines2)
    
    return df