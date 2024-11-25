# -*- coding: utf-8 -*-

import pandas as pd
from tqdm import tqdm
import numpy as np
import math
import os

'''
usage:
    analyze_lncrna_expression.py <counts>.tsv <output_path>
'''

def is_specific(sample_count, other_samples, counts):
    for other_sample in other_samples:
        if counts[other_sample] > 0:
            if sample_count / counts[other_sample] < 5:
                return False
    return True

def calc_tissue_specificity(counts, groups, sample_names):
    higher_tpm = sample_names[0]
    lower_tpm = sample_names[0]
    for sample in sample_names:
        if counts[sample] > counts[higher_tpm]:
            higher_tpm = sample
        if counts[sample] < counts[lower_tpm]:
            lower_tpm = sample
    lower_tpm = counts[lower_tpm]

    male_mean = 0.0
    males = 0
    female_mean = 0.0
    females = 0

    relevant = False
    tissues_expressed = []
    sample_counts = {}
    for sample in sample_names:
        if counts[sample] >= 1.0:
            relevant = True
            tissues_expressed.append(sample)
            sample_counts[sample] = counts[sample]
        if "Male" in sample:
            male_mean += counts[sample]
            males += 1
        elif "Female" in sample:
            female_mean += counts[sample]
            females += 1

    sex_specific = np.nan
    if male_mean > 0 and female_mean == 0.0:
        sex_specific = "Male"
    elif female_mean > 0 and male_mean == 0.0:
        sex_specific = "Female"

    male_mean = male_mean / males if males > 0 else 0
    female_mean = female_mean / females if females > 0 else 0
    raw_fc = np.inf
    if male_mean > 0 and female_mean > 0:
        raw_fc = female_mean/male_mean
    elif (not male_mean > 0) and (not female_mean > 0):
        raw_fc = 1.0
    fold_change = abs(math.log(raw_fc,2))

    if not relevant:
        return 'Not Expressed', np.nan, lower_tpm, male_mean, female_mean, fold_change, len(tissues_expressed), sex_specific

    other_samples = sample_names[:]
    
    #other_samples.remove(higher_tpm)
    #specific = is_specific(counts[higher_tpm], other_samples, counts)
    group_name = higher_tpm.split('_')[0]
    for s in groups[group_name]:
        other_samples.remove(s)
    
    specifics_bool = [is_specific(counts[s], other_samples, counts)
                       for s in groups[group_name]]
    n_specifics = sum(specifics_bool)
    tissue_specific = group_name if n_specifics > 0 else None
    
    if tissue_specific != None:
        return 'Tissue Specific', group_name, lower_tpm, male_mean, female_mean, fold_change, len(tissues_expressed), sex_specific
    else:
        expressed_in_all = True
        for sample in sample_names:
            if counts[sample] < 1.0:
                expressed_in_all = False
                break
        if expressed_in_all:
            return 'Housekeeping', np.nan, lower_tpm, male_mean, female_mean, fold_change, len(tissues_expressed), sex_specific
        else:
            return 'Mixed Expression', np.nan, lower_tpm, male_mean, female_mean, fold_change, len(tissues_expressed), sex_specific

'''def get_max_expression(name):
    lines = df[df.index == name]
    maxes = lines['max_expression'].tolist()[0]
    return maxes'''

def find_lncrna_classes(df_path, output_path):
    #df_path = sys.argv[1]
    #output_path = sys.argv[2]

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    print("Reading counts")
    df = pd.read_csv(df_path, sep=',', header=0, index_col='gene')
    
    groups = {'Heart': ['Heart_Male'], 
            'Brain': ['Brain_Female', 'Brain_Male'], 
            'Liver': ['Liver_Female', 'Liver_Male'], 
            'Skin': ['Skin_Female', 'Skin_Male'],
            'Muscle':['Muscle_Female'],
            'Gonad': ['Gonad_Female', 'Gonad_Male'], 
            'Lung': ['Lung_Female', 'Lung_Male'],
            'Kidney': ['Kidney_Female', 'Kidney_Male']}

    sample_names = list(df.columns)

    print(groups)
    print(sample_names)
    print(df)

    df['max_expression'] = df.max(axis=1)

    print("Analyzing expression of each lncrna")
    bar = tqdm(total=len(df))
    colnames_vec = ['Classification', 'Specific_Tissue', 'Lowest_Expression_Count', 'Male_Mean_Expression', 
        'Female_Mean_Expression', 'Log_FC_Expression', 'Samples_With_Expression', 'Expressed_Only_In_This_Sex']
    results = []
    for index, row in tqdm(df.iterrows()):
        new_row = dict(row)
        new_row['gene'] = index
        classes = calc_tissue_specificity(row, groups, sample_names)
        for i, name in enumerate(colnames_vec):
            new_row[name] = classes[i]
        results.append(new_row)
        bar.update(1)
    bar.close()

    results_df = pd.DataFrame(results)
    results_df_path = output_path + '/ncrna_classification.tsv'
    results_df.to_csv(results_df_path, sep= '\t', index=False)

    expressed = results_df[results_df['Classification'] != 'Not Expressed']
    expressed.set_index('gene', inplace=True)
    for key in colnames_vec:
        del expressed[key]
    del expressed['max_expression']
    expressed_df_path = output_path + '/ncrna_expressed.tsv'
    expressed.to_csv(expressed_df_path, sep=',')

    return results_df_path, expressed_df_path
