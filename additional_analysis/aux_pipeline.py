import subprocess
import os
from sys import getsizeof
from tqdm import tqdm
import sys
import pandas as pd
import obonet

from analyze_lncrna_expression import find_lncrna_classes

def runCommand(cmd, print_cmd=True):
    if print_cmd:
        print("\t> " + cmd)
    process = subprocess.call(cmd, shell=True)
    return process

def get_go_list(p):
    return set([l.rstrip("\n").split()[0] for l in open(p,'r').readlines()])

def write_names_to_file(names, f):
    with open(f,'w') as stream:
        for name in names:
            stream.write(name+"\n")
    return names

def read_id2gos(id2gos_file, id2gos = {}):
    with open(id2gos_file, 'r') as stream:
        print("Reading " + id2gos_file)
        for line in stream:
            cells = line.rstrip("\n").split("\t")
            transcript_name = cells[0]
            go_names = cells[1].split(";")
            if not transcript_name in id2gos:
                id2gos[transcript_name] = set()
            for go_name in go_names:
                id2gos[transcript_name].add(go_name)
    return id2gos

def genes_annotations_in_set(gene_id, associations, interest_list, descriptions, correct_id):
    if gene_id in associations:
        ann = associations[gene_id]
        interest_ann = ann.intersection(interest_list)
        return [go + ": " + descriptions[correct_id[go]] for go in interest_ann]
        #return list(interest_ann)
    else:
        return []

def write_gene_lists(df2, output_path, descriptions, correct_id):
    print("Writing gene lists")
    housekeeping_df = df2[df2['Classification'] == 'Housekeeping'][['Name', 'Lowest_Expression_Count', 'Samples_With_Expression']]
    housekeeping_df.sort_values(by=['Samples_With_Expression', 'Lowest_Expression_Count'], 
                                ascending = [False, False], inplace=True)
    housekeeping_df.to_csv(output_path +"/housekeeping.tsv", sep="\t")
    write_names_to_file(housekeeping_df['Name'].tolist(), 
                        output_path +"/housekeeping.txt")

    sex_diff_df = df2[df2['Diff_Sex'] == True][['Name', 'Classification', 'Specific_Tissue', 'Male_Mean_Expression', 'Female_Mean_Expression', 'Log_FC_Expression', 'Number_Of_DE_Packages']]
    sex_diff_df.sort_values(by=['Number_Of_DE_Packages', 'Log_FC_Expression'], 
                                ascending = [False, False], inplace=True)
    sex_diff_df.to_csv(output_path +"/sex_diff.tsv", sep="\t")
    write_names_to_file(sex_diff_df['Name'].tolist(), 
                        output_path +"/sex_diff.txt")

    growth_df = df2[df2['Involved_in_Growth'] == True][['Name', 'Classification', 'Specific_Tissue', 'Growth_Functions', 'Growth_Functions_Percent']]
    growth_df['Functions'] = growth_df.apply(lambda row: str(genes_annotations_in_set(row['Name'], 
                                                            associations, growth_names, descriptions, correct_id)),
                                            axis=1) 
    growth_df.sort_values(by=['Growth_Functions', 'Growth_Functions_Percent'], 
                                ascending = [False, False], inplace=True)
    growth_df.to_csv(output_path +"/involved_in_growth.tsv", sep="\t")
    write_names_to_file(growth_df['Name'].tolist(), 
                        output_path +"/involved_in_growth.txt")
    growth_hk_df = growth_df[growth_df['Classification'] == 'Housekeeping']
    growth_hk_df.to_csv(output_path +"/involved_in_growth-housekeeping.tsv", sep="\t")
    write_names_to_file(growth_hk_df['Name'].tolist(), 
                        output_path +"/involved_in_growth-housekeeping.txt")

    maturation_df = df2[df2['Involved_in_Maturation'] == True][['Name', 'Classification', 'Specific_Tissue', 'Maturation_Functions', 'Maturation_Functions_Percent', 'Diff_Sex', 'Male_Mean_Expression', 'Female_Mean_Expression']]
    maturation_df['Functions'] = maturation_df.apply(lambda row: str(genes_annotations_in_set(row['Name'], 
                                                            associations, maturation_names, descriptions, correct_id)),
                                            axis=1)
    maturation_df.sort_values(by=['Maturation_Functions', 'Maturation_Functions_Percent'], 
                                ascending = [False, False], inplace=True)
    maturation_df.to_csv(output_path +"/involved_in_maturation.tsv", sep="\t")
    write_names_to_file(maturation_df['Name'].tolist(), 
                        output_path +"/involved_in_maturation.txt")

def make_tissue_summary(df2_path, output_path):
    print("Reading results to print summary")
    df2 = pd.read_csv(df2_path, sep="\t", header=0)
    tissue_data = []
    summary_path = output_path+"/tissue_sumarry.tsv"
    name_lists = {}
    for tissue_name, tissue_df in tqdm(
            list(df2[df2['Classification'] == 'Tissue Specific'].groupby('Specific_Tissue'))):
        print(tissue_name)
        names = set(tissue_df['Name'].tolist())
        sex_diff = set(tissue_df[tissue_df['Diff_Sex'] == True]['Name'].tolist())
        growth = set(tissue_df[tissue_df['Involved_in_Growth'] == True]['Name'].tolist())
        maturation = set(tissue_df[tissue_df['Involved_in_Maturation'] == True]['Name'].tolist())
        tissue_data.append(
            {'Tissue Name': tissue_name, 
            'Tissue Specific': len(tissue_df),
            'Differentially Expressed by Sex': len(sex_diff),
            'Growth Genes': len(growth),
            'Maturation Genes': len(maturation)
            })
        
        name_list_prefix = output_path +"/"+tissue_name+"."
        print("\tWriting genes")
        names = write_names_to_file(names, name_list_prefix+"tissue.txt")

        if "Skin" in tissue_name:
            male_df = tissue_df[tissue_df['Male_Mean_Expression'] > tissue_df['Female_Mean_Expression']]
            female_df = tissue_df[tissue_df['Male_Mean_Expression'] < tissue_df['Female_Mean_Expression']]
            male_expressed = male_df['Name'].tolist()
            female_expressed = female_df['Name'].tolist()
            female_skin_de = sex_diff.intersection(female_expressed)
            male_skin_de = sex_diff.intersection(male_expressed)
            male_specific = write_names_to_file(sorted(male_skin_de), name_list_prefix+"male_expressed_de.txt")
            female_specific = write_names_to_file(sorted(female_skin_de), name_list_prefix+"female_expressed_de.txt")
        print("\tDone")
        #name_lists[tissue_name] = names
    print("Calculating for Mixed Expression")
    print("\Mixed Expression")
    tissue_data.append({'Tissue Name': "Mixed Expression", 
            'Tissue Specific': len(df2[df2['Classification'] == 'Mixed Expression']),
            'Differentially Expressed by Sex': 
                len(df2[df2['Classification'] == 'Mixed Expression'][df2['Diff_Sex'] == True]['Name'].tolist()),
            'Growth Genes': 
                len(df2[df2['Classification'] == 'Mixed Expression'][df2['Involved_in_Growth'] == True]['Name'].tolist()),
            'Maturation Genes': 
                len(df2[df2['Classification'] == 'Mixed Expression'][df2['Involved_in_Maturation'] == True]['Name'].tolist())})
    print("\tHouse keeping")
    tissue_data.append({'Tissue Name': "Housekeeping", 
            'Tissue Specific': len(df2[df2['Classification'] == 'Housekeeping']),
            'Differentially Expressed by Sex': 
                len(df2[df2['Classification'] == 'Housekeeping'][df2['Diff_Sex'] == True]['Name'].tolist()),
            'Growth Genes': 
                len(df2[df2['Classification'] == 'Housekeeping'][df2['Involved_in_Growth'] == True]['Name'].tolist()),
            'Maturation Genes': 
                len(df2[df2['Classification'] == 'Housekeeping'][df2['Involved_in_Maturation'] == True]['Name'].tolist())})
    print("\tNot expressed")
    tissue_data.append({'Tissue Name': "Not Expressed", 
            'Tissue Specific': len(df2[df2['Classification'] == 'Not Expressed']),
            'Differentially Expressed by Sex': 
                len(df2[df2['Classification'] == 'Not Expressed'][df2['Diff_Sex'] == True]['Name'].tolist()),
            'Growth Genes': 
                len(df2[df2['Classification'] == 'Not Expressed'][df2['Involved_in_Growth'] == True]['Name'].tolist()),
            'Maturation Genes': 
                len(df2[df2['Classification'] == 'Not Expressed'][df2['Involved_in_Maturation'] == True]['Name'].tolist())})
    print("Writing dataframe")
    summary_df = pd.DataFrame(data=tissue_data,columns=tissue_data[0].keys())
    summary_df.to_csv(summary_path, sep='\t', index=False, header=True)
    print(summary_df)

if __name__ == '__main__':
    gigas_dir = sys.argv[1]
    analysis_dir = gigas_dir + '/analysis_for_paper'
    samples_tpm = gigas_dir + '/go_predict/samples_tpm.tsv'
    results_df_path = analysis_dir + '/ncrna_classification.tsv'
    expressed_df_path = analysis_dir + '/ncrna_expressed.tsv'
    '''filter_ncrna_cmd = ['/usr/bin/Rscript', '--vanilla', 'filter_ncrna_counts.R', samples_tpm]
    runCommand(' '.join(filter_ncrna_cmd))

    results_df_path, expressed_df_path = find_lncrna_classes('counts_noduplicates_noprot.csv', analysis_dir)'''

    print('Performing DE')
    de_file_path = analysis_dir + '/sex_de_filtered.csv'
    de_cmd = ['/usr/bin/Rscript', '--vanilla', 'edge.R', expressed_df_path, de_file_path]
    #runCommand(' '.join(de_cmd))

    print('Saving DE results to classification DF')
    de_df = pd.read_csv(de_file_path)
    de_genes = de_df['gene'].tolist()
    results_df = pd.read_csv(results_df_path, sep='\t')
    results_df['Diff_Sex'] = results_df.apply(lambda row: row['gene'] in de_genes, axis=1)
    results_df.to_csv(results_df_path, sep= '\t', index=False)

    inferred_functions_1_path = gigas_dir + '/go_predict/ALL/SPR.c0.92.pval0.05.fdr0.05.tsv'
    associations_path = gigas_dir + '/annotation/step_26-make_id2go/id2go.tsv'
    associations = read_id2gos(associations_path)
    read_id2gos(inferred_functions_1_path, id2gos=associations)

    obo_path = analysis_dir + '/go.obo'
    if not os.path.exists(obo_path):
        runCommand(' '.join(['wget', '-O', obo_path, 'http://purl.obolibrary.org/obo/go.obo']))

    graph = obonet.read_obo(obo_path)
    roots = ['GO:0005575', 'GO:0003674', 'GO:0008150']
    print("Reading graph")
    obo_nodes = graph.nodes(data=True)
    go_ids = [id_ for id_, data in obo_nodes]
    correct_id = {}
    descriptions={}
    print("Solving redundant GO ids")
    for ID in tqdm(go_ids):
        if "alt_id" in obo_nodes[ID]:
            for alt_id in obo_nodes[ID]["alt_id"]:
                correct_id[alt_id] = ID
        correct_id[ID] = ID
        descriptions[ID] = obo_nodes[ID]['name']
    print("Solved " + str(len(correct_id.keys())))

    data_dir = '../data'
    growth_names_path = data_dir + '/function_sets/growth_functions.txt'
    maturation_names_path = data_dir + '/function_sets/sexual_maturation_funcs.txt'
    print("Reading GO lists")
    growth_names = get_go_list(growth_names_path)
    maturation_names = get_go_list(maturation_names_path)

    print("\tGrowth columns")
    results_df['Growth_Functions'] = results_df.apply(
        lambda row: len(genes_annotations_in_set(row['gene'], associations, growth_names, descriptions, correct_id)), axis=1)
    results_df['Growth_Functions_Percent'] = results_df.apply(
        lambda row: (row['Growth_Functions'] / len(associations[row['gene']])) if row['gene'] in associations else 0.0,
        axis=1)
    results_df['Involved_in_Growth'] = results_df.apply(
        lambda row: row['Growth_Functions'] > 0, axis=1)

    print("\tMaturation Columns")
    results_df['Maturation_Functions'] = results_df.apply(
        lambda row: len(genes_annotations_in_set(row['gene'],associations, maturation_names, descriptions, correct_id)), axis=1)
    results_df['Maturation_Functions_Percent'] = results_df.apply(
        lambda row: (row['Maturation_Functions'] / len(associations[row['gene']])) if row['gene'] in associations else 0.0,
        axis=1)
    results_df['Involved_in_Maturation'] = results_df.apply(
        lambda row: row['Maturation_Functions'] > 0, axis=1)

    df2 = results_df.copy(deep=True)
    for sample_name in list(results_df.columns):
        if sample_name.endswith('_Male') or sample_name.endswith('_Female'):
            del df2[sample_name]
    df2['Name'] = df2['gene']
    del df2['gene']
    df2['Number_Of_DE_Packages'] = 2
    df2_path = analysis_dir+"/tissue_analysis.tsv"
    df2.to_csv(df2_path, sep='\t', index=False, header=True)

    write_gene_lists(df2, analysis_dir, descriptions, correct_id)
    make_tissue_summary(df2_path, analysis_dir)