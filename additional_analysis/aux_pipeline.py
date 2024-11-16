import json
import subprocess
import os
from sys import getsizeof
from tqdm import tqdm
import sys
import pandas as pd
import obonet
import numpy as np
from ete3 import NCBITaxa

from analyze_lncrna_expression import find_lncrna_classes

ncbi = NCBITaxa()

data_dir = '../data'
growth_names_path = data_dir + '/function_sets/growth_functions.txt'
maturation_names_path = data_dir + '/function_sets/sexual_maturation_funcs.txt'

def runCommand(cmd, print_cmd=True):
    if print_cmd:
        print("\t> " + cmd)
    process = subprocess.call(cmd, shell=True)
    return process

def get_tax_name(id):
    ids = ncbi.translate_to_names([id])
    return str(ids[0])

def evol_sim(taxid1, taxid2):
    try:
        lineage1 = ncbi.get_lineage(taxid1)
    except ValueError as err:
        print(err)
        return 0
    
    l1 = set(lineage1)
    l2 = set(ncbi.get_lineage(taxid2))
    common_taxid = l1.intersection(l2)
    return len(common_taxid)

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

def taxid_from_dbid(row):
    if 'sseqid' in row:
        sseqid = row['sseqid']
        if '_' in sseqid:
            taxid_str = sseqid.split("_")[-1]
            try:
                taxid = int(taxid_str)
                return taxid
            except Exception as ex:
                pass
    
    return None

def find_sequence_homologs(transcripts_fasta, genomes_path, analysis_dir):
    all_rna_names = []
    for rawline in open(transcripts_fasta, 'r').readlines():
        if rawline.startswith('>'):
            all_rna_names.append(rawline.lstrip('>'))
    n_rnas = len(all_rna_names)

    align_results = []
    genomes = []
    for rawline in open(genomes_path, 'r'):
        family, name, taxid, url = rawline.rstrip('\n').split(',')
        print(family, name, url)
        if url.startswith('http'):
            genome_dir = analysis_dir + '/genome_'+name.replace(' ', '_')
            genome_path = genome_dir + '/' + url.split('/')[-1]
            if not os.path.exists(genome_dir):
                os.mkdir(genome_dir)
            if not os.path.exists(genome_path):
                runCommand("cd "+genome_dir+" && wget "+url)
        else:
            genome_path = url
        assert os.path.exists(genome_path)
        genomes.append([family, name, taxid, genome_path])
    
    for family, name, taxid, fasta_db in genomes:
        result_paf = analysis_dir + '/align_to_'+name.replace(' ', '_')+'.paf'
        print('Aligning to', fasta_db)
        created = os.path.exists(result_paf)
        
        if not created:
            index_1 = fasta_db + '.mmi'
            index_2 = os.path.dirname(fasta_db) + '/index.mmi'
            correct_indexes = [i for i in [index_1, index_2] if os.path.exists(i)]
            if len(correct_indexes) == 0:
                print('No index found for', fasta_db)
                index_cmd = 'minimap2 -d ' + index_2 + ' ' + fasta_db
                runCommand(index_cmd)
                correct_indexes.append(index_2)
            
            index_path = correct_indexes[0]
            cmd = " ".join(["minimap2",
                "-x sr -t", str(8),
                index_path, transcripts_fasta,
                ">", result_paf])
            code = runCommand(cmd)
        else:
            print(result_paf, 'already calculated')
        align_results.append([family, name, taxid, result_paf])
    dfs = []
    for family, name, taxid, result_paf in align_results:
        print('Loading', result_paf)
        minimap_df = pd.read_csv(result_paf, sep='\t', header=None, index_col=False,
                names=["qseqid","qseq_len","qstart","qend","strand",
                    "sseqid","sseq_len","sstart","send","matchs",
                    "block_len","quality","13th","14th","15th","16th","17th","18th"])
        minimap_df = minimap_df.astype({"qstart": 'int32', "qend": 'int32', "qseq_len": "int32",
                    "sstart": 'int32', "send": 'int32', "sseq_len": "int32",
                    "quality": 'int32', "block_len": "int32", "matchs": "int32"})
        minimap_df['taxid'] = int(taxid)
        minimap_df['name'] = name
        minimap_df['family'] = family
            
        dfs.append(minimap_df)

    minimap_df = pd.concat(dfs)

    print('Getting species names:')
    all_taxids = [n for n in minimap_df['taxid'].unique().tolist() if type(n) == int]
    print(len(all_taxids), 'taxon ids')
    species_names_vec = ncbi.translate_to_names(all_taxids)
    taxid_to_name = {all_taxids[i]: species_names_vec[i] for i in range(len(all_taxids)) if type(species_names_vec[i]) == str}
    tax_closeness = {t: evol_sim(t, 113544) if t in taxid_to_name else 0 for t, name in taxid_to_name.items()}

    minimap_df["species"] = minimap_df.apply(
        lambda row: taxid_to_name[row['taxid']] if row['taxid'] in taxid_to_name else None, axis=1)

    minimap_df = minimap_df[minimap_df['species'] != "Arapaima gigas"]

    minimap_df["id"] = np.arange(len(minimap_df))

    print("Calculating taxonomic closeness")
    minimap_df["common_taxid"] = minimap_df.apply(
        lambda row: tax_closeness[row['taxid']] if row['taxid'] in tax_closeness else None, axis=1)

    print("Filtering...")

    minimap_df["qcovs"] = minimap_df.apply(
        lambda row: (row["qend"]-row["qstart"]) / row["qseq_len"], axis=1)
    minimap_df["identity"] = minimap_df.apply(
        lambda row: row["matchs"] / row["block_len"], axis=1)

    print(str(minimap_df.head()))
    print(str(len(minimap_df)) + " alignments")

    high_th = 0.80
    #minimap_df = minimap_df.reindex().copy(deep=True)
    minimap_df['cov_id_min'] = minimap_df[['identity','qcovs']].min(axis=1)
    minimap_df['is_similar'] = minimap_df['cov_id_min'] >= high_th

    print("Finding best hits")
    
    print("Getting ncbi common names")
    common_names = ncbi.get_common_names(set(minimap_df['taxid'].unique()))
    common_names[41665] = "bony fishes"
    species_dfs = []
    min_hits_to_list_species = 19

    print('Finding best hits by species')
    taxid_bar =tqdm(total=len(taxid_to_name.keys()))
    for taxid, species_hits in minimap_df.groupby(["taxid"], sort=True):
        if len(species_hits) > 2000:
            print(taxid, len(species_hits), 'hits')
        species_hits2 = species_hits.copy(deep=True)
        best_hits = set()
        rna_groups = tqdm(species_hits2.groupby(['qseqid'])) if len(species_hits) > 2000 else species_hits2.groupby(['qseqid'])
        for gigas_rna, rna_hits in rna_groups:
            sorted = rna_hits.sort_values(["identity","qcovs","common_taxid"],
                ascending=[False,False,False])
            hit = sorted.iloc[0]
            best_hits.add(hit['id'])
        species_best_hits = species_hits[species_hits['id'].isin(best_hits)].copy(deep=True)
        if taxid in common_names:
            common_name = common_names[taxid]
        else:
            common_name = None
        species_best_hits['common_name'] = common_name
        n_similar = len(species_best_hits[species_best_hits['is_similar']])
        if n_similar < min_hits_to_list_species:
            species_best_hits['aligned_to'] = 'others'
        else:
            species_best_hits['aligned_to'] = common_name
        if len(species_hits) > 2000:
            print('\tbest hits:', len(species_best_hits))
        if len(species_hits) > 2000:
            print('\tsimilar:', n_similar)
        species_dfs.append(species_best_hits)
        taxid_bar.update(1)
    taxid_bar.close()
    
    minimap_bests = pd.concat(species_dfs)

    print(str(len(minimap_bests)) + " total ncRNA with homologs.")

    homologs_path = analysis_dir + '/homologs.tsv'
    minimap_bests.to_csv(homologs_path, sep='\t')

    homology_json = {}
    for name, species_hits in minimap_bests.groupby(["aligned_to"], sort=True):
        similar_names = species_hits[species_hits['is_similar']]['qseqid'].unique().tolist()
        n_similar = len(similar_names)
        perc_aligned = n_similar / n_rnas
        print(name, ' best hits are', n_similar)
        hits = []
        taxon_id = []
        scientific_name = []
        for _, row in tqdm(species_hits.iterrows()):
            taxid = row['taxid']
            if taxid == taxid and taxid != None:
                taxon_id.append(taxid)
            name = row['species']
            if name == name and name != None:
                scientific_name.append(name)

            rna_name = row['qseqid']
            cov = row['qcovs']
            identity = row['identity']
            score = row['quality']

            hits.append({'rna': rna_name, 'cov': cov, 'id': identity, 'score': score})
        
        homology_json[name] = {
            'name': name,
            'perc_gigas_mapped': n_similar,
            'perc_aligned': perc_aligned,
            'taxid': list(set(taxon_id)),
            'scientific_name': list(set(scientific_name)),
            'similar_rnas': similar_names,
            'hits': hits
        }
    homologs_json_path = analysis_dir + '/homologs.json'
    json.dump(homology_json, open(homologs_json_path, 'w'), indent=4)
    return minimap_bests, homology_json

def expand_types_review(gigas_dir, homolog_df_path, min_for_hit=0.8):
    print('Loading homolog information')
    homologs_df = pd.read_csv(homolog_df_path, sep='\t')
    homologs_df = homologs_df[homologs_df['cov_id_min'] >= min_for_hit]
    print('Loading rna types json')
    types_json = json.load(open(gigas_dir + '/annotation/step_24-review_annotations/type_review.json'))

    print("Calculating Relevant alignments")
    has_relevant_hit = set()
    relevant_hits_at_distant_species = set()
    by_species = []
    for species_name, species_df in homologs_df.groupby(['species']):
        common_names = species_df['common_name'].unique().tolist()
        common_names = [c for c in common_names if type(c) == str and not(c in ['False'])]
        if len(common_names) > 0:
            species_name = common_names[0]
        relevant_hits = species_df['qseqid'].unique().tolist()
        has_relevant_hit.update(relevant_hits)
        
        if len(relevant_hits) < 19:
            relevant_hits_at_distant_species.update(relevant_hits)
        else:
            by_species.append({
                'Species': species_name.title(),
                'A. gigas ncRNAs With Sequence Similarity': len(relevant_hits)
            })
    by_species.append({'Species': 'All', 
        'A. gigas ncRNAs With Sequence Similarity': len(has_relevant_hit)})
    by_species.sort(key=lambda b: b['A. gigas ncRNAs With Sequence Similarity'], reverse=True)
    by_species.append({'Species': 'Others', 
        'A. gigas ncRNAs With Sequence Similarity': len(relevant_hits_at_distant_species)})
    
    print("Calculating Perc With Sequence Homology")
    homology_by_type = {}
    for type_dict in types_json:
        type_name = type_dict['ncRNA Type']
        ncRNAs = type_dict['ncRNA IDs']
        relevant_ncrnas = has_relevant_hit.intersection(ncRNAs)
        homology_by_type[type_name] = (len(relevant_ncrnas) / len(ncRNAs)) * 100
    
    print('Loading rna types df')
    types_df_path = gigas_dir + '/annotation/step_24-review_annotations/type_review.tsv'
    types_df = pd.read_csv(types_df_path, sep='\t')
    types_df["Found With High Similarity in Other Species"] = types_df.apply(
        lambda row: homology_by_type[row['ncRNA Type']], axis=1)
    
    types_df.to_csv(os.path.dirname(homolog_df_path)+'/types_review.tsv', sep='\t', index=False, decimal=',')

    species_df = pd.DataFrame(by_species)
    species_df.to_csv(os.path.dirname(homolog_df_path)+'/similar_species.tsv', sep='\t', index=False, decimal=',')

def make_id2gos(out_file, id2go_file):
    id2gos = {}
    with open(id2go_file, 'r') as stream:
        print("Reading " + id2go_file)
        for line in stream:
            cells = line.rstrip("\n").split("\t")
            transcript_name = cells[0]
            go_name = cells[1]
            if not transcript_name in id2gos:
                id2gos[transcript_name] = set()
            id2gos[transcript_name].add(go_name)
    with open(out_file, 'w') as stream:
        print("\tWriting in id2gos format")
        for name, gos in id2gos.items():
            stream.write(name+"\t"+";".join(gos)+"\n")
    return id2gos

if __name__ == '__main__':
    gigas_dir = sys.argv[1]
    #niloticus_genome_path = sys.argv[2]
    #arowana_genome_path = sys.argv[3]
    #rnacentral_db_path = sys.argv[4]
    analysis_dir = gigas_dir + '/analysis_for_paper'
    samples_tpm = gigas_dir + '/go_predict/samples_tpm.tsv'
    results_df_path = analysis_dir + '/ncrna_classification.tsv'
    expressed_df_path = analysis_dir + '/ncrna_expressed.tsv'
    transcriptome_fasta = gigas_dir + '/annotation/step_25-write_transcriptome/transcriptome.fasta'
    genomes_df_path = 'other_genomes.csv'
    if not os.path.exists(analysis_dir):
        os.mkdir(analysis_dir)

    homolog_df_path = analysis_dir+'/homologs.tsv'
    homolog_json_path = analysis_dir+'/homologs.json'
    #if not os.path.exists(homolog_json_path):
    homolog_df, homolog_json = find_sequence_homologs(transcriptome_fasta, genomes_df_path, analysis_dir)
    expand_types_review(gigas_dir, homolog_df_path)
    
    de_file_path = analysis_dir + '/sex_de_filtered.csv'
    '''if not os.path.exists('counts_noduplicates_noprot.csv'):
        filter_ncrna_cmd = ['/usr/bin/Rscript', '--vanilla', 'filter_ncrna_counts.R', samples_tpm]
        runCommand(' '.join(filter_ncrna_cmd))
    results_df_path, expressed_df_path = find_lncrna_classes('counts_noduplicates_noprot.csv', analysis_dir)

    print('Performing DE')
    
    if not os.path.exists(de_file_path):
        de_cmd = ['/usr/bin/Rscript', '--vanilla', 'edge.R', expressed_df_path, de_file_path]
        runCommand(' '.join(de_cmd))'''

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


    runCommand("mkdir " + outdir)
    tissue_names = ['Brain', 'Gonad', 'Heart', 
                    'Kidney', 'Liver', 'Lung',
                    'Muscle', 'Skin']
    print("Parsing gene lists")
    lists_to_enrich = []
    for tissue in tissue_names:
        lists_to_enrich.append((tissue,gene_list_dir+"/"+tissue+".tissue.txt"))

    for sex in ["female", "male"]:
        lists_to_enrich.append(("Skin-"+sex+'-Sex_DE', gene_list_dir+"/Skin."+sex+"_expressed.txt"))
        lists_to_enrich.append(("Skin-"+sex+'-Sex_DE', gene_list_dir+"/Skin."+sex+"_expressed_de.txt"))
    lists_to_enrich.append(("Sex_DE", gene_list_dir+"/sex_diff.txt"))

    lists_to_enrich.append(("housekeeping",gene_list_dir+"/housekeeping.txt"))
    #lists_to_enrich.append(("growth",gene_list_dir+"/involved_in_growth.txt"))
    #lists_to_enrich.append(("growth_housekeeping",gene_list_dir+"/involved_in_growth-housekeeping.txt"))
    #lists_to_enrich.append(("maturation",gene_list_dir+"/involved_in_maturation.txt"))
    associations_file_path = outdir + "/associations.tsv"
    all_ncgene_ids_file = outdir + "/ncrna_pop.txt"
    associations = make_id2gos(associations_file_path, predictions_file)
    # = make_population_from_associations(outdir, associations)
    enrichments_dir = outdir + "/enrichments"
    runCommand("mkdir " + enrichments_dir)

    # %%
    for name, list_file in lists_to_enrich:
        if not os.path.exists(list_file):
            print("Could not find " + list_file)

    cmds = {name: " ".join(["find_enrichment.py",
                            "--pval="+str(max_pval)+" --indent",
                            "--obo", obo_path,
                            "--outfile", 
                            enrichments_dir+"/"+name+".tsv",
                            list_file, all_ncgene_ids_file,
                            associations_file_path])
            for name, list_file in lists_to_enrich}

    for tissue_name, cmd in tqdm(cmds.items()):
        runCommand(cmd)