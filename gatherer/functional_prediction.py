import os
import threading
import numpy as np
import pandas as pd
from tqdm import tqdm
from gatherer.util import *
import math
from scipy.stats.stats import pearsonr, spearmanr
from scipy.special import comb
from scipy.stats import hypergeom
import multiprocessing
import networkx
import obonet
import statsmodels.stats.multitest as multitest
import dcor
from minepy import MINE
from sys import getsizeof
import warnings
from gatherer.confidence_levels import *

def is_constant(vec):
    first = vec[0]
    for i in range(len(vec)):
        if vec[i] != first:
            return False
    return True

def calc_sobolev(reads1,reads2):
    z1 = ((np.square(reads1) / np.sum(np.square(reads1))) -
        (np.square(reads2) / np.sum(np.square(reads2))))
    fourier_transform = np.fft.fft(z1)
    w = (2*math.pi
        *(np.array([x for x in range(0,len(fourier_transform))], dtype=np.float32))
        /(len(fourier_transform)))
    abs_fourier = np.abs(fourier_transform)
    squared_fourier = np.square(abs_fourier)
    w_plus_one = 1+w
    weighted_fourier = w_plus_one*squared_fourier
    fourier_sum = np.sum(weighted_fourier)
    sum_power = np.power(fourier_sum,0.5)
    return sum_power

def calc_sobolev_norm(reads1, reads2):
    sob = calc_sobolev(reads1,reads2)
    reduced_sum = sob / 14.675
    corrected_corr = 1.0 - reduced_sum
    if corrected_corr < 0.0:
            corrected_corr = 0.0
    return corrected_corr

def calc_fisher_information(reads1,reads2):
    square1 = np.square(reads1)
    square2 = np.square(reads2)
    square_sum1 = np.sum(square1)
    square_sum2 = np.sum(square2)
    try:
        div1 = square1 / square_sum1
        div2 = square2 / square_sum2
        t = div1 * div2
        total = np.sum(np.sqrt(t))
        raw_corr = np.arccos(total)
        return raw_corr
    except Warning:
        #print("Numpy invalid division1:\n\t"+str(square1)+"\n\t"+str(square2))
        return 0.0

def calc_fisher_information_norm(reads1,reads2):
    fisher = calc_fisher_information(reads1,reads2)
    norm_fisher = fisher / 1.57
    corrected_corr = 1.0 - fisher
    if corrected_corr < 0.0:
        corrected_corr = 0.0
    return corrected_corr

def calc_mic(reads1,reads2,mine):
    mine.compute_score(reads1, reads2)
    return mine.mic()

def prs(reads1,reads2):
    pearson_corr, p = pearsonr(reads1, reads2)
    return pearson_corr

def spr(reads1,reads2):
    spearman_corr, p = spearmanr(reads1, reads2)
    return spearman_corr

def abs_func(func):
    return lambda x, y: abs(func(x,y))

def dc(reads1,reads2):
    return dcor.distance_correlation(reads1, reads2)

def filter_by_min(val, min_value):
    if val >= min_value:
        return val
    else:
        return None

def metric_with_filter(metric, th_value, metric_name):
    if metric_name == "SOB" or metric_name == "FSH":
        return lambda a,b: geometric_filter(metric(a,b), th_value)
    else:
        return lambda a,b: normal_filter(metric(a,b), th_value)

def get_mic():
    mine = MINE()
    mic = lambda reads1,reads2: calc_mic(reads1,reads2,mine)
    return mic

method_names = {"MIC": get_mic(), "PRS": abs_func(prs), "SPR": abs_func(spr), "DC": dc, 
                "FSH": calc_fisher_information, "SOB": calc_sobolev}

'''
def ai_find_coexpression_process(pid, coding_rows, regulators, predictor, return_dict):
    warnings.filterwarnings('error')
    coding_noncoding_pairs = []
    methods = []
    method_ids = ['DC','MIC','SPR','PRS','SOB','FSH']
    method_names["SOB"] = calc_sobolev
    method_names["FSH"] = calc_fisher_information
    method_names["SPR"] = spr
    method_names["PRS"] = prs
    for method_name in method_ids:
        methods.append(method_names[method_name])

    iterator = range(len(regulators))
    if pid == 0:
        iterator = tqdm(range(len(regulators)))
    for i in iterator:
        row2 = regulators.iloc[i]
        reads2 = np.array(row2.values[1:],dtype=np.float32)
        rows = coding_rows[coding_rows[coding_rows.columns[0]] != row2[coding_rows.columns[0]]]
        correlation_lists = []
        name_pairs = []
        for j in range(len(rows)):
            row1 = coding_rows.iloc[j]
            reads1 = np.array(row1.values[1:],dtype=np.float32)
            corrs = [methods[method_index](reads1,reads2) 
                    for method_index in range(len(methods))]
            corrs_nan = [x if x != None else np.nan for x in corrs]
            corrs_valid = [x if x != np.nan else 0.0 for x in corrs_nan]
            correlation_lists.append(corrs_valid)
            name_pairs.append((row1[0], row2[0]))
        interacting = predictor.predict_many(correlation_lists)
        predicted = 0
        for index in range(len(interacting)):
            if interacting[index]:
                coding_noncoding_pairs.append(name_pairs[index])
                predicted += 1
        if predicted > 0:
            print("\n" + str(predicted) + " predictions added to " + str(len(coding_noncoding_pairs)))
    return_dict[pid] = coding_noncoding_pairs
'''

def ai_find_coexpression_process(pid, coding_rows, nc_rows, predictor, return_dict):
    warnings.filterwarnings('error')
    coding_noncoding_pairs = []
    methods = []
    method_ids = ['DC','MIC','SPR','PRS','SOB','FSH']
    method_names["SOB"] = calc_sobolev
    method_names["FSH"] = calc_fisher_information
    method_names["SPR"] = spr
    method_names["PRS"] = prs
    for method_name in method_ids:
        methods.append(method_names[method_name])

    iterator = range(len(coding_rows))
    if pid == 0:
        iterator = tqdm(range(len(coding_rows)))
    for i in iterator:
        row1 = coding_rows.iloc[i]
        reads1 = np.array(row1.values[1:],dtype=np.float32)
        correlation_lists = []
        name_pairs = []
        for name, row2 in nc_rows.iterrows():
            reads2 = np.array(row2.values[1:],dtype=np.float32)
            corrs = [methods[method_index](reads1,reads2) 
                    for method_index in range(len(methods))]
            corrs_nan = [x if x != None else np.nan for x in corrs]
            corrs_valid = [x if x != np.nan else 0.0 for x in corrs_nan]
            correlation_lists.append(corrs_valid)
            name_pairs.append((row1[0], row2[0]))

        interacting = predictor.predict_many(correlation_lists)
        predicted = 0
        for index in range(len(interacting)):
            if interacting[index]:
                coding_noncoding_pairs.append(name_pairs[index])
                predicted += 1
        if predicted > 0:
            print("\n" + str(predicted) + " predictions added to " + str(len(coding_noncoding_pairs)))
    print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
    return_dict[pid] = coding_noncoding_pairs

def calc_all(pid, coding_rows, regulators, metrics_used, return_dict):
    warnings.filterwarnings('error')
    coding_noncoding_pairs = []

    methods = []
    method_ids = metrics_used
    method_names["SOB"] = calc_sobolev
    method_names["FSH"] = calc_fisher_information

    for method_name in method_ids:
        methods.append(method_names[method_name])
    iterator = range(len(regulators))
    if pid == 0:
        iterator = tqdm(range(len(regulators)))
    for i in iterator:
        row2 = regulators.iloc[i]
        reads2 = np.array(row2.values[1:],dtype=np.float32)
        rows = coding_rows[coding_rows[coding_rows.columns[0]] != row2[coding_rows.columns[0]]]
        for i in range(len(rows)):
            row1 = coding_rows.iloc[i]
            reads1 = np.array(row1.values[1:],dtype=np.float32)
            for i in range(len(method_ids)):
                corr = methods[i](reads1,reads2)
                if corr != None:
                    coding_noncoding_pairs.append((row1[0], row2[0], corr, method_ids[i]))
    return_dict[pid] = coding_noncoding_pairs

def leave_one_out(pid, coding_rows, regulators, method_ids, min_thresholds, return_dict):
    minimum_coefs = min_thresholds
    #print("Turning warnings into errors.")
    warnings.filterwarnings('error')
    valid_metrics = []
    for metric_name in method_ids:
        if minimum_coefs[metric_name] != None:
            valid_metrics.append(metric_name)
    method_ids = valid_metrics

    coding_noncoding_pairs = []

    methods = []
    names = []
    for method_name in method_ids:
        methods.append(
            metric_with_filter(
                method_names[method_name],minimum_coefs[method_name], method_name))
    '''if show:
        print("Regulators list:\n\t"+str(regulators.head())+"\nLen="+str(len(regulators)))
        print(str(coding_rows.head())+"\nLen="+str(len(coding_rows)))'''
    for i in range(len(regulators)):
        row2 = regulators.iloc[i]
        reads2 = np.array(row2.values[1:],dtype=np.float32)
        rows = coding_rows[coding_rows[coding_rows.columns[0]] != row2[coding_rows.columns[0]]]
        for i in range(len(rows)):
            row1 = coding_rows.iloc[i]
            reads1 = np.array(row1.values[1:],dtype=np.float32)
            for i in range(len(method_ids)):
                corr = methods[i](reads1,reads2)
                if corr != None:
                    coding_noncoding_pairs.append((row1[0], row2[0], 
                                                    corr, method_ids[i]))
            
    #print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
    return_dict[pid] = coding_noncoding_pairs

def try_find_coexpression_process(pid, coding_rows, nc_rows, method_ids, 
    min_thresholds, return_dict):
    minimum_coefs = min_thresholds

    valid_metrics = []
    for metric_name in method_ids:
        if minimum_coefs[metric_name] != None:
            valid_metrics.append(metric_name)
    method_ids = valid_metrics

    coding_noncoding_pairs = []

    methods = []
    names = []
    for method_name in method_ids:
        methods.append(
            metric_with_filter(
                method_names[method_name],minimum_coefs[method_name],method_name))
        #methods.append(method_names[method_name])
    #print("Methods="+str(methods))
    iterator = range(len(coding_rows))
    if pid == 0:
        iterator = tqdm(range(len(coding_rows)))
    for i in iterator:
        row1 = coding_rows.iloc[i]
        reads1 = np.array(row1.values[1:],dtype=np.float32)
        for name, row2 in nc_rows.iterrows():
            reads2 = np.array(row2.values[1:],dtype=np.float32)
            for i in range(len(method_ids)):
                corr = methods[i](reads1,reads2)
                if corr != None:
                    coding_noncoding_pairs.append((row1[0], row2[0], corr, method_ids[i]))
    #print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
    return_dict[pid] = coding_noncoding_pairs

def get_ancestors(graph, parent_id):
    ancestors = networkx.ancestors(graph, parent_id)
    return ancestors

def get_descendants(graph, parent_id):
    descendants = networkx.descendants(graph, parent_id)
    return descendants

def pvalue(m, N, M, n):
    '''
    M = number of successes in the population
    m = number of drawn successes
    N = sample size
    n = n
    '''
    '''comb_sum = 0
    for i in range(m, min(n,M)+1):
        comb_sum += comb(M, i, exact=True)*comb(N-M, n-i, exact=True)
    p = comb_sum / comb(N, n, exact=True)
    return p'''
    #return hypergeom.sf(m, N, M, n)
    return hypergeom.sf(m, N, M, n)

def pvalue_process(N, params, pid, return_dict):
    return_dict[pid] = [(i,pvalue(m, N, M, n)) for i, n, M, m in params]

def get_valid_associations(genes_coexpressed_with_ncRNA, genes_annotated_with_term, possible_gene_term,
                            min_n, min_M, min_m):
    n_lens = [len(genes_coexpressed_with_ncRNA[gene]) for gene, term in possible_gene_term]
    M_lens = [len(genes_annotated_with_term[term]) for gene, term in possible_gene_term]
    m_lens = [len(genes_annotated_with_term[term]
                .intersection(genes_coexpressed_with_ncRNA[gene]))
                for gene, term in possible_gene_term]
    valid = [n_lens[i] >= min_n and M_lens[i] >= min_M and m_lens[i] >= min_m for i in range(len(possible_gene_term))]
    return valid, n_lens, M_lens, m_lens


def parallel_pvalues(N, possible_gene_term, valid_gene_term, 
                    n_lens, M_lens, m_lens, 
                    threads, available_memory):
    min_for_chunks = threads * 24 * 1024
    if available_memory < min_for_chunks:
        print("Setting available memory to minimum")
        available_memory = min_for_chunks
    print("Memory available for pvalue calculation: " + str(available_memory/1024))
    print("Dividing chunks for parallel p-value calculation")
    print("Possible associations: " + str(len(possible_gene_term)))
    params = []
    for i in range(len(possible_gene_term)):
        if valid_gene_term[i]:
            #gene_terms.append(possible_gene_term[i])
            params.append((i,n_lens[i],M_lens[i],m_lens[i]))
    print("Calculating pvalues for " + str(len(params)) + " associations.")
    one_size = getsizeof(params[0])
    print("Size of one element: " + str(one_size/1024))
    params_per_run = int(available_memory / one_size)
    run_chunks = list(chunks(params, params_per_run))

    index_pvalue = []
    print("Pvalues per chunk: " + str(len(run_chunks[0])))
    for run_chunk in tqdm(run_chunks):
        calcs_per_thread = int(len(run_chunk) / threads)+1
        calc_chunks = list(chunks(run_chunk, calcs_per_thread))
        #print("Calculating p-values")

        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        processes = []
        showing = False
        i = 0
        for chunk in calc_chunks:
            p = multiprocessing.Process(target=pvalue_process,
                    args=(N, chunk, i, return_dict, ))
            processes.append(p)
            #print("Spawned process")
            p.start()
            if showing:
                showing = False
            i += 1

        for p in processes:
            p.join()
            #print("Finished process " + str(p.pid))

        #print("Merging results")
        for value in return_dict.values():
            index_pvalue += value
        #print("Saved results from processes")
        manager._process.terminate()
        manager.shutdown()
    print("Aggregating pvalue information to (gene,term) pairs")
    gene_term_pvalue = []
    for i in tqdm(range(len(index_pvalue))):
        index,pvalue = index_pvalue[i]
        if pvalue != np.nan:
            gene,term = possible_gene_term[index]
            gene_term_pvalue.append((gene,term,float(pvalue)))
    return gene_term_pvalue
