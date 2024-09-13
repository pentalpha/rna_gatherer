"""Function predictor for lncRNA that uses gene coexpression.

This script is used to apply statistics in order to figure out
possible Gene Ontology terms for lncRNAs, based on gene 
expression counts and an annotation for coding genes. The gene
expression is read from a counts table (see test_data/counts)
and the expression of lncRNA and coding genes is compared
by calculating correlation coefficients through several
differente metrics. The pairs of (lncRNA,gene) with 
coefficients passing the required minimum value to be 
considered inside the confidence level are taken as
coexpressed.
The associations between coding genes and GO terms in the 
functional annotation (see test_data/annotation) of the
coding genes are passed as possible associations for their
lncRNA coexpressed pairs. Statistical tests are used to filter
which of these possible associations are statistically
relevant, calculating P-Values and FDRs.
The associations passing the filtering are grouped together in
a output file, the functional prediction.

Author: Pitágoras Alves (github.com/pentalpha)
"""

import sys
import obonet
import networkx
import numpy as np
import argparse
import multiprocessing
import json

from gatherer.functional_prediction import *
from gatherer.util import *
from gatherer.confidence_levels import *
from gatherer.bioinfo import short_ontology_name, long_ontology_name
'''
Loading configurations and command line arguments
'''
from config import configs, require_files

class GoPredictor:

    #all_methods = ["MIC","DC","PRS","SPR","SOB","FSH"]
    all_ontologies = ["molecular_function","cellular_component","biological_process"]
    default_threads = max(2, multiprocessing.cpu_count()-1)
    min_correlation = 0.9

    """
    count_reads_path:       Count reads table path. TSV file with where the first column must be the gene names and the
                                following columns are the FPKM normalized counts of each sample.
    regulators:             a list of regulators, where each line contains the name of one gene.
    coding_gene_ontology_path:   
                            Path to functional annotation file of the genes.
    tempDir:                Output directory.
    processes:              CPUs to use.
    threshold:              Sets an unique correlation coefficient threshold for all methods
    k_min_coexpressions:    The minimum number of ncRNAs a Coding Gene must be coexpressed with.
                                Increasing the value improves accuracy of functional assignments, but
                                may restrict the results. Default: 1.
    pvalue:                 Maximum pvalue. Default: 0.05
    fdr:                    Maximum FDR. Default: 0.05
    min_m:                  Minimum m value. Default: 1
    min_M:                  Minimum M value. Default: 1
    min_n:                  Minimum n value. Default: 1
    cache_usage:            Portion of the cache memory to use for storing the counts table.
    cache_size:             Sets the size of the cache memory. Default: auto-detection of CPU cache size.
    ontology_type:          One of the following: molecular_function (default), cellular_component, 
                                biological_process or ALL.
    """
    def __init__ (self, count_reads_path, regulators, coding_gene_ontology_path, tempDir,
            processes = default_threads, threshold = min_correlation,
            k_min_coexpressions = 1, pvalue = 0.05, fdr = 0.05, min_m = 1, 
            min_M = 1, min_n = 1, cache_usage = 0.6, cache_size = None,
            ontology_type = 'molecular_function'):
        self.coding_gene_ontology_path = coding_gene_ontology_path
        mandatory_files = ["go_obo"]
        require_files(mandatory_files)
        #set configuration values
        self.confs = {}
        for conf in configs:
            self.confs[conf] = configs[conf]

        self.display_cache = cache_size

        self.available_cache = get_cache(usage=cache_usage)
        if cache_size :
            self.available_cache = int(int( cache_size ) * cache_usage)
        print("Available cache memory: " + str(int(self.available_cache/1024)) + "KB")
        self.go_path = self.confs["go_obo"]

        self.ontology_types_arg = ontology_type.split(",")
        if self.ontology_types_arg[0] == "ALL":
            self.ontology_types_arg = GoPredictor.all_ontologies
        
        self.threads = processes
        self.min_threshold_SPR = threshold

        self.pval = pvalue
        self.fdr = fdr
        self.K = k_min_coexpressions

        self.regulators_max_portion = 0.4

        '''
        Creating output directory
        '''

        if not os.path.exists(tempDir):
            os.mkdir(tempDir)

        self.correlations_dir = tempDir + "/correlations"
        if not os.path.exists(self.correlations_dir):
            os.mkdir(self.correlations_dir)

        '''
        Pre-processing of the count-reads table
        '''

        reads = pd.read_csv(count_reads_path, sep='\t')
        print(str(len(reads)) + " raw rows.")
        reads["constant"] = reads.drop([reads.columns[0]], axis=1).apply(
                lambda row: is_constant(np.array(row.values,dtype=np.float32)),axis=1
            )
        mask = reads["constant"] == False
        self.reads = reads[mask]
        del self.reads["constant"]
        print(str(len(self.reads)) + " rows after removing constant rows.")
        print(self.reads.head())

        print("Reading regulators")
        self.regulators = regulators
        '''with open(regulators_path,'r') as stream:
            for line in stream.readlines():
                self.regulators.append(line.rstrip("\n").lstrip(">"))'''

        '''
        Looking for metrics not calculated yet.
        '''

        self.correlation_file = self.correlations_dir+"/SPR.tsv"
        delete_if_empty(self.correlation_file, min_cells=3, sep="\t")
        self.tempDir = tempDir
    
    def get_metric_file(self, metric_name):
        return open(self.correlations_dir + "/" + metric_name + ".tsv", 'a+')
    
    def find_correlated(self, reads, regulators, tempDir, method_stream, threads):
        """Find coexpressed pairs using a set of metrics."""
        if len(reads) < threads*2:
            threads = len(reads)/2
        coding_noncoding_pairs = []
        func = try_find_coexpression_process
        genes_per_process = int(len(reads) / threads)
        limit = len(reads)-1
        end = 0
        last = -1
        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        processes = []
        last_pid = 0
        print('Spawning', threads, 'threads')
        for i in range(threads-1):
            start = last+1
            end = start + genes_per_process
            if end >= limit:
                end = limit
            parcial_df = reads.iloc[start:end]
            #pid, coding_rows, nc_rows, min_threshold, return_dict
            p = multiprocessing.Process(target=func, 
                args=(i, parcial_df, regulators, self.min_threshold_SPR, return_dict, ))
            processes.append(p)
            #print("Spawned process from gene " + str(start) + " to " + str(end))
            p.start()
            last = end

            last_pid = i
            if end == limit:
                break
        if end < limit:
            parcial_df = reads.iloc[end:limit]
            p = multiprocessing.Process(target=func, 
                args=(last_pid+1, parcial_df, regulators, self.min_threshold_SPR, return_dict, ))
            processes.append(p)
            #print("Spawned process from gene " + str(end) + " to " + str(limit))
            p.start()
        
        for p in processes:
            p.join()

        #print("Merging results")
        for value in return_dict.values():
            coding_noncoding_pairs += value

        #print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
        self.output = tempDir+"/correlated.tsv"
        for coding_name, noncoding_name, corr, method_name in coding_noncoding_pairs:
            method_stream.write("\t".join([coding_name,noncoding_name,str(corr)]) + "\n")
        manager._process.terminate()
        manager.shutdown()
        del manager

    def calc_correlation_file(self):
        print("Calculating correlation coefficients")
        print("Separating regulators from regulated.")
        print("\tRegulator IDs: " + str(len(self.regulators)))
        mask = self.reads[self.reads.columns[0]].isin(self.regulators)
        regulators_reads = self.reads.loc[mask]
        print("\tRegulator IDs in dataframe: " 
            + str(len(regulators_reads[regulators_reads.columns[0]].tolist())))
        non_regulators_reads = self.reads.loc[~mask]

        print(str(len(non_regulators_reads)) + " regulated.")
        print(str(len(regulators_reads)) + " regulators.")
        
        available_size = self.available_cache

        '''
        Split the table into cache-sized smaller parts.
        '''
        max_for_regulators = available_size*self.regulators_max_portion
        #print("Available for regulators: " + str(int(max_for_regulators/1024)) + "KB")
        #regs_size = getsizeof(regulators_reads)
        #print("Regulators size: " + str(int(regs_size/1024)) + "KB")
        #regulator_dfs = [regulators_reads]
        print("Dividing regulator rows:")
        regulator_dfs = split_df_to_max_mem(regulators_reads, max_for_regulators)
        available_size -= getsizeof(regulator_dfs[0])
        print("Dividing non-regulator rows:")
        dfs = split_df_to_max_mem(non_regulators_reads, available_size)
        
        '''print("Chunks for regulated: " + str(len(dfs)) 
        + "\nChunks for regulators: " + str(len(regulator_dfs)))'''

        df_pairs = []
        for df in dfs:
            for regulator_df in regulator_dfs:
                df_pairs.append((df,regulator_df))

        
        method_stream = self.get_metric_file('SPR')

        '''
        Calculate the correlations for each part of the table
        '''
        i = 1
        for df,regulator_df in tqdm(df_pairs):
            self.find_correlated(df, regulators_reads, self.tempDir, 
                        method_stream, self.threads)
            i += 1
        
        method_stream.close()
    
    def predict(self, ontology_type="molecular_function"):
        """Predict functions based on the loaded correlation coefficients."""
        
        K = self.K
        min_n = self.min_n
        min_M = self.min_M
        min_m = self.min_m
        ontology_type_mini = short_ontology_name(ontology_type)
        threshold = self.min_threshold_SPR
        pval_threshold = self.pval
        fdr_threshold = self.fdr
        
        print("Current method = SPR, ontology type = " + ontology_type
                + ", pvalue = " + str(pval_threshold) + ", fdr = " + str(fdr_threshold))
        print("Current thresholds = " + str(threshold))
        
        out_dir = self.tempDir+"/"+ontology_type_mini
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        longer = []
        if K != 1 or min_n != 1 or min_M != 1 or min_m != 1:
            longer = ["K"+str(K),"n"+str(min_n),"M"+str(min_M),"m"+str(min_m)]
        run_mode = "LNC"
        
        output_file = (out_dir + "/" +".".join(['SPR',"c"+str(threshold),
                run_mode,"pval"+str(pval_threshold),"fdr"+str(fdr_threshold)]
                + longer + ["tsv"])).replace(".LNC.",".")
        output_file = output_file.replace(".thNone.",".")

        if os.path.exists(output_file):
            print(output_file + " already created, skiping.")
            return output_file, False

        print("Selecting significant genes for processing")
        valid_coding_genes = {}
        valid_genes_coexpressed_with_ncRNA = {}
        total = 0
        correct_method = 0
        valid_corr = 0
        for rna in self.genes_coexpressed_with_ncRNA.keys():
            for gene in self.genes_coexpressed_with_ncRNA[rna]:
                metric = 'SPR'
                key = (rna,gene,metric)
                if key in self.correlation_values.keys():
                    corr = float(self.correlation_values[key])
                    correct_method += 1
                    if compare_to_th(corr, threshold, metric):
                    #if corr >= thresholds[metric] or corr <= -thresholds[metric]:
                        valid_corr += 1
                        if not gene in valid_coding_genes:
                            valid_coding_genes[gene] = 0
                        valid_coding_genes[gene] += 1
                        if not rna in valid_genes_coexpressed_with_ncRNA:
                            valid_genes_coexpressed_with_ncRNA[rna] = set()
                        valid_genes_coexpressed_with_ncRNA[rna].add(gene)
                total += 1
        
        print("len(valid_genes_coexpressed_with_ncRNA)=" + str(len(valid_genes_coexpressed_with_ncRNA)))
        print("valid correlations loaded:", str(valid_corr))
        #print("Discarding coding genes with too little correlations with regulators.")
        genes_to_discard = set()
        for coding_gene in valid_coding_genes.keys():
            if valid_coding_genes[coding_gene] < K:
                genes_to_discard.add(coding_gene)

        for gene in genes_to_discard:
            del valid_coding_genes[gene]
            for rna in valid_genes_coexpressed_with_ncRNA.keys():
                if gene in valid_genes_coexpressed_with_ncRNA[rna]:
                    valid_genes_coexpressed_with_ncRNA[rna].remove(gene)
        #print("len(valid_genes_coexpressed_with_ncRNA)=" + str(len(valid_genes_coexpressed_with_ncRNA)))

        valid_id2gos = {}
        valid_genes_annotated_with_term = {}

        print("Reading annotation of coding genes from " + self.coding_gene_ontology_path)
        id2gos = self.onto_id2gos[ontology_type]
        for _id in id2gos.keys():
            if _id in valid_coding_genes.keys():
                valid_id2gos[_id] = id2gos[_id]
        
        genes_annotated_with_term = self.onto_genes_annotated_with_term[ontology_type]
        for term in genes_annotated_with_term.keys():
            for _id in genes_annotated_with_term[term]:
                if _id in valid_coding_genes.keys():
                    if not term in valid_genes_annotated_with_term.keys():
                        valid_genes_annotated_with_term[term] = set()
                    valid_genes_annotated_with_term[term].add(_id)

        genes_annotated_with_term2 = {}

        #print("Extending associations of terms to genes by including children")
        found = 0
        for go in valid_genes_annotated_with_term.keys():
            genes = set()
            genes.update(valid_genes_annotated_with_term[go])
            before = len(genes)
            if go in self.graph:
                found += 1
                childrens = get_ancestors(self.graph, go)
                for children_go in childrens:
                    if children_go in valid_genes_annotated_with_term:
                        genes.update(valid_genes_annotated_with_term[children_go])
            genes_annotated_with_term2[go] = genes
        #print(str((float(found)/len(valid_genes_annotated_with_term.keys()))*100) + "% of the GO terms found in network.")
        valid_genes_annotated_with_term = genes_annotated_with_term2

        print("Listing possible associations between rnas and GOs")
        possible_gene_term = []
        for rna in valid_genes_coexpressed_with_ncRNA.keys():
            for go in valid_genes_annotated_with_term.keys():
                possible_gene_term.append((rna, go))
        if len(possible_gene_term) == 0:
            print("No possible association to make, under current parameters and data."
                + " Suggestions: Try a different correlation threshold or a different method.")
            return "", False
        #print("len(valid_genes_coexpressed_with_ncRNA)=" + str(len(valid_genes_coexpressed_with_ncRNA)))
        #print("len(valid_genes_annotated_with_term)= " + str(len(valid_genes_annotated_with_term)))
        #print("Possible gene,term = " + str(len(possible_gene_term)))
        valid_gene_term, n_lens, M_lens, m_lens = get_valid_associations(valid_genes_coexpressed_with_ncRNA,
                                                valid_genes_annotated_with_term,
                                                possible_gene_term,
                                                min_n=min_n, min_M=min_M, min_m=min_m)

        print("Calculating p-values")
        gene_term_pvalue = parallel_pvalues(self.N, possible_gene_term, 
                                            valid_gene_term, n_lens, M_lens, m_lens, 
                                            self.threads, self.available_cache)
        print("Calculating corrected p-value (FDR)")
        pvalues = [pval for gene, term, pval in gene_term_pvalue]
        reject, fdrs, alphacSidak, alphacBonf = multitest.multipletests(pvalues, alpha=0.05, method='fdr_by')
        #print("Finished calculating pvalues, saving now")
        with open(self.tempDir + "/association_pvalue.tsv", 'w') as stream:
            for i in range(len(gene_term_pvalue)):
                if valid_gene_term[i]:
                    rna, term, pvalue = gene_term_pvalue[i]
                    stream.write(rna+"\t"+term+"\t"+str(pvalue)+"\t" + str(fdrs[i]) + "\n")

        print("Selecting relevant pvalues and fdr")
        relevant_pvals = []
        rna_id2gos = {}
        pval_passed = 0
        fdr_passed = 0
        for i in tqdm(range(len(gene_term_pvalue))):
            rna, term, pvalue = gene_term_pvalue[i]
            fdr = fdrs[i]
            if pvalue <= pval_threshold:
                pval_passed += 1
                if fdr <= fdr_threshold:
                    fdr_passed += 1
                    relevant_pvals.append((rna, term, pvalue, fdr))
                    if not rna in rna_id2gos:
                        rna_id2gos[rna] = set()
                    rna_id2gos[rna].add(term)
        print(str(pval_passed) + " rna->go associations passed p-value threshold ("
                + str((pval_passed/len(gene_term_pvalue))*100) + "%)")
        print(str(fdr_passed) + " rna->go associations passed fdr threshold ("
                + str((fdr_passed/len(gene_term_pvalue))*100) + "%)")
        
        print("Writing results")
        print("Output annotation is " + output_file)
        with open(output_file, 'w') as stream:
            for rna, term, pvalue, fdr in relevant_pvals:
                stream.write("\t".join([rna,term,ontology_type,str(pvalue),str(fdr)])+"\n")
        return output_file, True

    def run_analysis(self):
        '''
        Calculate any missing metrics
        '''
        if not os.path.exists(self.correlation_file):
            self.calc_correlation_file()

        self.coding_genes = {}
        self.genes_coexpressed_with_ncRNA = {}
        self.correlation_values = {}

        '''
        Load correlation coefficients from the files where they are stored.
        '''
        print("Loading correlations from " + self.correlation_file + ".")
        min_value = self.min_threshold_SPR
        load_condition = lambda x: x >= min_value
        
        with open(self.correlation_file,'r') as stream:
            lines = 0
            invalid_lines = 0
            loaded = 0
            for raw_line in stream.readlines():
                cells = raw_line.rstrip("\n").split("\t")
                if len(cells) == 3:
                    gene = cells[0]
                    rna = cells[1]
                    corr = cells[2]
                    corr_val = float(corr)
                    
                    #if corr_val >= min_value:
                    if load_condition(corr_val):
                        if not gene in self.coding_genes:
                            self.coding_genes[gene] = 0
                        self.coding_genes[gene] += 1
                        
                        if not rna in self.genes_coexpressed_with_ncRNA:
                            self.genes_coexpressed_with_ncRNA[rna] = set()
                        self.genes_coexpressed_with_ncRNA[rna].add(gene)
                        self.correlation_values[(rna,gene,'SPR')] = corr
                        loaded += 1
                else:
                    invalid_lines += 1
                lines += 1
            
            if lines == 0:
                print("Fatal error, no correlations could be loaded from "
                        + self.correlation_file + "\n(The file may be "
                        + "corrupted or just empty)")
                quit()
            else: 
                print(str(float(invalid_lines)/lines)
                    + " lines without proper number of columns (4 columns)")
                print(str(loaded), "total correlations loaded from file")
        print("correlation_values = "+str(len(self.correlation_values.keys())))
        print("genes_coexpressed_with_ncRNA = "+str(len(self.genes_coexpressed_with_ncRNA.keys())))
        print("coding_genes = "+str(len(self.coding_genes.keys())))
        self.N = len(self.reads)

        print("Loading GO network.")
        self.graph = obonet.read_obo(self.go_path)

        self.onto_id2gos = {"biological_process":{},"molecular_function":{},"cellular_component":{}}
        self.onto_genes_annotated_with_term = {"biological_process":{},"molecular_function":{},"cellular_component":{}}

        print("Reading annotation of coding genes from " + self.coding_gene_ontology_path)
        with open(self.coding_gene_ontology_path,'r') as stream:
            lines = 0
            invalid_lines = 0
            associations = 0
            for raw_line in stream.readlines():
                cells = raw_line.rstrip("\n").split("\t")
                if len(cells) == 3 or len(cells) == 4:
                    gene = cells[0]
                    go = cells[1]
                    onto = cells[2]
                    if onto in self.onto_id2gos.keys():
                        id2gos = self.onto_id2gos[onto]
                        genes_annotated_with_term = self.onto_genes_annotated_with_term[onto]
                        #coding_genes.add(gene)
                        if not (gene in self.coding_genes):
                            gene = gene.upper()
                        if gene in self.coding_genes:
                            if not gene in id2gos:
                                id2gos[gene] = set()
                            id2gos[gene].add(go)
                            if not go in genes_annotated_with_term:
                                genes_annotated_with_term[go] = set()
                            genes_annotated_with_term[go].add(gene)
                            associations += 1
                        else:
                            invalid_lines += 1
                            print("Invalid coding gene" + gene)
                    else:
                        invalid_lines += 1
                else:
                    invalid_lines += 1
                lines += 1
            print(str(float(invalid_lines)/lines)
                    + " lines without proper number of columns (3 or 4 columns)")
            print(str(associations) + " valid associations loaded.")
        
        output_files = []
        created = 0
        for onto in self.ontology_types_arg:
            out_file, made = self.predict(ontology_type=onto)
            if out_file != "":
                output_files.append((out_file,onto))
            if made:
                created += 1
        
        print("Writing annotation file with all ontologies")
    
        if len(output_files) > 1 and created > 1:
            lines = []
            ontos = set()
            for output_file,onto_value in output_files:
                with open(output_file,'r') as stream:
                    new_lines = [line for line in stream.readlines()]
                    lines += new_lines
                ontos.add(onto_value)
            ontos = list(ontos)
            ontos.sort()
            ontos_str = "_".join([short_ontology_name(str(onto)) 
                                for onto in ontos])
            if len(ontos) == 3:
                ontos_str = "ALL"
            onto_dir = self.tempDir + "/" + ontos_str
            if not os.path.exists(onto_dir):
                os.mkdir(onto_dir)
            output_file = (onto_dir + "/" + ".".join(['SPR',
                                "c"+str(self.min_threshold_SPR),
                                "pval"+str(self.pval),"fdr"+str(self.fdr),"tsv"]
                            ))
            with open(output_file,'w') as stream:
                for line in lines:
                    stream.write(line)