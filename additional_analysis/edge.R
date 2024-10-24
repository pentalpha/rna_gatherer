#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args
library(edgeR)
library(dplyr)
library(stringr)
#
#input_path <- '/home/pita/main/experiments/rna-gatherer_gigas_2024/analysis_for_paper/ncrna_expressed.tsv'
input_path <- args[1]
#output_path <- '/home/pita/main/experiments/rna-gatherer_gigas_2024/analysis_for_paper/sex_de.tsv'
output_path <- args[2]

x <- read.delim(input_path, sep=',', row.names = 'gene')
nrow(x)
#x$X <- NULL
is_male <- grepl('_Male', colnames(x))
group <- factor(as.numeric(is_male))

y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
#keep <- filterByExpr(y, )
#y <- y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y)

fit_exac <- exactTest(y)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)

tops_exac = topTags(fit_exac, n=nrow(fit_exac$table))
tops = topTags(qlf, n=nrow(qlf$table))
nrow(tops_exac$table); nrow(tops$table)
df_exac <- as_tibble(tops_exac$table)
df_exac$gene <- row.names(tops_exac$table)
df <- as_tibble(tops$table)
df$gene <- row.names(tops$table)
nrow(df_exac); nrow(df)
relevant_exac <- filter(df_exac, PValue < 0.01 & abs(logFC) >= 1.5)
relevant <- filter(df, PValue < 0.01)
nrow(relevant_exac); nrow(relevant)

relevant_exac$up_regulated <- relevant_exac$logFC > 0
relevant_exac$FDR <- NULL
relevant$up_regulated <- relevant$logFC > 0
write.csv(relevant_exac, output_path)
write.csv(relevant, paste(output_path, '.1', sep=''))
