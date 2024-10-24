#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args
library(edgeR)
library(dplyr)
library(stringr)
input_path <- args[1]
output_path <- args[2]
x0 <- read.csv(input_path, sep='\t')
numeric_cols <- colnames(x0)[2:length(colnames(x0))]
for (a in numeric_cols){
  x0[a] <- as.numeric(unlist(x0[a]))
}
names <- x0$Name
x0$Name <- NULL
x1 <- rowsum(x0, names)
#write.csv(x1, 'counts_noduplicates.csv')
x1$gene <- row.names(x1)
x2 <- as_tibble(x1)

non_protein <- filter(x2, str_detect(gene, 'tr\\|', negate=TRUE))
non_protein <- filter(non_protein, str_detect(gene, 'sp\\|', negate=TRUE))
nrow(non_protein)
write.csv(non_protein, 'counts_noduplicates_noprot.csv', row.names = FALSE)

x <- read.delim('counts_noduplicates_noprot.csv', sep=',', row.names = 'gene')
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
relevant_exac <- filter(df_exac, FDR < 0.05 & PValue < 0.05)
relevant <- filter(df, FDR < 0.05 & PValue < 0.05)
nrow(relevant_exac); nrow(relevant)

relevant_exac$up_regulated <- relevant_exac$logFC > 0
relevant$up_regulated <- relevant$logFC > 0
write.csv(relevant_exac, output_path)
write.csv(relevant, paste(output_path, '.1', sep=''))
