#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args
#library(edgeR)
library(dplyr)
library(stringr)
input_path <- args[1]

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