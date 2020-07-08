#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)                
library(tximport)         
library(readr)
library(rlist) # required to append to lists

# data
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)

# performance
library(BiocParallel)
library(tictoc)

tic()

# 0. user defined variables
register(MulticoreParam(8))
setwd("~/scratch/")
salmon_dir = "/Users/alomana/backups/corona/cdi/results/salmon"
metadata_file = '/Volumes/omics4tb2/alomana/projects/cdi/data/metadata/metadata.txt'
results_dir = '/Volumes/omics4tb2/alomana/projects/cdi/results/deseq2/unfiltered/'

# 1. build annotation reference
k = keys(EnsDb.Hsapiens.v86, keytype = "TXNAME")
tx2gene = select(EnsDb.Hsapiens.v86, k, "GENEID", "TXNAME")

# 2. read metadata
metadata = read.table(metadata_file, header = TRUE)
metadata

# 2. store abundance for all samples
files = file.path(salmon_dir, metadata$sample, 'filtered_recalculated_quant.sf')
files

txi = tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps=TRUE)
head(txi$counts)

tpm = txi$abundance
colnames(tpm) = metadata$sample
store = paste(results_dir, 'DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)