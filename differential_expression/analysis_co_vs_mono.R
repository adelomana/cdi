#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tictoc")

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
register(MulticoreParam(4))
setwd("~/scratch/")
salmon_dir = "/Users/alopez/projects_isb/cdi/results/expression/salmon"
metadata_file = '/Users/alopez/projects_isb/cdi/data/metadata/metadata.txt'
results_dir = '/Users/alopez/projects_isb/cdi/results/deseq2/unfiltered/'

# 1. build annotation reference
k = keys(EnsDb.Hsapiens.v86, keytype = "TXNAME")
tx2gene = select(EnsDb.Hsapiens.v86, k, "GENEID", "TXNAME")

# 2. read metadata
metadata = read.table(metadata_file, header = TRUE)
metadata

# 2. store abundance for all samples
files = file.path(salmon_dir, metadata$sample, 'filtered_recalculated_quant.sf')
files

# 4. define analysis function
tag = 'co_vs_mono_time_zero'
reference = c(1,2,3)
testing = c(13,14,15)
  
print(tag)
print('reference:')
print(reference)
print('testing:')
print(testing)
  
### f.1. define hypothesis metadata
selected_samples = c(reference, testing)
hypothesis_metadata = metadata[selected_samples, ] 
print(hypothesis_metadata)
  
### f.2. define and read kallisto files
files = file.path(salmon_dir, hypothesis_metadata$sample, 'filtered_recalculated_quant.sf')
names(files) = hypothesis_metadata$sample
print(files)
txi = tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps=TRUE)
  
### f.3. import into DESeq object
dds = DESeqDataSetFromTximport(txi, colData=hypothesis_metadata, design=~condition)
dds$condition = relevel(dds$condition, ref = "mono")
  
### f.4. filtering
threshold = 10
keep = rowSums(counts(dds)) >= threshold  
dds = dds[keep, ]
print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))
  
### f.5 analysis
print('analysis')
dds = DESeq(dds, parallel=TRUE)

# f.6. filter, annotate, format and store
print('filter')
res = results(dds, lfcThreshold=1, parallel=TRUE)
filt1 = res[which(res$pvalue < 0.05), ]
filt2 = filt1[which(filt1$padj < 0.1), ]
print(paste('DEGs found', dim(filt2)[1], sep=' '))
write.table(filt2, file='testing.tsv', quote=FALSE, sep='\t')

print('annotate')
df = as.data.frame(filt2)
df['common'] = rownames(df)
selected = rownames(df)
info = select(EnsDb.Hsapiens.v86, selected, c("GENEBIOTYPE", "GENENAME"), "GENEID")
info['common'] = info$GENEID

descriptions = tryCatch({
  descriptions = select(org.Hs.eg.db, keys=selected, columns=c("GENENAME"), keytype="ENSEMBL")
}, error = function(e) {
  print('Warning: no description found for ENSEMBL IDs')
  descriptions = data.frame('ENSEMBL'=selected, 'GENENAME'=rep('Not found',each=length(selected)))
})

names(descriptions)[names(descriptions) == "GENENAME"] <- "DESCRIPTION" # arrow is needed here!
descriptions['common'] = descriptions$ENSEMBL
dh = merge(df, info, by='common')
di = merge(dh, descriptions, by='common')

print('format')
formatted = di[ , c(8,10,9,12,2,3,6,7)]
up = formatted[formatted$log2FoldChange > 0, ]
down = formatted[formatted$log2FoldChange < 0, ]
sorted_up = up[order(up$log2FoldChange, decreasing=TRUE), ]
sorted_down = down[order(down$log2FoldChange), ]

print('store')
store = paste(results_dir, tag, '_up', '.tsv', sep='')
write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
store = paste(results_dir, tag, '_down', '.tsv', sep='')
write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)

print('---')
  
toc()

txi$abundance['ENSG00000142173', ]
