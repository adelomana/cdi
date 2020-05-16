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
head(metadata)
metadata

# 3. arrange metadata for hypothesis testing
hypos = list()

# 3.1. mono trajectory
reference_samples = c(1:3)
hypo = list('reference'=reference_samples, 'tag'='mono_time_one_vs_time_zero', 'testing' = c(4,5,6))
hypos = list.append(hypos, hypo)

hypo = list('reference'=reference_samples, 'tag'='mono_time_four_vs_time_zero', 'testing' = c(7,8,9))
hypos = list.append(hypos, hypo)

hypo = list('reference'=reference_samples, 'tag'='mono_time_twentyfour_vs_time_zero', 'testing' = c(10,11,12))
hypos = list.append(hypos, hypo)

# 3.1. co trajectory
reference_samples = c(13,14,15)
hypo = list('reference'=reference_samples, 'tag'='co_time_one_vs_time_zero', 'testing' = c(16,17,18))
hypos = list.append(hypos, hypo)

hypo = list('reference'=reference_samples, 'tag'='co_time_four_vs_time_zero', 'testing' = c(19,20,21))
hypos = list.append(hypos, hypo)

hypo = list('reference'=reference_samples, 'tag'='co_time_twentyfour_vs_time_zero', 'testing' = c(22,23,24))
hypos = list.append(hypos, hypo)

# 4. define analysis function
compare = function(hypo) {
  
  tag = hypo$tag
  reference = hypo$reference
  testing = hypo$testing
  
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
  dds = DESeqDataSetFromTximport(txi, colData=hypothesis_metadata, design=~time)
  dds$time = relevel(dds$time, ref = "zero")
  
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
  #print(descriptions)
  
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
  
  print(txi$abundance['ENSG00000142173', ])
  
  print('---')
  
}

# 5. iterate function
for (hypo in hypos) {
  tempo = compare(hypo)
}

toc()
