#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install(c("DESeq2", "tximport", "EnsDb.Hsapiens.v86"))
#BiocManager::install('readr')

library('DESeq2')
library("tximport")
library("EnsDb.Hsapiens.v86")
library('readr')
 
# 0. user defined variables
setwd("~/scratch/")
salmon_dir = "/Users/alomana/corona/cdi/results/salmon"
metadata_file = '/Volumes/omics4tb2/alomana/projects/cdi/data/metadata/metadata.txt'
selection = c(1,2,3,13,14,15)

# 1. build annotation reference
txdb = EnsDb.Hsapiens.v86
k = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")

# 2. read metadata
cdi_metadata = read.table(metadata_file, header = TRUE)
rownames(cdi_metadata) <- cdi_metadata$sample
cdi_metadata = cdi_metadata[selection, ]
cdi_metadata

# 2. define kallisto files
files = file.path(salmon_dir, cdi_metadata$sample, 'quant.sf')
files

# 3. read salmon files 
txi = tximport(files, type="salmon", tx2gene=tx2gene)
head(txi$counts)

# 4. create object for DESeq2
dds = DESeqDataSetFromTximport(txi, colData=cdi_metadata, design=~condition)

# 5. filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 6. differential expression

# 6.1. DEGs between mono- and co-culture at time zero
dds$condition = relevel(dds$condition, ref = "mono")
dds = DESeq(dds)
res = results(dds, lfcThreshold=log2(3/2))
summary(res)

# 6.2. filtering
filt1 = res[which(res$pvalue < 0.05), ]
summary(filt1)

filt2 = filt1[which(filt1$padj < 0.1), ]
summary(filt2)

resOrdered <- filt2[order(filt2$log2FoldChange), ]

# 7. write filtered results
write.table(as.data.frame(resOrdered), file="hypothesis.one.csv", sep='\t')