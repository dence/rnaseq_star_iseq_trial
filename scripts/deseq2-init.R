#modification for working in hipergator cluster environment
#Daniel ence
#October 22, 2020
args=commandArgs(trailingOnly=TRUE)
counts_file=args[1]
output_file=args[2]
sample_file=args[3]
log_file=args[4]
threads=args[5]



log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("tidyverse")

parallel <- FALSE
if (threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(counts_file, header=TRUE, row.names="gene", check.names=FALSE)

#need to sort cts columns
cts <- cts[ , order(names(cts))]

coldata <- read.table(sample_file, header=TRUE, row.names="sample", check.names=FALSE)

#need to sort coldata by value of "sample" columns
row.names(coldata) -> coldata$tmp_sample_names
coldata <- arrange(coldata, tmp_sample_names) %>% select(-tmp_sample_names)




dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~ condition)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=output_file)
