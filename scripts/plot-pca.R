#modification for working in hipergator cluster environment
#Daniel Ence
#December 20, 2020
args=commandArgs(trailingOnly=TRUE)
rds_file=args[1]
out_file=args[2]
labels=args[3]
log_file=args[4]


log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

# load deseq2 data
dds <- readRDS(rds_file)

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
svg(out_file)
plotPCA(counts, intgroup=labels)
dev.off()
