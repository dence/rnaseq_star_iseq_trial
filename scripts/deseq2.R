#modification for working in hipergator cluster environment
#Daniel Ence
#December 20, 2020
args=commandArgs(trailingOnly=TRUE)
rds_file=args[1]
out_table=args[2]
ma_plot=args[3]
log_file=args[4]
threads=args[5]
curr_contrast=args[6:7]

log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(threads))
    parallel <- TRUE
}

dds <- readRDS(rds_file)

contrast <- c("condition", curr_contrast)
res <- results(dds, contrast=contrast, parallel=parallel)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)
# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage


# store results
svg(ma_plot)
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file=out_table)
