message ("Loading the DESeq2 library...")

library(DESeq2)

# Log versions of R and packages
session_info <- capture.output(sessionInfo())
session_info_file <- paste("R_session_info.log",sep="")
write.table(session_info, file=session_info_file, quote=FALSE,sep="\t",eol="\n")

args<-commandArgs(TRUE)

count_file <- snakemake@input[["counts"]]
design_file <- snakemake@input[["matrix"]]
results_path <- snakemake@output[["pvals"]]
hist_path <- snakemake@output[["hist"]]
normcounts_path <- snakemake@output[["normcounts"]]
is_paired <- FALSE
do_filter <- TRUE

message ("Reading the count and design files into memory")
counts <- read.table(count_file, header=TRUE, row.names=1, sep="\t")

# Filtering rows of 0s and then rounding since DESeq doesn't work with fractional counts
#counts <- round(counts[ rowSums(counts) > 0,])

# new filtering (kick out if more than half of samples has zero counts)
counts <- counts[rowSums(counts > 0) > ncol(counts)/2,]

head(counts)

# rounding counts to nearest integers (otherwise deseq2 complains)
counts <- round(counts)

design <- read.table(design_file, header=TRUE, row.names=1, sep="\t")

# Paired or unpaired analysis
if (is_paired) {
	message ("A paired analysis will be performed")
	samples <- data.frame(row.names=colnames(counts), condition=design)
	deseq2obj <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design=~condition.factor1+condition.factor2)
} else {
	message ("An unpaired analysis will be performed")
	samples <- data.frame(row.names=colnames(counts), condition=design$factor)
	deseq2obj <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design=~condition)
}

message ("Performing differential expression analysis using DESeq2...")
deseq2objde <- DESeq(deseq2obj)

results <- results(deseq2objde, cooksCutoff=Inf)


if (do_filter) {
    message ("Filtering genes/KOs/reactions with very low mean expression values...")
    use<-results$baseMean>metadata(results)$filterThreshold
    results <- results [ use, ]
}

head(results)

# store results
gene_ranks <- order(results$pvalue)
sorted_genes <- rownames(results)[gene_ranks]
results <- results[gene_ranks,]
write.table(results, file=results_path, quote=FALSE,sep="\t",eol="\n")

# tablename2 <- paste(output_prefix, "_pval_fc.txt", sep="")
# #pval_and_fc <- data.frame(row.names=paste("\"", row.names(results), "\"", sep=""), results$pvalue, results$log2FoldChange)
# pval_and_fc <- data.frame(row.names=row.names(results), results$pvalue, results$log2FoldChange)
# write.table(pval_and_fc, tablename2, quote=FALSE,sep="\t")

# store pvalue histogram
pdf(hist_path)
pvals <- results$pvalue
hist(pvals,breaks=100)
dev.off()


# store normalized counts
deseq2obj <- estimateSizeFactors(deseq2obj)
normcounts <- counts(deseq2obj, normalized=TRUE)
normcounts <- normcounts[sorted_genes,]
write.table(normcounts, file=normcounts_path, sep="\t", quote=FALSE)
