message ("Loading BioNet...")

library(BioNet)

pval_file <- snakemake@input[["deseq_out"]]
output_prefix <- paste("scores", snakemake@wildcards[["experiment"]], sep="/")
examine_prefix <- paste("eXamine", snakemake@wildcards[["experiment"]], sep="/")

deseq_tab <- read.table(pval_file, header=TRUE, row.names=1, sep="\t")
pvals <- deseq_tab[,"pvalue", drop=FALSE]
head(pvals)

#head(pvals)

# BioNet does not like p-values with p = 0
# set 0 pvals to next smallest value
pvals$pvalue[pvals$pvalue == 0] <- unique(sort(pvals$pvalue))[2]

head(pvals)

p_table <- paste(output_prefix, "_pvals.txt", sep="")
write.table(data.frame("geneID"=rownames(pvals), pvals), p_table, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

message("Fitting the bum model to the p-value distribution...")
bum <- fitBumModel(pvals, plot=FALSE)


mat <- c(bum$lambda, bum$a)
bumtablename <- paste(output_prefix, "_bum_fit.txt", sep="")
write.table(x=mat, file=bumtablename,sep="",quote=FALSE, row.names=FALSE, col.names=FALSE)
message ("done.")
