#!/usr/local/bin/Rscript
# written by Mohammed El-Kebir
# adapted by Gunnar W. Klau

if(!require(topGO)){
  suppressMessages(suppressWarnings(source("https://bioconductor.org/biocLite.R")))
  suppressMessages(suppressWarnings(biocLite("topGO")))
}

suppressMessages(library(topGO))

args <- commandArgs(trailingOnly = TRUE)
mappingFilename <- args[1]
geneListFilename <- args[2]
go_domain <- args[3]
output <- args[4]

geneID2GO <- readMappings(mappingFilename)
geneNames <- names(geneID2GO)

myInterestingGenes <- as.character(read.delim(geneListFilename, header=FALSE)$V1)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

GOdata <- new("topGOdata", ontology = go_domain, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

results <- runTest(GOdata,algorithm="classic",statistic="fisher")
table <- GenTable(GOdata, classicFisher=results,topNodes=20)


head(table,n=10)

#table2 <- cbind(table$GO.ID, table$classicFisher)
#write.table(table2, file=output, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

table2 <- cbind(table$GO.ID, table$classicFisher, table$Term) #, table$Annotated, table$Significant, table$Expected
write.table(table2, file=output, col.names = c("ID", "P-value", "Label"), row.names = FALSE, quote = FALSE, sep = '\t')
