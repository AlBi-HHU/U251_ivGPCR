
#Um das aufuehren zu koennen muessen folgende Packte installiert sein
# g++, gfortran, gcc, libblas-dev, libgfortran, lipopenblas
if(!require(biomaRt)){
  suppressMessages(suppressWarnings(source("https://bioconductor.org/biocLite.R")))
  suppressMessages(suppressWarnings(biocLite("biomaRt")))
}

suppressMessages(library(biomaRt))
mart_hsa <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl", version = snakemake@config[["ENSEMBL"]][["version"]])
go = getBM(attributes=c('hgnc_symbol','go_id','go_linkage_type','namespace_1003','name_1006','definition_1006'),filters='', mart=mart_hsa)

#go = getBM(attributes=c('hgnc_symbol','go_id'),filters='',values = '', mart=mart_hsa)

write.table(go, sep = "\t", na = "NA", dec= ".", col.names = FALSE, row.names = FALSE, quote = FALSE)
