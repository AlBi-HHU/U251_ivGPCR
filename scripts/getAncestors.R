#!/usr/local/bin/Rscript

if (!require(GO.db)) {
  suppressMessages(suppressWarnings(source(
    "https://bioconductor.org/biocLite.R"
  )))
  suppressMessages(suppressWarnings(biocLite("GO.db")))
}
suppressMessages(library("GO.db"))
args <- commandArgs(trailingOnly = TRUE)

bp <- as.list(GOBPANCESTOR)
mf <- as.list(GOMFANCESTOR)
cc <- as.list(GOCCANCESTOR)

#fc <- file(args[1])
#a <- readLines(fc)
b <- read.table(args[1], sep = "\t", fill = T)

for (i in 1:nrow(b)) {
  if (b[i, "V2"] == "biological_process") {
    bp_list <- paste(unlist(bp[b[i, "V1"]]), collapse = "\t")
    cat(sprintf("%s\t%s\t%s\n", b[i, "V1"], b[i, "V2"], bp_list))
    
  }
  if (b[i, "V2"] == "molecular_function") {
    mf_list <- paste(unlist(mf[b[i, "V1"]]), collapse = "\t")
    cat(sprintf("%s\t%s\t%s\n", b[i, "V1"], b[i, "V2"], mf_list))
    
  }
  if (b[i, "V2"] == "cellular_component") {
    cc_list <- paste(unlist(cc[b[i, "V1"]]), collapse = "\t")
    cat(sprintf("%s\t%s\t%s\n", b[i, "V1"], b[i, "V2"], cc_list))
    
  }
}

#-----original------
# 
# for (x in a)
# {
#   if (x == "") break
#   b <- strsplit(x, split="\t")[[1]]
#   cat(b[1],sep="")
#   cat("\t",sep="")
#   cat(b[2],sep="")
# 
#   if (b[2] == "biological_process")
#   {
#     c <- bp[b[1]]
#     for (x in c) cat(paste0("\t",x),sep="")
#     cat("\n",sep="")
#   }
#   if (b[2] == "molecular_function")
#   {
#     c <- mf[b[1]]
#     for (x in c) cat(paste0("\t",x),sep="")
#     cat("\n",sep="")
#   }
#   if (b[2] == "cellular_component")
#   {
#     c <- cc[b[1]]
#     for (x in c) cat(paste0("\t",x),sep="")
#     cat("\n",sep="")
#   }
# }
