library("BiocManager")
library("devtools")
library("rBLAST")
library("biomartr")
library("Biostrings")
library("GenomicRanges")
library("IRanges")
library("S4Vectors")
library("seqinr")
library("utils")
library("metagenomeFeatures")
library("qiime2R")
library("rentrez")
library("DBI")
library("XML")

all_seqs <- list()
all_names <- c()

for (f in list.files("/Users/mattracz/Projects/Microbial_real/AA_seqs/kept", full.names=TRUE)){
  #seqs <- read.fasta(f, as.string = TRUE)
  seqname <- basename(f)
  #print(seqname)
  aa <- toupper(read.fasta(f, as.string = TRUE))
  names(aa) <- rep(seqname, length(aa))
  all_seqs <- c(all_seqs, aa)           # Add sequences to the list
  all_names <- c(all_names, names(aa))
}

write.fasta(sequences = all_seqs, names = all_names, file.out = "/Users/mattracz/Projects/Microbial_real/FINALSEQS")

