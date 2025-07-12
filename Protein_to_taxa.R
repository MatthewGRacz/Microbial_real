library("metagenomeFeatures")
library("qiime2R")
library("Biostrings")
library("rBLAST")
library("BiocManager")
library("devtools")
library("biomartr")
library("GenomicRanges")
library("IRanges")
library("S4Vectors")
library("seqinr")
library("utils")
Sys.setenv(PATH = paste("/Users/mattracz/Projects/Microbial_real/ncbi-blast-2.16.0+/bin", Sys.getenv("PATH"), sep=":"))
#allows BLAST+ package to run

#system("which blastn")
#system("which blastp")

blast_setup <- blast(db = "/Users/mattracz/Projects/Microbial_real/BLAST+databases/nr/nr.000", remote = FALSE, type = "blastp")
#sets up blast database for 16S RNA

get_order <- function(){
  
  #files in AA_seqs that are not directories (actual files)
  files <- list.files("/Users/mattracz/Projects/Microbial_real/AA_seqs/", full.names=TRUE)[file.info(list.files("/Users/mattracz/Projects/Microbial_real/AA_seqs/", full.names=TRUE))$isdir == FALSE]  # Filter out directories
  
  for(s in files){
    #for each AA sequence just created and stored in AA_seqs
    current_seq <- read.fasta(s) #AA sequence used
    current_seq <- paste(current_seq[[1]], collapse = "") #get AA sequence as an element
    
    AAseq <- AAStringSet(x=current_seq, start=1, end=nchar(current_seq), width=NA, use.names=TRUE)
    #make the AA sequence an AAStringSet so that BLAST can operate on it
    
    blast_data <- predict(blast_setup, AAseq) #gets data of protein
    #print(blast_data) #print results of blast search on protein
    
    if(nrow(blast_data)==0){ #if no match
      print("nothing")
      individual_file_name <- strsplit(basename(s), "/Users/mattracz/Projects/Microbial_real/AA_seqs/")[[1]]
      file.rename(from=s, to=paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/rejected/", individual_file_name,  sep=""))
      next
    }
    
    
    taxsearch <- entrez_search(db="nuccore", term=blast_data$sseqid[1])
    tax_id <- entrez_summary(db="nuccore", id=taxsearch$ids[1])$taxid
    #print(tax_id)#print that taxonomic ID
    
    taxonomic_info <- entrez_fetch(db="taxonomy", id=tax_id, rettype="xml", parsed=TRUE)
    #then fetch the full taxonomic lineage from that record
    
    #print(taxonomic_info) #View full taxonomic lineage
    
    order <- xpathSApply(taxonomic_info, "//Taxon[Rank='order']/ScientificName", xmlValue)
    #get the family of its taxonomic lineage and convert it to text content (not necessarily a string, might be a character code)
    
    #print(family) #print the family's name
    
    print(order)
    
    individual_file_name <- strsplit(basename(s), "/Users/mattracz/Projects/Microbial_real/AA_seqs/")[[1]]
    file.rename(from=s, to=paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/kept/", individual_file_name,  sep=""))
    
    
  } 
  
}

get_order()


