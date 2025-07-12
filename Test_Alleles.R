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
packages.install()



Sys.setenv(PATH = paste("/Users/mattracz/Projects/Microbial_real/ncbi-blast-2.16.0+/bin", Sys.getenv("PATH"), sep = ":"))
#allows BLAST+ package to run
setwd("/Users/mattracz/Projects/Microbial_real") #sets scope as Microbial_real

#system("which blastn")
#system("which blastp")

blast_setup <- blast(db = "/Users/mattracz/Projects/Microbial_real/BLAST+databases/nr/nr.000", remote = FALSE, type = "blastp")
#initializes blast nr database for proteins

alphabetadata <- function(ABdata){ #analyzes AB Allele Protein data
  
  cut_point = 82
  proteindata <- read.csv(ABdata, header=TRUE)
  num_sequences <- nrow(proteindata)/2
  #head(proteindata)

  for(i in 1:num_sequences){ #for each unique allele
    
    alpha_one <- substr(proteindata[(2*i)-1, 4], start=1, stop=cut_point)
    #print(alpha_one)
  
    beta_one <- substr(proteindata[(2*i)-1, 4], start=cut_point+1, stop=nchar(proteindata[(2*i)-1, 4]))
    #print(beta_one)
  
    alpha_two <- substr(proteindata[(2*i), 4], start=1, stop=cut_point)
    #print(alpha_two)
  
    beta_two <- substr(proteindata[(2*i), 4], start=cut_point+1, stop=nchar(proteindata[2*i, 4]))
    #print(beta_two)
    
    a1b1 <- paste(alpha_one, beta_one, sep="") #alpha one + beta one
    #a1b1 <- gsub(pattern="\\*", replacement="", a1b1) #replace any missing/corrupted proteins with a blank, BLAST can still recognize them (else it throws an error due to incompatibility of * with any amino acid)
    filepath <- paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/", "a1b1_allele", i, sep="")
    #saves to DNA seqs under name a1b1_allele + the number allele it's for
    
    write.fasta(
      sequences = list(strsplit(a1b1, "")[[1]]), #returns the single element of what is in the first row, but the only data in it is just the amino acid sequence, so that is returned (""+AA, AA as an element is returned, rather than as a string or with any formatting problems)
      names="", #if you add a name, it'll come before the AA sequence in the fasta file and cause problems
      file.out = filepath #Plus, the files are already named directly via their file paths
    )
  
    a2b1 <- paste(alpha_two, beta_one, sep="") #repeat for other three combined AA sequences for the allele
    #a2b1 <- gsub(pattern="\\*", replacement="", a2b1)
    filepath <- paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/", "a2b1_allele", i, sep="")
  
    write.fasta(
      sequences = list(strsplit(a2b1, "")[[1]]),
      names="",
      file.out = filepath
    )
  
    a1b2 <- paste(alpha_two, beta_one, sep="")
    #a1b2 <- gsub(pattern="\\*", replacement="", a1b2)
    filepath <- paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/", "a1b2_allele", i, sep="")
    
    write.fasta(
      sequences = list(strsplit(a1b2, "")[[1]]),
      names="",
      file.out = filepath
    )
    
    a2b2 <- paste(alpha_two, beta_one, sep="")
    #a2b2 <- gsub(pattern="\\*", replacement="", a2b2)
    filepath <- paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/", "a2b2_allele", i, sep="")
    
    write.fasta(
      sequences = list(strsplit(a2b2, "")[[1]]),
      names="",
      file.out = filepath
    )
  } 
}

#alphabetadata("AB Alleles Protein.csv")

syngnathidae_checker <- function(){
  
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
    
    tax_id <- entrez_summary(db="protein", id=entrez_search(db="protein", term=blast_data$sseqid[[1]])[[1]])$taxid
    #print(tax_id)#print that taxonomic ID
    
    taxonomic_info <- entrez_fetch(db="taxonomy", id=tax_id, rettype="xml", parsed=TRUE)
    #then fetch the full taxonomic lineage from that record
    
    #print(taxonomic_info) #View full taxonomic lineage
    
    family <- xpathSApply(taxonomic_info, "//Taxon[Rank='family']/ScientificName", xmlValue)
    #get the family of its taxonomic lineage and convert it to text content (not necessarily a string, might be a character code)
    
    #print(family) #print the family's name
    
    if(!("Syngnathidae" %in% family)){ #if family is not a Syngnathidae
      individual_file_name <- strsplit(basename(s), "/Users/mattracz/Projects/Microbial_real/AA_seqs/")[[1]]
      file.rename(from=s, to=paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/rejected/", individual_file_name,  sep=""))
      #get individual file's name, then move it to rejected directory
      
      print("Not a Syngnathid!")
      
    } else{ #if family is a Syngnathidae
      individual_file_name <- strsplit(basename(s), "/Users/mattracz/Projects/Microbial_real/AA_seqs/")[[1]]
      file.rename(from=s, to=paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/kept/", individual_file_name,  sep=""))
      #get individual file's name, then move it to kept directory
      
      print("Behold a Syngnathid!")
    }
    
    
  } 
  #Syngnathidae checker ends, with all Syngnathidae in ./kept, and all else in ./rejected
  
}

syngnathidae_checker()


