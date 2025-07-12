library("seqinr")
nfiles = 10

rna_generate <-function(rand_sequence, seqlength){
  for(i in 1:seqlength) {

    n=sample(1:4, 1, replace=TRUE)
    
    if(n ==1) {
      rand_sequence <- c(rand_sequence, "A")  # Add "A" to the sequence
    } else if(n ==2) {
      rand_sequence <- c(rand_sequence, "G")  # Add "G" to the sequence
    } else if(n ==3) {
      rand_sequence <- c(rand_sequence, "U")  # Add "U" to the sequence
    } else if(n ==4) {
      rand_sequence <- c(rand_sequence, "C")  # Add "C" to the sequence
    }
  }
  return(rand_sequence)
}

dna_generate <-function(rand_sequence, seqlength){
  for(i in 1:seqlength) {

    n=sample(1:4, 1, replace=TRUE)
    
    if(n ==1) {
      rand_sequence <- c(rand_sequence, "A")  # Add "A" to the sequence
    } else if(n ==2) {
      rand_sequence <- c(rand_sequence, "G")  # Add "G" to the sequence
    } else if(n ==3) {
      rand_sequence <- c(rand_sequence, "T")  # Add "T" to the sequence
    } else if(n ==4) {
      rand_sequence <- c(rand_sequence, "C")  # Add "C" to the sequence
    }
  }
  
  return(rand_sequence)
}

aa_generate <-function(rand_sequence, seqlength){
  for(i in 1:seqlength) {

    n=sample(1:20, 1, replace=TRUE)
    
    if(n ==1) {
      rand_sequence <- c(rand_sequence, "A")  # Add "A" to the sequence
    } else if(n ==2) {
      rand_sequence <- c(rand_sequence, "C")  # Add "C" to the sequence
    } else if(n ==3) {
      rand_sequence <- c(rand_sequence, "D")  # Add "D" to the sequence
    } else if(n ==4) {
      rand_sequence <- c(rand_sequence, "E")  # Add "E" to the sequence
    }else if(n ==5) {
      rand_sequence <- c(rand_sequence, "F")  # Add "F" to the sequence
    } else if(n ==6) {
      rand_sequence <- c(rand_sequence, "G")  # Add "G" to the sequence
    } else if(n ==7) {
      rand_sequence <- c(rand_sequence, "H")  # Add "H" to the sequence
    }else if(n ==8) {
      rand_sequence <- c(rand_sequence, "I")  # Add "I" to the sequence
    } else if(n ==9) {
      rand_sequence <- c(rand_sequence, "K")  # Add "K" to the sequence
    } else if(n ==10) {
      rand_sequence <- c(rand_sequence, "L")  # Add "L" to the sequence
    }else if(n ==11) {
      rand_sequence <- c(rand_sequence, "M")  # Add "M" to the sequence
    } else if(n ==12) {
      rand_sequence <- c(rand_sequence, "N")  # Add "N" to the sequence
    } else if(n ==13) {
      rand_sequence <- c(rand_sequence, "P")  # Add "P" to the sequence
    }else if(n ==14) {
      rand_sequence <- c(rand_sequence, "Q")  # Add "Q" to the sequence
    } else if(n ==15) {
      rand_sequence <- c(rand_sequence, "R")  # Add "R" to the sequence
    } else if(n ==16) {
      rand_sequence <- c(rand_sequence, "S")  # Add "S" to the sequence
    }else if(n ==17) {
      rand_sequence <- c(rand_sequence, "T")  # Add "T" to the sequence
    } else if(n ==18) {
      rand_sequence <- c(rand_sequence, "V")  # Add "V" to the sequence
    } else if(n ==19) {
      rand_sequence <- c(rand_sequence, "W")  # Add "W" to the sequence
    }else if(n ==20) {
      rand_sequence <- c(rand_sequence, "Y")  # Add "Y" to the sequence
    }
  }
  
  return(rand_sequence)
}


for(filenum in 1:nfiles){
  rand_sequence <- character()  # Create an empty character vector for the RNA sequence
  seqlength = 10000000

  #rand_sequence <- aa_generate(rand_sequence, seqlength)
  rand_sequence <- rna_generate(rand_sequence, seqlength)
  #rand_sequence <- dna_generate(rand_sequence, seqlength)
  
  # Convert the character vector to a single string
  rand_string <- paste(rand_sequence, collapse = "")
  
  filepath <- paste("/Users/mattracz/Projects/Microbial_real/RNA_seqs/", "random_gen_sequence", filenum, ".fasta", sep="")
  #filepath <- paste("/Users/mattracz/Projects/Microbial_realAA_seqs/", "random_gen_sequence", filenum, ".fasta", sep="")
  
  #print(filepath)
  
  # Write the RNA sequence to a FASTA file
  write.fasta(sequence=rand_string, names=NULL, file.out=filepath, open="w")

}



