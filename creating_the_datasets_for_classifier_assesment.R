################################################################################
#
# creating the datasets for classifier assessment 
# LIDo rotation 2
# Cameron Ferguson, 03-04-2023
#
################################################################################

############################ setting up the session ############################

#importing the libraries 
library(Biostrings)
library(data.table)
library(parallel)

#setting the working directories 
setwd('~/Documents/LIDo rotation 2/NCBI data/sequnce_data/Data_for_classifier_assessment/')

#defining the functions needed 
mutate_nucleotide <- function(DNA_seq, diff = 0.1, pos = 3, prevent_stop_codons = T, random_position = F) {
  #defining the nucleotides 
  code <- c('A', 'T', 'G', 'C')
  
  #splitting the sequence into triplets
  if ((floor(length(DNA_seq[[1]])/3)*3) < length(DNA_seq[[1]])) {
    split_seq <- lapply(seq(1, (floor(length(DNA_seq[[1]])/3)*3), 3), function(x){
      DNA_seq[[1]][x:(x+2)]
    })
   split_seq <- append(split_seq, DNA_seq[[1]][-c(1:(floor(length(DNA_seq[[1]])/3)*3))])
    
  } else {
    split_seq <- lapply(seq(1, length(DNA_seq[[1]]), 3), function(x){
      DNA_seq[[1]][x:(x+2)]
    }) 
  }
  
  #using sample to randomly sample triplets 
  sites <- sample(1:((length(DNA_seq[[1]])/3)-1), 
                  size = ((length(DNA_seq[[1]])/3)-1)*diff,
                  replace = F)
  
  #randomly changing the nucleotides at those positions 
  for (x in sites) {
    if (random_position == T) {
      pos <- sample(c(1,2,3), 
                    size = 1, 
                    replace = F)
    }
    
    nucleotide <- as.character(split_seq[[x]][pos]) 
    
    if (nucleotide != 'A' & nucleotide != 'C' & nucleotide != 'G' & nucleotide != 'C') {
      split_seq[[x]][pos] <- sample(code,
                                    size = 1, 
                                    replace = F) 
      
    } else {
      split_seq[[x]][pos] <- sample(code[-grep(pattern = split_seq[[x]][pos], code)],
                                    size = 1, 
                                    replace = F) 
      
    }
    
    if (prevent_stop_codons == T) {
      if (as.character(split_seq[[x]]) == 'TAA' | as.character(split_seq[[x]]) == 'TAG' | as.character(split_seq[[x]]) == 'TGA') {
        if (pos == 1) {
          code2 <- c('A', 'G', 'C')
          split_seq[[x]][pos] <- sample(code2,
                                        size = 1, 
                                        replace = F) 
          
        } else {
          code2 <- c('T', 'C')
          split_seq[[x]][pos] <- sample(code2,
                                        size = 1, 
                                        replace = F) 
        }
      }
    }
  }
  
  #collapsing the nested DNAString object into 1 object 
  mut_DNA_seq <- do.call(c, split_seq)
  
  return(mut_DNA_seq)
  
}

mutate_homopolymer <- function(DNA_seq, diff = 0.002, min_length = 4) {
  #converting the dnaString object to a list of characters 
  DNA_list <- unlist(lapply(as.list(DNA_seq[[1]]), function(x){as.character(x)}))
  
  #using the run length encoding function to compute the lengths and values of 
  #runs of equal values in the converted dnaString vector
  rle_data <- rle(DNA_list)
  
  #converting the data to a data.table so its easier to work with 
  rle_data_table <- data.table(end_positions = cumsum(rle_data$lengths), #using the cumulative Sums function to add together the lengths
                               lengths = rle_data$lengths,
                               values = rle_data$values)
  
  #filtering the table for homopolymers equal too or grater then the given length
  rle_data_filt <- rle_data_table[lengths >= min_length, ,]
  
  #using sample to randomly sample homopolymers 
  sites <- sample(1:nrow(rle_data_filt), 
                  size = (nrow(rle_data_filt))*diff,
                  replace = F)
  
  #subsetting just the sites we wont to mutate
  rle_data_filt <- rle_data_filt[sites, ,]
  
  #sorting the table into ascending order based on the positions
  rle_data_filt <- rle_data_filt[order(end_positions), ,]
  
  #adjusting the positions for the frame shift 
  rle_data_filt$adjusted_pos <- unlist(lapply(1:nrow(rle_data_filt),
                                              function(x){
                                                rle_data_filt$end_positions[x] + (x-1)
                                              }))
  DNA_seq_mut <- DNA_seq[[1]]
  
  if (nrow(rle_data_filt) == 0) {
    return(DNA_seq_mut)
    
  } else {
    for (x in 1:nrow(rle_data_filt)) {
      temp <- rle_data_filt[x, ,]
      DNA_seq_mut <- append(DNA_seq_mut, 
                            temp$values, 
                            after = temp$adjusted_pos)
    }
    
    return(DNA_seq_mut)
    
  }
}


########################## mutating single nucleotides #########################

#looking at the distribution of the data it appears that ~98% of the data are 
#between 0-30% different from the refseq sequence. With ~92% of the data between
#0-20% different 

#creating a list of % differences 
diff_list <- seq(0.05, 0.3, 0.05)

#grabbing the names of the folders and files we need 
folders1 <- list.files()[2:6]
folders2 <- list.files('sym_data_1/SNPs/')[2:7]
files <- list.files('sym_data_1/SNPs/0_percent_different/')

#use forking here to run the loop in parallel
mclapply(folders1, function(Sym_set){
  
  for (i in 1:length(diff_list)) {
    x <- diff_list[i]
    y <- folders2[i]
    print(c(x, y))
    for (z in files) {
      DNA_sequence <- readDNAStringSet(paste(Sym_set,
                                             '/',
                                             'SNPs/0_percent_different/',
                                             z,
                                             sep = ''),
                                       format = 'fasta')
      
      mut_DNA_seq <- mutate_nucleotide(
        DNA_sequence,
        diff = x,
        prevent_stop_codons = T,
        random_position = T
      )
      
      header <- paste('> ',
                      names(DNA_sequence),
                      ' ',
                      x,
                      ' different')
      
      new_fasta <- paste(header,
                         mut_DNA_seq,
                         sep = '\n')
      
      write.table(
        new_fasta,
        paste(Sym_set,
              '/',
              'SNPs/',
              y,
              '/',
              z,
              sep = ''),
        quote = F,
        col.names = F,
        row.names = F,
        sep = '\n'
      )
      
    }
  }
  
}, mc.cores = 6, mc.preschedule = F)


###################### adding nucleotides to homopolymers ######################

#grabbing the names of the folders and files we need 
folders1 <- list.files()[2:6]
folders2 <- list.files('sym_data_1/SNPs/')
files <- list.files('sym_data_1/SNPs/0_percent_different/')

mclapply(folders1, function(Sym_set){
  
  for (x in 1:length(folders2)) {
    print(folders2[x])
    for (y in files) {
      DNA_sequence <- readDNAStringSet(paste(Sym_set,
                                             '/',
                                             'SNPs/',
                                             folders2[x],
                                             '/',
                                             y,
                                             sep = ''),
                                       format = 'fasta')
      
      mut_homopol_DNA_seq <- mutate_homopolymer(DNA_sequence)
      
      header <- paste('> ',
                      names(DNA_sequence),
                      ' with homopolymer errors',
                      sep = '')
      
      new_fasta <- paste(header,
                         mut_homopol_DNA_seq,
                         sep = '\n')
      
      write.table(
        new_fasta,
        paste(Sym_set,
              '/',
              'SNPS_and_homopolymers/',
              folders2[x],
              '/',
              y,
              sep = ''),
        quote = F,
        col.names = F,
        row.names = F,
        sep = '\n'
      )
    }
  }
  
}, mc.cores = 6, mc.preschedule = F)

