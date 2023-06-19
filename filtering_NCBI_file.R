################################################################################
#
# filtering the NCBI sequences
# Rotation 2
# Cameron Ferguson, 20-03-2023
#
################################################################################

############################ setting up the session ############################

#importing the libraries required 
library(data.table)
library(genbankr)
library(stringi)
library(rentrez)
library(ggplot2)

#setting the working directory 
setwd('~/Documents/LIDo rotation 2/NCBI data/')

#importing the dataset 
NCBI_data <- fread('NCBI_data.tsv', header = T, sep = '\t')


############################ exploring the dataset #############################

#what are the proportion of the genuses:
#creating a summary table
genus_summary <- data.table(Genus = names(summary(as.factor(NCBI_data$Genus))),
                            Number = summary(as.factor(NCBI_data$Genus)))

genus_summary$Genus[1] <- 'Unknown'
genus_summary$labels <- paste(round((genus_summary$Number/4360)*100, digits = 2), '%', sep = '')

#visualizing the data
ggplot(genus_summary, aes(x = '', y = Number, fill = Genus)) + geom_col() +
  coord_polar(theta = "y") + 
  theme_void() +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.475)) 


########################## downloading the sequences ###########################

lapply(1:length(NCBI_data$Accession), function(x){
  
  seq_data <- entrez_fetch(db = "nucleotide", 
                           id = NCBI_data$Accession[x],
                           rettype = "fasta")
  write.table(seq_data, 
              paste('sequnce_data/Full_sequences/', 
                    NCBI_data$Accession[x], 
                    '.fasta', 
                    sep = ''), 
              
              sep = '\n', 
              quote = F, 
              row.names = F, 
              col.names = F)
  
  print(paste('Acession', 
              x, 
              ':',
              NCBI_data$Accession[x],
              sep = ' '))
  
  print(paste((x/4360)*100, 
              '%', 
              'complete', 
              sep = ' '))
  
})


######################### downloading the annotations ##########################

lapply(1:length(NCBI_data$Accession), function(x){
  
  seq_data <- entrez_fetch(db = "nucleotide", 
                           id = NCBI_data$Accession[x],
                           rettype = "gbwithparts")
  
  if (NCBI_data$Genus[x] == 'Alphacoronavirus') {
    write.table(seq_data, 
                paste('sequnce_data/Annotation/Alphacoronavirus/', 
                      NCBI_data$Accession[x], 
                      '.gbk', 
                      sep = ''), 
                
                sep = '\n', 
                quote = F, 
                row.names = F, 
                col.names = F)
    
  } else if (NCBI_data$Genus[x] == 'Alphapironavirus') {
    write.table(seq_data, 
                paste('sequnce_data/Annotation/Alphapironavirus/', 
                      NCBI_data$Accession[x], 
                      '.gbk', 
                      sep = ''), 
                
                sep = '\n', 
                quote = F, 
                row.names = F, 
                col.names = F)
    
  } else if (NCBI_data$Genus[x] == 'Betacoronavirus') {
    write.table(seq_data, 
                paste('sequnce_data/Annotation/Betacoronavirus/', 
                      NCBI_data$Accession[x], 
                      '.gbk', 
                      sep = ''), 
                
                sep = '\n', 
                quote = F, 
                row.names = F, 
                col.names = F)
    
  } else if (NCBI_data$Genus[x] == 'Deltacoronavirus') {
    write.table(seq_data, 
                paste('sequnce_data/Annotation/Deltacoronavirus/', 
                      NCBI_data$Accession[x], 
                      '.gbk', 
                      sep = ''), 
                
                sep = '\n', 
                quote = F, 
                row.names = F, 
                col.names = F)
    
  } else if (NCBI_data$Genus[x] == 'Gammacoronavirus') {
    write.table(seq_data, 
                paste('sequnce_data/Annotation/Gammacoronavirus/', 
                      NCBI_data$Accession[x], 
                      '.gbk', 
                      sep = ''), 
                
                sep = '\n', 
                quote = F, 
                row.names = F, 
                col.names = F)
    
  } else {
    write.table(seq_data, 
                paste('sequnce_data/Annotation/Unkown/', 
                      NCBI_data$Accession[x], 
                      '.gbk', 
                      sep = ''), 
                
                sep = '\n', 
                quote = F, 
                row.names = F, 
                col.names = F)
    
  } 
  
  print(paste('Acession', 
              x, 
              ':',
              NCBI_data$Accession[x],
              sep = ' '))
  
  print(paste((x/4360)*100, 
              '%', 
              'complete', 
              sep = ' '))
  
})


################ pulling out the proteins we are interested in #################

for (file in list.files(paste('sequnce_data/Annotation', '/', sep = ''))) {
  
  lapply(list.files(paste('sequnce_data/Annotation/', file, '/', sep = '')), function(x){
    
    gbk_file <- parseGenBank(paste('sequnce_data/Annotation/', file, '/', x, sep = ''))
    
    gene_of_intrest <- unlist(gbk_file$FEATURES[grep('\\bORF1ab\\b',
                                                     gbk_file$FEATURES)])
    
    if (is.null(gene_of_intrest) == T) {
      print(paste(file,
                  ',',
                  x,
                  'does not contain this protein',
                  sep = ' '))
      
    } else {
      
      gene_of_intrest <- matrix(c(names(gene_of_intrest), gene_of_intrest), ncol = 2)
      
      amino_acid_seq <- gene_of_intrest[grep(pattern = 'translation', gene_of_intrest), 2]
      
      seqID <- stri_split(x, fixed = '.', simplify = T)[,1]
      
      fasta_header <- paste(paste('>', seqID), 
                            gene_of_intrest[grep(pattern = 'seqnames', gene_of_intrest), 2][1], 
                            gene_of_intrest[grep(pattern = 'product', gene_of_intrest), 2][1],
                            sep = ', ')
      
      fasta_data <- paste(fasta_header, amino_acid_seq[1], sep = '\n')
      
      write.table(fasta_data, 
                  paste('sequnce_data/protein_seq/',
                        file,
                        '/ORF1ab/', 
                        seqID, 
                        '.fasta', 
                        sep = ''), 
                  
                  sep = '\n', 
                  quote = F, 
                  row.names = F, 
                  col.names = F)
      
    }
    
  })
}





for (file in list.files(paste('sequnce_data/Annotation', '/', sep = ''))) {
  
  test <- lapply(list.files(paste('sequnce_data/Annotation/', file, '/', sep = '')), function(x){
    
    gbk_file <- parseGenBank(paste('sequnce_data/Annotation/', 
                                   file, '/', 
                                   x, 
                                   sep = ''))
      
      taxa_list <- gbk_file$SOURCE$lineage[10:length(gbk_file$SOURCE$lineage)]
      
      #for whether the strand is refseq or not look at gbk_file$KEYWORDS
      refseq <- gbk_file$KEYWORDS
      
      taxa_list[length(taxa_list)] <- paste(taxa_list[length(taxa_list)],
                                            refseq,
                                            x,
                                            sep = '£')
      
      temp_matrix <- matrix(taxa_list,
                            ncol = length(taxa_list),
                            byrow = T)
      
      write.table(temp_matrix, 
                  'sequnce_data/protein_seq/test.csv', 
                  append = T,
                  sep = ',', 
                  quote = F, 
                  row.names = F, 
                  col.names = F)
      
      
    
    
  })
  
}

test <- fread('sequnce_data/protein_seq/test.csv', header = F, sep = '\n')
test2 <- stri_split(test$V1, regex = ',', simplify = T)

extra_info <- data.table(accession = stri_split(stri_split(test$V1, 
                                                           regex = '£', 
                                                           simplify = T)[, 3], 
                                                fixed = '.',
                                                simplify = T)[, 1],
                         
                         genus = stri_split(test2[, 1], 
                                            fixed = '.', 
                                            simplify = T)[, 1], 
                         
                         subgenus = stri_split(test2[, 2], 
                                               fixed = '.', 
                                               simplify = T)[, 1],
                         
                         keywords = stri_split(test$V1, 
                                               regex = '£', 
                                               simplify = T)[, 2])

fwrite(extra_info, 
       file = 'NCBI_extra_info.tsv', 
       quote = F, 
       sep = '\t')

