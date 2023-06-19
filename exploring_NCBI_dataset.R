################################################################################
#
# Exploring the NCBI sequence data
# Rotation 2
# Cameron Ferguson, 23-03-2023
#
################################################################################

############################ setting up the session ############################

#importing the libraries required 
library(data.table)
library(genbankr)
library(stringi)
library(ggplot2)

#setting the working directory 
setwd('~/Documents/LIDo rotation 2/NCBI data/')

#importing the dataset(s) 
NCBI_data <- fread('NCBI_data.tsv', header = T, sep = '\t')
extra_data <- fread('NCBI_extra_info.tsv', header = T, sep = '\t')


####################### exploring the inital NCBI dataset ######################

#what are the proportion of the genuses:
#creating a summary table
genus_summary <- data.table(Genus = names(summary(as.factor(NCBI_data$Genus))),
                            Number = summary(as.factor(NCBI_data$Genus)))

genus_summary$Genus[1] <- 'Unknown'
genus_summary$labels <- paste(round((genus_summary$Number/4360)*100, digits = 2), 
                              '%', 
                              sep = '')

#visualizing the data
ggplot(genus_summary, aes(x = '', y = Number, fill = Genus)) + geom_col() +
  coord_polar(theta = "y") + 
  theme_void() +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.475)) 


################## exploring the extra data from the gbk files #################

#creating a summary table looking at subgenus
subgenus_summary <- data.table(Subgenus = names(summary(as.factor(extra_data$subgenus))),
                               Number = summary(as.factor(extra_data$subgenus)))

#cleaning up the dataset 
subgenus_summary$Subgenus[1] <- 'Unknown'

#combining rows with the same type of info but different names 
subgenus_summary$Number[1] <- sum(subgenus_summary$Number[c(1:4, 23, 30)])
subgenus_summary <- subgenus_summary[!c(2:4, 23, 30), ,]

#visualizing the data
ggplot(subgenus_summary, aes(x = Subgenus, y = (subgenus_summary$Number/4360)*100)) +
  geom_bar(stat = "identity", fill = 'dark gray') + 
  theme_classic() + 
  ylab('Percentage of sequences (%)') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#creating a summary table looking at whether the sequences is refseq or not 
keyword_summary <- data.table(Keyword = names(summary(as.factor(extra_data$keywords))),
                              Number = summary(as.factor(extra_data$keywords)))

#based off inspection of the summary table only 64 samples are refseq with the vast majority not being
#specified at all.
#listing the accessions of the refseq files:
extra_data[grep('RefSeq', extra_data$keywords), ,]


######## looking at the subgenus breakdown of the spike protein samples ########

#creating a list of the spike protein accessions 
files <- unlist(lapply(c('sequnce_data/protein_seq/Alphacoronavirus/spike/',
                         'sequnce_data/protein_seq/Alphapironavirus/spike/',
                         'sequnce_data/protein_seq/Betacoronavirus/spike/',
                         'sequnce_data/protein_seq/Deltacoronavirus/spike/',
                         'sequnce_data/protein_seq/Gammacoronavirus/spike/',
                         'sequnce_data/protein_seq/Unkown/spike/'), 
                       function(x){list.files(x)}))

file_acc <- stri_split(files, 
                       fixed = '.', 
                       simplify = T)[, 1]

#subsetting the extra data 
sub_extra_data <- extra_data[unlist(lapply(file_acc, function(x){
  
  grep(x, extra_data$accession)
  
})), ,]

#creating a summary table 
spike_subgenus_summary <- data.table(Subgenus = names(summary(as.factor(sub_extra_data$subgenus))),
                               Number = summary(as.factor(sub_extra_data$subgenus)))

#cleaning up the dataset 
spike_subgenus_summary$Subgenus[1] <- 'Unknown'

#combining rows with the same type of info but different names 
spike_subgenus_summary$Number[1] <- sum(spike_subgenus_summary$Number[c(1:4, 23, 30)])
spike_subgenus_summary <- spike_subgenus_summary[!c(2:4, 23, 30), ,]


######### looking at the subgenus breakdown of the orf protein samples #########

#creating a list of the spike protein accessions 
files <- unlist(lapply(c('sequnce_data/protein_seq/Alphacoronavirus/ORF1ab/',
                         'sequnce_data/protein_seq/Alphapironavirus/ORF1ab/',
                         'sequnce_data/protein_seq/Betacoronavirus/ORF1ab/',
                         'sequnce_data/protein_seq/Deltacoronavirus/ORF1ab/',
                         'sequnce_data/protein_seq/Gammacoronavirus/ORF1ab/',
                         'sequnce_data/protein_seq/Unkown/ORF1ab/'), 
                       function(x){list.files(x)}))

file_acc <- stri_split(files, 
                       fixed = '.', 
                       simplify = T)[, 1]

#subsetting the extra data 
sub_extra_data <- extra_data[unlist(lapply(file_acc, function(x){
  
  grep(x, extra_data$accession)
  
})), ,]

#creating a summary table 
ORF1ab_subgenus_summary <- data.table(Subgenus = names(summary(as.factor(sub_extra_data$subgenus))),
                                     Number = summary(as.factor(sub_extra_data$subgenus)))

#cleaning up the dataset 
ORF1ab_subgenus_summary$Subgenus[1] <- 'Unknown'


#combining the fasta files 
first_alignment <- sub_extra_data[subgenus == 'Sarbecovirus', ,]

for (x in first_alignment$accession) {
  fasta <- Biostrings::readAAStringSet(paste('sequnce_data/protein_seq/Betacoronavirus/spike/',
                                             x,
                                             '.fasta',
                                             sep = ''), 
                                       format = 'fasta')
  
  Biostrings::writeXStringSet(fasta,
                              'alignment/combined_Betacoronavirus_Sarbecovirus.fasta',
                              format = 'fasta',
                              append = T)
}






