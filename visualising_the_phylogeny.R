################################################################################
#
# visualising the phylogeny 
# Rotation 2
# Cameron Ferguson, 07-04-2023
#
################################################################################

############################ setting up the session ############################

#importing the libraries required 
library(RColorBrewer)
library(data.table)
library(stringi)
library(ggplot2)
library(ggtree)
library(ape)

#setting the working directory 
setwd('~/Documents/LIDo rotation 2/NCBI data/')

#importing the dataset(s) 
taxonomy_data <- fread('NCBI_extra_info.tsv', header = T, sep = '\t')


######################## data cleaning and preparation #########################

#importing the IQ tree
tree <- read.tree('phylogeny/spike/test_tree')

#changing the metadata labels to match the tree tip labels 
taxonomy_data2 <- data.table(taxonomy_data[, 1,],
                            
                            tree_labs = unlist(lapply(taxonomy_data[, 1,], 
                                                      function(x) {
                                                        paste(x, '_', sep = '')
                                                      })),
                            
                            genus = unlist(lapply(taxonomy_data[, 2,], 
                                                  function(x) {
                                                    ifelse(is.na(x) == T,
                                                           yes = 'Unknown',
                                                           no = x)
                                                  })),
                              
                            subgenus = unlist(lapply(taxonomy_data[, 3,], 
                                                     function(x) {
                                                       ifelse(is.na(x) == T,
                                                              yes = 'Unknown',
                                                              no = x)
                                                     })),
                            
                            taxonomy_data[, 4,])


#matching tips of tree to metadata and use to construct annotation dataframe
metaMatch <- taxonomy_data[match(tree$tip.label, taxonomy_data$tree_labs), , ] 

#using data.frame over data.table as data.table doesn't suport rownames which we need for the heatmap
dfMatch <- data.frame(Name = metaMatch$accession,
                      genus = metaMatch$genus, 
                      subgenus = metaMatch$subgenus)

rownames(dfMatch) <- tree$tip.label

tree$tip.label <- dfMatch$Name

#cleaning up the keywords column
x <- stri_split(stri_split(metaMatch$keywords, 
                           fixed = ';',
                           simplify = T)[,1], 
                fixed = '.', 
                simplify = T)[,1]

x <- ifelse(x != 'RefSeq',
            yes = 'None-RefSeq',
            no = x)

metaMatch$RefSeq <- ifelse(is.na(x) == T,
                           yes = 'None-RefSeq',
                           no = x)

#creating an additional dataframe to be used by gheatmap for annotating host species
heat_Map_Anotation <- data.frame(RefSeq = as.factor(metaMatch$RefSeq))
row.names(heat_Map_Anotation) <- tree$tip.label

#removing a sample that im not sure about 
tree.root <- drop.tip(tree, 'MK611985')


############################### Plotting the tree ##############################

#creating the initial plot 
p <- ggtree(tree.root, layout = 'circular')

#creating the colour pallet 
col <- brewer.pal(5, "Dark2")

#adding in the tip labels 
p <- p %<+% dfMatch + geom_tippoint(aes(color = genus), size = 1.4) +
  scale_color_manual(values = c(col[1:2], col[4], col[5], col[3]))

#annotating the tree with a heatmap 
gheatmap(p, heat_Map_Anotation, width = 0.05, colnames_angle = 90, hjust = 1, font.size = 3) + 
  guides(fill = guide_legend(title = "RefSeq")) + 
  theme(plot.margin = margin(2, 2, 30, 2)) +
  scale_fill_manual(values = c("light grey", "red")) 


######################### looking at pairwise distance #########################

#importing the alignment
alignment <- read.FASTA('alignment/spike/DNA/combined_genus_aligned_filtered.fasta', 
                        type = 'DNA')

#simplifying the alignment names 
names(alignment) <- stri_split(names(alignment),
                               fixed = ' ',
                               simplify = T)[,1]

#creating a summary table of the refseq sequences per genus
metaMatch[, .(names(summary(as.factor(RefSeq))), 
              summary(as.factor(RefSeq))), 
          by = genus]

#calculating the pairwise distance between sequences and closest refseq file
alignment_dist <- dist.dna(alignment,
                           model = 'raw', 
                           as.matrix = T, 
                           pairwise.deletion = T)


#pulling out the accessions for the refseq samples 
RefSeq_sequences <- metaMatch[RefSeq == 'RefSeq', accession,]

#subsetting the distance matrix for the RefSeq samples 
alignment_RefSeq_dist <- alignment_dist[,unlist(lapply(RefSeq_sequences, 
                                                       function(x){
                                                         grep(pattern = x, 
                                                              colnames(alignment_dist))
                                                         }))]

#creating a table of samples and there closest RefSeq samples with distances 
closest_RefSeq <- data.table(accession = unlist(lapply(1:nrow(alignment_RefSeq_dist),
                                                       function(x){
                                                         rownames(alignment_RefSeq_dist)[x]
                                                       })),
                             
                             genus = unlist(lapply(1:nrow(alignment_RefSeq_dist), 
                                                   function(x){
                                                     metaMatch[accession == stri_split(rownames(alignment_RefSeq_dist)[x], 
                                                                                       fixed = '.', 
                                                                                       simplify = T)[,1], genus,]
                                                   })),
                             
                             RefSeq_samples = lapply(1:nrow(alignment_RefSeq_dist), 
                                                     function(x){
                                                       names(which(alignment_RefSeq_dist[x,] == min(alignment_RefSeq_dist[x,])))
                                                       }),
                             
                             Distance = unlist(lapply(1:nrow(alignment_RefSeq_dist), 
                                                      function(x){
                                                        min(alignment_RefSeq_dist[x,])
                                                      })))

#replacing the RefSeq samples closest RefSeq sample with NA's
positions <- unlist(lapply(RefSeq_sequences, 
                           function(x){
                             grep(pattern = x, 
                                  closest_RefSeq$accession)
                             }))
for (x in positions) {
  closest_RefSeq$RefSeq_samples[x] <- NA
  closest_RefSeq$Distance[x] <- NA
}

#saving the table 
fwrite(closest_RefSeq, 
       'closest_RefSeq.tsv',
       sep = '\t',
       quote = F,
       row.names = F)


closest_RefSeq$RefSeq_samples <- as.character(closest_RefSeq$RefSeq_samples)
closest_RefSeq$RefSeq_samples <- as.factor(closest_RefSeq$RefSeq_samples)

#creating a summary of the distances of the different RefSeq sequences 
RefSeq_summary <- closest_RefSeq[, .(min_distance = min(Distance), 
                                     max_distance = max(Distance)), 
                                 by = RefSeq_samples][ !c(2,42), ,]
