################################################################################
#
# analysing the results of the classifier assessment 
# rotation 2
# Cameron Ferguson, 03-05-2023
#
################################################################################

############################ setting up the session ############################

#importing the libraries 
library(equatiomatic)
library(data.table)
library(ggplot2)
library(emmeans)
library(ggpubr)
library(MASS)

#setting the working directory
setwd('~/Documents/LIDo rotation 2/NCBI data/sequnce_data/classifier_assessment_results/')

#importing the dataset
classifier_results <- fread('classification_results.csv', header = T, sep = ',')


######################### data cleaning and preparation ########################

#checking data types in the table
str(classifier_results)

#converting the columns to the correct data types
classifier_results$accession <- as.factor(classifier_results$accession)
classifier_results$classifier <- as.factor(classifier_results$classifier)
classifier_results$genus <- as.factor(classifier_results$genus)
classifier_results$subgenus <- as.factor(classifier_results$subgenus)
classifier_results$level_of_difference <- as.factor(classifier_results$level_of_difference)
classifier_results$sym_set <- as.factor(classifier_results$sym_set)
classifier_results$homopolymers <- as.logical(classifier_results$homopolymers)
classifier_results$correct_genus <- as.logical(classifier_results$correct_genus)
classifier_results$correct_species <- as.logical(classifier_results$correct_species)
classifier_results$correct_species_num <- ifelse(classifier_results$correct_species == TRUE,
                                                 yes = 1,
                                                 no = 0)

#converting the time from s to ms so we are dealing with whole numbers 
classifier_results$`time_taken(ms)` <- classifier_results$`time_taken(s)`*1000

#creating a new dataset and averaging across the sym_sets
level <- c(0 ,0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
homopol <- c(T, F)
sym_levels <- c(1, 2, 3, 4, 5)

#creating an empty table for the data 
new_data <- data.table(sym_set = NA, 
                       homopolymer_errors = NA, 
                       level_of_difference = NA, 
                       classifier = NA, 
                       `mean_time(ms)` = NA, 
                       proportion_identified_species = NA)[.0]

for (z in sym_levels) {
  for (y in homopol) {
    for (x in level) {
      temp2 <- data.table( rep(as.factor(z), 3),
                           
                           rep(as.logical(y), 3),
                           
                           rep(as.factor(x), 3),
                           
                           classifier_results[homopolymers == y & sym_set == z & level_of_difference == x,
                                              c(round(mean(`time_taken(ms)`), digits = 0)), 
                                              by = classifier],
                           
                           classifier_results[homopolymers == y & sym_set == z & level_of_difference == x,
                                              c((sum(correct_species_num)/24)), 
                                              by = classifier][, 2,]) 
      
      names(temp2) <- c('sym_set', 'homopolymer_errors', 'level_of_difference', 'classifier', 
                        'mean_time(ms)', 'proportion_identified_species')
      
      new_data <- rbind(new_data, temp2)
    }
  }
}


############################# assessing the dataset ############################

#exploring the dataset 
ggplot(new_data, aes(`mean_time(ms)`)) + geom_density() + facet_wrap(~homopolymer_errors)

ggplot(new_data[classifier == 'kaiju', ,], aes(`mean_time(ms)`)) + geom_density() + 
  facet_wrap(~~homopolymer_errors)

ggplot(new_data[classifier == 'kraken2', ,], aes(`mean_time(ms)`)) + geom_density() + 
  facet_wrap(~~homopolymer_errors)

ggplot(new_data[classifier == 'centrifuge', ,], aes(`mean_time(ms)`)) + geom_density() + 
  facet_wrap(~~homopolymer_errors)

#looking at the datasets it appears time is following a poisson distribution and given that 
#whether or not the correct prediction was made or not is being measured by a proportion which 
#is a bounded form of count data it will be modelled by a binomial/logistic distribution 

#worth noting that the visual inspection of the distributions indicate that the presence of 
#homopolymer errors shifts the distribution to the right 


############################# analysing the dataset ############################

###### modelling the time taken

#creating a saturated model
satmod <- with(new_data, glm(`mean_time(ms)` ~ classifier + level_of_difference + homopolymer_errors + proportion_identified_species, 
                             family = poisson()))
summary(satmod)

#checking the model for overdispersion by calculating the Residual deviance/degrees of freedom
1622.1/199
#the dispersion parameter is much greater then 1 so we should switch the distribution to a 
#negative binomial

satmodNB <- with(new_data, glm.nb(`mean_time(ms)` ~ classifier + level_of_difference + homopolymer_errors + proportion_identified_species, 
                                  link = 'log'))
summary(satmodNB)

#model simplification
mod1 <- with(new_data, glm.nb(`mean_time(ms)` ~ classifier + homopolymer_errors + proportion_identified_species, 
                              link = 'log'))
AIC(satmodNB, mod1)
BIC(satmodNB, mod1)
anova(satmodNB, mod1)

summary(mod1)

mod2 <- with(new_data, glm.nb(`mean_time(ms)` ~ classifier + homopolymer_errors, 
                              link = 'log'))

AIC(mod1, mod2)
BIC(mod1, mod2)
anova(mod1, mod2)

summary(mod2)

mod3 <- with(new_data, glm.nb(`mean_time(ms)` ~ classifier, 
                              link = 'log'))

AIC(mod2, mod3) #indicates model2 is better
BIC(mod2, mod3) #indicates model2 is better
anova(mod2, mod3) 

#this suggests that model2 is the minimally adequate model

#using emmeans to look at whats driving the classifiers significance
emmeans(mod2, pairwise~classifier, type = 'response')
emmeans(mod2, pairwise~homopolymer_errors, type = 'response')


###### modelling whether the correct outcome was predicted 

satmodB <- with(new_data, glm(proportion_identified_species ~ classifier + level_of_difference + homopolymer_errors + `mean_time(ms)`, 
                              family = binomial(link = 'logit')))

summary(satmodB)

#checking for over dispersion
3.3971/199
#model appears to be massively under dispersed

satmodQB <- with(new_data, glm(proportion_identified_species ~ classifier * level_of_difference + homopolymer_errors + `mean_time(ms)`, 
                               family = quasibinomial(link = 'logit')))
anova(satmodQB, test = 'F')

#model simplification

modQB <- with(new_data, glm(proportion_identified_species ~ classifier + level_of_difference + homopolymer_errors + `mean_time(ms)`, 
                            family = quasibinomial(link = 'logit')))
anova(satmodQB, modQB, test = 'F') #can't remove interaction term


modQB1 <- with(new_data, glm(proportion_identified_species ~ classifier * level_of_difference + homopolymer_errors, 
                             family = quasibinomial(link = 'logit')))

anova(satmodQB, modQB1, test = 'F')
anova(modQB1, test = 'F')

modQB2 <- with(new_data, glm(proportion_identified_species ~ classifier * level_of_difference, 
                             family = quasibinomial(link = 'logit')))

anova(modQB1, modQB2, test = 'F')
anova(modQB2, test = 'F') #minimal model 

#converting model to LaTex 
use_coefs = TRUE
extract_eq(modQB2, wrap = T)

############################### plotting the data ##############################

#without homopolymer errors 

time_plot <- ggplot(new_data[homopolymer_errors == F, ,], 
                    aes(level_of_difference, `mean_time(ms)`)) +
  geom_jitter(aes(col = classifier)) +
  
  geom_point(data = new_data[homopolymer_errors == F, 
                             c(median(`mean_time(ms)`)), 
                             by = c('classifier', 'level_of_difference')],
             
             aes(level_of_difference, V1)) +
  stat_boxplot(data = new_data[homopolymer_errors == F & classifier == 'kraken2', ,], 
               geom ='errorbar', 
               width = 0.25) +
  stat_boxplot(data = new_data[homopolymer_errors == F & classifier == 'kaiju', ,], 
               geom ='errorbar', 
               width = 0.25) +
  stat_boxplot(data = new_data[homopolymer_errors == F & classifier == 'centrifuge', ,], 
               geom ='errorbar', 
               width = 0.25) +
  scale_y_log10() + 
  ylab('log 10 of mean time taken (ms)') +
  xlab('Proportion of distance from RefSeq sequence') +
  theme_classic()

prediction_plot <- ggplot(new_data[homopolymer_errors == F, ,], 
                          aes(as.numeric(as.character(level_of_difference)), 
                              proportion_identified_species)) +
  geom_jitter(aes(col = classifier)) +
  geom_smooth(aes(group = classifier, col = classifier),
              method = glm, 
              formula = y~x, 
              method.args = list(family = quasibinomial(link = 'logit')), 
              se = F) +
  scale_x_continuous(labels = as.factor(seq(0, 0.3, 0.05)), 
                     breaks = seq(0, 0.3, 0.05)) +
  ylab('Proportion of species identified correctly') +
  xlab('Proportion of distance from RefSeq sequence') +
  theme_classic()



#combining the plots into 1 large plot 
ggarrange(time_plot, prediction_plot, ncol = 2, nrow = 1, labels = c('A', 'B'))


#with homopolymer errors 

time_plot2 <- ggplot(new_data[homopolymer_errors == T, ,], 
                     aes(level_of_difference, `mean_time(ms)`)) +
  geom_jitter(aes(col = classifier)) +
  
  geom_point(data = new_data[homopolymer_errors == T, 
                             c(median(`mean_time(ms)`)), 
                             by = c('classifier', 'level_of_difference')],
             
             aes(level_of_difference, V1)) +
  stat_boxplot(data = new_data[homopolymer_errors == T & classifier == 'kraken2', ,], 
               geom ='errorbar', 
               width = 0.25) +
  stat_boxplot(data = new_data[homopolymer_errors == T & classifier == 'kaiju', ,], 
               geom ='errorbar', 
               width = 0.25) +
  stat_boxplot(data = new_data[homopolymer_errors == T & classifier == 'centrifuge', ,], 
               geom ='errorbar', 
               width = 0.25) +
  scale_y_log10() + 
  ylab('log 10 of mean time taken (ms)') +
  xlab('Proportion of distance from RefSeq sequence') +
  theme_classic()

prediction_plot2 <- ggplot(new_data[homopolymer_errors == T, ,], 
                           aes(as.numeric(as.character(level_of_difference)), 
                              proportion_identified_species)) +
  geom_jitter(aes(col = classifier)) +
  geom_smooth(aes(group = classifier, col = classifier),
              method = glm, 
              formula = y~x, 
              method.args = list(family = quasibinomial(link = 'logit')), 
              se = F) +
  scale_x_continuous(labels = as.factor(seq(0, 0.3, 0.05)), 
                     breaks = seq(0, 0.3, 0.05)) +
  ylab('Proportion of species identified correctly') +
  xlab('Proportion of distance from RefSeq sequence') +
  theme_classic()



#combining the plots into 1 large plot 
ggarrange(time_plot2, prediction_plot2, ncol = 2, nrow = 1, labels = c('A', 'B'))

#creating a combined figure
ggarrange(time_plot, prediction_plot, time_plot2, prediction_plot2, ncol = 2, nrow = 2, labels = c('A', 'B', 'C', 'D'))
