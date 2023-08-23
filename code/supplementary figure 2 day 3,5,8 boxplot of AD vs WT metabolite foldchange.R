######################### day 3 ##########################################################
library(tidyverse)

#METABOLITES

#Load Counts
metabolite_counts <- read.csv("./input/metabolomics/metabolite_counts.csv")
rownames(metabolite_counts) <- metabolite_counts$ï..PARENT_SAMPLE_NAME
metabolite_counts$ï..PARENT_SAMPLE_NAME <- NULL

metabolite_features <- read.csv("./input/metabolomics/metabolite_features.csv")

colnames(metabolite_counts) <- metabolite_features$CHEMICAL_NAME


#Load metadata
metabolite_metadata <- read.csv("./input/metabolomics/metabolite_metadata.csv")


#metabolite_degs is the strain comparisons, so AD as numerator and N2 as denominator to 
#see what is diff from N2 and has all the metabolites and their fold change (FC) and q and p value
#Load Comparisons
metabolite_degs_day3 <- read.csv("./input/angelina_code/metabolite_degs_day3.csv")


#FOLD CHANGE COMPARISONS----------------------------------
sig_metabolites_day3 <- filter(metabolite_degs_day3, denominator == "N2")
#sig_metabolites_day3 <- filter(metabolite_degs_day3, q.value < 0.05)
sig_metabolites_day3 <- sig_metabolites_day3
#sig_metabolites_data <- filter(metabolite_degs, metabolite %in% sig_metabolites) 


ggplot(sig_metabolites_day3, aes(x = numerator, y = abs(log(FC, 2)))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.1, color = "black", stroke = 1) +
  xlab("") + ylab("Absolute log2 fold change") +
  theme_light()
#add shape = 1 after stroke to make it the dot black outline with no fill

#geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "black") #geom points with jitter
#geom_boxplot(outlier.shape = NA) # remove outlier from boxplot




#################################################################################
########################### day 5 ###############################################


#metabolite_degs is the strain comparisons, so AD as numerator and N2 as denominator to 
#see what is diff from N2 and has all the metabolites and their fold change (FC) and q and p value
#Load Comparisons
metabolite_degs_day5 <- read.csv("./input/angelina_code/metabolite_degs_day5.csv")


#FOLD CHANGE COMPARISONS----------------------------------
sig_metabolites_day5 <- filter(metabolite_degs_day5, denominator == "N2")
#sig_metabolites_day5 <- filter(metabolite_degs_day5, q.value < 0.05)
sig_metabolites_day5 <- sig_metabolites_day5
#sig_metabolites_data <- filter(metabolite_degs, metabolite %in% sig_metabolites) 


ggplot(sig_metabolites_day5, aes(x = numerator, y = abs(log(FC, 2)))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.1, color = "black", stroke = 1) +
  xlab("") + ylab("Absolute log2 fold change") +
  theme_light()
#add shape = 1 after stroke to make it the dot black outline with no fill

#geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "black") #geom points with jitter
#geom_boxplot(outlier.shape = NA) # remove outlier from boxplot


#################################################################################
########################### day 8 ###############################################


#metabolite_degs is the strain comparisons, so AD as numerator and N2 as denominator to 
#see what is diff from N2 and has all the metabolites and their fold change (FC) and q and p value
#Load Comparisons
metabolite_degs_day8 <- read.csv("./input/angelina_code/metabolite_degs_day8.csv")


#FOLD CHANGE COMPARISONS----------------------------------
sig_metabolites_day8 <- filter(metabolite_degs_day8, denominator == "N2")
#sig_metabolites_day8 <- filter(metabolite_degs_day8, q.value < 0.05)
sig_metabolites_day8 <- sig_metabolites_day8
#sig_metabolites_data <- filter(metabolite_degs, metabolite %in% sig_metabolites) 


ggplot(sig_metabolites_day8, aes(x = numerator, y = abs(log(FC, 2)))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.1, color = "black", stroke = 1) +
  xlab("") + ylab("Absolute log2 fold change") +
  theme_light()
#add shape = 1 after stroke to make it the dot black outline with no fill

#geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "black") #geom points with jitter
#geom_boxplot(outlier.shape = NA) # remove outlier from boxplot

