######################################################
############ DIFFERENTIAL EXPRESSED GENES ############ 
######################################################

#############
## LIBRARY ##
#############

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") #DESeq2 is now built under Bioconductor and needs to be installed differently. 

BiocManager::install("DESeq2") #Downloads DESeq2 from Bioconductor
BiocManager::install("Rsamtools") #Downloads SAMTools from Bioconductor
library(pheatmap)
library(dplyr)
library(Rsamtools)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

#######################
## DIRECTORY SETTING ##
#######################

setwd="~/BIT107-A2-LPT-SAVE/BIT107-A2/DEGdata/" #Site for all DEG data produced to be inputted to (Change LPT for PC use)

####################
## IMPORTING DATA ##
####################

counts_data <- read.csv("DEGdata/count_matrix.csv", #Counts
                        header = T,
                        row.names = 1)

colnames(counts_data) #Checking samples are column names
head(counts_data) #Checking the data imported correctly

target_data <- read.csv("DEGdata/design.csv", #Target meta data
                        header = T, 
                        row.names = 1)

colnames(target_data) #Checking samples are column names
head(target_data) #Checking the data imported correctly

###################
## FACTOR LEVELS ## - Infected vs Mock
###################

target_data$Treatment <- factor(target_data$Treatment) #Setting Treatments as a factor argument
target_data$Sequencing <- factor(target_data$Sequencing) #Setting Sequencing method as a factor argument

#############################
## CREATING A DESeq OBJECT ##
#############################

DEG <- DESeq2::DESeqDataSetFromMatrix(countData = counts_data, #Adding Counts Data
                                      colData = target_data, #Adding targets data
                                      design = ~Sequencing + Treatment) #Factors for Comparison
#If we're looking at a multi-factor analysis, we want to input our primary factor should be inputted last.

#######################
## SETTING REFERENCE ##
#######################

DEG$Treatment <- factor(DEG$Treatment, levels = c("untreated", "treated")) 
#The first level in the list is what becomes the base level, which the other levels will be compared to. We want this to be the control/mock group

## GENE FILTERING ## - double check if this is done with Sarah's
keep <- rowSums(counts(DEG)) >= 1 #Give an expression with readcounts more than 1 will be stored here. 
#This number should be justified in writing!
DEG <- DEG[keep,] #Selecting genes with more than 1 read. 

##############
## ANALYSIS ##
##############

DEG <- DESeq2::DESeq(DEG) #Perform the Analysis
results_DEG <- DESeq2::results(DEG) #Coalating the results into a dataframe
results_DEG #Printing the results into the console

#Need to change the DESeq object into a dataframe
results_DEG <- as.data.frame(results_DEG) # Produces and R dataframe

#Changing table order by acending padj
results_DEG_acending <- results_DEG[order(results_DEG$padj),]
head(results_DEG_acending)


