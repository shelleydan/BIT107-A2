############ DIFFERENTIAL EXPRESSED GENES ############ 

### LIBRARY ###
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

## DIRECTORY SETTING ##
setwd="~/BIT107-A2-LPT-SAVE/BIT107-A2/DEGdata/" #Site for all DEG data produced to be inputted to (Change LPT for PC use)

## IMPORTING DATA ##
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

## FACTOR LEVELS ## - Infected vs Mock

target_data$Treatment <- factor(target_data$Treatment) #Setting Treatments as a factor argument
target_data$Sequencing <- factor(target_data$Sequencing) #Setting Sequencing method as a factor argument




