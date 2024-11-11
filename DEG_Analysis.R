############ DIFFERENTIAL EXPRESSED GENES ############ 

### LIBRARY ###
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") #DESeq2 is now built under Bioconductor and needs to be installed differently. 

BiocManager::install("DESeq2") #Downloads DESeq2 from Bioconductor
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

## DIRECTORY SETTING ##
setwd="~/BIT107-A2-LPT-SAVE/BIT107-A2/DEGdata/" #Site for all DEG data produced to be inputted to (Change LPT for PC use)
