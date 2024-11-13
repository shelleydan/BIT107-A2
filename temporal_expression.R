
####################################################
############ TEMPORTAL RNA-seq ANALYSIS ############
####################################################

## REFERENCES ##
#Varoquaux, N. and Purdom, E. 2020. A pipeline to analyse time-course gene expression data. 
#Available at: https://f1000research.com/articles/9-1447 [Accessed: 12 November 2024].

#KEGG Analysis
# https://www.youtube.com/watch?v=SMBF4DyRiuo


#############
## LIBRARY ##
#############

install.packages("BiocManager") #Some Packages are now built under Bioconductor and needs to be installed differently. 
library(BiocManager)

BiocManager::install("timecoursedata")
BiocManager::install("moanin")
BiocManager::install("topGO")
BiocManager::install("biomaRt")
BiocManager::install("KEGGprofile") #This was removed from Bioconductor - Another method will be used
BiocManager::install("BiocWorkflowTools")
install.packages("NMF")
install.packages("ggfortify")

library(timecoursedata)
library(moanin)
library(topGO)
library(biomaRt)
library(KEGGprofile)
library(BiocWorkflowTools)
library(NMF)
library(ggfortify)

###########################
## DATA INPUT & CLEANING ##
###########################

#NOTE!: Start with a Fresh Environment

counts_data <- read.csv("rawdata/GSE217504_host_counts_matrix.csv", #Input Counts Data
                        header = T, 
                        row.names = 1)

colnames(counts_data) #Checking samples are column names
head(counts_data) #Checking the data imported correctly

target_data <- read.delim("rawdata/targets.txt", sep = "", header = T) #Inputting the targets
target_data <- target_data[,8:10] #Seperating out the only columns we want
target_data$Condition <- as.factor(target_data$Condition)
target_data$Test <- as.factor(target_data$Test)

colnames(target_data) #Checking samples are column names
head(target_data) #Checking the data imported correctly

# Remove the data which have no controls
vals <- c(4, 12, 48) # Taking only the timepoints that have controls and tests
target_data <- target_data[target_data$Condition %in% vals, ] #Removing Bad Data from Targets
target_data <- data.frame(target_data, row.names = 3) #Setting sample names as the rownames
counts_data <- counts_data[, colnames(counts_data) %in% row.names(target_data)] #Removing bad data from Counts

#DESeq2 doesn't like column names not in the same order as the targets, here I reorder the column names alpha-numerically
counts_data_order = sort(colnames(counts_data)) #Re-order column names to alpha-numerical
counts_data<- counts_data[, counts_data_order] #Apply the reordering

#Checking Data
summary(target_data)
summary(counts_data)






















