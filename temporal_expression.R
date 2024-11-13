
################################################################################
########################## TEMPORTAL RNA-seq ANALYSIS ##########################
################################################################################

## REFERENCES ##
#Varoquaux, N. and Purdom, E. 2020. A pipeline to analyse time-course gene expression data. 
#Available at: https://f1000research.com/articles/9-1447 [Accessed: 12 November 2024].

#KEGG Analysis
# https://www.youtube.com/watch?v=SMBF4DyRiuo


################################################################################
############################# LIBRARY INSTALL/LOAD #############################
################################################################################

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

################################################################################
############################# DATA INPUT & CLEANING ############################
################################################################################

#NOTE!: Start with a Fresh Environment

data <- read.csv("rawdata/GSE217504_host_counts_matrix.csv", #Input Counts Data
                        header = T, 
                        row.names = 1)

meta <- read.delim("rawdata/targets.txt", 
                          sep = "", 
                          header = T) #Inputting the targets
meta <- meta[,8:10] #Seperating out the only columns we want
meta$Test <- as.factor(meta$Test)

# Remove the data which have no controls
vals <- c(4, 12, 48) # Taking only the timepoints that have controls and tests
meta <- meta[meta$Condition %in% vals, ] #Removing Bad Data from Targets
meta <- data.frame(meta, row.names = 3) #Setting sample names as the rownames
meta$Replicate <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3) #Manually adding replicate numbers
data <- data[, colnames(data) %in% row.names(meta)] #Removing bad data from Counts

#DESeq2 doesn't like column names not in the same order as the targets, here I reorder the column names alpha-numerically
data_order = sort(colnames(data)) #Re-order column names to alpha-numerical
data <- data[, data_order] #Apply the reordering

#Checking Data
summary(data)
summary(testData)
summary(counts_data)

################################################################################
################ COLOUR CODING FOR WHOLE SCRIPT STANDARDISATION ################
################################################################################

#Defining Colours
test_colours <- c("mock" = "Green",
                  "infected" = "Red")

condition_colours <- c("4" = "Blue",
                       "12" = "Purple",
                       "48" = "Black")

rep_marker <- c(15, 17, 19)
names(rep_marker) = c(1,2,3)

#Assigning Colours
ann_colours <- list(
  Condition=condition_colours,
  Test=test_colours
)

ann_markers <- list(
  Replicate = rep_marker
)

################################################################################
############################## TEMPORAL MODELLING ##############################
################################################################################

moaninModel <- create_moanin_model(data=data, 
                                    meta=meta,
                                    group_variable = "Test",
                                    time_variable = "Condition")

show(moaninModel)




