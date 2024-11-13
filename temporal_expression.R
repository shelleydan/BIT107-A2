
####################################################
############ TEMPORTAL RNA-seq ANALYSIS ############
####################################################

## REFERENCES ##
#Varoquaux, N. and Purdom, E. 2020. A pipeline to analyse time-course gene expression data. 
#Available at: https://f1000research.com/articles/9-1447 [Accessed: 12 November 2024].

################
## DATA INPUT ##
################

targets <- read.delim("rawdata/targets.txt", sep = "", header = T) #Inputting the targets
targets <- targets[,8:10] #Seperating out the only columns we want

counts <- read.csv("rawdata/GSE217504_host_counts_matrix.csv",
                   header = T,
                   row.names = 1)




















