
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
meta$Test[meta$Test == 'infected'] <- 'I'
meta$Test[meta$Test == 'mock'] <- 'M'
meta$Replicate <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3) #Manually adding replicate numbers
str(meta)

data <- data[, colnames(data) %in% row.names(meta)] #Removing bad data from Counts

#Moanin Model doesn't like column names not in the same order as the targets, here I reorder the column names alpha-numerically
data_order = sort(colnames(data)) #Re-order column names to alpha-numerical
data <- data[, data_order] #Apply the reordering

#Checking Data
summary(data)

################################################################################
################ COLOUR CODING FOR WHOLE SCRIPT STANDARDISATION ################
################################################################################

#Defining Colours
test_colours <- c("M" = "Blue",
                  "I" = "Red")

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

meta$Test <- factor(meta$Test, levels(meta$Test)[c(2,1)])
#2 is Mock, 1 is infected

################################################################################
############################## TEMPORAL MODELLING ##############################
################################################################################
moaninModel <- create_moanin_model(data=data, 
                                    meta=meta,
                                    group_variable = "Test",
                                    time_variable = "Condition")

show(moaninModel) # Gives a summary of the produced model

#Defining groups to contrasts for analysis
contrasts <- create_timepoints_contrasts(moaninModel, "M", "I") 
#These contrasts output all the groups for comparison for all the timepoints, e.g. mock.4hrs_vs_infected.4hrs

hour_de_analysis <- DE_timepoints(moaninModel, contrasts, use_voom_weights = T) #Set vooms to F for microarray data
#This provides a comparative table of pvals, logfoldchange(lfc) and adjusted p (qval)


#PLOT: Histogram to highlight DEGs at different timepoints
time = names(condition_colours)
perWeek_barplot(hour_de_analysis, labels=time, main="Mock vs Infected", las=3) #More genes are differentially expressed at 12 hours

#Code for counting the number of unique combinations of time points that are DE
getTimepoint<-function(x){sapply(strsplit(gsub("_qval","",x),"\\."),.subset2,3)}
qval_colnames = colnames(hour_de_analysis)[
  grepl("qval", colnames(hour_de_analysis))]
signifCombos<-apply(hour_de_analysis[, qval_colnames], 1, 
                    function(x){paste(getTimepoint(qval_colnames[which(x<0.05)]),collapse=",")})
signifCombos<-signifCombos[signifCombos!=""]
tabCombos<-table(signifCombos)


#PLOT!: Gene Expression of top10 overtime. 
head(hour_de_analysis)

#top10DEGs are the top 10 DEGs in the 4hr analysis determined on DEG_Analysis.R

exampleGenes<-names(signifCombos[signifCombos=="4,12,48"][row.names(top10DEGs)])
exampleGenes <- row.names(top10DEGs)
plot_splines_data(moaninModel, 
                  subset_data=exampleGenes,
                  colors=ann_colours$Test,
                  smooth=TRUE,
                  ylim = c(-10000, 10000))

head(hour_de_analysis)

################################################################################
#The paper continues to do analysis between different conditions, as this data #
#only contains infected or mock, this isn't possible.                          #
################################################################################


################################################################################
########################### KEGG ANALYSIS OVER-TIME ############################
################################################################################

#1. Filtering genes
# Then rank by fisher's p-value and take max the number of genes of interest
# Filter out q-values for the pvalues table
fishers_pval = pvalues_fisher_method(pvalues)
qvalues = apply(pvalues, 2, p.adjust)
fishers_qval = p.adjust(fishers_pval)

genes_to_keep = row.names(
  log_fold_change_max[
    (rowSums(log_fold_change_max > 2) > 0) &
      (fishers_qval < 0.05), ])
# Keep the data corresponding to the genes of interest in another variable.
# by subsetting the `moanin_model`, which contains the data.
de_moanin_model = moanin_model[genes_to_keep,]






