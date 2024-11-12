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
BiocManager::install("apeglm")

library(pheatmap)
library(dplyr)
library(tidyverse)
library(Rsamtools)
library(apeglm)
library(DESeq2)
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

#Changing table order by ascending padj
results_DEG_acending <- results_DEG[order(results_DEG$padj),]
head(results_DEG_acending)

###########################
## FILTERING THE RESULTS ##
###########################

#Setting a column for the Volcano plot
results_DEG$diffexpressed <- "NO"
results_DEG$diffexpressed[results_DEG$log2FoldChange > 1 & results_DEG$padj < 0.05] <- "UP"
results_DEG$diffexpressed[results_DEG$log2FoldChange < -1 & results_DEG$padj < 0.05] <- "DOWN"

########################################
## SAVING THE RAW-COUNTS AND FILTERED ##
########################################

#COPY FROM POSIT

write.csv(results_DEG, "output_data/raw_DEG_results.csv") #Output of the non-filtered DESeq2 results

############################
###### PRETTY FIGURES ######
############################

#####################
## DISPERSION PLOT ##
#####################

#We expect that when a gene's read count increases the dispersion of that same gene decreases
plotDispEsts(DEG)

##################################
## PRINCIPLE COMPONENT ANALYSIS ##
##################################

# Variance stabilisation transformation
vst <- DESeq2::vst(DEG, blind = F)

# Generating the PCA Plot
DESeq2::plotPCA(vst, 
                intgroup=c("Sequencing", "Treatment")) #Applying layers like found above.

#LOOK AT PCA WITH GGPLOT

##############
## HEATMAPS ##
##############

#This used pheatmaps to analyse some gene expression clusting.

# Sample-to-sample distance matrix (normalised counts)

sampleDist <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDist) # Generate a matrix
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255) #Setting the colours and range of those. 

#Generating the Heatmap
pheatmap(sampleDistMatrix, #Matrix input
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         color = colours) #as defined above

#This heatmap colours show the different between the samples. So the darkest blue shows no difference (e.g. when the same samples are plotted against eachother)

# Log transformed normalised counts (using top 10 genes)

DEG_padj_orders <- results_DEG[order(results_DEG$padj),] #Ascending order of padj values

top10DEGs <- results_DEG[order(results_DEG$padj), ][1:10,] #Selecting the top 10
top10DEGs_names <- row.names(top10DEGs) #Extracting names of the top 10 genes

rld <- rlog(DEG, blind = F) 
pheatmap(assay(rld)[top10DEGs_names,], #Subset by labels extracted
         cluster_rows = T, #Adds column tree-clustering
         show_rownames = T, 
         cluster_cols = T, #Adds rows tree-clustering
         annotation_col = anno_info) 

anno_info <- as.data.frame(colData(DEG)[,c("Sequencing", "Treatment")])



# Z-Scores with top 10 genes

z_calc <- function(x) { #Function for calculating z-scores
  ((x - mean(x)) / sd(x))
}

normal_counts <- counts(DEG, normalized = T)

zscore_all <- t(apply(normal_counts, 1, z_calc)) #Z-scores for all genes

zscore_sub <- zscore_all[top10DEGs_names,] #Subsetting for the top 10 genes

pheatmap(zscore_sub)

#############
## MA PLOT ## - Gene Expression vs log2FoldChange
#############

#These plots are used to see the distribution of gene expressions

plotMA(DEG, ylim=c(-2,2))

#Removing Noise

resLFC <- lfcShrink(DEG, 
                    coef = "Treatment_treated_vs_untreated", 
                    type = "apeglm") #This gives a ref

plotMA(resLFC, ylim=c(-2,2))


######################################
## VOLCANO PLOT OF DIFFERENTIALTION ##
######################################

#Setting the theme for standardization
theme_set(theme_classic(base_size = 15) +
            theme(
              axis.title.y = element_text(face = 'bold', margin = margin(0,20,0,0), size = rel(1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

#Volcano Plot with ggplot2
ggplot(data = results_DEG, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = deflabel)) +
  geom_vline(xintercept = c(1, -1), col = "gray", linetype = "dashed") + #Adding vertical lines to show fold change cut off
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed") + #Adding Horizontal line to show p-value cut off
               geom_point() + #Setting Point size
               scale_color_manual(values = c("#00AFBB", "gray", "#bb0c00"), #Changing plot colours
                                  labels = c("Downregulated","Not Significant", "Upregulated")) + #Changing label titles
               coord_cartesian(ylim = c(0, 250), xlim = c(-5,5)) + #Applying figure axis limits
               scale_x_continuous(breaks = seq(-10, 10, 2)) + # Setting continuous breaks (min, max, step)
               labs(color = "Regulation", # Changing the colour legend title
                    x = expression("log"[2]*" Fold Change"), #Changing the x axis
                    y = expression("-log"[10]*" p-adjusted")) + #Changing the y axis
               ggtitle("DEGs")  # Setting a Figure Title
               geom_text_repel(max.overlaps = Inf) #Adding labels which we express in ggplot line 1
             










