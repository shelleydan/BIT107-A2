######################################################
############ DIFFERENTIAL EXPRESSED GENES ############ 
######################################################

## SCRIPT REFERENCES ##
# Batut et al - https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html
# Cristofides - https://github.com/ecologysarah?tab=repositories

#############
## LIBRARY ##
#############

install.packages("BiocManager") #Some Packages are now built under Bioconductor and needs to be installed differently. 

#Bioconductor Package Installs
BiocManager::install("DESeq2") #Downloads DESeq2 from Bioconductor
BiocManager::install("Rsamtools") #Downloads SAMTools from Bioconductor
BiocManager::install("apeglm")

#Library Loading
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

setwd="~/BIT107-A2-LPT-SAVE/BIT107-A2/DEGdata/" 
#Site for all DEG data (inputs and outputs)

####################
## IMPORTING DATA ##
####################

counts_data <- read.csv("rawdata/GSE217504_host_counts_matrix.csv", #Input Counts Data
                        header = T, 
                        row.names = 1)

colnames(counts_data) #Checking samples are column names
head(counts_data) #Checking the data imported correctly

target_data <- read.delim("rawdata/targets.txt", sep = "", header = T) #Inputting the targets
target_data <- target_data[,8:10] #Seperating out the only columns we want
target_data$Condition <- as.factor(target_data$Condition)

colnames(target_data) #Checking samples are column names
head(target_data) #Checking the data imported correctly

# Remove the data which have no controls
vals <- c(4,12,48) #List of values to subset from
target_data <- target_data[target_data$Condition %in% vals, ] #Removing Bad Data from Targets
counts_data <- counts_data[, colnames(counts_data) %in% target_data$Sample] #Removing bad data from Counts


###################
## FACTOR LEVELS ## - Infected vs Mock
###################

target_data$Test <- factor(target_data$Test) #Setting Treatments as a factor argument
target_data$Condition <- factor(target_data$Condition) #Setting Sequencing method as a factor argument

#############################
## CREATING A DESeq OBJECT ##
#############################

DEG <- DESeq2::DESeqDataSetFromMatrix(countData = counts_data, 
                                      #Adding Counts Data, 2:8 removes column 1 with row names in. 
                                      colData = target_data, #Adding targets data
                                      design = ~Condition + Test) #Factors for Comparison
#If we're looking at a multi-factor analysis, we want to input our primary factor should be inputted last.

#######################
## SETTING REFERENCE ##
#######################

DEG$Test <- factor(DEG$Test, levels = c("mock", "infected")) 
#The first level in the list is what becomes the base level, which the other levels will be compared to. We want this to be the control/mock group

## GENE FILTERING ## - double check if this is done with Sarah's
#keep <- rowSums(counts(DEG)) >= 1 #Give an expression with readcounts more than 1 will be stored here. 
#This number should be justified in writing!
#DEG <- DEG[keep,] #Selecting genes with more than 1 read. 

##############
## ANALYSIS ##
##############

DEG <- DESeq2::DESeq(DEG) #Perform the Analysis
results_DEG <- DESeq2::results(DEG) #Coalating the results into a dataframe
results_DEG #Printing the results into the console

#Need to change the DESeq object into a dataframe
results_DEG <- as.data.frame(results_DEG) # Produces and R dataframe

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
plotDispEsts(DEG, 
             ylab = "Dispersion",
             xlab = "Mean of Normalised Counts") #This is run on the Analysis Data

##################################
## PRINCIPLE COMPONENT ANALYSIS ##
##################################

# Variance stabilisation transformation
vst <- DESeq2::vst(DEG, blind = F)

# Generating the PCA Plot
DESeq2::plotPCA(vst, 
                intgroup= "Test") #Applying layers like found above.

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

topDEGs <- results_DEG[order(results_DEG$padj), ][1:10,] #Selecting the top 10
topDEGs_names <- row.names(topDEGs) #Extracting names of the top 10 genes

rld <- rlog(DEG, blind = F) #Performling a log transformation

anno_info <- as.data.frame(colData(DEG)[, c("Condition", "Test")]) #Setting Annotation Levels
anno_info$Condition <- as.character(anno_info$Condition) #Changing Condition from continuous to ordinal

pheatmap(assay(rld)[topDEGs_names,], #Subset by labels extracted
         cluster_rows = T, #Adds column tree-clustering
         show_rownames = T, 
         cluster_cols = T, #Adds rows tree-clustering
         annotation_col = anno_info) 

# Z-Scores with top 10 genes

z_calc <- function(x) { #Function for calculating z-scores
  ((x - mean(x)) / sd(x))
}

top10DEGs <- results_DEG[order(results_DEG$padj), ][1:10,]
top10DEGs_names <- row.names(top10DEGs)

normal_counts <- counts(DEG, normalized = T) #Normalising the counts
zscore_all <- t(apply(normal_counts, 1, z_calc)) #Z-scores for all genes
zscore_sub <- zscore_all[top10DEGs_names,] #Subsetting for the top 10 genes

pheatmap(zscore_sub)

#############
## MA PLOT ## - Gene Expression vs log2FoldChange
#############

#These plots are used to see the distribution of gene expressions
#The default alpha for MA plots is 0.1

plotMA(DEG, ylim=c(-2,2)) #Setting the alpha value to 0.05

#Removing Noise

resLFC <- lfcShrink(DEG, 
                    coef = "Test_infected_vs_mock", 
                    type = "apeglm") #This gives a ref

plotMA(resLFC, ylim=c(-2,2), alpha = 0.05)

## MA Plot with ggplot 2
colnames(results_DEG)


###########################
## FILTERING THE RESULTS ## - For Volcano Plot
###########################

#Setting a column for the Volcano plot
results_DEG$diffexpressed <- "NO"
results_DEG$diffexpressed[results_DEG$log2FoldChange > 1 & results_DEG$pvalue < 0.05] <- "UP"
results_DEG$diffexpressed[results_DEG$log2FoldChange < -1 & results_DEG$pvalue < 0.05] <- "DOWN"

#####################################
## VOLCANO PLOT OF DIFFERENTIATION ##
#####################################

#Extracting GeneIDs from the row.names without needing for a new GeneID column
top10DEGs <- results_DEG[order(results_DEG$pvalue), ][1:10,]
results_DEG$difflabel <- ifelse(row.names(results_DEG) %in% row.names(top10DEGs), row.names(results_DEG), NA)
summary(results_DEG$difflabel)

#Setting the theme for standardization
theme_set(theme_classic(base_size = 15) +
            theme(
              axis.title.y = element_text(face = 'bold', 
                                          margin = margin(0,20,0,0), 
                                          size = rel(1), 
                                          color = 'black'),
              axis.title.x = element_text(hjust = 0.5, 
                                          face = "bold", 
                                          margin = margin(20,0,0,0), 
                                          size = rel(1), 
                                          color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

#Volcano Plot with ggplot2
ggplot(data = results_DEG, aes(x = log2FoldChange, 
                               y = -log10(pvalue), 
                               col = diffexpressed, 
                               label = difflabel)) +
  geom_vline(xintercept = c(1, -1), #Manually adding cutoff lines 
             col = "gray", 
             linetype = "dashed") + #Adding vertical lines to show fold change cut off
  geom_hline(yintercept = c(-log10(0.05)), 
             col = "gray", 
             linetype = "dashed") + #Adding Horizontal line to show p-value cut off
  geom_point() + #Setting Point size
  scale_shape_manual(values = 6) +
  scale_color_manual(values = c("#00AF99", "gray", "#bb0c00"), #Changing plot colours
                     labels = c("Downregulated","Not Significant", "Upregulated")) + #Changing label titles
  coord_cartesian(ylim = c(0, 5), xlim = c(-5,5)) + #Applying figure axis limits
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # Setting continuous breaks (min, max, step)
  labs(color = "Regulation", # Changing the colour legend title
       x = expression("log"[2]*" Fold Change"), #Changing the x axis
       y = expression("-log"[10]*" p-adjusted")) + #Changing the y axis
  ggtitle("DEGs") + # Setting a Figure Title
  geom_text_repel(max.overlaps = Inf) #Adding labels which we express in ggplot line 1
             








