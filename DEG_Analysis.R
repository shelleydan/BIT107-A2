################################################################################
############ DIFFERENTIAL EXPRESSED GENES ############ 
################################################################################

## SCRIPT REFERENCES ##
# Batut et al - https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html
# Cristofides - https://github.com/ecologysarah?tab=repositories

################################################################################
## LIBRARY ##
################################################################################

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
library(ggpubr)


######################## SETTING GGPLOT2 DOCUMENT THEME ########################


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
              plot.title = element_text(hjust = 0,
                                        size = 15),
              legend.title=element_blank()
            ))

################################################################################
## DIRECTORY SETTING ##
#######################

setwd="~/BIT107-A2-LPT-SAVE/BIT107-A2/DEGdata/" 
#Site for all DEG data (inputs and outputs)

####################
## IMPORTING DATA ##
################################################################################

#READ IN COUNTS DATA
counts_data <- read.csv("rawdata/GSE217504_host_counts_matrix.csv",
                        header = T, 
                        row.names = 1)

#READ IN AND EXTRACT TARGETS DATA
target_data <- read.delim("rawdata/targets.txt", 
                          sep = "", 
                          header = T) #Inputting the targets.
target_data <- target_data[,8:10] #Separating out the necessary columns.
target_data$Condition <- as.factor(target_data$Condition) #Setting 'Condition' (AKA Time) to factors.

#SUBSET THE TARGETS DATA (4,12,48 hrs which all have controls)
target_data <- data.frame(target_data, row.names = 3) #Setting sample names as the rownames

target_data_4 <- target_data[target_data$Condition %in% 4, ] #Targets for 4hrs
target_data_12 <- target_data[target_data$Condition %in% 12, ] #Targets for 12hrs
target_data_48 <- target_data[target_data$Condition %in% 48, ] #Tagrets for 48hrs

#ORDER THE COUNT DATA (DESeq2 requires).
counts_data_order = sort(colnames(counts_data)) #Re-order column names to alpha-numerical
counts_data<- counts_data[, counts_data_order]

#SUBSET THE COUNT DATA (4, 12, 48hrs).
counts_data_4 <- counts_data[, colnames(counts_data) %in% row.names(target_data_4)] #Counts for 4hrs
counts_data_12 <- counts_data[, colnames(counts_data) %in% row.names(target_data_12)] #Counts for 12hrs
counts_data_48 <- counts_data[, colnames(counts_data) %in% row.names(target_data_48)] #Counts for 48hrs

#TIDYING THE ENVIRONMENT
rm(counts_data)
rm(target_data)
rm(counts_data_order)

################################################################################
## FACTOR LEVELS ##
################################################################################

#FACTOR LEVELS FOR 4HRS
target_data_4$Test <- factor(target_data_4$Test) 
target_data_4$Condition <- factor(target_data_4$Condition) 

#FACTOR LEVELS FOR 12HRS
target_data_12$Test <- factor(target_data_12$Test) 
target_data_12$Condition <- factor(target_data_12$Condition) 

#FACTOR LEVELS FOR 48HRS
target_data_48$Test <- factor(target_data_48$Test) 
target_data_48$Condition <- factor(target_data_48$Condition) 

################################################################################
## CREATING A DESeq MATRIX ##
################################################################################

#If we're looking at a multi-factor analysis, we want to input our primary factor should be inputted last.

#The first level in the list is what becomes the base level, which the other levels will be compared to. We want this to be the control/mock group

#4HRS 
DEG_4 <- DESeqDataSetFromMatrix(countData = counts_data_4, 
                                colData = target_data_4, #Adding targets data
                                design = ~Test) #Factors for Comparison
DEG_4$Test <- factor(DEG_4$Test, levels = c("mock", "infected")) 

#12HRS
DEG_12 <- DESeqDataSetFromMatrix(countData = counts_data_12, 
                                 colData = target_data_12, #Adding targets data
                                 design = ~Test) #Factors for Comparison
DEG_12$Test <- factor(DEG_12$Test, levels = c("mock", "infected")) 

#48HRS
DEG_48 <- DESeqDataSetFromMatrix(countData = counts_data_48, 
                                 colData = target_data_48, #Adding targets data
                                 design = ~Test) #Factors for Comparison
DEG_48$Test <- factor(DEG_48$Test, levels = c("mock", "infected")) 

## GENE FILTERING ## - double check if this is done with Sarah's
#keep <- rowSums(counts(DEG)) >= 1 #Give an expression with readcounts more than 1 will be stored here. 
#This number should be justified in writing!
#DEG <- DEG[keep,] #Selecting genes with more than 1 read. 

################################################################################
## PERFORMING THE ANALYSIS ##
################################################################################

#4HRS
DEG_4 <- DESeq2::DESeq(DEG_4) #Perform the Analysis
results_DEG_4 <- DESeq2::results(DEG_4) #Coalating the results into a dataframe
summary(results_DEG_4$padj)
results_DEG_4 <- as.data.frame(results_DEG_4) # Produces and R dataframe

#12HRS
DEG_12 <- DESeq2::DESeq(DEG_12) #Perform the Analysis
results_DEG_12 <- DESeq2::results(DEG_12) #Coalating the results into a dataframe
summary(results_DEG_12$padj)
results_DEG_12 <- as.data.frame(results_DEG_12) # Produces and R dataframe

#48HRS
DEG_48 <- DESeq2::DESeq(DEG_48) #Perform the Analysis
results_DEG_48 <- DESeq2::results(DEG_48) #Coalating the results into a dataframe
summary(results_DEG_48$padj)
results_DEG_48 <- as.data.frame(results_DEG_48) # Produces and R dataframe

################################################################################
## SAVING THE RAW-COUNTS AND FILTERED ##
################################################################################

#Output of the non-filtered DESeq2 results
write.csv(results_DEG_4, "output_data/raw_DEG_results_4.csv") 
write.csv(results_DEG_12, "output_data/raw_DEG_results_12.csv")
write.csv(results_DEG_48, "output_data/raw_DEG_results_48.csv")

################################################################################
###### PRETTY FIGURES ######
################################################################################

################################################################################
## DISPERSION PLOT ##
################################################################################

#We expect that when a gene's read count increases the dispersion of that same gene decreases
par(mfrow=c(2,2))
#4HRS
plotDispEsts(DEG_4, 
             ylab = "Dispersion",
             xlab = "Mean of Normalised Counts") #This is run on the Analysis Data

#12HRS
plotDispEsts(DEG_12, 
             ylab = "Dispersion",
             xlab = "Mean of Normalised Counts")

#48HRS
plotDispEsts(DEG_48, 
             ylab = "Dispersion",
             xlab = "Mean of Normalised Counts")

par(mfrow=c(1,1))
################################################################################
## PRINCIPLE COMPONENT ANALYSIS ##
################################################################################

#4HRS
# Variance stabilisation transformation
vst_4 <- DESeq2::vst(DEG_4, blind = F)

# Generating the PCA Plot
DESeq2::plotPCA(vst_4, 
                intgroup= "Test") #Applying layers like found above.

#12HRS
# Variance stabilisation transformation
vst_12 <- DESeq2::vst(DEG_12, blind = F)

# Generating the PCA Plot
DESeq2::plotPCA(vst_12, 
                intgroup= "Test") #Applying layers like found above.

#48HRS
# Variance stabilisation transformation
vst_48 <- DESeq2::vst(DEG_48, blind = F)

# Generating the PCA Plot
DESeq2::plotPCA(vst_48, 
                intgroup= "Test") #Applying layers like found above.

################################################################################
## HEATMAPS ##
################################################################################

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

################################################################################
## MA PLOT ## - Gene Expression vs log2FoldChange
################################################################################

#These plots are used to see the distribution of gene expressions
#The default alpha for MA plots is 0.1

plotMA(DEG, ylim=c(-2,2), alpha = 0.05) #Setting the alpha value to 0.05

#Removing Noise

resLFC <- lfcShrink(DEG, 
                    coef = "Test_mock_vs_infected", 
                    type = "apeglm") #This gives a ref

plotMA(resLFC, ylim=c(-2,2), alpha = 0.05)

## MA Plot with ggplot 2
colnames(results_DEG)

par(mfrow=c(2,2))

#################### VOLCANO PLOT OF DIFFERENTIATION AT 4HRS ###################

## FILTERING THE RESULTS ## - For Volcano Plot
#Setting a column for the Volcano plot
results_DEG_4$diffexpressed <- "NO"
results_DEG_4$diffexpressed[results_DEG_4$log2FoldChange > 1 & results_DEG_4$pvalue < 0.05] <- "UP"
results_DEG_4$diffexpressed[results_DEG_4$log2FoldChange < -1 & results_DEG_4$pvalue < 0.05] <- "DOWN"

## VOLCANO PLOT OF DIFFERENTIATION ##

#Extracting GeneIDs from the row.names without needing for a new GeneID column
top10DEGs_4 <- results_DEG_4[order(results_DEG_4$pvalue), ][1:10,]
results_DEG_4$difflabel <- ifelse(row.names(results_DEG_4) %in% row.names(top10DEGs_4), row.names(results_DEG_4), NA)
summary(results_DEG_4$difflabel)

#Volcano Plot with ggplot2
VP4 <- ggplot(data = results_DEG_4, aes(x = log2FoldChange, 
                               y = -log10(pvalue), 
                               col = diffexpressed, 
                               label = difflabel)) +
  geom_vline(xintercept = c(1, -1), #Manually adding cutoff lines 
             col = "gray", 
             linetype = "dashed") + #Adding vertical lines to show fold change cut off
  geom_hline(yintercept = c(-log10(0.05)), 
             col = "gray", 
             linetype = "dashed") + #Adding Horizontal line to show p-value cut off
  geom_point(shape=19) + #Setting Point size
  guides(color = guide_legend(override.aes = list(size = 3))) +
  #scale_shape_manual(values = 6) +
  scale_color_manual(values = c("#00A", "gray", "#bb0c00"), #Changing plot colours
                     labels = c("Downregulated","Not Significant", "Upregulated")) + #Changing label titles
  coord_cartesian(ylim = c(0, 30), xlim = c(-10,10)) + #Applying figure axis limits
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # Setting continuous breaks (min, max, step)
  labs(x = expression("log"[2]*" Fold Change"), #Changing the x axis
       y = expression("-log"[10]*" p-adjusted")) + #Changing the y axis
  ggtitle("Timestamp: 4hrs") + # Setting a Figure Title
  geom_text_repel(max.overlaps = 2000) #Adding labels which we express in ggplot line 1

#################### VOLCANO PLOT OF DIFFERENTIATION AT 12HRS ###################

## FILTERING THE RESULTS ## - For Volcano Plot
#Setting a column for the Volcano plot
results_DEG_12$diffexpressed <- "NO"
results_DEG_12$diffexpressed[results_DEG_12$log2FoldChange > 1 & results_DEG_12$pvalue < 0.05] <- "UP"
results_DEG_12$diffexpressed[results_DEG_12$log2FoldChange < -1 & results_DEG_12$pvalue < 0.05] <- "DOWN"

## VOLCANO PLOT OF DIFFERENTIATION ##

#Extracting GeneIDs from the row.names without needing for a new GeneID column
top10DEGs_12 <- results_DEG_12[order(results_DEG_12$pvalue), ][1:10,]
results_DEG_12$difflabel <- ifelse(row.names(results_DEG_12) %in% row.names(top10DEGs_12), row.names(results_DEG_12), NA)
summary(results_DEG_12$difflabel)

#Volcano Plot with ggplot2
VP12 <- ggplot(data = results_DEG_12, aes(x = log2FoldChange, 
                               y = -log10(pvalue), 
                               col = diffexpressed, 
                               label = difflabel)) +
  geom_vline(xintercept = c(1, -1), #Manually adding cutoff lines 
             col = "gray", 
             linetype = "dashed") + #Adding vertical lines to show fold change cut off
  geom_hline(yintercept = c(-log10(0.05)), 
             col = "gray", 
             linetype = "dashed") + #Adding Horizontal line to show p-value cut off
  geom_point(shape=19) + #Setting Point size
  guides(color = guide_legend(override.aes = list(size = 3))) +
  #scale_shape_manual(values = 6) +
  scale_color_manual(values = c("#00A", "gray", "#bb0c00"), #Changing plot colours
                     labels = c("Downregulated","Not Significant", "Upregulated")) + #Changing label titles
  coord_cartesian(ylim = c(0, 30), xlim = c(-10,10)) + #Applying figure axis limits
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # Setting continuous breaks (min, max, step)
  labs(x = expression("log"[2]*" Fold Change"), #Changing the x axis
       y = expression("-log"[10]*" p-adjusted")) + #Changing the y axis
  ggtitle("Timestamp: 12hrs") + # Setting a Figure Title
  geom_text_repel(max.overlaps = 2000) #Adding labels which we express in ggplot line 1

#################### VOLCANO PLOT OF DIFFERENTIATION AT 48HRS ###################

## FILTERING THE RESULTS ## - For Volcano Plot
#Setting a column for the Volcano plot
results_DEG_48$diffexpressed <- "NO"
results_DEG_48$diffexpressed[results_DEG_48$log2FoldChange > 1 & results_DEG_48$pvalue < 0.05] <- "UP"
results_DEG_48$diffexpressed[results_DEG_48$log2FoldChange < -1 & results_DEG_48$pvalue < 0.05] <- "DOWN"

## VOLCANO PLOT OF DIFFERENTIATION ##

#Extracting GeneIDs from the row.names without needing for a new GeneID column
top10DEGs_48 <- results_DEG_48[order(results_DEG_48$pvalue), ][1:10,]
results_DEG_48$difflabel <- ifelse(row.names(results_DEG_48) %in% row.names(top10DEGs_48), row.names(results_DEG_48), NA)
summary(results_DEG_48$difflabel)

#Volcano Plot with ggplot2
VP48 <- ggplot(data = results_DEG_48, aes(x = log2FoldChange, 
                               y = -log10(pvalue), 
                               col = diffexpressed, 
                               label = difflabel)) +
  geom_vline(xintercept = c(1, -1), #Manually adding cutoff lines 
             col = "gray", 
             linetype = "dashed") + #Adding vertical lines to show fold change cut off
  geom_hline(yintercept = c(-log10(0.05)), 
             col = "gray", 
             linetype = "dashed") + #Adding Horizontal line to show p-value cut off
  geom_point(shape=19) + #Setting Point size
  guides(color = guide_legend(override.aes = list(size = 3))) +
  #scale_shape_manual(values = 5) +
  scale_color_manual(values = c("#00A", "gray", "#bb0c00")) + #Changing plot colours
  coord_cartesian(ylim = c(0, 150), xlim = c(-10,10)) + #Applying figure axis limits
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # Setting continuous breaks (min, max, step)
  labs( # Changing the colour legend title
       x = expression("log"[2]*" Fold Change"), #Changing the x axis
       y = expression("-log"[10]*" p-adjusted")) + #Changing the y axis
  ggtitle("Timestamp: 48hrs") + # Setting a Figure Title
  geom_text_repel(max.overlaps = 2000) #Adding labels which we express in ggplot line 1


#PLOT 3X VOLCANO PLOTS IN ONE WINDOW
ggarrange(VP4, VP12, VP48,
          labels = c("A", "B", "C"),
          ncol = 2, 
          nrow = 2,
          common.legend = TRUE, 
          legend = "bottom")


 
              
