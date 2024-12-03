################################################################################
######################### DIFFERENTIAL EXPRESSED GENES ######################### 
################################################################################

## REFERENCES ##
# Batut et al - https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html
# Cristofides - https://github.com/ecologysarah?tab=repositories

## LIBRARY ---------------------------------------------------------------------

install.packages("BiocManager") #Some Packages are now built under Bioconductor and needs to be installed differently. 

#Bioconductor Package Installs
BiocManager::install("DESeq2") #Downloads DESeq2 from Bioconductor
BiocManager::install("Rsamtools") #Downloads SAMTools from Bioconductor
BiocManager::install("apeglm") #Downloads apeglm from Bioconductor

#Library Loading
library(pheatmap)
library(dplyr)
library(tidyverse)
library(Rsamtools)
library(apeglm)
library(DESeq2)
library(dendextend)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(ggdendro)
library(plotly)
library(ggplotify)

## SETTING GGPLOT2 DOCUMENT THEME -----------------------------------------------

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


## IMPORTING DATA (COUNTS AND TARGETS) -----------------------------------------

#READ IN COUNTS DATA
counts_data <- read.csv("2_rawdata/GSE217504_host_counts_matrix.csv",
                        header = T, 
                        row.names = 1)

#READ IN AND EXTRACT TARGETS DATA
target_data <- read.delim("2_rawdata/targets.txt", 
                          sep = "", 
                          header = T) #Inputting the targets.

target_data <- target_data[,8:10] #Separating out the necessary columns.
target_data$Condition <- as.factor(target_data$Condition) #Setting 'Condition' (AKA Time) to factors.

names(target_data)[names(target_data) == 'Test'] <- 'group'
names(target_data)[names(target_data) == 'Condition'] <- 'timepoint'
names(target_data)[names(target_data) == 'Sample'] <- 'sample'

#SUBSET THE TARGETS DATA (4,12,48 hrs which all have controls)
target_data <- data.frame(target_data, row.names = 3) #Setting sample names as the rownames

target_data_4 <- target_data[target_data$timepoint %in% 4, ] #Targets for 4hrs
target_data_12 <- target_data[target_data$timepoint %in% 12, ] #Targets for 12hrs
target_data_48 <- target_data[target_data$timepoint %in% 48, ] #Tagrets for 48hrs
target_data_TEMP <- target_data[target_data$timepoint %in% c(4,12,48), ] #Tagrets for TEMPORAL

#ORDER THE COUNT DATA (DESeq2 requires).
counts_data_order = sort(colnames(counts_data)) #Re-order column names to alpha-numerical
counts_data<- counts_data[, counts_data_order]

#SUBSET THE COUNT DATA (4, 12, 48hrs).
counts_data_4 <- counts_data[, colnames(counts_data) %in% row.names(target_data_4)] #Counts for 4hrs
counts_data_12 <- counts_data[, colnames(counts_data) %in% row.names(target_data_12)] #Counts for 12hrs
counts_data_48 <- counts_data[, colnames(counts_data) %in% row.names(target_data_48)] #Counts for 48hrs
counts_data_TEMP <- counts_data[, colnames(counts_data) %in% row.names(target_data_TEMP)] #Counts for 48hrs

#TIDYING THE ENVIRONMENT
rm(counts_data)
rm(target_data)
rm(counts_data_order)

#SAVING COUNTS FOR EXTERNAL ANALYSIS
write.csv(counts_data_TEMP, "2_rawdata/1_counts/TEMP_counts.txt")
write.csv(counts_data_4, "2_rawdata/1_counts/4HRS_counts.txt")
write.csv(counts_data_12, "2_rawdata/1_counts/12HRS_counts.txt")
write.csv(counts_data_48, "2_rawdata/1_counts/48HRS_counts.txt")

#SAVING TARGETS FOR EXTERNAL ANALYSIS
write.csv(target_data_TEMP, "2_rawdata/2_targets/TEMP_targets.txt")
write.csv(target_data_4, "2_rawdata/2_targets/4HRS_targets.txt")
write.csv(target_data_12, "2_rawdata/2_targets/12HRS_targets.txt")
write.csv(target_data_48, "2_rawdata/2_targets/48HRS_targets.txt")

## ASSIGNING FACTOR LEVELS -----------------------------------------------------

#FACTOR LEVELS FOR 4HRS
target_data_4$group <- factor(target_data_4$group) 
target_data_4$timepoint <- factor(target_data_4$timepoint) 

#FACTOR LEVELS FOR 12HRS
target_data_12$group <- factor(target_data_12$group) 
target_data_12$timepoint <- factor(target_data_12$timepoint) 

#FACTOR LEVELS FOR 48HRS
target_data_48$group <- factor(target_data_48$group) 
target_data_48$timepoint <- factor(target_data_48$timepoint) 

#FACTOR LEVELS FOR TEMPORAL
target_data_TEMP$group <- factor(target_data_TEMP$group) 
target_data_TEMP$timepoint <- factor(target_data_TEMP$timepoint)

## CREATING DESEQ2 MATRICIES ---------------------------------------------------

#If we're looking at a multi-factor analysis, we want to input our primary factor should be inputted last.

#The first level in the list is what becomes the base level, which the other levels will be compared to. We want this to be the control/mock group

#4HRS 
DEG_4 <- DESeqDataSetFromMatrix(countData = counts_data_4, 
                                colData = target_data_4, #Adding targets data
                                design = ~group) #Factors for Comparison
DEG_4$group <- factor(DEG_4$group, levels = c("mock", "infected")) 

#12HRS
DEG_12 <- DESeqDataSetFromMatrix(countData = counts_data_12, 
                                 colData = target_data_12, #Adding targets data
                                 design = ~group) #Factors for Comparison
DEG_12$group <- factor(DEG_12$group, levels = c("mock", "infected")) 

#48HRS
DEG_48 <- DESeqDataSetFromMatrix(countData = counts_data_48, 
                                 colData = target_data_48, #Adding targets data
                                 design = ~group) #Factors for Comparison
DEG_48$group <- factor(DEG_48$group, levels = c("mock", "infected")) 

#TEMPORAL
DEG_TEMP <- DESeqDataSetFromMatrix(countData = counts_data_TEMP, 
                                   colData = target_data_TEMP, #Adding targets data
                                   design = ~ group + timepoint + group:timepoint) #Factors for Comparison
DEG_TEMP$group <- factor(DEG_TEMP$group, levels = c("mock", "infected"))
DEG_TEMP$timepoint <- factor(DEG_TEMP$timepoint, levels = c(4,12,48))

## PERFORMING DESEQ2 ANALYSIS --------------------------------------------------

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

#ACROSS TIME
DEG_TEMP <- DESeq2::DESeq(DEG_TEMP) #Perform the Analysis
results_DEG_TEMP <- DESeq2::results(DEG_TEMP) #Coalating the results into a dataframe
summary(results_DEG_TEMP$padj)
results_DEG_TEMP <- as.data.frame(results_DEG_TEMP) # Produces and R dataframe

## SAVING THE RAW-COUNTS AND FILTERED ------------------------------------------

#Output of the non-filtered DESeq2 results
write.csv(results_DEG_4, "3_output_data/raw_DEG_results_4.csv") 
write.csv(results_DEG_12, "3_output_data/raw_DEG_results_12.csv")
write.csv(results_DEG_48, "3_output_data/raw_DEG_results_48.csv")
write.csv(results_DEG_TEMP, "3_output_data/raw_DEG_results_TEMPORAL.csv")

## COMPARATIVE SAMPLE HCA ANALYSIS ---------------------------------------------
counts_TEMP <- t(counts_data_TEMP) # Needs to be transposed for HCA Analysis 

HCA_TEMP <- scale(counts_TEMP)

# Dissimilarity matrix - There is little difference between using 10000 or 20000 genes. 
dist_TEMP <- dist(HCA_TEMP[,1:20000], method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(dist_TEMP, method = "complete" )
hc1 <- as.dendrogram(hc1)
hc1 <- dendro_data(hc1) # Ready to Plot

order <- hc1$labels$label
hca_targets <- target_data_TEMP[order,]
hc1$labels$timepoint <- hca_targets$timepoint
hc1$labels$group <- hca_targets$group

cols <- c("4.infected"="#f04546","12.infected"="#3591d1","48.infected"="#62c76b","4.mock"="#f04546","12.mock"="#3591d1","48.mock"="#62c76b")

cols <- interaction(hc1$labels$timepoint, hc1$labels$group)

ggHCA <- ggplot(segment(hc1)) + 
  geom_text(data = label(hc1),
            show.legend = T,
            aes(label = label, 
                x = x, 
                y = 0, 
                color = cols, 
                angle = 90,
                vjust = 1.5,
                hjust = 0,
                fontface = 2)) +
  guides(colour = guide_legend(title="Condition", 
                               override.aes=list(alpha=1, shape=19, size=5))) +
  geom_segment(show.legend = F,
               colour = "gray",
               aes(x = x, 
                   y = y,
                   xend = xend, 
                   yend = yend)) +
  labs(x = "Samples",
       y = "Distance") +
  theme(aspect.ratio=1,
        legend.position="bottom",
        legend.box.background = element_rect(colour = "black"))

ggHCA # This will be arranges later with the PCA in one plot

## PCA of All Samples ----------------------------------------------------------

#TEMPORAL
# Variance stabilisation transformation
vst_TEMP <- DESeq2::vst(DEG_TEMP, blind = F)

# Generating the PCA Plot
TEMPData <- plotPCA(vst_TEMP, 
                    intgroup = c("group", "timepoint"),
                    returnData = T)

TEMPData$group.1 #Some reason the infected/mock is stored as group.1 instead and group is a combination of both factor levels. 

percentVar <- round(100 * attr(TEMPData, "percentVar"))
ggPCA <- ggplot(TEMPData, aes(PC1, PC2, color=group.1, shape=timepoint)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  theme(aspect.ratio=1,
        legend.position="top",
        legend.box.background = element_rect(colour = "black"))

ggPCA # This will be arranged/viewd later

## MOCK DEG COMPARISON 4hrs vs 12hrs -------------------------------------------
MOCK_DEG_TARGET <- target_data_TEMP[target_data_TEMP$timepoint %in% c(4,12), ] #Targets for 4hrs
MOCK_DEG_TARGET <- MOCK_DEG_TARGET[MOCK_DEG_TARGET$group %in% "mock", ] #Targets for 12hrs

MOCK_DATA <- counts_data_TEMP[, colnames(counts_data_TEMP) %in% row.names(MOCK_DEG_TARGET)]

#Performing Analysis w/ DESeq

DEG_MOCK <- DESeqDataSetFromMatrix(countData = MOCK_DATA, 
                                colData = MOCK_DEG_TARGET, #Adding targets data
                                design = ~timepoint) #Factors for Comparison
DEG_MOCK$group <- factor(DEG_MOCK$group, levels = c(4,12)) 

DEG_MOCK <- DESeq2::DESeq(DEG_MOCK) #Perform the Analysis
results_DEG_MOCK <- DESeq2::results(DEG_MOCK) #Coalating the results into a dataframe
summary(results_DEG_MOCK$padj)
results_DEG_MOCK <- as.data.frame(results_DEG_MOCK) # Produces and R dataframe

#Volcano 4hrs vs 12hrs MOCK
#Setting a column for the Volcano plot
results_DEG_MOCK$diffexpressed <- "NO"
results_DEG_MOCK$diffexpressed[results_DEG_MOCK$log2FoldChange > 1 & results_DEG_MOCK$pvalue < 0.05] <- "UP"
results_DEG_MOCK$diffexpressed[results_DEG_MOCK$log2FoldChange < -1 & results_DEG_MOCK$pvalue < 0.05] <- "DOWN"

## VOLCANO PLOT OF DIFFERENTIATION ##

#Extracting GeneIDs from the row.names without needing for a new GeneID column
top10DEGs_4 <- results_DEG_MOCK[order(results_DEG_MOCK$padj), ][1:10,]
results_DEG_MOCK$difflabel <- ifelse(row.names(results_DEG_MOCK) %in% row.names(top10DEGs_4), row.names(results_DEG_MOCK), NA)
summary(results_DEG_MOCK$difflabel)

#Volcano Plot with ggplot2
plot_4_12 <- ggplot(data = results_DEG_MOCK, aes(x = log2FoldChange, 
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
  scale_color_manual(values = c("cornflowerblue", "gray", "#bb0c00"), #Changing plot colours
                     labels = c("Downregulated","Not Significant", "Upregulated")) + #Changing label titles
  coord_cartesian(ylim = c(0, 250), xlim = c(-10,10)) + #Applying figure axis limits
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # Setting continuous breaks (min, max, step)
  labs(x = expression("log"[2]*" Fold Change"), #Changing the x axis
       y = expression("-log"[10]*" p-adjusted")) + #Changing the y axis
  ggtitle("4Hrs vs 12Hrs") + # Setting a Figure Title
  geom_text_repel(max.overlaps = 2000, size=5, show.legend = F) #Adding labels which we express in ggplot line 1

## MOCK DEG COMPARISON 12hrs vs 48hrs ------------------------------------------
MOCK_DEG_TARGET1 <- target_data_TEMP[target_data_TEMP$timepoint %in% c(12,48), ] #Targets for 4hrs
MOCK_DEG_TARGET1 <- MOCK_DEG_TARGET1[MOCK_DEG_TARGET1$group %in% "mock", ] #Targets for 12hrs

MOCK_DATA1 <- counts_data_TEMP[, colnames(counts_data_TEMP) %in% row.names(MOCK_DEG_TARGET1)]

#Performing Analysis w/ DESeq
DEG_MOCK1 <- DESeqDataSetFromMatrix(countData = MOCK_DATA1, 
                                   colData = MOCK_DEG_TARGET1, #Adding targets data
                                   design = ~timepoint) #Factors for Comparison
DEG_MOCK1$timepoint <- factor(DEG_MOCK1$timepoint, levels = c(12,48)) 

DEG_MOCK1 <- DESeq2::DESeq(DEG_MOCK1) #Perform the Analysis
results_DEG_MOCK1 <- DESeq2::results(DEG_MOCK1) #Coalating the results into a dataframe
summary(results_DEG_MOCK1$padj)
results_DEG_MOCK1 <- as.data.frame(results_DEG_MOCK1) # Produces and R dataframe

#Volcano 4hrs vs 12hrs MOCK
#Setting a column for the Volcano plot
results_DEG_MOCK1$diffexpressed <- "NO"
results_DEG_MOCK1$diffexpressed[results_DEG_MOCK1$log2FoldChange > 1 & results_DEG_MOCK1$padj < 0.05] <- "UP"
results_DEG_MOCK1$diffexpressed[results_DEG_MOCK1$log2FoldChange < -1 & results_DEG_MOCK1$padj < 0.05] <- "DOWN"

## VOLCANO PLOT OF DIFFERENTIATION ##

#Extracting GeneIDs from the row.names without needing for a new GeneID column
top10DEGs_mock1 <- results_DEG_MOCK1[order(results_DEG_MOCK1$padj), ][1:10,]
results_DEG_MOCK1$difflabel <- ifelse(row.names(results_DEG_MOCK1) %in% row.names(top10DEGs_mock1), row.names(results_DEG_MOCK), NA)
summary(results_DEG_MOCK1$difflabel)

#Volcano Plot with ggplot2
plot_12_48 <- ggplot(data = results_DEG_MOCK1, aes(x = log2FoldChange, 
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
  scale_color_manual(values = c("cornflowerblue", "gray", "#bb0c00"), #Changing plot colours
                     labels = c("Downregulated","Not Significant", "Upregulated")) + #Changing label titles
  coord_cartesian(ylim = c(0, 250), xlim = c(-10,10)) + #Applying figure axis limits
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # Setting continuous breaks (min, max, step)
  labs(x = expression("log"[2]*" Fold Change"), #Changing the x axis
       y = expression("-log"[10]*" p-adjusted")) + #Changing the y axis
  ggtitle("12Hrs vs 48Hrs") + # Setting a Figure Title
  geom_text_repel(max.overlaps = 2000, size=5, show.legend = F) #Adding labels which we express in ggplot line 1

plot_12_48

ggarrange(ggHCA, ggPCA, plot_4_12, plot_12_48,
          labels = c("A", "B", "C", "D"),
          ncol = 2, 
          nrow = 2,
          common.legend = F, 
          legend = "bottom") + 
  bgcolor("White") +
  border("White")

#Save the above figure
ggsave("4_figures/SAMPLE COMPARISONS/Figure_1_HCA_PCA_VP.png",
       height = 40,
       width = 40,
       units = "cm",
       dpi = 500)

## DISPERSION PLOT -------------------------------------------------------------

#We expect that when a gene's read count increases the dispersion of that same gene decreases

# These plots are for review only and are not included in the final text. 

par(mfrow=c(2,2))
#4HRS
plotDispEsts(DEG_4, 
             ylab = "Dispersion",
             xlab = "Mean of Normalised Counts",
             title.main = "Test") #This is run on the Analysis Data
?plotDispEsts

#12HRS
plotDispEsts(DEG_12, 
             ylab = "Dispersion",
             xlab = "Mean of Normalised Counts")

#48HRS
plotDispEsts(DEG_48, 
             ylab = "Dispersion",
             xlab = "Mean of Normalised Counts")

par(mfrow=c(1,1))

## HEATMAPS --------------------------------------------------------------------

#This used pheatmaps to analyse some gene expression clusting.

# Log transformed normalised counts (using top 10 genes)

DEG_padj_orders <- results_DEG[order(results_DEG$padj),] #Ascending order of padj values

topDEGs <- results_DEG_TEMP[order(results_DEG_TEMP$padj), ][1:10,] #Selecting the top 10
topDEGs_names <- row.names(topDEGs) #Extracting names of the top 10 genes

rld <- rlog(DEG_TEMP, blind = F) #Performling a log transformation

anno_info <- as.data.frame(colData(DEG_TEMP)[, c("timepoint", "group")]) #Setting Annotation Levels
anno_info$timepoint <- as.character(anno_info$timepoint) #Changing Condition from continuous to ordinal

pheatmap(assay(rld)[topDEGs_names,], #Subset by labels extracted
         cluster_rows = T, #Adds column tree-clustering
         show_rownames = T, 
         cluster_cols = T, #Adds rows tree-clustering
         annotation_col = anno_info) 

# Sample-to-sample distance matrix (normalised counts)

sampleDist <- dist(t(assay(vst_TEMP)))

sampleDistMatrix <- as.matrix(sampleDist) # Generate a matrix
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255) #Setting the colours and range of those. 

#Generating the Heatmap
pheatmap(sampleDistMatrix, #Matrix input
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         color = colours,
         annotation_col = anno_info) #as defined above

#This heatmap colours show the different between the samples. So the darkest blue shows no difference (e.g. when the same samples are plotted against eachother)

## MA PLOT ---------------------------------------------------------------------

#These plots are used to see the distribution of gene expressions
#The default alpha for MA plots is 0.1

plotMA(DEG_TEMP, ylim=c(-2,2), alpha = 0.05) #Setting the alpha value to 0.05

#Removing Noise

DEG_TEMP$group

resLFC <- lfcShrink(DEG_TEMP, 
                    coef = "group_infected_vs_mock", 
                    type = "apeglm") #This gives a ref

plotMA(resLFC, ylim=c(-3,3), alpha = 0.05)

## VOLCANO PLOT 4HRS -----------------------------------------------------------

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
  scale_color_manual(values = c("cornflowerblue", "gray", "#bb0c00"), #Changing plot colours
                     labels = c("Downregulated","Not Significant", "Upregulated")) + #Changing label titles
  coord_cartesian(ylim = c(0, 30), xlim = c(-10,10)) + #Applying figure axis limits
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # Setting continuous breaks (min, max, step)
  labs(x = expression("log"[2]*" Fold Change"), #Changing the x axis
       y = expression("-log"[10]*" p-adjusted")) + #Changing the y axis
  ggtitle("Timestamp: 4hrs") + # Setting a Figure Title
  geom_text_repel(max.overlaps = 2000, size=5) #Adding labels which we express in ggplot line 1

## VOLCANO PLOT 12HRS ----------------------------------------------------------

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
  scale_color_manual(values = c("cornflowerblue", "gray", "#bb0c00"), #Changing plot colours
                     labels = c("Downregulated","Not Significant", "Upregulated")) + #Changing label titles
  coord_cartesian(ylim = c(0, 30), xlim = c(-10,10)) + #Applying figure axis limits
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # Setting continuous breaks (min, max, step)
  labs(x = expression("log"[2]*" Fold Change"), #Changing the x axis
       y = expression("-log"[10]*" p-adjusted")) + #Changing the y axis
  ggtitle("Timestamp: 12hrs") + # Setting a Figure Title
  geom_text_repel(max.overlaps = 2000, size=5) #Adding labels which we express in ggplot line 1

## VOLCANO PLOT 48HRS ----------------------------------------------------------

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
  scale_color_manual(values = c("cornflowerblue", "gray", "#bb0c00")) + #Changing plot colours
  coord_cartesian(ylim = c(0, 150), xlim = c(-10,10)) + #Applying figure axis limits
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # Setting continuous breaks (min, max, step)
  labs( # Changing the colour legend title
    x = expression("log"[2]*" Fold Change"), #Changing the x axis
    y = expression("-log"[10]*" p-adjusted")) + #Changing the y axis
  ggtitle("Timestamp: 48hrs") + # Setting a Figure Title
  geom_text_repel(max.overlaps = 2000, size=5) #Adding labels which we express in ggplot line 

## PLOTTING ALL VOLCANO PLOTS IN ONE WINDOW ------------------------------------

legendVP <- get_legend(VP12)
legendVP <- as.ggplot(legendVP)

ggarrange(VP4, VP12, VP48, legendVP,
          labels = c("A", "B", "C"),
          ncol = 2, 
          nrow = 2,
          common.legend = T, 
          legend = "none") + 
  bgcolor("White") +
  border("White")

#Save the above figure
ggsave("4_figures/DEG/Figure_2_VP_4_12_48.png",
       height = 30,
       width = 30,
       units = "cm",
       dpi = 500)
