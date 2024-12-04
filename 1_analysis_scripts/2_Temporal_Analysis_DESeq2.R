
# TEMPORTAL RNA-SEQ ANALYSIS: DESeq Time Course

# DESeq2 can perform some time course analysis but fails for individual genes

# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

## Library ---------------------------------------------------------------------

library("DESeq2")
library("pheatmap")

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

## READING IN DATA -------------------------------------------------------------

rm(list = ls(all.names = T))

counts_data_TEMP <- read.csv('2_rawdata/1_counts/TEMP_counts.txt', header = T, row.names = 1)
target_data_TEMP <- read.csv('2_rawdata/2_targets/TEMP_targets.txt', header = T, row.names = 1, stringsAsFactors = T)

target_data_TEMP$timepoint <- as.factor(target_data_TEMP$timepoint)

## Overall Gene Count Trajectory -----------------------------------------------

DEG_TEMP <- DESeqDataSetFromMatrix(countData = counts_data_TEMP, 
                                   colData = target_data_TEMP, #Adding targets data
                                   design = ~ group + timepoint + group:timepoint) #Factors for Comparison

DEG_TEMP$group <- factor(DEG_TEMP$group, levels = c("mock", "infected"))
DEG_TEMP$timepoint <- factor(DEG_TEMP$timepoint, levels = c(4,12,48))

ddsTC <- DESeqDataSet(DEG_TEMP, ~ group + timepoint + group:timepoint)

ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ group + timepoint)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),], 4)

fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                  intgroup = c("group","timepoint"), 
                  returnData = TRUE)
fiss$minute <- as.numeric(as.character(fiss$timepoint))

ggplot(fiss,
       aes(x = minute, 
           y = count, 
           color = group, 
           group = group)) + 
  geom_point() + 
  stat_summary(fun.y=mean, 
               geom="line") +
  labs(x = "Time (hrs)",
       y = "Counts",
       color = "Group")
  scale_y_log10()

ggsave("4_figures/Figure_3_Temp_DESeq_Counts.png",
       height = 20,
       width = 20,
       units = "cm",
       dpi = 500)
