
# TEMPORTAL RNA-SEQ ANALYSIS: DESeq Time Course

# DESeq2 can perform some time course analysis but fails for individual genes

# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

## Library ---------------------------------------------------------------------

BiocManager::install("fission")
library("fission")
library("DESeq2")
library("pheatmap")
data("fission")

head(fission)

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

## PHEATMAP SAVE FUNCTION ------------------------------------------------------

#Credit: https://gist.github.com/mathzero/a2070a24a6b418740c44a5c023f5c01e

save_pheatmap <- function(x, filename, width=12, height=12){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  if(grepl(".png",filename)){
    png(filename, width=width, height=height, units = "in", res=300)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  else if(grepl(".pdf",filename)){
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  else{
    print("Filename did not contain '.png' or '.pdf'")
  }
}


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

ggsave("4_figures/Temp_DESeq_Counts.png",
       height = 20,
       width = 20,
       units = "cm",
       dpi = 500)

## Statistical Testing

resultsNames(ddsTC)
res12 <- results(ddsTC, name="groupinfected.timepoint12", test="Wald")
res12[which.min(res12$padj),] # Greatest DEG at 12hrs

res48 <- results(ddsTC, name="groupinfected.timepoint48", test="Wald")
res48[which.min(res48$padj),] # Greatest DEG at 48hrs

## Further Analysis

betas <- coef(ddsTC)
colnames(betas)

#Note for this comparison, 4hrs is being used as the baseline.

topGenes <- head(order(resTC$padj),30)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
map <- pheatmap(mat, breaks=seq(from=-thr, 
                         to=thr, 
                         length=101),
         cluster_col=FALSE)

#Heatmap width and height best at 10
save_pheatmap(map, "4_figures/Temp_Heatmap.png", 10, 10)




