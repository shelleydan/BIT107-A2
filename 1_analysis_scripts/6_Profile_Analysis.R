
######################### KEGG Profile Analysis of DEGs ########################

## Reference -------------------------------------------------------------------

# https://github.com/ACSoupir/Bioinformatics_YouTube.git


## Library --------------------------------------------------------------------

library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggupset)

# Functions -------------------------------------------------------------------

# Function: Adjacency matrix to list
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Using the Enricher Function in Cluster Profiler ----------------------------

# Loading in our data (4HRS)
DEG4 <- read.csv("3_output_data/raw_DEG_results_4.csv", 
                 header = T)
names(DEG4)[names(DEG4) == 'X'] <- 'gene_symbol'
DEG4_Sig <- DEG4[DEG4$diffexpressed == 'UP',] #| DEG4$diffexpressed == 'DOWN',]

# Loading in our data (12HRS)
DEG12 <- read.csv("3_output_data/raw_DEG_results_12.csv", 
                 header = T)
names(DEG12)[names(DEG12) == 'X'] <- 'gene_symbol'
DEG12_Sig <- DEG12[DEG12$diffexpressed == 'UP'| DEG12$diffexpressed == 'DOWN',]

# Loading in our data (48HRS)
DEG48 <- read.csv("3_output_data/raw_DEG_results_48.csv", 
                 header = T)
names(DEG48)[names(DEG48) == 'X'] <- 'gene_symbol'
DEG48_Sig <- DEG48[DEG48$diffexpressed == 'UP'| DEG48$diffexpressed == 'DOWN',]

# Note: In the tutorial they use this step to add values of 'UP' and 'DOWN', but for plotting in DESeq2 I already performed this. 

# Using these values of 'UP' & 'DOWN', we subset the data.

DEG4_res_list <- split(DEG4_Sig, DEG4_Sig$diffexpressed)
DEG12_res_list <- split(DEG12_Sig, DEG12_Sig$diffexpressed)
DEG48_res_list <- split(DEG48_Sig, DEG48_Sig$diffexpressed)

## Aquiring Reference Genes ---------------------------------------------------

# Gene sets can be downloaded here: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# Gene Sets used in this analysis are found in '5_Gene_Sets'
DEG_Genes <- DEG4$gene_symbol #All time points have the same total genes

gmt_files <- list.files(path = '5_Gene_Lists/', 
                        pattern = '.gmt', 
                        full.names = T)

for (file in gmt_files) {
  pwl2 <- read.gmt(file)
  pwl2 <- pwl2[pwl2$gene %in% DEG_Genes,]
  filename_gs <- paste(gsub('c.\\.', '', gsub('.v2024.1.*$', '', file)), '.RDS', sep = '')
  saveRDS(pwl2, filename_gs)
}

## Environment Housekeeping ----------------------------------------------------
rm(DEG4_Sig)
rm(DEG12_Sig)
rm(DEG48_Sig)

## Filtering Gene Sets (REACTOME, KEGG, GO) ------------------------------------
# We look for profiles enriched in both our up regulated and down regulated. 

# Settings
gene_path <- '5_Gene_Lists/'
output <- '3_output_data/2_pathway_analysis/'

name_of_comparison <- 'mockvsinfected' # for our filename
background_genes <- 'reactome' # for our filename
bg_genes <- readRDS(paste0(gene_path, 'reactome.RDS')) # background genes.
padj_cutoff <- 0.05 # p-adjusted threshold.
genecount_cutoff <- 5 # minimum number of genes in the pathway.
filename <- paste0(output, 'clusterProfiler/', name_of_comparison, '_', background_genes) # filename of our PEA results

if(background_genes == 'KEGG'){
  bg_genes <- readRDS(paste0(gene_path, 'kegg_legacy.RDS'))
} else if(background_genes == 'reactome'){
  bg_genes <- readRDS(paste0(gene_path, 'reactome.RDS'))
} else if(background_genes == 'go.bp'){
  bg_genes <- readRDS(paste0(gene_path, 'go.bp.RDS'))
} else {
  stop('Invalid background genes. Select one of the following: KEGG, Reactome, GO, or add new pwl to function')
}

## Performing Cluster Analysis -------------------------------------------------

DEG4_res <- lapply(names(DEG4_res_list),
                   function(x) enricher(gene = DEG4_res_list[[x]]$gene_symbol,
                   TERM2GENE = bg_genes))

names(DEG4_res) <- names(DEG4_res_list) # Apply the UP and DOWN tags

DEG4_df_res <- lapply(names(DEG4_res), function(x) rbind(DEG4_res[[x]]@result))
names(DEG4_df_res) <- names(DEG4_res)
DEG4_df_res <- do.call(rbind,DEG4_df_res)

# Subset the pathways by (i) padj value, (ii) Gene count

target_pws <- unique(DEG4_df_res$ID[DEG4_df_res$p.adjust < padj_cutoff & DEG4_df_res$Count > genecount_cutoff])

DEG4_df_res <- DEG4_df_res[DEG4_df_res$ID %in% target_pws,]

## ENRICHMENT ------------------------------------------------------------------

enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = DEG4_df_res,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 organism = "human",
                 ontology = "UNKNOWN",
                 gene = DEG4$gene_symbol,
                 keytype = "UNKNOWN",
                 universe = unique(bg_genes$gene),
                 gene2Symbol = character(0),
                 geneSets = bg_genes)
class(enrichres)

## VISUALISATION ---------------------------------------------------------------

p1 <- barplot(enrichres, showCategory = 20)
p2 <- mutate(enrichres, qscore = -log(p.adjust, base = 10)) %>% 
  barplot(x = "qscore")
p3 <- dotplot(enrichres, showCategory = 15) + ggtitle("Mock vs Infected")

enrichres_eid <- setReadable(enrichres, 'org.Hs.eg.db', 'SYMBOL')
p4 <- cnetplot(enrichres)

?cnetplot


enrichres2 <- pairwise_termsim(enrichres)

treeplot(enrichres2)

emapplot(enrichres2)

upsetplot(enrichres)





















#OLD METHOD


## Library Prep ----------------------------------------------------------------
BiocManager::install("pathview")
BiocManager::install("gage")
BiocManager::install("gageData")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("GOplot")
library(pathview)
library(gage)
library(gageData)


library(biomartr)

library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)

## Importing GO Database -------------------------------------------------------

go<-read.csv("2_rawdata/5_GO_Data/go_map.csv")
fullgo<-c(as.list(GOTERM),as.list(GOOBSOLETE))

## Importing DEG Data ----------------------------------------------------------

DEG4 <- read.csv("3_output_data/raw_DEG_results_4.csv", 
                 header = T, row.names = 1)

DEG12 <- read.csv("3_output_data/raw_DEG_results_12.csv", 
                  header = T, row.names = 1)

DEG48 <- read.csv("3_output_data/raw_DEG_results_48.csv", 
                  header = T, row.names = 1)

## Filtering Data for Only Significant -----------------------------------------

DEG4_Sig_up <- rownames(DEG4[DEG4$diffexpressed == 'UP',])
DEG4_Sig_down <- rownames(DEG4[DEG4$diffexpressed == 'DOWN',])
DEG12_Sig_up <- rownames(DEG12[DEG12$diffexpressed == 'UP',])
DEG12_Sig_down <- rownames(DEG12[DEG12$diffexpressed == 'DOWN',])
DEG48_Sig <- DEG48[DEG48$diffexpressed == 'UP'| DEG48$diffexpressed == 'DOWN',]
DEG48_Sig_d<- DEG48[DEG48$diffexpressed == 'DOWN',]

rm(DEG4, DEG12, DEG48) # Remove un-filtered data frame

## GO Analysis 4HRS ------------------------------------------------------------

GO_4UP_res <- enrichGO(gene = DEG4_Sig_up, 
                       OrgDb = "org.Hs.eg.db",
                       keyType = 'SYMBOL',
                       ont = "BP")
GO4UP <- plot(barplot(GO_4UP_res, showCategory = 20))

GO_4DOWN_res <- enrichGO(gene = DEG4_Sig_down, 
                         OrgDb = "org.Hs.eg.db",
                         keyType = 'SYMBOL',
                         ont = "BP")
GO4DOWN <- plot(barplot(GO_4DOWN_res, showCategory = 20))

## GO Analysis 12HRS -----------------------------------------------------------

GO_12UP_res <- enrichGO(gene = DEG12_Sig_up, 
                        OrgDb = "org.Hs.eg.db",
                        keyType = 'SYMBOL',
                        ont = "BP")
GO12UP <- plot(barplot(GO_12UP_res, showCategory = 20))

GO_12DOWN_res <- enrichGO(gene = DEG12_Sig_down, 
                          OrgDb = "org.Hs.eg.db",
                          keyType = 'SYMBOL',
                          ont = "BP")
GO12DOWN <- plot(barplot(GO_12DOWN_res, showCategory = 20))

## GO Analysis 48HRS -----------------------------------------------------------

GO_48UP_res <- enrichGO(gene = DEG48_Sig_up, 
                        OrgDb = "org.Hs.eg.db",
                        keyType = 'SYMBOL',
                        ont = "BP")
GO48UP <- plot(barplot(GO_48UP_res, showCategory = 20))

GO_48DOWN_res <- enrichGO(gene = DEG48_Sig_down, 
                          OrgDb = "org.Hs.eg.db",
                          keyType = 'SYMBOL',
                          ont = "BP")
head(GO48DOWN)

GO48DOWN <- plot(barplot(GO_48DOWN_res, showCategory = 20))
GO48DOWN

data(EC)
head(EC$david)

testdata <- circle_dat(GO48DOWN, DEG48_Sig_d)

list <- as.data.frame(rownames(DEG48_Sig))
print(list)

EC$genelist

GOBubble(testdata, labels = 3)










## Producing an Array ----------------------------------------------------------
foldChange_4 <- DEG4_Sig$log2FoldChange
names(foldChange_4) <- symbols4$ENTREZID
head(foldChange_4)

foldChange_12 <- DEG12_Sig$log2FoldChange
names(foldChange_12) <- symbols12$ENTREZID
head(foldChange_12)

foldChange_48 <- DEG48_Sig$log2FoldChange
names(foldChange_48) <- symbols48$ENTREZID
head(foldChange_48)


## Bringing Gage Database ------------------------------------------------------

#The data for Homo sapiens.
data(go.sets.hs)
data(go.subs.hs)

## Analyzing only Biological Properties First ----------------------------------

#Pulling only BP
gobpsets <- go.sets.hs[go.subs.hs$BP]

#4HRS
gobpDEG4 <- gage(exprs = foldChange_4,
                 gsets =  gobpsets, 
                 same.dir = T)
view(gobpDEG4)

#12HRS
gobpDEG12 <- gage(exprs = foldChange_12,
                  gsets =  gobpsets, 
                  same.dir = T)
view(gobpDEG12)

#48HRS
gobpDEG48 <- gage(exprs = foldChange_48,
                  gsets =  gobpsets, 
                  same.dir = T)
view(gobpDEG4)

## KEGG Analysis ---------------------------------------------------------------

#Preparing for Analysis
data(kegg.sets.hs)

rm(keggRes4)
head(foldChange_4)
keggRes4 <- gage(exprs = foldChange_4, 
                 gsets = kegg.sets.hs, 
                 same.dir = T)
view(keggRes4$greater) # Note there's a lot of NAs but I think that might make sense

## Plotting Pathways -----------------------------------------------------------
keggRes4Pathways <- data.frame(id = rownames(keggRes4$greater), keggRes4$greater) %>%
  tibble::as_tibble() %>%
  filter(row_number() <= 2) %>% #Selecting the top 2 Pathways
  .$id %>%
  as.character()
keggRes4Pathways

keggRes4IDs <- substr(keggRes4Pathways, start = 1, stop = 8) #Extracting IDs
keggRes4IDs

## Plotting Pathways

tmp = sapply(keggRes4IDs, function(pid) pathview(gene.data = foldChange_4, pathway.id = pid, species = 'hsa'))





