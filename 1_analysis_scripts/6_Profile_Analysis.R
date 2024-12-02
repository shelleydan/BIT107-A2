
######################### KEGG Profile Analysis of DEGs ########################

## Reference -------------------------------------------------------------------

rm(list = ls(all.names = T))
gc()
## Library --------------------------------------------------------------------

library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggupset)
library(ggpubr)

## Setting Path Variables ------------------------------------------------------
gene_path <- '5_Gene_Lists/'
output <- '3_output_data/2_pathway_analysis/'



# Functions -------------------------------------------------------------------

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
DEG4df <- as.data.frame(cbind(DEG4$gene_symbol, DEG4$pvalue, DEG4$padj, DEG4$log2FoldChange, DEG4$diffexpressed))
names(DEG4df)[names(DEG4df) == 'V1'] <- 'gene_symbol'
names(DEG4df)[names(DEG4df) == 'V2'] <- 'pval'
names(DEG4df)[names(DEG4df) == 'V3'] <- 'padj'
names(DEG4df)[names(DEG4df) == 'V4'] <- 'log2fc'
names(DEG4df)[names(DEG4df) == 'V5'] <- 'diffexpressed'

DEG4df <- DEG4df[DEG4df$diffexpressed == 'UP'| DEG4df$diffexpressed == 'DOWN',]
deg_results_list <- split(DEG4df, DEG4df$diffexpressed)

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
# These need to be the whole number of genes used in the experiement - the 46863 in this assingment
DEG_Genes <- DEG4$gene_symbol #All time points have the same total genes

gmt_files <- list.files(path = '5_Gene_Lists/', 
                        pattern = '.gmt', 
                        full.names = T)

# Reads in the gmt file and filters for our genes and makes a new subset
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
rm(DEG12)
rm(DEG48)
rm(file)
rm(filename_gs)
rm(gmt_files)
rm(pwl2)

## Filtering Gene Sets (REACTOME, KEGG, GO) ------------------------------------
# We look for profiles enriched in both our up regulated and down regulated. 

# Settings
name_of_comparison <- 'mock_vs_infected_4hrs' # for our filename

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

res <- lapply(names(deg_results_list),
                    function(x) enricher(gene = deg_results_list[[x]]$gene_symbol,
                   TERM2GENE = bg_genes))

names(res) <- names(deg_results_list) # Apply the UP and DOWN tags

res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res) # Apply the UP and DOWN tags
res_df <- do.call(rbind, res_df)

res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', 
                                                 rownames(res_df)))


# Subset the pathways by (i) padj value, (ii) Gene count
target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff])

res_df <- res_df[res_df$ID %in% target_pws,]

# Moving Onward from here, use res_df

# Select only upregulated genes in Severe
res_df <- res_df %>% filter(diffexpressed == 'UP') %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df) <- res_df$ID

# For visualisation purposes, let's shorten the pathway names
res_df$Description <- gsub('(H|h)iv', 'HIV', 
                           gsub('pd 1', 'PD-1',
                                gsub('ecm', 'ECM', 
                                     gsub('(I|i)nterleukin', 'IL', 
                                          gsub('(R|r)na', 'RNA', 
                                               gsub('(D|d)na', 'DNA',
                                                    gsub(' i ', ' I ', 
                                                         gsub('(A|a)tp ', 'ATP ', 
                                                              gsub('(N|n)adh ', 'NADH ', 
                                                                   gsub('(N|n)ad ', 'NAD ',
                                                                        gsub('t cell', 'T cell',
                                                                             gsub('b cell', 'B cell',
                                                                                  gsub('built from .*', ' (...)',
                                                                                       gsub('mhc', 'MHC',
                                                                                            gsub('mhc class i', 'MHC I', 
                                                                                                 gsub('mhc class ii', 'MHC II', 
                                                                                                      stringr::str_to_sentence(
                                                                                                        gsub('_', ' ',  
                                                                                                             gsub('GOBP_|KEGG_|REACTOME_', '', res_df$Description)))))))))))))))))))


## ENRICHMENT ------------------------------------------------------------------

enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = res_df,
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

p3 <- dotplot(enrichres, showCategory = 15)

p4 <- enrichplot::cnetplot(enrichres)

heatplot(enrichres, showCategory = 5)

enrichres2 <- pairwise_termsim(enrichres)

treeplot(enrichres2)

emapplot(enrichres2)

upsetplot(enrichres)

ggarrange(p1, p3,
          labels = c("A", "B"),
          ncol = 2, 
          nrow = 1,
          common.legend = TRUE, 
          legend = "bottom") + 
  bgcolor("White") +
  border("White")

ggsave("4_figures/emaplot4hr.png",
       height = 15,
       width = 30,
       units = "cm",
       dpi = 500)
