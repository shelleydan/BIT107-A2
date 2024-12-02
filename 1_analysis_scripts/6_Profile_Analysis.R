
######################### KEGG Profile Analysis of DEGs ########################

## Reference -------------------------------------------------------------------

# Start Fresh
rm(list = ls(all.names = T))

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

## Functions -------------------------------------------------------------------

#DID WE USE THIS?
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Loading in our data (4HRS) --------------------------------------------------

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
deg_results_list_4 <- split(DEG4df, DEG4df$diffexpressed)

## Loading in our data (12HRS) -------------------------------------------------

DEG12 <- read.csv("3_output_data/raw_DEG_results_12.csv", 
                  header = T)
names(DEG12)[names(DEG12) == 'X'] <- 'gene_symbol'
DEG12df <- as.data.frame(cbind(DEG12$gene_symbol, DEG12$pvalue, DEG12$padj, DEG12$log2FoldChange, DEG12$diffexpressed))
names(DEG12df)[names(DEG12df) == 'V1'] <- 'gene_symbol'
names(DEG12df)[names(DEG12df) == 'V2'] <- 'pval'
names(DEG12df)[names(DEG12df) == 'V3'] <- 'padj'
names(DEG12df)[names(DEG12df) == 'V4'] <- 'log2fc'
names(DEG12df)[names(DEG12df) == 'V5'] <- 'diffexpressed'

DEG12df <- DEG12df[DEG12df$diffexpressed == 'UP'| DEG12df$diffexpressed == 'DOWN',]
deg_results_list_12 <- split(DEG12df, DEG12df$diffexpressed)

## Loading in our data (48HRS) -------------------------------------------------

DEG48 <- read.csv("3_output_data/raw_DEG_results_48.csv", 
                  header = T)
names(DEG48)[names(DEG48) == 'X'] <- 'gene_symbol'
DEG48df <- as.data.frame(cbind(DEG48$gene_symbol, DEG48$pvalue, DEG48$padj, DEG48$log2FoldChange, DEG48$diffexpressed))
names(DEG48df)[names(DEG48df) == 'V1'] <- 'gene_symbol'
names(DEG48df)[names(DEG48df) == 'V2'] <- 'pval'
names(DEG48df)[names(DEG48df) == 'V3'] <- 'padj'
names(DEG48df)[names(DEG48df) == 'V4'] <- 'log2fc'
names(DEG48df)[names(DEG48df) == 'V5'] <- 'diffexpressed'

DEG48df <- DEG48df[DEG48df$diffexpressed == 'UP'| DEG48df$diffexpressed == 'DOWN',]
deg_results_list_48 <- split(DEG48df, DEG48df$diffexpressed)

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
rm(DEG12)
rm(DEG48)
rm(file)
rm(filename_gs)
rm(gmt_files)
rm(pwl2)

## Filtering Gene Sets (REACTOME, KEGG, GO) ------------------------------------
# We look for profiles enriched in both our up regulated and down regulated. 

# Settings
name_of_comparison <- 'mock_vs_infected' # for our filename

background_genes <- 'reactome' # for our filename

bg_genes <- readRDS(paste0(gene_path, 'reactome.RDS')) # background genes.

up_or_down <- 'DOWN' # Set this to UP or DOWN

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

## Performing Cluster Analysis (AUTO) ------------------------------------------

res4 <- lapply(names(deg_results_list_4),
               function(x) enricher(gene = deg_results_list_4[[x]]$gene_symbol,
                                    TERM2GENE = bg_genes))
res12 <- lapply(names(deg_results_list_12),
               function(x) enricher(gene = deg_results_list_12[[x]]$gene_symbol,
                                    TERM2GENE = bg_genes))
res48 <- lapply(names(deg_results_list_48),
               function(x) enricher(gene = deg_results_list_48[[x]]$gene_symbol,
                                    TERM2GENE = bg_genes))

names(res4) <- names(deg_results_list_4) # Apply the UP and DOWN tags
names(res12) <- names(deg_results_list_12) # Apply the UP and DOWN tags
names(res48) <- names(deg_results_list_48) # Apply the UP and DOWN tags

res_df4 <- lapply(names(res4), function(x) rbind(res4[[x]]@result))
names(res_df4) <- names(res4) # Apply the UP and DOWN tags
res_df4 <- do.call(rbind, res_df4)

res_df12 <- lapply(names(res12), function(x) rbind(res12[[x]]@result))
names(res_df12) <- names(res12) # Apply the UP and DOWN tags
res_df12 <- do.call(rbind, res_df12)

res_df48 <- lapply(names(res48), function(x) rbind(res48[[x]]@result))
names(res_df48) <- names(res48) # Apply the UP and DOWN tags
res_df48 <- do.call(rbind, res_df48)

res_df4 <- res_df4 %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', 
                                                 rownames(res_df4)))

res_df12 <- res_df12 %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', 
                                                 rownames(res_df12)))

res_df48 <- res_df48 %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', 
                                                 rownames(res_df48)))


# Subset the pathways by (i) padj value, (ii) Gene count
target_pws4 <- unique(res_df4$ID[res_df4$p.adjust < padj_cutoff & res_df4$Count > genecount_cutoff])
res_df4 <- res_df4[res_df4$ID %in% target_pws4,]

target_pws12 <- unique(res_df12$ID[res_df12$p.adjust < padj_cutoff & res_df12$Count > genecount_cutoff])
res_df12 <- res_df12[res_df12$ID %in% target_pws12,]

target_pws48 <- unique(res_df48$ID[res_df48$p.adjust < padj_cutoff & res_df48$Count > genecount_cutoff])
res_df48 <- res_df48[res_df48$ID %in% target_pws48,]

# HERE FILTER BASED ON UP OR DOWN REGULATED AS SPECIFIED IN SETTINGS

res_df4 <- res_df4 %>% filter(diffexpressed == up_or_down) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df4) <- res_df4$ID

# For visualisation purposes, let's shorten the pathway names
res_df4$Description <- gsub('(H|h)iv', 'HIV', 
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
                                                                                                             gsub('GOBP_|KEGG_|REACTOME_', '', res_df4$Description)))))))))))))))))))

res_df12 <- res_df12 %>% filter(diffexpressed == up_or_down) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df12) <- res_df12$ID

# For visualisation purposes, let's shorten the pathway names
res_df12$Description <- gsub('(H|h)iv', 'HIV', 
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
                                                                                                             gsub('GOBP_|KEGG_|REACTOME_', '', res_df12$Description)))))))))))))))))))

res_df48 <- res_df48 %>% filter(diffexpressed == up_or_down) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df48) <- res_df48$ID

# For visualisation purposes, let's shorten the pathway names
res_df48$Description <- gsub('(H|h)iv', 'HIV', 
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
                                                                                                             gsub('GOBP_|KEGG_|REACTOME_', '', res_df48$Description)))))))))))))))))))


## ENRICHMENT ANALYSIS

enrichres4 <- new("enrichResult",
                 readable = FALSE,
                 result = res_df4,
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
class(enrichres4) # CHECK THIS HAS BECOME AN ENRICHMENT CLASS

enrichres12 <- new("enrichResult",
                 readable = FALSE,
                 result = res_df12,
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
class(enrichres12) # CHECK THIS HAS BECOME AN ENRICHMENT CLASS

enrichres48 <- new("enrichResult",
                 readable = FALSE,
                 result = res_df48,
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
class(enrichres48) # CHECK THIS HAS BECOME AN ENRICHMENT CLASS

## VISUALISATION ---------------------------------------------------------------

# This process will automaticall produce figures (arranges) for all options above

## 4HRS PLOTS ------------------------------------------------------------------

p4.1 <- barplot(enrichres4, showCategory = 20)

p4.2 <- dotplot(enrichres4, showCategory = 15)

p4.3 <- enrichplot::cnetplot(enrichres4)

ggarrange(
  p4.3,
  ggarrange(p4.1, p4.2, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A") + 
  bgcolor("White") +
  border("White")

file4.1 <- paste0('4_figures/PATHWAY_ENRICHMENT/3_plot_4hrs_', up_or_down, '.png')

ggsave(file4.1,
       height = 30,
       width = 30,
       units = "cm",
       dpi = 500)

enrichres4.2 <- pairwise_termsim(enrichres4)

p4.4 <- emapplot(enrichres4.2)

p4.5 <- upsetplot(enrichres4)

ggarrange(p4.4, p4.5,
          labels = c("A", "B"),
          ncol = 1, 
          nrow = 2,
          common.legend = F, 
          legend = "bottom") + 
  bgcolor("White") +
  border("White")

file4.2 <- paste0('4_figures/PATHWAY_ENRICHMENT/2_plot_4hrs_', up_or_down, '.png')

ggsave(file4.2,
       height = 30,
       width = 30,
       units = "cm",
       dpi = 500)

## 12HRS PLOTS -----------------------------------------------------------------

p12.1 <- barplot(enrichres12, showCategory = 20)

p12.2 <- dotplot(enrichres12, showCategory = 15)

p12.3 <- enrichplot::cnetplot(enrichres12)

ggarrange(
  p12.3,
  ggarrange(p12.1, p12.2, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A") + 
  bgcolor("White") +
  border("White")

file12.1 <- paste0('4_figures/PATHWAY_ENRICHMENT/3_plot_12hrs_', up_or_down, '.png')

ggsave(file12.1,
       height = 30,
       width = 30,
       units = "cm",
       dpi = 500)

enrichres12.2 <- pairwise_termsim(enrichres12)

p12.4 <- emapplot(enrichres12.2)

p12.5 <- upsetplot(enrichres12)

ggarrange(p12.4, p12.5,
          labels = c("A", "B"),
          ncol = 1, 
          nrow = 2,
          common.legend = F, 
          legend = "bottom") + 
  bgcolor("White") +
  border("White")

file12.2 <- paste0('4_figures/PATHWAY_ENRICHMENT/2_plot_12hrs_', up_or_down, '.png')

ggsave(file12.2,
       height = 30,
       width = 30,
       units = "cm",
       dpi = 500)

## 48HRS PLOTS -----------------------------------------------------------------

p48.1 <- barplot(enrichres48, showCategory = 10)

p48.2 <- dotplot(enrichres48, showCategory = 10)

p48.3 <- enrichplot::cnetplot(enrichres48)

ggarrange(
  p48.3,
  ggarrange(p48.1, p48.2, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A") + 
  bgcolor("White") +
  border("White")

file48.1 <- paste0('4_figures/PATHWAY_ENRICHMENT/3_plot_48hrs_', up_or_down, '.png')

ggsave(file48.1,
       height = 40,
       width = 40,
       units = "cm",
       dpi = 500)

enrichres48.2 <- pairwise_termsim(enrichres48)

p48.4 <- emapplot(enrichres48.2)

p48.5 <- upsetplot(enrichres48)

ggarrange(p48.4, p48.5,
          labels = c("A", "B"),
          ncol = 1, 
          nrow = 2,
          common.legend = F, 
          legend = "bottom") + 
  bgcolor("White") +
  border("White")

file48.2 <- paste0('4_figures/PATHWAY_ENRICHMENT/2_plot_48hrs_', up_or_down, '.png')

ggsave(file48.2,
       height = 30,
       width = 30,
       units = "cm",
       dpi = 500)
