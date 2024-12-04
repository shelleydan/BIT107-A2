
## THIS IS FOR SUPPLEMENTARY ONLY, USE 6_PROFILE_ANALYSIS.R ##

######################### KEGG Profile Analysis of DEGs ########################

## Reference -------------------------------------------------------------------

# Biostatsquid.com [Accessed: 02/12/2024]

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

## Loading in our data (4HRS) --------------------------------------------------

DEGTEMP <- read.csv("3_output_data/raw_DEG_results_TEMPORAL.csv", 
                    header = T)
names(DEGTEMP)[names(DEGTEMP) == 'X'] <- 'gene_symbol'
DEGTEMPdf <- as.data.frame(cbind(DEGTEMP$gene_symbol, 
                                 DEGTEMP$pvalue, 
                                 DEGTEMP$padj, 
                                 DEGTEMP$log2FoldChange, 
                                 DEGTEMP$diffexpressed))
names(DEGTEMPdf)[names(DEGTEMPdf) == 'V1'] <- 'gene_symbol'
names(DEGTEMPdf)[names(DEGTEMPdf) == 'V2'] <- 'pval'
names(DEGTEMPdf)[names(DEGTEMPdf) == 'V3'] <- 'padj'
names(DEGTEMPdf)[names(DEGTEMPdf) == 'V4'] <- 'log2fc'
names(DEGTEMPdf)[names(DEGTEMPdf) == 'V5'] <- 'diffexpressed'

DEGTEMPdf <- DEGTEMPdf[DEGTEMPdf$diffexpressed == 'UP'| DEGTEMPdf$diffexpressed == 'DOWN',]
deg_results_list_TEMP <- split(DEGTEMPdf, DEGTEMPdf$diffexpressed)

## Aquiring Reference Genes ---------------------------------------------------

# Gene sets can be downloaded here: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# Gene Sets used in this analysis are found in '5_Gene_Sets'
# These need to be the whole number of genes used in the experiement - the 46863 in this assingment
DEG_Genes <- DEGTEMP$gene_symbol #All time points have the same total genes

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
rm(file)
rm(filename_gs)
rm(gmt_files)
rm(pwl2)

## Filtering Gene Sets (REACTOME, KEGG, GO) ------------------------------------
# We look for profiles enriched in both our up regulated and down regulated. 

################################### SETTINGS ################################### 
name_of_comparison <- 'mock_vs_infected' # for our filename

background_genes <- 'GO' # Choose GO, KEGG or REACTOME

bg_genes <- readRDS(paste0(gene_path, 'reactome.RDS')) # background genes.

up_or_down <- 'UP' # Set this to UP or DOWN

padj_cutoff <- 0.05 # p-adjusted threshold.

genecount_cutoff <- 5 # minimum number of genes in the pathway.

filename <- paste0(output, 'clusterProfiler/', name_of_comparison, '_', background_genes) # filename of our PEA results

if(background_genes == 'KEGG'){
  bg_genes <- readRDS(paste0(gene_path, 'kegg_legacy.RDS'))
} else if(background_genes == 'REACTOME'){
  bg_genes <- readRDS(paste0(gene_path, 'reactome.RDS'))
} else if(background_genes == 'GO'){
  bg_genes <- readRDS(paste0(gene_path, 'go.bp.RDS'))
} else {
  stop('Invalid background genes. Select one of the following: KEGG, Reactome, GO, or add new pwl to function')
}

## Performing Cluster Analysis (AUTO) ------------------------------------------

resTEMP <- lapply(names(deg_results_list_TEMP),
                  function(x) enricher(gene = deg_results_list_TEMP[[x]]$gene_symbol,
                                       TERM2GENE = bg_genes))

names(resTEMP) <- names(deg_results_list_TEMP) # Apply the UP and DOWN tags

res_dfTEMP <- lapply(names(resTEMP), function(x) rbind(resTEMP[[x]]@result))
names(res_dfTEMP) <- names(resTEMP) # Apply the UP and DOWN tags
res_dfTEMP <- do.call(rbind, res_dfTEMP)

res_dfTEMP <- res_dfTEMP %>% mutate(minuslog10padj = -log10(p.adjust),
                                    diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', 
                                                         rownames(res_dfTEMP)))

# Subset the pathways by (i) padj value, (ii) Gene count
target_pws <- unique(res_dfTEMP$ID[res_dfTEMP$p.adjust < padj_cutoff & res_dfTEMP$Count > genecount_cutoff])
res_dfTEMP <- res_dfTEMP[res_dfTEMP$ID %in% target_pws,]

# HERE FILTER BASED ON UP OR DOWN REGULATED AS SPECIFIED IN SETTINGS

res_dfTEMP <- res_dfTEMP %>% filter(diffexpressed == up_or_down) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_dfTEMP) <- res_dfTEMP$ID

# For visualisation purposes, let's shorten the pathway names
res_dfTEMP$Description <- gsub('(H|h)iv', 'HIV', 
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
                                                                                                                 gsub('GOBP_|KEGG_|REACTOME_', '', res_dfTEMP$Description)))))))))))))))))))

## ENRICHMENT ANALYSIS

enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = res_dfTEMP,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 organism = "human",
                 ontology = "UNKNOWN",
                 gene = DEGTEMP$gene_symbol,
                 keytype = "UNKNOWN",
                 universe = unique(bg_genes$gene),
                 gene2Symbol = character(0),
                 geneSets = bg_genes)
class(enrichres) # CHECK THIS HAS BECOME AN ENRICHMENT CLASS
?new

## VISUALISATION ---------------------------------------------------------------

# This process will automatically produce figures (arranges) for all options above

## 4HRS PLOTS ------------------------------------------------------------------

#Do not run if looking at downregulated reactome
p4.1 <- barplot(enrichres, showCategory = 20)

p4.2 <- dotplot(enrichres, showCategory = 15)

p4.3 <- enrichplot::cnetplot(enrichres)

ggarrange(
  p4.3,
  ggarrange(p4.1, p4.2, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A") + 
  bgcolor("White") +
  border("White")

file4.1 <- paste0('4_figures/Figure_X_3_plot_TEMP_', up_or_down,'_',background_genes, '.png')

ggsave(file4.1,
       height = 55,
       width = 35,
       units = "cm",
       dpi = 500)


# The following two plots are good for overview, but due to the number of labels it makes it impossible to read. Therefore, these will not be included. 
# 
# enrichres2 <- pairwise_termsim(enrichres)
#
# p4.4 <- emapplot(enrichres2)
# 
# p4.5 <- upsetplot(enrichres)
# 
# ggarrange(p4.4, p4.5,
#           labels = c("A", "B"),
#           ncol = 1, 
#           nrow = 2,
#           common.legend = F, 
#           legend = "bottom") + 
#   bgcolor("White") +
#   border("White")
# 
# file4.2 <- paste0('4_figures/Supplementary/Figure_X_2_plot_TEMP_', up_or_down,'_',background_genes, '.png')
# 
# ggsave(file4.2,
#        height = 30,
#        width = 30,
#        units = "cm",
#        dpi = 500)
