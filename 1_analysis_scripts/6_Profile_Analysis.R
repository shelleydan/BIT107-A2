
######################### KEGG Profile Analysis of DEGs ########################

## Reference -------------------------------------------------------------------

# https://github.com/ACSoupir/Bioinformatics_YouTube.git

## Library Prep ----------------------------------------------------------------
BiocManager::install("pathview")
BiocManager::install("gage")
BiocManager::install("gageData")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("biomartr")
library(pathview)
library(gage)
library(gageData)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomartr)

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

DEG4_Sig <- DEG4[DEG4$diffexpressed == 'UP' | DEG4$diffexpressed == 'DOWN',]
DEG12_Sig <- DEG12[DEG12$diffexpressed == 'UP' | DEG12$diffexpressed == 'DOWN',]
DEG48_Sig <- DEG48[DEG48$diffexpressed == 'UP' | DEG48$diffexpressed == 'DOWN',]

rm(DEG4, DEG12, DEG48) # Remove un-filtered data frame

## Assigning ENSAMBL IDs -------------------------------------------------------

symbols4 <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = row.names(DEG4_Sig),
                                  columns = c('ENTREZID'), 
                                  keytype = 'SYMBOL')

symbols12 <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = row.names(DEG12_Sig),
                                  columns = c('ENTREZID'), 
                                  keytype = 'SYMBOL')

symbols48 <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = row.names(DEG48_Sig),
                                  columns = c('ENTREZID'), 
                                  keytype = 'SYMBOL')

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





