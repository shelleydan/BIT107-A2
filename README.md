<h1>Time Course Analysis of Differential Gene Expression.</h1>
<h3>Working analysis for differentially expressed genes (DEGs) at individual timepoints and across time using three seperate analytical techniques (TiSA, Moanin, Temporal DESeq2).</h3>

<h2>Introduction</h2>

The analysis provided in this git repository supports the assignments of BIT101 & BIT107 of MSc Big Data Biology with the focus to:

* Download RNAseq metadata from their GEO IDs,
* Analyse DEGs,
* Investigate the temporal nature of DEGs,
* Investigate pathways effected. 

>[!TIP]
> The analysis provided here, and written in the BIT101 Scientific Communication Style Report uses data under the GEOID GSE217504.

<h2>Visualisation of the Workflow</h2>

![workflow](https://github.com/user-attachments/assets/39bd6a87-b149-4c77-8e9e-ba4f5e78228f)

<h2>RMarkdown Method</h2>

A comprehensive RMarkdown can be found in `0_RMarkdown_Summary/Assessment_Summary.Rmd`. To recreate all elements of the BIT101 report submission knit this markdown. Note, it will take approximately 5 minutes to render.

Figure 5 & Supplementary Figure 7 were created using the TiSA method (mentioned below) and have been included in the RMarkdown as images. Manual editing of Figure 5 was needed to arrange similarly to Figure 6. 

<h2>Overview of Scripts</h2>

All scripts can be found within the `1_analysis_scripts` directory, with the exception of `1_analysis_scripts/TiSA_rmarkdown_method/5_TS_analysis.rmd1`.

<table style="width:100%">
  <tr>
    <th style="width:20%">Script ID</th>
    <th style="width:50%">Function</th>
    <th style="width:10%">Language</th>
    <th style="width:20%">Reference</th>
  </tr>
  <tr>
    <td>0_matrix_download.sh </td>
    <td>This is a **bash** script and needs to be run using the command line, this will download the matrix for any given GEOID provided.</td>
    <td>Bash</td>
    <td>N/A</td>
  </tr>
  <tr>
    <td>1_DEG_Analysis.R</td>
    <td>DEGs are identified across time and at each individual timepoint. This script produces figures used for differential analysis.</td>
    <td>R</td>
    <td><a href=https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html>Batut et al., 2016 [1], </a><br>  <a href=https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010752>Hiltemann et al., 2023 [3], </a><br> <a href=https://pubmed.ncbi.nlm.nih.gov/29953864/>Batut et al., 2018 [2],</a></td>
  </tr>
  <tr>
    <td>2_Temporal_Analysis_<br>DESeq2.R</td>
    <td>The first script of the time series analysis, conventional DESeq2 provides basic temporal information for the change of counts over time.</td>
    <td>R</td>
    <td><a href=https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html>Love et al., 2019 [6],</a></td>
  </tr>
  <tr>
    <td> 3_Temporal_Analysis_<br>Moanin.R</td>
    <td>This model takes genes of interest(*) and fits a spline curve to track their change over time.</td>
    <td>R</td>
    <td><a href=https://f1000research.com/articles/9-1447 >Varoquaux & Purdom., 2020 [7], </a></td>
  </tr>
  <tr>
    <td>4_TS_Analysis_Prep.R</td>
    <td>This is a preperation script for the following time series analysis which requires the input data to be properly formatted; metadata needs to be comma seperated, and each sample to have its own count matrix which is tab seperated. </td>
    <td>R</td>
    <td>N/A</td>
  </tr>
  <tr>
    <td>5_TS_analysis.rmd</td>
    <td>This is a sample taken from the larger TiSA pipeline, this was chosen due to its automation. This runs a time series analysis for genes of interest(*).</td>
    <td>Rmd</td>
    <td><a href=https://academic.oup.com/nargab/article/5/1/lqad020/7069286>Lefol et al., 2023 [5],</a></td>
  </tr>
  <tr>
    <td>6.1_Profile_Analysis.R</td>
    <td>Profile analysis makes use of the KEGG, REACTOME and GO:BP pathway enrichment databases and used locally in R, this script is able to accept different settings to use a user specified database.</td>
    <td>R</td>
    <td><a href=https://github.com/biostatsquid>BioStatsSquid [4]</a></td>
  </tr>
  <tr>
    <td>6.2_SUPPLEMENTARY<br>_Profile_Analysis.R</td>
    <td>Same functionality as `6.1_Profile_Analysis.R` but looks specifically at the individual timepoints (4, 12, 48hrs)</td>
    <td>R</td>
    <td><a href=https://github.com/biostatsquid>BioStatsSquid [4]</a></td>
  </tr>
</table>

> [!TIP]
> (*) - The genes of interest used in this analysis have been kept constant between Moanin and TiSA to compare both methods. These genes include: LOC112268133, CLIC4, NDNF, MGLL, CTH, ACE2, ACE, & TMPRSS2. <br>
> <br>
> ACE2, ACE and TMPRSS2 were specifically selected due to their importance in viral entry and the remaining selected at random from a list of the top 10 p-adjusted DEGs.  


<h2>Installation</h2>

This workflow is designed to be able to run as soon as it is downloaded. All of the scripts that need to be run can be found in 1_analysis_scripts and are numbered based on the order of analysis (also highlighted in the visualisation of the workflow above). 

> [!TIP]
> The script to download the metadata (0_matrix_download.sh) is a bash script and will need to be run on the command line.

**Installation** of this workflow is simple, **clone this repository into RStudio**; this workflow was prepared locally with the github extension in R and **not** on [Posit](https://posit.sponsa.bios.cf.ac.uk/).

To successfully do this, copy **https://github.com/shelleydan/BIT107-A2.git** and open a **new project** in RStudio. Select **Version Control** and **Git**, there paste the copied github repository URL. This will clone the repository. 

<h2>Running the Analysis.</h2>

1. Run `0_matrix_download.sh`, you will be prompted to input a GEOID, to follow the analysis in this repository use `GSE217504`.
2. Next, a targets.txt will be generated from the bash script; manually exported into RStudio, saved into `2_rawdata`. Unfortunately, the count matrices all have unique HTMLs and must be downloaded manually.
3. The first script to use is `1_analysis_scripts/1_DEG_Analysis.R` - this will clean the rawdata and apply the DESeq2 Analysis to individual timepoints.
4. The beginning of the temporal analysis starts with `1_analysis_scripts/2_Temporal_Analysis_DESeq2.R` to apply some basic analysis.
5. The Moanin model (`1_analysis_scripts/3_Temporal_Analysis_Moanin.R`) is the secomnd temporal analysis, this will take the genes of intesest and fit a spline to their expression overtime.
6. The final temporal model is the Time Series Analysis, where we use an RMarkdown method (further instructions for this method can be found below - it's a bit different).
7. The next stage of analysis uses a sum of the DEGs at all timepoints and applies a pathway analysis. 

**Time Series Analysis**
1. The RMarkdown to use for this analysis is found here: `1_analysis_scripts/TiSA_rmarkdown_method/5_TS_analysis.rmd` - this pipeline was designed by Lefol et al. (2023)
2. From the initial DESeq2 analysis, we save the cleaned targets and count data to `2_rawdata/1_counts` and `2_rawdata/2_targets`.
3. Firstly, we need to prepare the data for analysis - this is done using `1_analysis_scripts/4_TS_Analysis_Prep.R`. **NB: The counts data needs to be tab seperated for the pipeline to function, do not alted the preperation script and all will be good.**
4. All parameters have been set within `1_analysis_scripts/TiSA_rmarkdown_method/5_TS_analysis.rmd`, including a list of genes used to compare the temporal methods. To run this analysis with the settings for this assingment simply **knit the RMarkdown**. 

**Note:** Performing this analysis on the assignment dataset may take some time - the initial run took ~3hrs, the .html is provided as an example of what the markdown would produce. 

<h2>Bug Reporting</h2>
Please understand that although this is for an assignment, this workflow may be continued to work on in the future - especially to reuse some of the code generated here. If anyone identifies any issues or errors with this analysis please submit an **issue.**

<h2>Known Issues</h2>
Some of the errors that can be produced are the following, but they do not impact the analysis:

*  Some errors producing the volcano plots with ggplot2 will be printed into the console - this is only due to the sheer volume of points to plot - the volcano plots will plot despite this. 
* The Moanin model will express that the fit of the spline is outside of the area - the area is maximised to balance resolution of the other genes, this is unfortunate but avoids sacrificing any other information.
* If you attempt to run `6.2_SUPPLEMENTARY_Profile_Analysis.R` be aware that at 4hrs for downregulated DEGs, and at 12hrs upregulated DEGs there are no significant pathways for REACTOME analysis - therefore this will flag an error and not work. 

If using the RMarkdown method found in `0_RMarkdown_Method`, none of this errors appear and the markdown is rendered. 

<h2>Referances</h2>

1 Batut, B. et al. 2016. Reference-based RNA-Seq data analysis. Available at: https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html [Accessed: 4 December 2024]. <br>
2 Batut, B. et al. 2018. Community-Driven Data Analysis Training for Biology. Cell Systems 6(6), pp. 752-758.e1. doi: 10.1016/j.cels.2018.05.012.<br>
3 Hiltemann, S. et al. 2023. Galaxy Training: A powerful framework for teaching! PLOS Computational Biology 19(1), p. e1010752. doi: 10.1371/journal.pcbi.1010752.<br>
4 Laura, L. BioStatSquid - Overview. Available at: https://github.com/biostatsquid [Accessed: 4 December 2024].<br>
5 Lefol, Y. et al. 2023. TiSA: TimeSeriesAnalysis—a pipeline for the analysis of longitudinal transcriptomics data. NAR Genomics and Bioinformatics 5(1), p. lqad020. doi: 10.1093/nargab/lqad020.<br>
6 Love, M.L., Anders, S., Kim, V. and Huber, W. 2019. RNA-seq workflow: gene-level exploratory analysis and differential expression. BioConductor. Available at: https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html.<br>
7 Varoquaux, N. and Purdom, E. 2020. A pipeline to analyse time-course gene expression data. Available at: https://f1000research.com/articles/9-1447 [Accessed: 4 December 2024].<br>

<h2>RStudio SessionInfo()</h2>

```
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8   
[3] LC_MONETARY=English_United Kingdom.utf8 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.utf8    

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggupset_0.4.0               enrichplot_1.22.0           DOSE_3.28.2                 org.Hs.eg.db_3.18.0        
 [5] clusterProfiler_4.10.1      ggfortify_0.4.17            NMF_0.28                    rngtools_1.5.2             
 [9] registry_0.5-1              BiocWorkflowTools_1.28.0    biomaRt_2.58.2              moanin_1.10.0              
[13] topGO_2.54.0                SparseM_1.84-2              GO.db_3.18.0                AnnotationDbi_1.64.1       
[17] graph_1.80.0                timecoursedata_1.12.0       BiocManager_1.30.25         ggplotify_0.1.2            
[21] plotly_4.10.4               ggdendro_0.2.0              factoextra_1.0.7            cluster_2.1.4              
[25] ggpubr_0.6.0                ggrepel_0.9.6               RColorBrewer_1.1-3          dendextend_1.19.0          
[29] DESeq2_1.42.1               SummarizedExperiment_1.32.0 Biobase_2.62.0              MatrixGenerics_1.14.0      
[33] matrixStats_1.4.1           apeglm_1.24.0               Rsamtools_2.18.0            Biostrings_2.70.3          
[37] XVector_0.42.0              GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0             
[41] S4Vectors_0.40.2            BiocGenerics_0.48.1         lubridate_1.9.3             forcats_1.0.0              
[45] stringr_1.5.1               purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                
[49] tibble_3.2.1                ggplot2_3.5.1               tidyverse_2.0.0             dplyr_1.1.4                
[53] pheatmap_1.0.12            

loaded via a namespace (and not attached):
  [1] splines_4.3.2           bitops_1.0-9            filelock_1.0.3          polyclip_1.10-7        
  [5] XML_3.99-0.17           lifecycle_1.0.4         rstatix_0.7.2           edgeR_4.0.16           
  [9] doParallel_1.0.17       lattice_0.21-9          MASS_7.3-60             backports_1.5.0        
 [13] magrittr_2.0.3          limma_3.58.1            rmarkdown_2.29          yaml_2.3.10            
 [17] cowplot_1.1.3           DBI_1.2.3               abind_1.4-8             zlibbioc_1.48.2        
 [21] ggraph_2.2.1            RCurl_1.98-1.16         yulab.utils_0.1.8       tweenr_2.0.3           
 [25] rappdirs_0.3.3          git2r_0.35.0            GenomeInfoDbData_1.2.11 tidytree_0.4.6         
 [29] codetools_0.2-19        DelayedArray_0.28.0     ggforce_0.4.2           xml2_1.3.6             
 [33] tidyselect_1.2.1        aplot_0.2.3             farver_2.1.2            gmp_0.7-5              
 [37] viridis_0.6.5           BiocFileCache_2.10.2    jsonlite_1.8.9          tidygraph_1.3.1        
 [41] Formula_1.2-5           iterators_1.0.14        bbmle_1.0.25.1          foreach_1.5.2          
 [45] tools_4.3.2             progress_1.2.3          treeio_1.26.0           Rcpp_1.0.13-1          
 [49] glue_1.8.0              gridExtra_2.3           SparseArray_1.2.4       xfun_0.49              
 [53] qvalue_2.34.0           usethis_3.0.0           withr_3.0.2             numDeriv_2016.8-1.1    
 [57] fastmap_1.2.0           fansi_1.0.6             digest_0.6.37           timechange_0.3.0       
 [61] R6_2.5.1                gridGraphics_0.5-1      colorspace_2.1-1        RSQLite_2.3.7          
 [65] utf8_1.2.4              generics_0.1.3          data.table_1.16.2       graphlayouts_1.2.0     
 [69] prettyunits_1.2.0       httr_1.4.7              htmlwidgets_1.6.4       S4Arrays_1.2.1         
 [73] scatterpie_0.2.4        pkgconfig_2.0.3         gtable_0.3.6            blob_1.2.4             
 [77] shadowtext_0.1.4        htmltools_0.5.8.1       carData_3.0-5           fgsea_1.28.0           
 [81] bookdown_0.41           scales_1.3.0            ClusterR_1.3.3          png_0.1-8              
 [85] ggfun_0.1.7             knitr_1.49              rstudioapi_0.17.1       tzdb_0.4.0             
 [89] reshape2_1.4.4          nlme_3.1-163            coda_0.19-4.1           curl_6.0.0             
 [93] bdsmatrix_1.3-7         cachem_1.1.0            parallel_4.3.2          HDO.db_0.99.1          
 [97] pillar_1.9.0            grid_4.3.2              vctrs_0.6.5             car_3.1-3              
[101] dbplyr_2.5.0            evaluate_1.0.1          mvtnorm_1.3-2           cli_3.6.3              
[105] locfit_1.5-9.10         compiler_4.3.2          rlang_1.1.4             crayon_1.5.3           
[109] ggsignif_0.6.4          emdbook_1.3.13          plyr_1.8.9              fs_1.6.5               
[113] stringi_1.8.4           viridisLite_0.4.2       gridBase_0.4-7          BiocParallel_1.36.0    
[117] munsell_0.5.1           lazyeval_0.2.2          GOSemSim_2.28.1         Matrix_1.6-1.1         
[121] patchwork_1.3.0         hms_1.1.3               bit64_4.5.2             KEGGREST_1.42.0        
[125] statmod_1.5.0           igraph_2.1.1            broom_1.0.7             memoise_2.0.1          
[129] ggtree_3.10.1           fastmatch_1.1-4         bit_4.5.0               gson_0.1.0             
[133] ape_5.8

```
