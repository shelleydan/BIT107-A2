<h1>Time Course Analysis of Differential Gene Expression.</h1>
<h3>Working analysis for differentially expressed genes (DEGs) at individual timepoints and across time using three seperate analytical techniques (TiSA, Moanin, Temporal DESeq2).</h3>

<h3 style="color: blue;">Introduction</h3>

The analysis provided in this git repository supports the assignments of BIT101 & BIT107 of MSc Big Data Biology with the focus to:

* Download RNAseq count matricies and metadata from their GEO IDs,
* Analyse DEGs,
* Investigate the temporal nature of DEGs,
* Investigate pathways and molecular systems modified. 

The analysis provided here, and written in the BIT101 Scientific Communication Style Report uses data under the $${\color{#0969DA}GEOID \space GSE217504.}$$

<h3>Visualisation of the Workflow</h3>

![workflow](https://github.com/user-attachments/assets/31711a91-9853-4fc1-9132-7c1dcf376f66)

<h3>Installation & Running Instructions</h3>

This workflow is designed to be able to run as soon as it is downloaded. All of the scripts that need to be run can be found in 1_analysis_scripts and are numbered based on the order of analysis (also highlighted in the visualisation of the workflow above). **NB:** the script to download the counts matrix and meta data (0_matrix_download.sh) is a bash script and will need to be run on the command line.

**Installation** of this workflow is simple, **clone this repository into RStudio**; this workflow was prepared locally with the github extension in R and **not** on [Posit](https://posit.sponsa.bios.cf.ac.uk/).

**Running the Analysis.**
1. Run 0_matrix_download.sh, you will be prompted to input a GEOID, to follow the analysis in this repository use `GSE217504`.
2. Next, the counts matrix and the targets.txt generated from the bash script can be manually exported and inported into RStudio, saved into `2_rawdata`.
3. The first script to use is `1_analysis_scripts/1_DEG_Analysis.R` - this will clean the rawdata and apply the DESeq2 Analysis to individual timepoints.
4. The beginning of the temporal analysis starts with `1_analysis_scripts/2_Temporal_Analysis_DESeq2.R` to apply some basic analysis.
5. The Moanin model (`1_analysis_scripts/3_Temporal_Analysis_Moanin.R`) is the secomnd temporal analysis, this will take the genes of intesest and fit a spline to their expression overtime.
6. The final temporal model is the Time Series Analysis, where we use an RMarkdown method (further instructions for this method can be found below - it's a bit different).

**Time Series Analysis**
1. The RMarkdown to use for this analysis is found here: `1_analysis_scripts/TiSA_rmarkdown_method/5_TS_analysis.rmd` - this pipeline was designed by Lefol et al. (2023)
2. From the initial DESeq2 analysis, we save the cleaned targets and count data to `2_rawdata/1_counts` and `2_rawdata/2_targets`.
3. Firstly, we need to prepare the data for analysis - this is done using `1_analysis_scripts/4_TS_Analysis_Prep.R`. **NB: The counts data needs to be tab seperated for the pipeline to function, do not alted the preperation script and all will be good.**
4. All parameters have been set withing `1_analysis_scripts/TiSA_rmarkdown_method/5_TS_analysis.rmd`, including a list of genes used to compare the temporal methods. To run this analysis with the settings for this assingment simply **knit the RMarkdown**. 

**Note:** Performing this analysis on the assignment dataset may take some time - the initial run took ~3hrs, the .html is provided as an example of what the markdown would produce. 

<h3>Bug Reporting</h3>
Please understand that although this is for an assignment, this workflow may be continued to work on in the future - especially to reuse some of the code generated here. If anyone identified any issues or errors with this analysis please submit an **issues report**.

<h3>Known Issues</h3>
Some of the errors that can be produced are the following, but they do not impact the analysis:

*  Some errors producing the volcano plots with ggplot2 will be printed into the console - this is only due to the sheer volume of points to plot - the volcano plots will plot despite this. 
* The Moanin model will express that the fit of the spline is outside of the area - the area is maximised to balance resolution of the other genes, this is unfortunate but avoids sacrificing any other information. 

<h3>Referances</h3>
This flow of RNAseq analysis was made using the following sources:

**_This section will be updated once the assignment is complete_**
