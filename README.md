# Tracking and Predicting Effectiveness of Adoptive Cell Transfer using T-cell Receptor Capture & Sequencing

## Contents
- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Demo](#demo)
- [Results](#results)

# Overview

Adoptive Cell Transfer (ACT) of in-vitro Tumor-Infiltrating Lymphocytes (TIL) is an effective treatment of solid tumours commonly resulting in sustainable responses with an objective clinical response rate of around 50% throughout multiple studies (Rosenberg & Restifo, 2015). The ability to sequence and track the T Cell Receptor (TCR) repertoire throughout ACT of in vitro TILs provides an enhanced analysis of a patient’s response. Previous studies have shown that an increased number of TILs correlates with patient response (Valpione et al., 2020), yet the diversity of the TCRs present in the ACT TIL-infusion product has been seen to have no impact on patient response (Mullinax et al., 2018). In this study, we sought to determine 1) whether T cell receptor (TCR) repertoire diversity of the ACT infusion product and pre-infusion baseline sample is associated with response and 2) whether ACT response and relapse can be accurately predicted and tracked through serial blood draws. In this study, 10 patients with cutaneous melanoma, mucosal melanoma, breast, or ovarian cancer were each given TIL ACT after chemotherapeutic depletion. In addition to the apheresis product, peripheral blood cells and cell-free DNA were collected before and after the transfer. The post-infusion samples were taken every 27 - 106 days (average 56 days) after the infusion. TCR repertoire was isolated using the hybrid-capture CapTCR-seq assay followed by comparison of TCR diversity and clonal tracking for each patient. Further analysis revealed that an increased clonality with a reduced richness and eveness of the baseline sample were present in higher response patients. Additionally, the expansive clones of the sample taken 4 weeks after the infusion were seen to originate from low frequency clonotypes present in the TIL infusion. These results demonstrate effective tracking methodologies for TIL ACT as well as a selective hyper-expansion of TIL clones in vivo associated with response to ACT.

# Repo Contents
This repo contains R scripts and figures used in the analysis for the CapTCR TIL Tracking project.

## Data Load and Organization
 - **DataLoad.R**: Loads specific chain, sample cohort, and clone fraction for specific patient's samples
 - **Organization.R**: Loads sample order for patients, required libraries, and a maximally distant color pallete
 - **ProductiveCloneFraction.R**: Readjusts clone fraction data based on spliced data frame
## Quality Control
 - **QCPlots.R**: Uses tabular log dataset from MiXCR to create and bind multiple quality control plots together showing different data points
 - **UnproductiveClones.R**: Analyzes the unproductive data consumption of the patient's clones and outputs excel sheet
## Plot Functions
- **CloneTracking.R**: Plots and tracks the redundant TIL clones through patient's timeline
- **InverseSimpsonDiversity.R**: Plots and tracks the inverse simpson diversity through patient's timeline
- **RelativeAbundance.R**: Plots and tracks the relative abundance of large and hyperexpansive clones through patient's timeline
- **PlotAlignment.R**: Aligns a group of patient's plots on a single x-axis
- **SampleCohortCorrelation.R**: Plots the correlative overlap between sample cohorts (gDNA, cDNA, cfDNA)
## Calculation Functions
 - **AlternateTILClones.R**: Calculating the amount of a sample's repertoire is taken up by TIL-related clones
 - **BaselineTILs.R**: Calculating the number of TIL-related clones in the baseline repertoire
 - **DeltaAnalysis.R**: Calculates the delta between 2 data points of different samples
 - **DiversityEval.R**: Calculates a patient's inverse simpson diversity throughout their timeline
 - **TILCharacteristicAnalysis.R**: Analyzes the amount of clones >5% in the TIL sample, total number of TCRs and the total amount of TIL clones
 - **Top10Comparison.R**: Calculates the clone fraction of the top 4 week clones in the baseline and TIL-infusion data
 - **Top50Analysis.R**: Finds the amount of a sample taken up by the top 50% of another sample

# System Requirements
## Hardware Requirements
All scripts in this repository were developed on a system with the following specs:

RAM: 8 GB
CPU: 1.7 GHz Intel Core i7
MAC OSX: 10.13.6

## Software Requirements
General software: R (version 3.6.1)

### R-packages
Visualization tools:
 - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
 - [ggalluvial](https://cran.r-project.org/web/packages/ggalluvial/index.html)
 - [cowplot](https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html#:~:text=The%20cowplot%20package%20is%20a,or%20mix%20plots%20with%20images.)

Analysis tools:
 - [immunarch](https://github.com/immunomind/immunarch)
 - [readxl](https://readxl.tidyverse.org/)

# Demo

## Data Load and Organization

Loading set-up functions

``` R
dir_code <- '/Users/JohnDoe/R/'
source(paste(dir_code, 'Data-Load-Functions/DataLoad.R', sep=""))
source(paste(dir_code, 'Data-Load-Functions/ProductiveCloneFraction.R', sep=""))
source(paste(dir_code, 'Data-Load-Functions/Organization.R', sep=""))
```
Example of loading data (Pre-installation of packages required)

``` R
dir_data <- '/Users/JohnDoe/Data/'
file_samplekeys <- 'TLML_samples_keys.xlsx'
Load_data("TLML_4_", "gDNA", "TRB", 0, dir_data, dir_data, file_samplekeys)
print(CDR3_fraction[1:10,])
```

## Quality Control

Using QC functions to print out QC plots and unproductive clones excel sheet

``` R
# Define quality control data variables and load in functions
dir_QCdata <- '/Users/JohnDoe/Data/Quality-Control/'
dir_output <- '/Users/JohnDoe/Output/'
source(paste(dir_code, 'Quality-Control/QCPlots.R', sep=""))
source(paste(dir_code, 'Quality-Control/UnproductiveClones.R, sep=""))

# Use the QC plots function and the unproductive clones analysis functions to print out plots and excel sheet
mixcrQC.fx('align_stats.csv', 'assemble_stats.csv', 'gDNA' , 'QC_Output_gDNA.png', dir_QCdata, dir_output)
dir_clones <- paste(dir_data, 'gDNA/', sep="")
clonestatfx('TRB', dir_clones, 'UnproductiveClones', dir_output)
```

## Plot Functions

Using the 3 plot functions

``` R
# Load in required functions for plot creation
source(paste(dir_code, 'Plot-Functions/CloneTracking.R', sep=""))
source(paste(dir_code, 'Plot-Functions/InverseSimpsonDiversity.R', sep=""))
source(paste(dir_code, 'Plot-Functions/RelativeAbundance.R', sep=""))

# Use functions to create 3 plots based on the data from gDNA and TRB from TLML_4_
ClonetrackPlot("TLML_4_", "gDNA", "TRB", 0, dir_data, dir_data, file_samplekeys, "Baseline")
print(myp)
DivPlot("TLML_4_", "gDNA", "TRB", 0, dir_data, dir_data, file_samplekeys, "Baseline", 500)
print(myp)
RelPlot("TLML_4_", "gDNA", "TRB", 0, dir_data, dir_data, file_samplekeys, "Baseline")
print(myp)
```
Printing out x-axis aligned plots

``` R
# Load in function for plot alignment
source(paste(dir_code, 'Plot-Functions/PlotAlignment.R', sep=""))
dir_output <- '/Users/JohnDoe/Output'
file_output <- 'HighCloneTrack'

# Use function to print out aligned version of the clone tracking plot for gDNA and TRB from the high expansion patients group
alignment_fig(high, "gDNA", "TRB", 0, dir_data, dir_data, file_samplekeys, clonetrack, "Baseline", dir_output, file_output)
```

## Calculation Functions

Example: Calculating Diversity Delta between baseline and the 4 week sample

```R
# Load the delta analysis function
source(paste(dir_code, 'Calculation-Functions/DeltaAnalysis.R', sep=""))
# Define the 2 desired samples for analysis
FW_sample <- TLML_4_gDNA[3]
BL_sample <- TLML_4_gDNA[1]

# Use function to print out the diversity delta between the 2 samples
DeltaA("TLML_4_", "gDNA", "TRB", 0, "Diversity", BL_sample, FW_sample, dir_data, dir_data, file_samplekeys)
```

# Results

## The level of expansion to the adoptive cell transfer of TILs can be tracked longitudinally using CapTCR with high accuracy
To demonstrate the ability to use CapTCR-seq for genomic DNA, RNA, and circulating free DNA, we analyzed the prevalence of all CDR3 sequences in the 4 week sample.
 - Genomic DNA correlated with RNA with an R2 value of 0.796 and cfDNA correlated with gDNA with an R2 value of 0.585. 

To demonstrate the effective monitoring of TIL ACT in peripheral blood samples, we tracked TIL-specific clonotypes over the patient timeline.
 -  The percentage of the sample 4 weeks after the infusion taken up by these TIL-related clones was calculated (median: 30.67%; range:55.31%).
    - High expansion patients had a median expansion of 45.77% (range: 22.66%)
    - Medium expanders with a median of 17.05% (range: 12.50%)
    - Low expansion patients had a median of 3.74% (range: 6.44%)
 - There was no overall correlation seen between the expansion of the patient’s TILs in the 4 week sample and their lesion evaluation. 
    - The 3 PR patients had the highest expansion as 51.05% of the 4 week sample was taken up by TIL-clones (range: 15.49%)
    - The 4 SD patients had a median expansion of 22.61% (range: 32.65%)
    - The 2 PD patients had a median expansion of 30.71% (range: 30.11%).

## Increased amount hyper-expansive clones post-infusion demonstrates high expansion and response
To demonstrate the effect of the clonality of the baseline, TIL infusion, and post-infusion samples on expansion, we analyzed the diversity and clonal structures of the patients.
 - A higher amount of the 4 week sample taken up by hyperexpanded (>5% of sample) clones in high expansion patients
     - High expansion patients showed a median of 16.62% (range: 31.61%)
     - Medium expansion patients showed a median of 15.27% (range: 19.75%)
     - Low expansion patients showed a median of 9.95% (range: 19.91%)
 - Higher response patients had a higher amount of the 4 week sample taken up by hyperexpanded clones
     - PR patients had a median of 16.62% (range: 7.55%) 
     - SD patients had a median of 15.31% (range: 37.98%)
     - PD patients had a median of 10.82% (range: 8.90%).

## Decreased Diversity of the TIL infusion product was seen in progressive response patients
To demonstrate the predictive role of the diversity of the TIL infusion product, the diversity and total clonotype metrics were analyzed for the TIL infusion product.
 - A decreased diversity of the TILs was seen in higher response patients.
     - PR patients had a median of 9.51 (range: 121.89)
     - SD patients had a median of 40.74 (range: 443.99)
     - PD patients had a median of 113.33 (range: 144.435)
 - The amount of hyperexpanded clones also positively correlated with the response of the patient. 
     - PR patients had 54.77% of (range: 72.27%)
     - SD patients had a median of 19.80% (range: 35.97%)
     - PD patients had a median of 10.94% (range: 21.87%).

## Small clonotypes undergo selected-expansion of TIL-infusion clones post-infusion
To demonstrate the selected-expansion of TIL-related clones we used the plots tracking TIL-specific clonotypes clone fraction to view specific clonal expansion. 
 - The top 50% of the TIL-infusion product took up a median of 3.28% (range: 32.66%) of the repertoire 4 weeks after the infusion
     - High expansion patients showed a median of 4.61% (range: 31.59%) 
     - Medium expansion patients showed a median of 2.53% (range: 11.22%) 
     - Low expansion patients both showed a median of 0.02% (range: 0.05%)
 - The top 50% of the 4 week sample took up a median of 0.51% (range: 79.50% of the TIL-infusion product
     - High expansion patients showed a median of 21.83% (range: 79.13%)
     - Medium expansion patients showed a median of 0.17% (range: 9.77%)
     - Low expansion patients showed a median of 0.03% (range: 0.06%).
To demonstrate the size of expanded clonotypes, the average percentage of the TIL infusion taken up by each of the top 10 clones in the 4 week sample was analyzed.
 - The top 10 clonotypes in the 4 week sample had a median clone percentage of 0.22% in the TIL infusion (range: 7.66%)
     - High expansion patients showed a higher average clone percentage with a median of 1.19% (range: 2.54%)
     - Medium expansion patients had a median of 0.02% (range: 0.02%) 
     - Low expansion patients had a median of 0.09% (range: 0.82%)

## HLA type HLA-A-301 shows decreased response to TIL therapy
To demonstrate the role of HLA types on TIL therapy response characteristics, HLA-type groups were analyzed for their different TIL characteristics and response rates.
 - 2 out of the 10 patients had a progressive disease lesion evaluation and those 2 patients were the only patients having the HLA type HLA-A301. One other patient showed the HLA-A-301 phenotype yet no lesion evaluation data was given. However, this patient was in our low TIL expansion group. 
