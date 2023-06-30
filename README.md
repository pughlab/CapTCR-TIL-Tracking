# Delineation and monitoring of T cell repertoire of adoptive cell transfer product during the treatment of advanced melanoma

# Abstract

Adoptive Cell Transfer (ACT) of in-vitro Tumour-Infiltrating Lymphocytes (TIL) is an effective treatment of solid tumours resulting in objective clinical responses in metastatic melanoma patients. The ability to sequence and track the T-cell repertoire throughout ACT of in vitro TILs provides a method to identify T-cell repertoire features associated with patients’ benefit from ACT. Identification of response biomarkers for patients receiving ACT of TILs has been limited. Conflicting evidence exists for biomarkers such as the quantity of TILs in the infusion product, while other potential biomarkers, such as diversity of the post-infusion peripheral repertoire, have not yet been studied. In this study, we sought to determine 1) the efficacy of using CapTCR-seq to track TILs in serial blood draws over the course of ACT immunotherapy 2) whether ACT response and relapse can be accurately predicted and tracked through serial blood draws. In this study, 9 patients with cutaneous (N = 7) or mucosal (N = 2) melanoma received TIL ACT after chemotherapeutic depletion. Hybrid-capture CapTCR-seq was conducted on apheresis product and pre-/post-transfer peripheral blood mononuclear cells (PBMC) and cell-free DNA. Comparison between PBMC DNA, PBMC RNA, and cfDNA repertoires demonstrated an increased presence of shared T cell clonotypes post-infusion when compared with baseline samples. Higher abundance of TIL clonotypes in the PBMC baseline and post-infusion DNA T-cell repertoires, presence of shared T cell clonotypes between timepoints, and a higher proportion of expanded clonotypes (above change point threshold) present in the PBMC baseline DNA T-cell repertoire was seen in responders when compared with non-responders according to RECIST criteria.  Additionally, presence of clonotypes in the baseline PBMC repertoire and TIL infusion product were predictors of expansion in the post-infusion PBMC repertoire. These results demonstrate effective tracking methodologies and suggest a predictive role for baseline repertoire statistics in response to the ACT of TILs.

---

# Organization and Quick Reproduction

## Running *Reproduce.R*
 - Running *Reproduce.R* reproduces the supplementary data, clone tracking, diversity, and overlap figures. Before running the script, ensure that the required libraries have been installed. Then, go into the repository parent directory and run the script. 

## Required Libraries
 - **bioseq**: 0.1.3
 - **changepoint**: 2.2.4
 - **cowplot**: 1.1.1
 - **changepoint**: 2.2.4
 - **ggalluvial**: 0.12.3
 - **ggplot2**: 3.3.6
 - **ggraph**: 2.0.6
 - **gridExtra**: 2.3
 - **magick**: 2.7.3
 - **openxlsx**: 4.2.5
 - **pheatmap**: 1.0.12
 - **readxl**: 1.4.1
 - **tidyverse**: 1.3.2

## Repo Organization
```bash
code
   |-- Calculation-Functions
   |   |-- changepoint.R
   |   |-- diversity.R
   |   |-- overlap.R
   |   |-- persistence.R
   |   |-- VJUsage_Step1.R
   |   |-- VJUsage_Step2.R
   |-- Data-Load-Functions
   |   |-- Initialize_Data.R
   |   |-- ProductiveCloneFraction.R
   |-- Plot-Functions
   |   |-- AbundanceDiversityOverlay.R
   |   |-- cfDNACloneTrack.R
   |   |-- cfDNACorrel.R
   |   |-- CloneTracking.R
   |   |-- InverseSimpsonDiversity.R
   |   |-- OverlapHeatmap.R
   |   |-- PlotAlignment.R
   |   |-- RelativeAbundance.R
   |   |-- Richnessboxplot.R
   |   |-- SampleCohortCorrelation.R
   |   |-- VJUsage_Step3.R
   |-- Quality-Control
   |   |-- QCPlots.R
   |   |-- UnproductiveClones.R
   |-- Reproduce.R
data
   |-- 10JUN2020-TILS project - clinical data.csv
   |-- CDR3_colors.csv
   |-- correlationmatrix.csv
   |-- mixcr_output.xlsx
   |-- overlapmatrix.csv
   |-- Overlay
   |   |-- clonetracking_overlay.png
   |   |-- diversity_overlay.png
   |   |-- heatmap_overlay.png
   |   |-- overlap_overlay.png
   |   |-- whitespace.png
   |-- Quality-Control
   |   |-- align_stats.csv
   |   |-- align_stats_merged.csv
   |   |-- align_stats_reseq.csv
   |   |-- assemble_stats.csv
   |   |-- assemble_stats_merged.csv
   |   |-- assemble_stats_reseq.csv
   |   |-- Unproductive-Clones
   |   |   |-- clonstats_TRBTLML_cDNAgenomic.csv
   |   |   |-- clonstats_TRBTLML_cfDNAgenomic.csv
   |   |   |-- clonstats_TRBTLML_gDNAgenomic.csv
   |-- SampleKeys.xlsx
   |-- VJTreemaps.csv
   |-- VJUsage.csv
results
   |-- changepoint.png
   |-- Circleplot.png
   |-- Clonetrack.png
   |-- Diversity.png
   |-- Overlap.png
   |-- SupplementaryData.xlsx
```
---

# Running the project

## Data

**1. mixcr_output.xlsx**: MiXCR output with additional 'patient' and 'cohort' columns
cloneno | filename | cloneCount | cloneFraction | allJHitsWithScore | allVHitsWithScore | nSeqCDR3 | aaSeqCDR3 | patient | cohort
 --- | --- | --- |  --- | --- | --- |  --- | --- | --- |  --- 
1 | CLONES_TRB.. | 3673 | 0.127 | TRBJ2-5 | TRBV2 | TGTCCAG.. | CASSPD.. | TLML_1_ | DNA

**2. SampleKeys.xslx**: Information regarding all samples analyzed
Project | Patient_ID | Initials | Timepoint | Sample_Year | Sample_Month | Sample_Day | Sample_Type | Sample_Cohort | Informatics_Name
 --- | --- | --- | --- | --- | --- |  --- | --- | --- | ---
TLML | 1 | AT | baseline | 2013 | 9 | 24 | blood | DNA | TLML-001-AT_A_B_DNA

**3. CDR3_colors.csv**: Desired colored clones for each patient
id | colored_clns | mycolors
 --- | --- | ---
TLML_1_ | CASSPANYEQYF | #1CE6FF

**4. VJUsage.csv**: Output from VJUsage_Step1.R
VJcombo | Number_aaCDR3 | H_aaCDR3 | Hnorm_aaCDR3 | Patient_id | Cycle | Cohort
 --- | --- | --- | --- | --- | --- | --- 
TRBV10 | 1 | 0 | 0 | TLML_1_ | baseline | DNA

**5. VJTreemaps.csv**: Output from VJUsage_Step2.R
cloneCount | cloneFraction | aaSeqCDR3 | nSeqCDR3 | VJcombo | Patient_id | Cycle | Cohort
 --- | --- | --- | --- | --- | --- | --- | ---
3203 | 0.076 | CASSPANY.. | TGTGCCAGC.. | TRBV14.. | TLML_1_ | baseline | DNA

## Loading Data and Organization

### **Initialize_Data.R**: Loads all data and saves into separate dataframes
**Variables Assigned**
 - *"Project" + "Patient_ID" + "Sample_Cohort" + "_samporder"*: List of *"Timepoint" + "_" + "Sample_Year" + "_" + "Sample_Month"* for a specific patient sample cohort
```R
TLML_1_DNA_samporder
[1] "baseline_2013_9" "infusion_2013_8" "FU_01_2014_1"    "FU_02_2014_4"    "FU_03_2014_7"
```
 - *"Project" + "Patient_ID" + "_" + "Sample_Cohort"*: Dataframe of mixcr_output for a specific patient sample cohort 
```R
TLML_1_DNA[1:5,]
filename        aaSeqCDR3 cloneFraction cloneCount
1 baseline_2013_9 CASSPDRGRYQETQYF    0.15566198       3673
2 baseline_2013_9 CATSDSGGLSNQPQHF    0.03979488        939
3 baseline_2013_9   CASRARELNTEAFF    0.03242075        765
4 baseline_2013_9 CASSRFAGGSGNTIYF    0.03085269        728
7 baseline_2013_9 CAGRSGKGAAYNEQFF    0.02212239        522

unique(TLML_1_DNA$filename)
[1] "baseline_2013_9" "FU_02_2014_4"    "infusion_2013_8" "FU_03_2014_7"    "FU_01_2014_1"
```
 - *"Project" + "Patient_ID" + "Sample_Cohort" + "_samporder" + "_" + "Timepoint" + "_" + "Sample_Year" + "_" + "Sample_Month"*: Dataframe of mixcr_output of a timepoint for a specific patient sample cohort
	 - **Note**: Dataframes remove singletons and unproductive clonotypes*
```R
TLML_1_DNA_baseline_2013_9[1:5,]
         filename        aaSeqCDR3 cloneFraction cloneCount
1 baseline_2013_9 CASSPDRGRYQETQYF    0.15566198       3673
2 baseline_2013_9 CATSDSGGLSNQPQHF    0.03979488        939
3 baseline_2013_9   CASRARELNTEAFF    0.03242075        765
4 baseline_2013_9 CASSRFAGGSGNTIYF    0.03085269        728
7 baseline_2013_9 CAGRSGKGAAYNEQFF    0.02212239        522

unique(TLML_1_DNA_baseline_2013_9$filename)
[1] "baseline_2013_9"
```
### *Reproduce.R*: Reproduce.R reproduces SupplementaryData.xlsx, VJUsage_Fig.png, cfDNA_fig.png, and Clone_RelDiv.png
**Reproducing Supplementary Data**
 - The directory *code/Calculation-Functions* consists of functions which run calculations which produce the numbers used in the paper and present in the supplementary materials. After the tables are created, they are exported as separate worksheets in an excel file.

**Reproducing clone tracking figure**
 - The directory *code/Plot-Functions/* contains the *CloneTracking.R* and *cfDNACorrel.R* scripts which produce the clone tracking figures for the DNA and cfDNA repertoire for a specific patient. These plots are used in the *PlotAlignment.R* script which aligns all patients within a certain group. This figure is exported as a png. The *Reproduce.R* script imports these images and overlays them onto the figure overlay image.
 
 **Reproducing diversity figure**
  - The directory *code/Plot-Functions/* contains the *AbundanceDiversityOverlay.R* function which produces the diversity plot for a specific patient and sample cohort. These plots are used in the *PlotAlignment.R* function which aligns all patients within a certain group. This figure is exported as a png. The *Reproduce.R* script imports these images and overlays them onto the figure overlay image.
  
**Reproduing overlap figure**
 - The directory *code/Plot-Functions/* contains the *OverlapHeatmap.R*, *Richnessboxplot.R*, and *cfDNACorrel.R* scripts. *OverlapHeatmap.R* exports a png of the heatmap created from the overlap matrix. *Richnessboxplot.R* exports a png of the boxplot richness comparison for the sample cohorts. *cfDNACorrel.R* produces the cfDNA correlation scatter plot for a specific patient which can be aligned and exported as a png using *PlotAlignment.R*. All exported pngs can be combined using the *Reproduce.R* script which overlays them onto the figure overlay image. 

---

# Results

Clone tracking                      |  Diversity                   | Overlap
:-------------------------:|:-------------------------:|:-------------------------:|
![](/results/Clonetrack.png)  |  ![](/results/Diversity.png) | ![](/results/Overlap.png)
