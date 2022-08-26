# Delineation and monitoring of T cell repertoire of adoptive cell transfer product during the treatment of advanced melanoma

# Abstract

Adoptive Cell Transfer (ACT) of in-vitro Tumor-Infiltrating Lymphocytes (TIL) is an effective treatment of solid tumours commonly resulting in sustainable responses with an objective clinical response rate of around 50% throughout multiple studies (Rosenberg & Restifo, 2015). The ability to sequence and track the T Cell Receptor (TCR) repertoire throughout ACT of in vitro TILs provides an enhanced analysis of a patient’s response. Previous studies have shown that an increased number of TILs correlates with patient response (Valpione et al., 2020), yet the diversity of the TCRs present in the ACT TIL-infusion product has been seen to have no impact on patient response (Mullinax et al., 2018). In this study, we sought to determine 1) whether T cell receptor (TCR) repertoire diversity of the ACT infusion product and pre-infusion baseline sample is associated with response and 2) whether ACT response and relapse can be accurately predicted and tracked through serial blood draws. In this study, 10 patients with cutaneous melanoma, mucosal melanoma, breast, or ovarian cancer were each given TIL ACT after chemotherapeutic depletion. In addition to the apheresis product, peripheral blood cells and cell-free DNA were collected before and after the transfer. The post-infusion samples were taken every 27 - 106 days (average 56 days) after the infusion. TCR repertoire was isolated using the hybrid-capture CapTCR-seq assay followed by comparison of TCR diversity and clonal tracking for each patient. Further analysis revealed that an increased clonality with a reduced richness and eveness of the baseline sample were present in higher response patients. Additionally, the expansive clones of the sample taken 4 weeks after the infusion were seen to originate from low frequency clonotypes present in the TIL infusion. These results demonstrate effective tracking methodologies for TIL ACT as well as a selective hyper-expansion of TIL clones in vivo associated with response to ACT.

---

# Organization and Quick Reproduction

## Running *Reproduce.R*
 - Running *Reproduce.R* reproduces the supplementary data, cfDNA figure, clone tracking / diversity figure, and the cfDNA figure. Before running the script, ensure that the required libraries have been installed. Then, go into the repository parent directory and run the script. 
 - Alternatively, it can be run on a cloud machine through [our Code Ocean repository] (https://codeocean.com/)

## Required Libraries
 - **bioseq**: 0.1.3
 - **cowplot**: 1.1.1
 - **ggalluvial**: 0.12.3
 - **ggplot2**: 3.3.6
 - **ggraph**: 2.0.6
 - **gridExtra**: 2.3
 - **magick**: 2.7.3
 - **openxlsx**: 4.2.5
 - **readxl**: 1.4.1
 - **tidyverse**: 1.3.2

## Repo Organization
```bash
code
   |-- Calculation-Functions
   |   |-- AlternateTILClones.R
   |   |-- BaselinePosition.R
   |   |-- BaselineTILs.R
   |   |-- CDR3perVJ.R
   |   |-- CDR3vsVJ_Expansion.R
   |   |-- Clonefraction_Correl.R
   |   |-- DeltaAnalysis.R
   |   |-- DiversityEval.R
   |   |-- DominantClones.R
   |   |-- Richness&EvennessCalculations.R
   |   |-- TCRConvergence.R
   |   |-- TimepointCharacteristic.R
   |   |-- Top10Comparison.R
   |   |-- Top50Analysis.R
   |   |-- VJUsage_Step1.R
   |   |-- VJUsage_Step2.R
   |   |-- cfDNACalculations.R
   |-- Data-Load-Functions
   |   |-- Initialize_Data.R
   |   |-- ProductiveCloneFraction.R
   |-- Plot-Functions
   |   |-- AbundanceDiversityOverlay.R
   |   |-- CloneTracking.R
   |   |-- InverseSimpsonDiversity.R
   |   |-- PlotAlignment.R
   |   |-- RelativeAbundance.R
   |   |-- SampleCohortCorrelation.R
   |   |-- VJUsage_Step3.R
   |   |-- cfDNACloneTrack.R
   |   |-- cfDNACorrel.R
   |-- Quality-Control
   |   |-- QCPlots.R
   |   |-- UnproductiveClones.R
   |-- Reproduce.R
data
   |-- 10JUN2020-TILS project - clinical data.csv
   |-- CDR3_colors.csv
   |-- Clone_RelDivFigOverlay.png
   |-- Quality-Control
   |   |-- Unproductive-Clones
   |   |   |-- clonstats_TRBTLML_cDNAgenomic.csv
   |   |   |-- clonstats_TRBTLML_cfDNAgenomic.csv
   |   |   |-- clonstats_TRBTLML_gDNAgenomic.csv
   |   |-- align_stats.csv
   |   |-- align_stats_merged.csv
   |   |-- align_stats_reseq.csv
   |   |-- assemble_stats.csv
   |   |-- assemble_stats_merged.csv
   |   |-- assemble_stats_reseq.csv
   |-- SampleKeys.xlsx
   |-- VJFigOverlay.png
   |-- VJTreemaps.csv
   |-- VJUsage.csv
   |-- cfDNAFigOverlay.png
   |-- mixcr_output.xlsx
result
   |-- Clone_RelDiv.png
   |-- Overview.png
   |-- Quality-Control
   |   |-- TLML1TLML26_QCDownsampling.png
   |   |-- TLML_Regular_cDNA_QC.png
   |   |-- TLML_Regular_cfDNA_QC.png
   |   |-- TLML_Regular_gDNA_QC.png
   |-- SampleCohortsCorrelation.png
   |-- SupplementaryData.xlsx
   |-- VJUsage_Fig.png
   |-- cfDNA_fig.png
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
 - *{Project}{Patient_ID}_{Sample_Cohort}_samporder*: List of *{Timepoint}_{Sample_Year}_{Sample_Month}* for a specific patient sample cohort
```R
TLML_1_DNA_samporder
[1] "baseline_2013_9" "infusion_2013_8" "FU_01_2014_1"    "FU_02_2014_4"    "FU_03_2014_7"
```
 - *{Project}{Patient_ID}_{Sample_Cohort}: Dataframe of mixcr_output for a specific patient sample cohort 
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
 - *{Project}{Patient_ID}_{Sample_Cohort}_{Timepoint}_{Sample_Year}_{Sample_Month}*: Dataframe of mixcr_output of a timepoint for a specific patient sample cohort
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
### **Reproduce.R**: Reproduce.R reproduces SupplementaryData.xlsx, VJUsage_Fig.png, cfDNA_fig.png, and Clone_RelDiv.png
**Reproducing Supplementary Data**
 - The directory *code/Calculation-Functions* consists of functions which run calculations which produce the numbers used in the paper and present in the supplementary materials. After the tables are created, they are exported as separate worksheets in an excel file.

**Reproducing clone tracking and diversity figure**
 - The directory *code/Plot-Functions* contains the *CloneTracking.R* function as well as the *AbundanceDiversityOverlay.R* function which create both plots for a specific patient sample cohort. These plots are used in the *PlotAlignment.R* function which aligns all patients within a certain group. This figure is exported as a png. The *Reproduce.R* script imports these images and overlays them onto the figure overlay image.

**Reproducing VJ usage**
 - The directory *code/Plot-Functions* contains the *VJUsage_Step3.R* function which creates the VJ circle plot for a specific patient sample cohort and exports it as a png. Th *Reproduce.R* script imports these images and overlays them onto the figure overlay image. 

**Reproducing cfDNA figure**
 - The directory *code/Plot-Functions* contains the *cfDNACorrel.R* function as well as the *cfDNACloneTrack.R* which create both plots for a specific patient sample cohort. These plots are used in the *PlotAlignment.R* function which aligns all patients within a certain group. This figure is exported as a png. The *Reproduce.R* script imports these images and overlays them onto the figure overlay image. 

---

# Results

## HLA type HLA-A-301 shows decrease response to TIL therapy
 - The role of general HLA-A3 phenotype was investigated. None of the partial response patients featured an HLA-A3 phenotype. However, 2 of the stable disease patients demonstrated an HLA-A3 phenotype (HLA-A32, HLA-A33) making up a total of 50% of the stable disease patients.
## Lower amount of CDR3s per VJ in the baseline repertoire is associated with a higher response
 - Responders demonstrated a lower amount of CDR3s per VJ than non-responders in the baseline repertoire(1.80, range: 1.78 - 2.04 vs 2.72, range: 2.23 - 4.06). This result was verified with a one-tailed t-test (p=0.023704).
## VJ usage analysis demonstrates CDR3 specific expansion
 - The null hypothesis was rejected for both TIL and background clones with p values of 5.16 x 10-6 and 2.96 x 10-7, respectively. The results demonstrate that if a VJ expands it is due to specific CDR3 expansion rather than a proportional expansion of all CDR3s within a VJ. 
## Higher TCR convergence in the baseline repertoire is associated with a higher response
 - t was found that responders had a higher TCR convergence at the baseline than non-responders (2.14%, range: 1.18 - 2.94% vs 0.81%, range: 0.20% - 2.40). This result was confirmed through a two-sample t-test (p=0.044226). 
## Reduced richness and evenness of the baselne TCR repertoire is associated with increased response
 - the various diversity-related statistics demonstrate an important role for the baseline and its relation to the post-infusion sampes in determining response. This association is in despite of chemodepletion of the baseline repertoire and a lack of association between the TIL infusion product and response. 
## Increased homeostatic abundance of large clonotypes in the baseline TCR repertoire is associated with increased response
 - Overall, the results demonstrate that although associations between homeostatic abundance of large clonotypes and the level of expansion, only the homeostatic abundance of large clonotypes in the baseline are associated with response.
## Baseline clonotype abundance is associated with abundance post-infusion
 - Overall, it is seen that the baseline abundance of a clonotype dictates its expansion post-infusion despite chemodepletion of the baseline repertoire. 
## Increased correlation between gDNA and cfDNA TCR repertoires post-infusion
 - the cfDNA repertoire doesn’t show the same hyperexpansion of specific TIL clones as shown in the gDNA repertoire yet an increased similarity is seen between the sample cohorts post-infusion. This increased correlation could represent the high TCR turnover of redundant TILs after adoptive cell-transfer.
## Change in the amount of a productive rearrangements in cfDNA samples did not correlate with TIL-expansion
 -  Responsive patients had a median of 0.47% (max: 2.72%, min: -1.85%), non-responsive patients showed a median of -1.43% (max: 11.71%, min: -2.98%). 
## Proportional expansion of the TIL Infusion Product is associated with response
 - It was seen that responders had a larger amount of the TIL Infusion product taken up by the top 50% of the 4 week sample (25.57% range: 21.83% - 79.50% vs 0.269% range: 0.055% - 21.36%). This result was confirmed using a one-tailed two-sample t-test (p=0.014).
## Propotional expansion of TIL clones in the baseline is associated with increased response
 - It was found that responders had a higher amount of the baseline repertoire taken up by TILs when compared with non-responders (24.23%, range: 19.73% - 37.77% vs 4.42%, range: 2.98% - 13.61%). The association between the amount of TIL clonotypes in the baseline repertoire and the response was verified by a one-tailed two-sample t-test (p=0.00095).
