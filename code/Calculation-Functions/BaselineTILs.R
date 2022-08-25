###################################
# Baseline TIL-clonotype counting #
###################################

# Calculating the amount of TIL-related clones in the baseline repertoire
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be DNA, RNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param char: desired characteristic to be analyzed ('clonotypes' or 'TCRs')

Base_TILclone <- function(patient, sampcohort, chain, char){
  
  # Setting the longitudinal order of the samples for patient
  samporder <- eval(as.name(paste0(patient, sampcohort, "_samporder")))
    
  # Pulls relative gDNA TIL Infusion product if the sample does not have an infusion
  if(length(grep("infusion", eval(as.name(paste0(patient, sampcohort, "_samporder"))))) == 0){
    DNA_infusion <- eval(as.name(paste0(patient, "DNA", "_samporder")))[2]
    CDR3_fraction <- rbind(CDR3_fraction, eval(as.name(paste0(patient, "DNA", "_", DNA_infusion))))
    samporder <- c(samporder[1], DNA_infusion,samporder[2:length(samporder)])
  }
    
    # Filtering the base data by the clones in the TIL data
    Base_filtered_data <- eval(as.name(paste0(patient, sampcohort, "_", samporder[1])))[eval(as.name(paste0(patient, sampcohort, "_", samporder[1])))$aaSeqCDR3 %in% eval(as.name(paste0(patient, sampcohort, "_", samporder[2])))$aaSeqCDR3,]

    if(char=="clonotypes"){
        # Find the number of TIL-clonotypes in the baseline sample
        result <<- nrow(Base_filtered_data)
    }
    
    if(char=="TCRs"){
        # Find the number of TIL-TCRs in the baseline sample
        result <<- sum(Base_filtered_data$cloneCount)
    }
}
