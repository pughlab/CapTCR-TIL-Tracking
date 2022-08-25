#####################################
# Timepoint Characteristic Analysis #
#####################################

# Analyzes the amount of clones >5%, total number of TCRs and the total amount of clones at a timepoint
# @param patient: specific patient code
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param char: The characteristic desired to be analyzed (Clonotypes, TCRCount, or Oligo)
# @param sample: the sample which is used in the analysis

Timepoint_char <- function(patient, sampcohort, chain, clnefrc, char, sample){
  
  # Loads longitudinal sample order
  samporder <- eval(as.name(paste0(patient, sampcohort, "_samporder")))

  #Loads sample data
  reference_data <- eval(as.name(paste0(patient, sampcohort, "_", sample)))
  
  # Prints the total amount of TIL clones if specified by the user
  if(char=="Clonotypes"){
    Total <<- nrow(reference_data)
  }
  
  if(char=="TCRCount"){
    TCRCount <<- sum(reference_data$cloneCount)
  }
  
  # Prints the fraction of clones >5% in the TIL sample if specified by the user
  if(char=="Oligo"){
    reference_data <- reference_data[which(reference_data$cloneFraction > clnefrc),]
    Sum_Oligo <<- sum(reference_data$cloneFraction)
  }
  if(char=="AvgCloneFrac"){
    AvgCloneFrac <<- sum(reference_data$cloneFraction)/nrow(reference_data)
  }
}
