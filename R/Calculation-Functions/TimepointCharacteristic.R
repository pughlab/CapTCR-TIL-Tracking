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
# @param dir_clones: directory where clone files are located
# @param dir_samplekeys: directory where the sample keys file is located
# @param file_samplekeys: name of the sample keys file

Timepoint_char <- function(patient, sampcohort, chain, clnefrc, char, sample, dir_clones, dir_samplekeys, file_samplekeys){
  
  # Loads longitudinal sample order
  samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
  options(scipen = 999)
  
  #Loads patient clone data
  Load_data(patient, sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys)
  
  reference_data <- CDR3_fraction[which(CDR3_fraction$filename==sample),]
  
  # Prints the total amount of TIL clones if specified by the user
  if(char=="Clonotypes"){
    Total <- nrow(reference_data)
    print(Total)
  }
  
  if(char=="TCRCount"){
    TCRCount <- sum(reference_data$cloneCount)
    print(TCRCount)
  }
  
  # Prints the fraction of clones >5% in the TIL sample if specified by the user
  if(char=="Oligo"){
    Sum_Oligo <- sum(reference_data$cloneFraction)
    print(Sum_Oligo)
  }
  if(char=="AvgCloneFrac"){
    AvgCloneFrac <- sum(reference_data$cloneFraction)/nrow(reference_data)
    print(AvgCloneFrac)
  }
}
