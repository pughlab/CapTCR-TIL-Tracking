persistence <- function(patient, sampcohort, timepoint) {
  # Computes the TIL persistence of a sample
  #
  # Args:
  #   patient: patient ID of sample
  #   sampcohort: sample cohort (DNA, RNA, cfDNA) of sample
  #   timepoint: full sample name
  # Returns: 
  #   sample TIL persistence
  
  # Load data -------------------
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  TIL_name <- eval(as.name(paste0(patient, "DNA_samporder")))[2]
  TIL_clones <- eval(as.name(paste0(patient, "DNA_", TIL_name)))$aaSeqCDR3
  # Calculate persistence -------------------
  persistent <- CDR3_df[which(CDR3_df$aaSeqCDR3 %in% TIL_clones),]
  sum(persistent$cloneFraction)
}
