########################
# VJ data calculations #
########################

# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param timepoint: Desired timepoint to be analyzed
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG

Convergence <- function(patient, sampcohort, chain, timepoint){
  # Returns the TCR convergence for patient and timepoint in given sampcohort
  # and chain. 
  CDR3_df <- clonotypes %>% 
    filter ( Patient_id %in% patient)%>%
    filter ( Cohort == cohort) %>%
    filter ( Cycle == timepoint) %>%
    filter ( cloneCount > 1) %>%
    mutate(VJcombo = gsub("[.]", "_", VJcombo))
  
  # Removing unproductive clones and recomputing the clone fraction
  CDR3_df <- CDR3_df[-c(grep("[*]", CDR3_df$aaSeqCDR3)),]
  CDR3_df <- CDR3_df[-c(grep("_", CDR3_df$aaSeqCDR3)),]
  for(clone in 1:length(CDR3_df$cloneFraction)){
    CDR3_df$cloneFraction[clone] <- CDR3_df$cloneCount[clone]/sum(CDR3_df[which(CDR3_df$Cycle==CDR3_df$Cycle[clone]),]$cloneCount)
  }
  
  CDR3_df$occurence <- NA
  for(CDR3 in 1:nrow(CDR3_df)){
    CDR3_df$occurence[CDR3] <- length(CDR3_df[which(CDR3_df$aaSeqCDR3==CDR3_df$aaSeqCDR3[CDR3]),]$aaSeqCDR3)
  }
  sum(CDR3_df[which(CDR3_df$occurence > 1),]$cloneFraction)
}
CDR3perVJ <- function(patient, sampcohort, chain, timepoint){
  # Returns the average CDR3 per VJ for patient at timepoint in given sampcohort
  # and chain
  
  CDR3_df <- clonotypes %>% 
    filter ( Patient_id %in% patient)%>%
    filter ( Cohort == cohort) %>%
    filter ( Cycle == timepoint) %>%
    mutate(VJcombo = gsub("[.]", "_", VJcombo))
  
  # Removing unproductive clones and recomputing the clone fraction
  CDR3_df <- CDR3_df[-c(grep("[*]", CDR3_df$aaSeqCDR3)),]
  CDR3_df <- CDR3_df[-c(grep("_", CDR3_df$aaSeqCDR3)),]
  for(clone in 1:length(CDR3_df$cloneFraction)){
    CDR3_df$cloneFraction[clone] <- CDR3_df$cloneCount[clone]/sum(CDR3_df[which(CDR3_df$Cycle==CDR3_df$Cycle[clone]),]$cloneCount)
  }
  length(CDR3_df$aaSeqCDR3)/length(unique(CDR3_df$VJcombo))
}