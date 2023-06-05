########################
# Diversity Evaluation #
########################

# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param timepoint: Desired timepoint to be analyzed
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG

SimpsonIndex <- function(patient, sampcohort, chain, timepoint){
  # Returns the Simpson index for patient at timepoint in given sampcohort
  # and chain. Note: Inverse simpson index is 1/simpson
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  sum(CDR3_df$cloneFraction^2)
}
ShannonIndex <- function(patient, sampcohort, chain, timepoint){
  # Returns the Shannon index for patient at timepoint in given sampcohort and
  # chain
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  -sum(CDR3_df$cloneFraction*log(CDR3_df$cloneFraction))
}
Richness <- function(patient, sampcohort, chain, timepoint){
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  nrow(CDR3_df)
}
NumberofCells <- function(patient, sampcohort, chain, timepoint){
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  sum(CDR3_df$cloneCount)
}
Eveness <- function(patient, sampcohort, chain, timepoint){
  # Returns the Eveness for patient at timepoint in given sampcohort and
  # chain. Note: Clonality is 1 - eveness
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  shannon <- ShannonIndex(patient, sampcohort, chain, timepoint)
  richness <- Richness(patient, sampcohort, chain, timepoint)
  shannon/log(richness)
}
MorisitaIndex <- function(patient, sampcohort, chain, timepoint1, timepoint2){
  # Returns the Morisita overlap index for patient between timepoint1 and
  # timepoint2 in given sampcohort and chain.
  
  # @param timepoint1/timepoint2: first and second timepoints to be compared
  CDR3_df1 <- eval(as.name(paste0(patient, sampcohort, "_", timepoint1)))
  CDR3_df2 <- eval(as.name(paste0(patient, sampcohort, "_", timepoint2)))
  simpson1 <- SimpsonIndex(patient, sampcohort, chain, timepoint1)
  simpson2 <- SimpsonIndex(patient, sampcohort, chain, timepoint2)
  
  total_unique <- unique(c(CDR3_df1$aaSeqCDR3, CDR3_df2$aaSeqCDR3))
  sum_numerator <- 0
  for(CDR3 in total_unique){
    frac1 <- CDR3_df1[which(CDR3_df1$aaSeqCDR3==CDR3),]$cloneFraction
    frac2 <- CDR3_df2[which(CDR3_df2$aaSeqCDR3==CDR3),]$cloneFraction
    if(length(frac1)!=0 & length(frac2)!=0){
      sum_numerator <- sum_numerator + frac1*frac2
    }
  }
  2*sum_numerator/(simpson1+simpson2)
}
HighExpandedClone <- function(patient, sampcohort, chain, timepoint, frac){
  # Returns the proportion of timepoint taken up by clones that have a clone 
  # fraction of atleast frac.
  
  # @param frac: the clone fraction which defines a highly expanded clone
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  sum(CDR3_df[which(CDR3_df$cloneFraction > frac),]$cloneFraction)
}