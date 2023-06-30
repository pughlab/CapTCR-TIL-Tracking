# Functions relating to calculating diversity indices
#
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param timepoint: Desired timepoint to be analyzed

SimpsonIndex <- function(patient, sampcohort, timepoint) {
  # Computes the Simpson's diversity index of a sample
  #
  # Returns: 
  #   sample simpson's diversity index
  # Note:
  #   inverse simpson diversity is 1 / simpson's index
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  sum(CDR3_df$cloneFraction^2)
}
ShannonIndex <- function(patient, sampcohort, timepoint) {
  # Computes the Shannon's diversity index of a sample
  #
  # Returns: 
  #   sample Shannon's diversity index
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  -sum(CDR3_df$cloneFraction*log(CDR3_df$cloneFraction))
}
Richness <- function(patient, sampcohort, timepoint) {
  # Computes the richness of a sample
  #
  # Returns: 
  #   sample richness
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  nrow(CDR3_df)
}
Eveness <- function(patient, sampcohort, timepoint) {
  # Computes the eveness of a sample
  #
  # Returns: 
  #   sample eveness
  # Note:
  #   clonality is 1 - eveness
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  shannon <- ShannonIndex(patient, sampcohort, timepoint)
  richness <- Richness(patient, sampcohort, timepoint)
  shannon/log(richness)
}
HighExpandedClone <- function(patient, sampcohort, timepoint, frac) {
  # Computes the abundance of clones greater than a specified threshold
  #
  # Args:
  #   frac: threshold for high expanded clone
  # Returns: 
  #   sample high expanded clone abundance
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  sum(CDR3_df[which(CDR3_df$cloneFraction > frac),]$cloneFraction)
}
