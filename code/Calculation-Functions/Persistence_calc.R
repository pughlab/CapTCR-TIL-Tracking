#############################
# Persistence and Dominance #
#############################

# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param timepoint: Desired timepoint to be analyzed
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG

Persistence <- function(patient, sampcohort, chain, timepoint){
  # Returns the persisistence of the TIL infusion product clonotypes in patient 
  # at timepoint for a given sampcohort and chain
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  TIL_name <- eval(as.name(paste0(patient, sampcohort, "_samporder")))[2]
  TIL_clones <- eval(as.name(paste0(patient,sampcohort, "_", TIL_name)))$aaSeqCDR3
  
  persistent <- CDR3_df[which(CDR3_df$aaSeqCDR3 %in% TIL_clones),]
  sum(persistent$cloneFraction)
}
BergerParkerIndex <- function(patient, sampcohort, chain, timepoint){
  # Returns the Berger Parker index in patient at timepoint for a given
  # sampcohort and chain. The index is the fraction that the most dominant 
  # clonotype takes up
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  max(CDR3_df$cloneFraction)
}
ChangePoint <- function(patient, sampcohort, chain, timepoint){
  # Returns the threshold clone count and plots the rank abudance curve for
  # the set of clonotypes and a horizontal line representing the threshold
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  
  cpt <- cpt.var(CDR3_df$cloneCount)
  threshindex <- cpt@cpts[1]
  threshold <- CDR3_df$cloneCount[threshindex]
  
  plot(CDR3_df$cloneCount, ylab = "Clone Count") + title(main = "Clone Count Rank Abundance Curve")
  for(i in threshold){
    abline(h=i)
  }
  threshold
}
GetAlluvialdf <- function(patient, sampcohort, chain, timepoint){
  # Returns the normal CDR3_df but with an added column containing whether that
  # clonotype is a part of the "Background" or "Expanded" group according to 
  # changepoint analysis
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  threshold <- ChangePoint(patient, sampcohort, chain, timepoint)
  CDR3_df$Group <- NA
  CDR3_df$Group[which(CDR3_df$cloneCount < threshold)] <- "Background"
  CDR3_df$Group[which(CDR3_df$cloneCount >= threshold)] <- "Expanded"
  CDR3_df
}
Timepointreference_changepoint <- function(patient, sampcohort, chain, timepointreference, group){
  # Returns the fraction of the expanded clones in the 4 week sample that are 
  # present in the section of timepointreference defined by group (i.e the 
  # background section of the baseline sample)  
  
  # @param timepointreference: timepoint to compare with the 4 week sample clones
  # @param group: Either 'Background' or 'Expanded'
  total_df <- data.frame()
  samporder <- eval(as.name(paste0(patient, sampcohort, "_samporder")))
  for(timepoint in samporder){
    CDR3_df <- GetAlluvialdf(patient, sampcohort, chain, timepoint)
    total_df <- rbind(total_df, CDR3_df)
  }
  FW_expanded_clones <- total_df$aaSeqCDR3[which(total_df$Group == "Expanded" & total_df$filename == samporder[3])]
  timepointreference_expanded_clones <- total_df$aaSeqCDR3[which(total_df$Group == group & total_df$filename == timepointreference)]
  length(which(FW_expanded_clones %in% timepointreference_expanded_clones))/length(FW_expanded_clones)
}
FractionExpandedClones <- function(patient, sampcohort, chain, timepoint){
  # Returns the fraction of the sample at timepoint which is made up of expanded
  # clones as given by changepoint analysis
  CDR3_df <- GetAlluvialdf(patient, sampcohort, chain, timepoint)
  sum(CDR3_df$cloneFraction[which(CDR3_df$Group == "Expanded")])
}