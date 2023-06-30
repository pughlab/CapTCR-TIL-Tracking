# Functions relating to change point analysis
#
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param timepoint: Desired timepoint to be analyzed

changepoint <- function(patient, sampcohort, timepoint) {
  # Computes the change point threshold
  #
  # Returns: 
  #   change point threshold
  
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  cpt <- cpt.var(CDR3_df$cloneCount)
  threshindex <- cpt@cpts[1]
  threshold <- CDR3_df$cloneCount[threshindex]
  threshold
}
get_cpdf <- function(patient, sampcohort, timepoint) {
  # Adds column describing whether clonotype is expanded to repertoire dataframe
  #
  # Returns: 
  #   dataframe with added column
  
  CDR3_df <- eval(as.name(paste0(patient, sampcohort, "_", timepoint)))
  threshold <- changepoint(patient, sampcohort, timepoint)
  CDR3_df$Group <- NA
  CDR3_df$Group[which(CDR3_df$cloneCount < threshold)] <- "Background"
  CDR3_df$Group[which(CDR3_df$cloneCount >= threshold)] <- "Expanded"
  CDR3_df
}
threshold_stats <- function(patient, sampcohort, timepoint, stat) {
  # Adds column describing whether clonotype is expanded to repertoire dataframe
  #
  # Args:
  #   stat: "percentage", "count", "abundance", "numberabove"
  # Returns: 
  #   dataframe with added column
  
  CDR3_df <- get_cpdf(patient, sampcohort, timepoint)
  if(stat == "percentage") {
    value <- CDR3_df$cloneFraction[tail(which(CDR3_df$Group=='Expanded', arr.ind=TRUE), n=1)]
  }
  else if(stat == "count") {
    value <- CDR3_df$cloneCount[tail(which(CDR3_df$Group=='Expanded', arr.ind=TRUE), n=1)]
  }
  else if(stat == "abundance") {
    value <- sum(CDR3_df$cloneFraction[which(CDR3_df$Group == 'Expanded')])
  }
  else if(stat == "numberabove") {
    value <- length(CDR3_df$cloneFraction[which(CDR3_df$Group == 'Expanded')])
  }
  value
}
expanded_reference <- function(patient1, sampcohort1, timepoint1, 
                               patient2, sampcohort2, timepoint2, stat) {
  # Computes the percentage of expanded clones in one sample present in another
  #
  # Args:
  #   patient1: patient ID of sample 1
  #   sampcohort1: sample cohort (DNA, RNA, cfDNA) of sample 1
  #   timepoint1: full name of sample 1
  #   patient2: patient ID of sample 2
  #   sampcohort2: sample cohort (DNA, RNA, cfDNA) of sample 2
  #   timepoint2: full name of sample 2
  #   stat: "Expanded", "Background", "total"
  # Returns: 
  #   Morisita's overlap index  
  
  # Load data -------------------
  CDR3_df1 <- get_cpdf(patient1, sampcohort1, timepoint1)
  CDR3_df1 <- CDR3_df1[which(CDR3_df1$Group == "Expanded"),]
  CDR3_df2 <- get_cpdf(patient2, sampcohort2, timepoint2)
  if(stat != "total"){
    CDR3_df2 <- CDR3_df2[which(CDR3_df2$Group == stat),]
  }
  # Calculate proportion -------------------
  length(which(CDR3_df1$aaSeqCDR3 %in% CDR3_df2$aaSeqCDR3))/nrow(CDR3_df1)
}
relativerisk <- function(patient, sampcohort, timepointreference) {
  # Computes the relative risk between timepointreference and the first post-
  # infusion sample
  #
  # Args:
  #   timepointreference: reference for relative risk analysis
  # Returns: 
  #   relative risk and confidence interval
  
  # Load data -------------------
  total_df <- data.frame()
  samporder <- eval(as.name(paste0(patient, sampcohort, "_samporder")))
  for(timepoint in samporder){
    CDR3_df <- get_cpdf(patient, sampcohort, timepoint)
    total_df <- rbind(total_df, CDR3_df)
  }
  FW_clones <- total_df[which(total_df$filename == samporder[3]),]
  FW_expanded <- FW_clones$aaSeqCDR3[which(FW_clones$Group == "Expanded")]
  FW_background <- FW_clones$aaSeqCDR3[which(FW_clones$Group == "Background")]
  timepoint_clones <-  total_df$aaSeqCDR3[which(total_df$filename == timepointreference)]
  
  total_fw = nrow(FW_clones)
  # Calcualte relative risk -------------------
  Expanded <- length(which(FW_expanded %in% timepoint_clones))
  Background <- length(which(FW_background %in% timepoint_clones))
  exp_ref <- length(which(FW_expanded %in% timepoint_clones))
  back_ref <- length(which(FW_background %in% timepoint_clones))
  exp_non <- length(FW_expanded) - exp_ref
  back_non <- length(FW_background) - back_ref
  total = exp_ref + back_ref + exp_non + back_non
  RR <- (exp_ref/(exp_ref+back_ref)) / (exp_non/(exp_non+back_non))
  SE_logRR <- sqrt(1/exp_ref + 1/back_ref - 1/(exp_ref+exp_non) - 1/(back_ref+back_non))
  conf_int <- c(exp(log(RR) - 1.96*SE_logRR), exp(log(RR) + 1.96*SE_logRR))
  # Return relative risk & confidence interval -------------------
  list("RelativeRisk" = RR, "ConfidenceInterval" = conf_int)
}
