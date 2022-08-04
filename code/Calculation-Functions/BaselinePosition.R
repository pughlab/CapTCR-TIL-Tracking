##############################
# Baseline Position Analysis #
##############################

# Analyzes the position in the baseline repertoire of the most abundant post-infusion clones
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param dir_clones: parent directory where clone files are located
# @param dir_samplekeys: directory where the sample keys excel file are located
# @param file_samplekeys: file name of the sample keys 

Base_Position <- function(patient,sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys){

  # Load data for the specific patient, sample cohort, and chain
  Load_data(patient,sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys)
  samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
  
  # Creates a dataframe of the post-infusion TIL clones
  PI_df <- CDR3_fraction[which(CDR3_fraction$filename %in% samporder[3:length(samporder)]),]
  PI_df <- PI_df[which(PI_df$aaSeqCDR3 %in% TIL_data$aaSeqCDR3==TRUE),]
  # Creates a new dataframe 'PI_clonefraction' consisting of a clonotype's average post-infusion clone fraction
  PI_clonefraction <- data.frame(matrix(0, ncol=2, nrow=length(unique(PI_df$aaSeqCDR3))))
  colnames(PI_clonefraction) <- c("aaSeqCDR3", "AvgPICloneFraction")
  PI_clonefraction$aaSeqCDR3 <- unique(PI_df$aaSeqCDR3)
  for(aaSeq in 1:length(unique(PI_df$aaSeqCDR3))){
    PI_clonefraction$AvgPICloneFraction[aaSeq] <- mean(PI_df[which(PI_df$aaSeqCDR3==unique(PI_df$aaSeqCDR3)[aaSeq]),]$cloneFraction)
  }
  
  # Creates a dataframe of the pre-infusion TIL clones
  TIL_BL <- Base_data[which(Base_data$aaSeqCDR3 %in% TIL_data$aaSeqCDR3==TRUE),]
  # Adds a column in the 'PI_clonefraction' dataframe consisting of the position of the clone in the baseline 
  PI_clonefraction$BL_position <- NA
  for(row in 1:length(PI_clonefraction$aaSeqCDR3)){
    if(PI_clonefraction$aaSeqCDR3[row]%in%TIL_BL$aaSeqCDR3==TRUE){
      PI_clonefraction$BL_position[row] <- which(TIL_BL$aaSeqCDR3 == PI_clonefraction$aaSeqCDR3[row])
    } else {
      PI_clonefraction$BL_position[row] <- NA
    }
  }
  # Removes post-infusion clones that are not in the baseline
  PI_clonefraction <- na.omit(PI_clonefraction[order(PI_clonefraction$AvgPICloneFraction, decreasing=TRUE),])
 # Add a column showing the ranking in terms of most abundant post-infusion TIL clone
  PI_clonefraction$row <- 1:nrow(PI_clonefraction)
  PI_clonefraction
}
