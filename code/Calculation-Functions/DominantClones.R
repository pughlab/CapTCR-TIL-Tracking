###########################
# Dominant Clone Analysis #
###########################

# Calculates various statistics about the dominant clones in a patient
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param dir_clones: parent directory where clone files are located
# @param dir_samplekeys: directory where the sample keys excel file are located
# @param file_samplekeys: file name of the sample keys 

DomClone <- function(patient, sampcohort, chain, clnefrc, 
         dir_clones, dir_samplekeys, file_samplekeys){
  
  # Loads data for the speicif patient, sample cohort, and chain
  Load_data(patient, sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys)
  samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
  
  # Creates a dataframe 'Hyperexpansive_Postinfusion' containing the post-infusion hyperexpansive clones
  Hyperexpansive <- CDR3_fraction[which(CDR3_fraction$cloneFraction > 0.05),]
  Hyperexpansive_Postinfusion <- Hyperexpansive[which(Hyperexpansive$filename %in% samporder[3:length(samporder)]==TRUE),]
  # Creates a dataframe 'Dominant_Clones' containing the hyperexpansive clones present in more than 1 post-infusion sample
  Redundant_Clones <- unique(Hyperexpansive_Postinfusion$aaSeqCDR3[which(duplicated(Hyperexpansive_Postinfusion$aaSeqCDR3)==TRUE)])
  Dominant_Clones <- data.frame()
  BL_fractions <- c()
  TIL_fractions <- c()
  for(clone in Redundant_Clones){
    Dominant_Clones <- rbind(Dominant_Clones, Hyperexpansive_Postinfusion[which(Hyperexpansive_Postinfusion$aaSeqCDR3==clone),])
    # Creates a list of the baseline and TIL product clone fractions for the specific clonotype
    BL_fractions <- c(BL_fractions, Base_data[which(Base_data$aaSeqCDR3==clone),]$cloneFraction)
    TIL_fractions <- c(TIL_fractions, TIL_data[which(TIL_data$aaSeqCDR3==clone),]$cloneFraction)
    # Prints whether the dominant clone is in the top 3 clones of the baseline repertoire
    print(clone %in% Base_data$aaSeqCDR3[1:3])
  }
  # Prints the amount of dominant clones present in the patient
  Total_Dominant <- length(Redundant_Clones)
  # Prints the average post-infusion clone fraction of the dominant clones
  Avg_CloneFraction_PI <- mean(Dominant_Clones$cloneFraction)
  # Prints the average baseline fraction of the dominant clones
  if(length(BL_fractions)==1){
    Avg_CloneFraction_BL <- BL_fractions
  } else {
    Avg_CloneFraction_BL <- mean(BL_fractions)
  }
  # Prints the average TIL product fraction of the dominant clones
  if(length(TIL_fractions)==1){
    Avg_CloneFraction_TIL <- TIL_fractions
  } else {
    Avg_CloneFraction_TIL <- mean(TIL_fractions)
  }
  # Prints the fraction of post-infusion samples in which the dominant-clone is dominant
  PercentLongevity <- length(unique(Dominant_Clones$filename))/(length(samporder)-2)
  # Creates and prints a dataframe 'DomClones_data' containing the statistics calcualted above
  DomClones_data <- data.frame(c("Total_Dominant", "Avg_CloneFraction_PI", "Avg_CloneFraction_BL", "Avg_CloneFraction_TIL", "PercentLongevity"), 
             c(Total_Dominant, Avg_CloneFraction_PI, Avg_CloneFraction_BL[[1]], Avg_CloneFraction_TIL[[1]], PercentLongevity))
  colnames(DomClones_data) <- c("Name", "Variable")
  DomClones_data
}