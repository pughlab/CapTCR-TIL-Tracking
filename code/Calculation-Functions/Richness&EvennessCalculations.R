#######################
# Richness & Evenness #
#######################

# Prints the richness, evenness, diversity and clonality in a patient
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param dir_clones: parent directory where clone files are located
# @param dir_samplekeys: directory where the sample keys excel file are located
# @param file_samplekeys: file name of the sample keys 
# @param primary: desired sample to appear first (Baseline or TIL)
# @param max: desired maximum diversity value, outliers will be shown as a triangle

RichEvCalc <- function(patient, sampcohort, chain, clnefrc, dir_clones, 
                    dir_samplekeys, file_samplekeys){
  
  # Loading in patient data for specific patient, sample cohort, and chain
  Load_data(patient, sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys)
  # Setting the longitudinal order of the samples for patient
  samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
  
  # Creating outline for diversity dataframe
  Div_df <- data.frame(matrix(NA, nrow=length(samporder), ncol=5))
  Div_df$X1 <- samporder
  colnames(Div_df) <- c("Filename", "Diversity", "Richness", "Evenness", "Clonality")
  # Calculates the diversity, clonality, richness, and evenness of a sample
  for(cycle in 1:length(samporder)){
    cycle_df <- CDR3_fraction[which(CDR3_fraction$filename==samporder[cycle]),]
    Richness <- nrow(cycle_df)
    Diversity <- sum(cycle_df$cloneFraction*log2(cycle_df$cloneFraction))
    Diversity_max <- 0
    for(i in 1:Richness){
      Diversity_max <- Diversity_max + (1/Richness)*log2(1/Richness)
    }
    Evenness <- Diversity/Diversity_max
    Clonality <- 1 - Evenness
    Div_df$Diversity[cycle] <- Diversity
    Div_df$Richness[cycle] <- Richness
    Div_df$Evenness[cycle] <- Evenness
    Div_df$Clonality[cycle] <- Clonality
  }
  # Orders filenames based on chronological sample order
  Div_df$Filename <- factor(Div_df$Filename, levels = c(samporder))
  levels(Div_df$Filename) <- c(samporder)
  Div_df <<- Div_df
  print(Div_df)
}
