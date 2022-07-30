###########################################
# Top 10 4 week clone fraction Comparison #
###########################################

# Calculates the clone fraction of the top 4 week clones in the baseline and TIL-infusion data
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param dir_clones: parent directory where clone files are located
# @param dir_samplekeys: directory where the sample keys excel file are located
# @param file_samplekeys: file name of the sample keys 

Top104WComp <- function(patient, sampcohort, chain, clnefrc, 
                        dir_clones, dir_samplekeys, file_samplekeys){
  
  # Loads in patient data
  Load_data(patient, sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys)
  
  # Loading patient sample order
  samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
  
  # Create Top104W dataframe out of the top 10 clones in the 4 week sample
  Top104W <- FW_data[1:10,]
  if(patient=="TLML_1_"){
    Top104W <- CDR3_fraction[which(CDR3_fraction$filename==TLML_1_gDNA[3]),][1:10,]
  }
  Top104W$TIL_cloneFraction <- replicate(10, 0)
  Top104W$Base_cloneFraction <- replicate(10, 0)
  
  # Fill in the TIL_cloneFraction and Base_cloneFraction with the clone fraction of the corresponding 4 week clone
  for(i in 1:10){
    if(length(TIL_data$cloneFraction[which(TIL_data$aaSeqCDR3 == Top104W$aaSeqCDR3[i])])>0){
      Top104W$TIL_cloneFraction[i] <- TIL_data$cloneFraction[which(TIL_data$aaSeqCDR3 == Top104W$aaSeqCDR3[i])]
    }
    if(length(Base_data$cloneFraction[which(Base_data$aaSeqCDR3 == Top104W$aaSeqCDR3[i])])>0){
      Top104W$Base_cloneFraction[i] <- Base_data$cloneFraction[which(Base_data$aaSeqCDR3 == Top104W$aaSeqCDR3[i])]
    }
  }
  
  # Calculates and prints average cloneFraction for the TIL and Base Data for the top 10 4 week clones
  TIL_cloneFraction <<- sum(Top104W$TIL_cloneFraction)/10
  Base_cloneFraction <<- sum(Top104W$Base_cloneFraction)/10
  Top104W_cloneFraction <<- sum(Top104W$cloneFraction)/10
  print(TIL_cloneFraction)
  print(Base_cloneFraction)
  print(sum(Top104W$cloneFraction)/10)
}
