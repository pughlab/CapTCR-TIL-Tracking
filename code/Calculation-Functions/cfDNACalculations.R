###############################
# cfDNA Specific Calculations #
###############################

# Calculates the gDNA overlap and total amount of TIL vs background clones for cfDNA samples
# @param patient: specific patient code 
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param analysis: Desired cfDNA analysis ("gDNAOverlap" or "TotalTIL&Background")

cfDNA_calc <- function(patient, chain, analysis){

  # Loading both gDNA and cfDNA data and assigning the 4 week and baseline data to differing data sets
  gDNA_samporder <- eval(as.name(paste0(patient, "DNA_samporder")))
  cfDNA_samporder <- eval(as.name(paste0(patient, "cfDNA_samporder")))
  gDNA_FW <- eval(as.name(paste0(patient, "DNA_", gDNA_samporder[3])))
  gDNA_Base <- eval(as.name(paste0(patient, "DNA_", gDNA_samporder[1])))
  gDNA_Total <- eval(as.name(paste0(patient, "DNA")))
  cfDNA_FW <- eval(as.name(paste0(patient, "cfDNA_", cfDNA_samporder[2])))
  cfDNA_Base <- eval(as.name(paste0(patient, "cfDNA_", cfDNA_samporder[1])))
  cfDNA_Total <- eval(as.name(paste0(patient, "cfDNA")))
  
  # Create outline for dataframe containing all CDR3s and their clonefraction in both sample cohorts and timepoints
  TotalCDR3 <- unique(c(cfDNA_Total$aaSeqCDR3, gDNA_Total$aaSeqCDR3))
  LongitudinalComparison <- data.frame(matrix(nrow=length(TotalCDR3), ncol=5, data=NA))
  colnames(LongitudinalComparison) <- c("aaSeqCDR3", "gDNA_FW_cloneFraction", "cfDNA_FW_cloneFraction", 
                                        "gDNA_Base_cloneFraction", "cfDNA_Base_cloneFraction")
 # Assign values in empty dataframe created previously
  LongitudinalComparison$aaSeqCDR3 <- TotalCDR3
  LongitudinalComparison$gDNA_FW_cloneFraction <- gDNA_FW$cloneFraction[match(LongitudinalComparison$aaSeqCDR3, gDNA_FW$aaSeqCDR3)]
  LongitudinalComparison$gDNA_Base_cloneFraction <- gDNA_Base$cloneFraction[match(LongitudinalComparison$aaSeqCDR3, gDNA_Base$aaSeqCDR3)]
  LongitudinalComparison$cfDNA_FW_cloneFraction <- cfDNA_FW$cloneFraction[match(LongitudinalComparison$aaSeqCDR3, cfDNA_FW$aaSeqCDR3)]
  LongitudinalComparison$cfDNA_Base_cloneFraction <- cfDNA_Base$cloneFraction[match(LongitudinalComparison$aaSeqCDR3, cfDNA_Base$aaSeqCDR3)]
  
  LongitudinalComparison[is.na(LongitudinalComparison)] <- 0
  
  # Creating a column consisting of the colors tracked in the figures (redundant TIL clones)
  CDR3_colors_id <- CDR3_colors[which(CDR3_colors$id==patient),]
  col <- as.character(CDR3_colors_id$mycolors)
  names(col) <- as.character(CDR3_colors_id$colored_clns)
  
  LongitudinalComparison$Color <- CDR3_colors_id$mycolors[match(LongitudinalComparison$aaSeqCDR3, CDR3_colors_id$colored_clns)]

  # Compute either the total amount of TIL and background clones or the overlap with the gDNA sample
  if (analysis == "TotalTIL&Background"){
    TIL_df <- na.omit(LongitudinalComparison)
    FW_TIL <- TIL_df[which(TIL_df$cfDNA_FW_cloneFraction > 0),]
    Base_TIL <- TIL_df[which(TIL_df$cfDNA_Base_cloneFraction > 0),]

    Background_df <- LongitudinalComparison[is.na(LongitudinalComparison$Color),]
    FW_Background <- Background_df[which(Background_df$cfDNA_FW_cloneFraction > 0),]
    Base_Background <- Background_df[which(Background_df$cfDNA_Base_cloneFraction > 0),]
    
    result <- data.frame(matrix(ncol=4, nrow=1, data=NA))
    colnames(result) <- c("Base_Background", "Base_TIL", "FW_Background", "FW_TIL")
    result$Base_Background[1] <- sum(Base_Background$cfDNA_Base_cloneFraction)
    result$Base_TIL[1] <- sum(Base_TIL$cfDNA_Base_cloneFraction)
    result$FW_Background[1] <- sum(FW_Background$cfDNA_FW_cloneFraction)
    result$FW_TIL[1] <- sum(FW_TIL$cfDNA_FW_cloneFraction)
  } else if(analysis == "gDNAOverlap"){
    result <- data.frame(matrix(ncol=2, nrow=1, data=NA))
    colnames(result) <- c("gDNA_Base", "gDNA_FW")
    Base_df <- LongitudinalComparison[which(LongitudinalComparison$cfDNA_Base_cloneFraction > 0),]
    FW_df <- LongitudinalComparison[which(LongitudinalComparison$cfDNA_FW_cloneFraction > 0),]
    
    result$gDNA_Base[1] <- sum(Base_df$gDNA_Base_cloneFraction)
    result$gDNA_FW[1] <- sum(FW_df$gDNA_FW_cloneFraction)
  }
  result <<- result
} 
