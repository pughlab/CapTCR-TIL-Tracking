############################################
# Amount of TIL clones in Alternate Sample #
############################################

# Calculating the amount of a sample's repertoire that is taken up by TIL-related clones
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be DNA, RNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param sample: the sample which is subjected to TIL clone calculation
# @param expanded: boolean parameter, if TRUE the clones that are larger in the reference 
#                  than the TIL-infusion will be selected. Default is FALSE

TIL_calc <- function(patient, sampcohort, chain, sample, expanded=FALSE){
    
  # Loading reference data
  reference_data <- eval(as.name(paste0(patient, sampcohort, "_", sample)))
  
  # Setting the longitudinal order of the samples for patient
  samporder <- eval(as.name(paste0(patient, sampcohort, "_samporder")))
  
  TIL_sample <- samporder[grep("infusion", eval(as.name(paste0(patient, sampcohort, "_samporder"))))]
  
  TIL_data <- eval(as.name(paste0(patient, sampcohort, "_", TIL_sample)))
    
  # Pulls relative gDNA TIL Infusion product if the sample does not have an infusion
  if(length(grep("infusion", eval(as.name(paste0(patient, sampcohort, "_samporder"))))) == 0){
      DNA_infusion <- eval(as.name(paste0(patient, "DNA", "_samporder")))[2]
      CDR3_fraction <- rbind(CDR3_fraction, eval(as.name(paste0(patient, "DNA", "_", DNA_infusion))))
      samporder <- c(samporder[1], DNA_infusion,samporder[2:length(samporder)])
    }
    # Creates dataframe including the clones which are prevalent in both the reference and TIL-infusion sample
    reference_data <- reference_data[reference_data$aaSeqCDR3 %in% TIL_data$aaSeqCDR3,]
    
   # Subsets the data frame to the expansive clones if desired by the user
    if(expanded==TRUE){
        reference_data$TIL_cloneFraction <- TIL_data$cloneFraction[match(reference_data$aaSeqCDR3, TIL_data$aaSeqCDR3)]
        reference_data <- reference_data[!(reference_data$cloneFraction < reference_data$TIL_cloneFraction),]
        reference_data[is.na(reference_data)] <- 0
    }
    
    # Calculates the amount of the sample is taken up by TIL-infusion clones and prints the result
    response <<- sum(reference_data$cloneFraction)
}
