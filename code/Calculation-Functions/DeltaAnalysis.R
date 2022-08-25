##################
# Delta Analysis #
##################

# Calculates the delta between 2 data points of different samples
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param datapoint: The data point wanted to be analyzed (Expansion, Productive, Expression, Diversity)
# @param first: The name of the first sample
# @param second: The name of the second sample

DeltaA <- function(patient, sampcohort, chain, datapoint, first, second){

    # Loads in patient samporder
    samporder <- eval(as.name(paste(patient, sampcohort, "_samporder", sep="")))
    
    # Loading CDR3_fraction & TIL_data
    CDR3_fraction <- eval(as.name(paste0(patient, sampcohort)))
    TIL_data <- eval(as.name(paste0(patient, sampcohort, "_", samporder[2])))
    
    # Pulls relative gDNA TIL Infusion product if the sample does not have an infusion
    if(length(grep("infusion", eval(as.name(paste0(patient, sampcohort, "_samporder"))))) == 0){
      DNA_infusion <- eval(as.name(paste0(patient, "DNA", "_samporder")))[2]
      CDR3_fraction <- rbind(CDR3_fraction, eval(as.name(paste0(patient, "DNA", "_", DNA_infusion))))
      samporder <- c(samporder[1], DNA_infusion,samporder[2:length(samporder)])
    }
    
    # Changes the name of the baseline to the first name in the sample order previously creates
    if(first=='BL'){
        first <- eval(as.name(paste(patient, sampcohort, "_samporder", sep="")))[1]
    }
    
    # Loads the clone data for the first and second identified samples
    df_first <- eval(as.name(paste0(patient, sampcohort, "_", first)))
    df_second <- eval(as.name(paste0(patient, sampcohort, "_", second)))

    # Finds the delta between the TIL clone fraction in the samples if user specified
    if(datapoint=="Expansion"){
        data_first <- sum(df_first[df_first$aaSeqCDR3 %in% TIL_data$aaSeqCDR3,]$cloneFraction)
        data_second <- sum(df_second[df_second$aaSeqCDR3 %in% TIL_data$aaSeqCDR3,]$cloneFraction)
        delta <- data_second - data_first
    }

    # Finds the delta between the total clone count in the samples if user specified
    if(datapoint=="Expression"){
        data_first <- sum(df_first$cloneCount)
        data_second <- sum(df_second$cloneCount)
        delta = data_second - data_first
    }
    
    # Finds the delta between the diversity in the samples if user specified
    if(datapoint=="Diversity"){
        N <- sum(df_first$cloneCount)
        sum_simpson <- 0
        for(n in df_first$cloneCount){
          sum_simpson <- sum_simpson + (n/N)^2
        }        
        data_first <- 1/sum_simpson

        N <- sum(df_second$cloneCount)
        sum_simpson <- 0
        for(n in df_second$cloneCount){
          sum_simpson <- sum_simpson + (n/N)^2
        }  
        data_second <- 1/sum_simpson
        
        delta = data_second - data_first
        
    }
    # Prints the resulting delta variable
}
