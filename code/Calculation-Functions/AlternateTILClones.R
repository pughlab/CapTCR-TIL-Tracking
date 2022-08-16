############################################
# Amount of TIL clones in Alternate Sample #
############################################

# Calculating the amount of a sample's repertoire that is taken up by TIL-related clones
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param sample: the sample which is subjected to TIL clone calculation
# @param expanded: boolean parameter, if TRUE the clones that are larger in the reference 
#                  than the TIL-infusion will be selected. Default is FALSE

TIL_calc <- function(patient, sampcohort, chain, clnefrc, sample, expanded=FALSE){
    
    #Loads longitudinal sample ordering
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))

    #Loads patient, sampcohort, chain specific data
    Load_data(patient, sampcohort, chain, clnefrc)
    
    # If a patient doesn't have a TIL file, it will use the TIL file from the related gDNA TIL data from the patient
    if(length(TIL_data$aaSeqCDR3) == 0){
      subset <- MiXCR_output[which(MiXCR_output$cohort=="gDNA" & MiXCR_output$patient==patient),]
      infusion <- unique(subset[grep("infusion|Infusion", subset$filename),]$filename)
      
      TIL_data <- as.data.frame(subset[which(subset$filename==infusion),])
        len <- length(samporder)
        TIL_data <- TIL_data[!duplicated(TIL_data$aaSeqCDR3),]
        TIL_data <- cbind(cloneno = row.names(TIL_data), 
                     filename = 'TIL Infusion Product', 
                   TIL_data)
        TIL_data <- TIL_data[, c("filename","aaSeqCDR3","cloneFraction", "cloneCount")]
        # Subset to include only clonotypes with more than specified clonal fraction    
        TIL_data <- TIL_data[TIL_data$cloneFraction > clnefrc,] 
        # Append the singletons
        TIL_data <- TIL_data[(TIL_data$cloneCount>1),]        
        # Removes unproductive clonotypes  
        TIL_data <- TIL_data[-c(grep("[*]", c(TIL_data$aaSeqCDR3, "CDR3*"))),]
        TIL_data <- TIL_data[-c(grep("_", c(TIL_data$aaSeqCDR3, "CDR3*"))),]

        # Readjusting clonal fraction values to spliced data
        ProdcloneFrac(TIL_data)
        TIL_data <- output_df
        # Generating dataframe with added gDNA TIL Infusion product and updated y-axis order for plots
        CDR3_fraction <- rbind(Base_data, TIL_data, FW_data)
        samporder <- c(samporder[1],'TIL Infusion Product',samporder[2:len])
    }
    # Creates dataframe including the clones which are prevalent in both the reference and TIL-infusion sample
    reference_data <- CDR3_fraction[which(CDR3_fraction$filename==sample),]
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
