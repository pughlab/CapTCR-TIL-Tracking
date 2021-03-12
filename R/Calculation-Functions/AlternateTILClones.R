############################################
# Amount of TIL clones in Alternate Sample #
############################################

# Calculating the amount of a sample's repertoire is taken up by TIL-related clones
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param sample: the sample which is subjected to TIL clone calculation
# @param expanded: boolean parameter, if TRUE the clones that are larger in the reference 
#                  than the TIL-infusion will be selected. Default is FALSE
# @param dir_clones: directory where clone files are located
# @param dir_samplekeys: directory where the sample keys file is located
# @param file_samplekeys: name of the sample keys file

TIL_calc <- function(patient, sampcohort, chain, clnefrc, sample, expanded=FALSE, dir_clones, dir_samplekeys, file_samplekeys){
    
    #Loads longitudinal sample ordering
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
    options(scipen = 999)
    
    #Loads patient, sampcohort, chain specific data
    Load_data(patient, sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys)
    
    # If a patient doesn't have a TIL file, it will use the TIL file from the related gDNA TIL data from the patient
    if(length(TIL_data$aaSeqCDR3) == 0){
        files <<- list.files(paste(dir_clones, "gDNA/CLONES_", chain, patient, "/", sep = ""))
        infusion <- files[grep('fusion', files)]
        len <<- length(samporder)
        TIL_data <- data.frame()
        TIL_data <- read.table(paste(dir_clones, "gDNA/CLONES_", chain, patient, "/", infusion, sep = ""), 
                   header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))
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
    }
    # Creates dataframe including the clones which are prevalent in both the reference and TIL-infusion sample
    reference_data <- CDR3_fraction[which(CDR3_fraction$filename==sample),]
    reference_data <- reference_data[reference_data$aaSeqCDR3 %in% TIL_data$aaSeqCDR3,]
    
   # Subsets the data frame to the expansive clones if desired by the user
    if(expanded==TRUE){
        reference_data$TIL_cloneFraction <- reference_data$cloneFraction[match(reference_data$aaSeqCDR3, TIL_data$aaSeqCDR3)]
        reference_data <- reference_data[!(reference_data$cloneFraction < reference_data$TIL_cloneFraction),]
        reference_data[is.na(reference_data)] <- 0
    }
    
    # Calculates the amount of the sample is taken up by TIL-infusion clones and prints the result
    response <- sum(reference_data$cloneFraction)
    print(response)
}