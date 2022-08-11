###################################
# Baseline TIL-clonotype counting #
###################################

# Calculating the amount of TIL-related clones in the baseline repertoire
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param char: desired characteristic to be analyzed ('clonotypes' or 'TCRs')

Base_TILclone <- function(patient, sampcohort, chain, clnefrc, char){
    
    # Loading specific locus and sampcohort 4W, TIL, and baseline samples for patient
    Load_data(patient, sampcohort, chain, clnefrc)
    
    if(length(TIL_data$aaSeqCDR3) == 0){
        files <<- list.files(paste(dir_clones, "gDNA/CLONES_", chain, patient, "/", sep = ""))
        infusion <- files[grep('fusion', files)]
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
        TIL_data <- TIL_data[-c(grep("[*]", TIL_data$aaSeqCDR3)),]
        TIL_data <- TIL_data[-c(grep("_", TIL_data$aaSeqCDR3)),]

        # Readjusting clonal fraction values to spliced data
        ProdcloneFrac(TIL_data)
        TIL_data <- output_df
    }
    
    # Filtering the base data by the clones in the TIL data
    Base_filtered_data <- Base_data[Base_data$aaSeqCDR3 %in% TIL_data$aaSeqCDR3,]

    if(char=="clonotypes"){
        # Find the number of TIL-clonotypes in the baseline sample
        result <- nrow(Base_filtered_data)
    }
    
    if(char=="TCRs"){
        # Find the number of TIL-TCRs in the baseline sample
        result <- sum(Base_filtered_data$cloneCount)
    }
    print(result)
}