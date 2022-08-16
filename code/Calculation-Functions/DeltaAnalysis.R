##################
# Delta Analysis #
##################

# Calculates the delta between 2 data points of different samples
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param datapoint: The data point wanted to be analyzed (Expansion, Productive, Expression, Diversity)
# @param first: The name of the first sample
# @param second: The name of the second sample

DeltaA <- function(patient, sampcohort, chain, clnefrc, datapoint, first, second){

    # Loads in patient samporder
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
    
    # Loads in patient data
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
    
    # Changes the name of the baseline to the first name in the sample order previously creates
    if(first=='BL'){
        first <- eval(as.name(paste(patient, sampcohort, sep="")))[1]
    }
    
    # Loads the clone data for the first and second identified samples
    df_first <- CDR3_fraction[which(grepl(first, CDR3_fraction$filename)==TRUE),]
    df_second <- CDR3_fraction[which(grepl(second, CDR3_fraction$filename)==TRUE),]

    # Finds the delta between the TIL clone fraction in the samples if user specified
    if(datapoint=="Expansion"){
        data_first <- sum(df_first[df_first$aaSeqCDR3 %in% TIL_data$aaSeqCDR3,]$cloneFraction)
        data_second <- sum(df_second[df_second$aaSeqCDR3 %in% TIL_data$aaSeqCDR3,]$cloneFraction)
        delta <- data_second - data_first
    }

    # Finds the delta between the productive clone fractions in the samples if user specified 
    if(datapoint=="Productive"){
        unproductive$filename <- as.character(unproductive$filename)
        # Changes the names of the file names in the 'CDR3_fraction' variable defined in the 'Load_data' function to readable format using the sample keys
        for(i in 1:nrow(unproductive)){
            unproductive$filename[i] <- sub("*genomic*", "DNA", unproductive$filename[i])
        }
        h <- 1
        for (f in unproductive$filename){
            j <- 1
            for(i in samplekeys$Informatics_Name){
                if(grepl(i, f) == TRUE){
                    file <- paste(samplekeys$Timepoint[j],
                                samplekeys$Sample_Year[j],
                                samplekeys$Sample_Month[j], sep="_")
                    test <- file
                    unproductive$filename[h] <- file
                }
                j <- j+1
            }
            h <- h + 1
        }
        
        df_first <- unproductive[which(grepl(first, unproductive$filename)==TRUE),]
        df_second <- unproductive[which(grepl(second, unproductive$filename)==TRUE),]
        
        data_first <- length(df_first$aaSeqCDR3[!grepl("[*]", df_first$aaSeqCDR3) & !grepl("_", df_first$aaSeqCDR3)])/length(df_first$aaSeqCDR3)
        data_second <- length(df_second$aaSeqCDR3[!grepl("[*]", df_second$aaSeqCDR3) & !grepl("_", df_second$aaSeqCDR3)])/length(df_second$aaSeqCDR3)
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