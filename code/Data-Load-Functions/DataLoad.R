################
# Loading Data #
################

# Loading specific chain, sample cohort, and clonefrac for a specific patient's samples
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param dir_data: parent directory where data is stored
dir_data <- "/Users/cameronkerr/CapTCRRepo/CapTCR-TIL-Tracking/data/"
GeneralDataLoad <- function(dir_data){
  # Loading in MiXCR output file
  MiXCR_output <<- read_xlsx(paste0(dir_data,"mixcr_output.xlsx"))
  
  # Loading data from the sample keys in order to convert file names into readable format
  samplekeys <- read_excel(paste(dir_data, "TLML_samples_keys.xlsx", sep=""))
  for(i in 1:nrow(samplekeys)){
    samplekeys$Informatics_Name[i] <- sub("*genomic*", "DNA", samplekeys$Informatics_Name[i]) 
  }
  samplekeys <<- samplekeys
  
  # Loading data from VJ usage analysis
  clonotypes <<- read_csv(file.path(dir_data , "VJTreemaps.csv"))
  cassette_stats <<- read_csv(file.path(dir_data , "VJUsage.csv"))
  
  # Loading in the CDR3 colors dataset
  CDR3_colors <- read_csv(file.path(CDR3_colors_datapath,CDR3_colors_filename))
}

Load_data <- function(patient, sampcohort, chain, clnefrc){
  compldfle <- as.data.frame(MiXCR_output[which(MiXCR_output$patient==patient & MiXCR_output$cohort==sampcohort),])
  
  # Longitudinally orders the samples from previously assigned variables
  samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
  
  # Subseting dataframe to necessary clonal and CDR3 data
  preCDR3_fraction <- compldfle[, c("filename","aaSeqCDR3","cloneFraction", "cloneCount")] 
  
  # Removes empty clones and clones with only 1 occurence
  preCDR3_fraction <- preCDR3_fraction[which(preCDR3_fraction$cloneCount>1),]
  unproductive <- preCDR3_fraction
  
  # Removes unproductive clonotypes (NOTE: An extra clone "CDR3*" is added to solve null set problem)
  preCDR3_fraction <- preCDR3_fraction[-c(grep("[*]", c(preCDR3_fraction$aaSeqCDR3, "CDR3*"))),]
  preCDR3_fraction <- preCDR3_fraction[-c(grep("_", c(preCDR3_fraction$aaSeqCDR3, "CDR3_"))),]
  
  # Readjusting clonal fraction values to spliced data
  ProdcloneFrac(preCDR3_fraction)
  CDR3_fraction <- output_df
  
  # Removing clonotypes above defined clone fraction
  CDR3_fraction <- CDR3_fraction[(CDR3_fraction$cloneFraction > clnefrc),]
  # Number of samples
  mysamples <- unique(CDR3_fraction$filename)
  
  # Changing the sample names to readable format through referencing sample keys dataframe 'samplekeys'
  CDR3_fraction$filename <- as.character(CDR3_fraction$filename)
  for(i in 1:nrow(CDR3_fraction)){
    CDR3_fraction$filename[i] <- sub("*genomic*", "DNA", CDR3_fraction$filename[i])
  }
  h <- 1
  for (f in CDR3_fraction$filename){
    j <- 1
    for(i in samplekeys$Informatics_Name){
      if(grepl(i, f) == TRUE){
        file = paste(samplekeys$Timepoint[j],
                     samplekeys$Sample_Year[j],
                     samplekeys$Sample_Month[j], sep="_")
        CDR3_fraction$filename[h] = file
      }
      j <- j+1
    }
    h <- h + 1
  }
  if(patient=="TLML_26"){
    CDR3_fraction <- CDR3_fraction[-c(which(CDR3_fraction$filename=='BL_2015_1')),]
  }
  
  CDR3_fraction <<- CDR3_fraction
  
  # Creating custom dataframes for the TIL Infusion product, 4W sample, and baseline sample
  TIL_data <<- dplyr::filter(CDR3_fraction, grepl('infusion', filename))
  FW_data <<- dplyr::filter(CDR3_fraction, grepl('4W', filename))
  Base_data <<- CDR3_fraction[which(CDR3_fraction$filename==samporder[1]),]
}
