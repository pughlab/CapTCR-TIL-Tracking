########################
# Diversity Evaluation #
########################

# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param dir_clones: parent directory where clone files are located
# @param dir_samplekeys: directory where the sample keys excel file are located
# @param file_samplekeys: file name of the sample keys 
# @param primary: desired sample to appear first (Baseline or TIL)
# @param max: desired maximum diversity value, outliers will be shown as a triangle

DivCalc <- function(patient, sampcohort, chain, clnefrc, dir_clones, 
                    dir_samplekeys, file_samplekeys){
    
    # Loading in patient data
    Load_data(patient, sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys)
    
    # Setting the longitudinal order of the samples for patient
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
    
    # Creating outline for diversity dataframe
    Div_df <- data.frame(matrix(NA, nrow=length(samporder), ncol=2))
    Div_df$X1 <- samporder
    colnames(Div_df) <- c("Filename", "Diversity")
    
    # Filling in dataframe values with previously loaded data
    for(i in 1:length(unique(CDR3_fraction$filename))){
        CDR3_df <- CDR3_fraction[which(CDR3_fraction$filename==samporder[i]),]
        CDR3_df <- CDR3_df[2:4]
        colnames(CDR3_df) <- c("CDR3.nt", "Proportion", "Clones")
        Div_df$Diversity[i] <- repDiversity(CDR3_df, .method="inv.simp")
    }
    
    # Orders filenames based on chronological sample order
    Div_df$Filename <- factor(Div_df$Filename, levels = c(samporder))
    levels(Div_df$Filename) <- c(samporder)
    
    print(Div_df)
}