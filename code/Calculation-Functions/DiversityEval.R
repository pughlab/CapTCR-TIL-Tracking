########################
# Diversity Evaluation #
########################

# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes

DivCalc <- function(patient, sampcohort, chain, clnefrc){
    
    # Loading in patient data
    Load_data(patient, sampcohort, chain, clnefrc)
    
    # Setting the longitudinal order of the samples for patient
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
    
    # Creating outline for diversity dataframe
    Div_df <- data.frame(matrix(NA, nrow=length(samporder), ncol=2))
    Div_df$X1 <- samporder
    colnames(Div_df) <- c("Filename", "Diversity")
    
    # Filling in dataframe values with previously loaded data
    for(i in 1:length(unique(CDR3_fraction$filename))){
        CDR3_df <- CDR3_fraction[which(CDR3_fraction$filename==samporder[i]),]
        N <- sum(CDR3_df$cloneCount)
        sum_simpson <- 0
        for(n in CDR3_df$cloneCount){
          sum_simpson <- sum_simpson + (n/N)^2
        }
        Div_df$Diversity[i] <- 1/sum_simpson
    }
    
    # Orders filenames based on chronological sample order
    Div_df$Filename <- factor(Div_df$Filename, levels = c(samporder))
    levels(Div_df$Filename) <- c(samporder)
    Div_df <<- Div_df
    print(Div_df)
}


