####################
# Top 50% Analysis #
####################

# Finds the amount of a sample taken up by the top 50% of another sample
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param Top50: sample which will have its top 50% compared
# @param Reference: sample which will be compared to 

Top50.fx <- function(patient, sampcohort, chain, clnefrc, Top50, Reference){
    
    # Loads patient sample order
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))

    # Loads patient clone data
    Load_data(patient, sampcohort, chain, clnefrc)

    #Creating dataframes for the clones of the Top50 and Reference samples
    Top50Samp <- CDR3_fraction[which(CDR3_fraction$filename==grep(Top50, samporder, value=TRUE)),]
    ReferenceSamp <- CDR3_fraction[which(CDR3_fraction$filename==grep(Reference, samporder, value=TRUE)),]
    
    # Creating dataframe showing both the clone fractions of the top 50 and reference sample datasets
    Reference_clonefraction <- NA
    Top50Data <- cbind(Top50Samp, Reference_clonefraction)
    Top50Data$Reference_clonefraction <- ReferenceSamp$cloneFraction[match(Top50Data$aaSeqCDR3, ReferenceSamp$aaSeqCDR3)]

    #Converting all NA values to 0
    Top50Data[is.na(Top50Data)] <- 0

    # Taking the top 50% of the Top50Data
    Top50Spliced <- data.frame()
    f <- 1
    Top50Frc <- 0
    for(i in Top50Data$cloneFraction){
        if(Top50Frc < 0.5){
            Top50Frc <- Top50Frc + i
            Top50Spliced <- rbind(Top50Spliced, Top50Data[f,])
            f <- f + 1
        }
    }

    # Calculating the fraction of the reference sample taken up by the clones in the top 50% of the Top50 sample
    ReferenceFrc <<- sum(Top50Spliced$Reference_clonefraction)
}