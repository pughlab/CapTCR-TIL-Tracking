#############################
# Inverse Simpson Diversity #
#############################

# Plots and tracks the inverse simpson diversity
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param dir_clones: parent directory where clone files are located
# @param dir_samplekeys: directory where the sample keys excel file are located
# @param file_samplekeys: file name of the sample keys 
# @param primary: desired sample to appear first (Baseline or TIL)
# @param max: desired maximum diversity value, outliers will be shown as a triangle

DivPlot <- function(patient, sampcohort, chain, clnefrc, dir_clones, 
                    dir_samplekeys, file_samplekeys, primary, max){
    
    # Loading in patient data
    Load_data(patient, sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys)
    
    # Setting the longitudinal order of the samples for patient
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
    
    # Creating outline for diversity dataframe
    Div_df <- data.frame(matrix(NA, nrow=length(samporder), ncol=3))
    Div_df$X1 <- samporder
    colnames(Div_df) <- c("Filename", "Diversity", "Shape")
    
    # Filling in dataframe values with previously loaded data
    for(i in 1:length(unique(CDR3_fraction$filename))){
        CDR3_df <- CDR3_fraction[which(CDR3_fraction$filename==samporder[i]),]
        CDR3_df <- CDR3_df[2:4]
        colnames(CDR3_df) <- c("CDR3.aa", "Proportion", "Clones")
        Div_df$Diversity[i] <- repDiversity(CDR3_df, .method="inv.simp")
        if(Div_df$Diversity[i]>max){
            Div_df$Shape[i] <- 17
            Div_df$Diversity[i] <- max
        }
        else{
            Div_df$Shape[i] <- 16
        }
    }

    # Reorders the order of samples based on desired primary sample
     if(primary=="TIL"){
        samporder <- samporder[c(which(grepl("infusion", samporder)==TRUE), which(grepl("infusion", samporder)==FALSE))]  
    }
    
    # Orders filenames based on chronological sample order
    Div_df$Filename <- factor(Div_df$Filename, levels = c(samporder))
    levels(Div_df$Filename) <- c(samporder)
    
    # Creates plot assigned to 'myp' global variable
    p <- ggplot(Div_df) + geom_line(aes(x=Filename, y=Diversity, group=1), stat="identity", color="#E7B800", size=1.5) +
                          geom_point(aes(x=Filename, y=Diversity, shape=factor(Shape)), colour="#E7B800", size = 4.0)+
                          scale_y_continuous(limits=c(0,max))
    myp <<- p + ylab(patient) +
        theme(axis.title.y = element_text(size = 13, angle=0, vjust=0.5),
              axis.title.x = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text = element_text(size = 13),
              axis.text.x = element_text(angle = 50, hjust = 1)) +
        theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "blank",
          plot.margin = unit(c(0.2,0,0,0),"cm"))
}