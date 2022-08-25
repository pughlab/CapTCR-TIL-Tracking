#############################
# Inverse Simpson Diversity #
#############################

# Plots and tracks the inverse simpson diversity
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be DNA, RNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param primary: desired sample to appear first (Baseline or TIL)
# @param max: desired maximum diversity value, outliers will be shown as a triangle

DivPlot <- function(patient, sampcohort, chain, primary, max){
    
    # Loading CDR3_fraction
    CDR3_fraction <- eval(as.name(paste0(patient, sampcohort)))
  
    # Setting the longitudinal order of the samples for patient
    samporder <- eval(as.name(paste0(patient, sampcohort, "_samporder")))
    
    # Creating outline for diversity dataframe
    Div_df <- data.frame(matrix(NA, nrow=length(samporder), ncol=3))
    Div_df$X1 <- samporder
    colnames(Div_df) <- c("Filename", "Diversity", "Shape")
    
    # Filling in dataframe values with previously loaded data
    for(i in 1:length(unique(CDR3_fraction$filename))){
        CDR3_df <- CDR3_fraction[which(CDR3_fraction$filename==samporder[i]),]
        N <- sum(CDR3_df$cloneCount)
        sum_simpson <- 0
        for(n in CDR3_df$cloneCount){
          sum_simpson <- sum_simpson + (n/N)^2
        }
        Div_df$Diversity[i] <- 1/sum_simpson
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
    
    Div_df <<- Div_df  
  
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
