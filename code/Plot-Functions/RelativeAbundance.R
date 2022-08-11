######################
# Relative Abundance #
######################

# Plots and tracks the relative abundance of large and hyperexpansive clones
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param primary: desired sample to appear first (Baseline or TIL)

RelPlot <- function(patient, sampcohort, chain, clnefrc, primary){

    # Loads patient data
    Load_data(patient, sampcohort, chain, clnefrc)
    
    # Creates longitudinally ordered sample list
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
      
    # Creates outline of relative abundance dataframe
    Rel_df <- data.frame(matrix(NA, nrow=length(samporder), ncol=3))
    Rel_df$X1 <- samporder
    colnames(Rel_df) <- c("Filename", "RelativeHyper", "RelativeExp")
    
    # Fills in relative abundance dataframe with amount of large and hyperexpansive clones
    for(i in 1:length(unique(CDR3_fraction$filename))){
        CDR3_df <- CDR3_fraction[which(CDR3_fraction$filename==samporder[i]),]
        Rel_df$RelativeHyper[i] <- sum(CDR3_df[which(CDR3_df$cloneFraction>=0.05),]$cloneFraction)
        Rel_df$RelativeExp[i] <- sum(CDR3_df[which(CDR3_df$cloneFraction<0.05 & CDR3_df$cloneFraction>0.005),]$cloneFraction)
    }
    
    # Reorders dataframe with desired primary sample
     if(primary=="TIL"){
        samporder <- samporder[c(which(grepl("infusion", samporder)==TRUE), which(grepl("infusion", samporder)==FALSE))]  
    }
    
    # Adjusts dataframe for plot production
    Rel_df <- rbind(
        data.frame(samporder, "count" = Rel_df$RelativeExp, "type" = "Exp"),
        data.frame(samporder, "count" = Rel_df$RelativeHyper, "type" = "Hyper")
    )
    
    
    # Orders filenames based on chronological sample order
    Rel_df$samporder <- factor(Rel_df$samporder, levels = c(samporder))
    levels(Rel_df$samporder) <- c(samporder)
  
    Rel_df <<- Rel_df
    
    
    # Creates relative abundance 'myp' plot global variable
    p <<- ggplot(Rel_df, aes(x=samporder, y=count, fill=type, alpha=0.5)) +
        geom_bar(stat="identity", width=1/3) + scale_fill_manual(values = c("steelblue","goldenrod")) +
        scale_y_continuous(limits=c(0,1.0)) 
    
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
