##################
# Clone Tracking #
##################

# Plots and tracks the redundant TIL clones
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be DNA, RNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param primary: desired sample to appear first (Baseline or TIL)

ClonetrackPlot <- function(patient, sampcohort, chain, primary){
  
  # Loading CDR3_fraction
  CDR3_fraction <- eval(as.name(paste0(patient, sampcohort)))
  
  # Setting the longitudinal order of the samples for patient
  samporder <- eval(as.name(paste0(patient, sampcohort, "_samporder")))
  
  # Pulls relative gDNA TIL Infusion product if the sample does not have an infusion
  if(length(grep("infusion", eval(as.name(paste0(patient, sampcohort, "_samporder"))))) == 0){
    DNA_infusion <- eval(as.name(paste0(patient, "DNA", "_samporder")))[2]
    CDR3_fraction <- rbind(CDR3_fraction, eval(as.name(paste0(patient, "DNA", "_", DNA_infusion))))
    samporder <- c(samporder[1], DNA_infusion,samporder[2:length(samporder)])
  }
  
  # Generates list of clones and assigns them to previously created distinct colour pallete
  CDR3_colors <- CDR3_colors[which(CDR3_colors$id==patient),]
  CDR3_presentcolors <- CDR3_colors[which(CDR3_colors$colored_clns %in% CDR3_fraction$aaSeqCDR3),]
  colored_clns <- CDR3_presentcolors$colored_clns
  non_clns <- CDR3_fraction$aaSeqCDR3[!CDR3_fraction$aaSeqCDR3 %in% colored_clns]
  mycolors <- CDR3_presentcolors$mycolors
  mycolors <- c(mycolors, rep("white", length(non_clns)))
  names(mycolors) <- c(colored_clns, non_clns)
  
  # Reorders the sample if the user specifies the TIL Infusion product to appear first
  if(primary=="TIL"){
    samporder <- samporder[c(which(grepl("infusion", samporder)==TRUE), which(grepl("infusion", samporder)==FALSE))]  
  }
  
  # Orders filenames based on chronological sample order
  CDR3_fraction$filename <- factor(CDR3_fraction$filename, levels = c(samporder))
  levels(CDR3_fraction$filename) <- c(samporder)
  CDR3_fraction <<- CDR3_fraction
  
  # Creates plot for the data frame and list of colors and associated CDR3 sequences
  if(patient == "TLML_1_" & sampcohort == "DNA"){
    p <- ggplot(data=CDR3_fraction, aes(x=factor(filename, level=c('baseline_2013_9', 'infusion_2013_8', 'FU_01_2014_1', 'FU_02_2014_4', 'FU_03_2014_7')), y=cloneFraction, fill=aaSeqCDR3, stratum=aaSeqCDR3, alluvium=aaSeqCDR3)) +
      geom_alluvium(decreasing=FALSE) +
      geom_stratum(decreasing = FALSE, stat = "alluvium") + 
      scale_fill_manual(breaks = names(mycolors[mycolors != "white"]), values=mycolors) +
      scale_x_discrete(limits=c('baseline_2013_9', 'infusion_2013_8', '4 week sample', 'FU_01_2014_1', 'FU_02_2014_4', 'FU_03_2014_7'), labels=c('baseline_2013_9', 'infusion_2013_8', '4 week sample', 'FU_01_2014_1', 'FU_02_2014_4', 'FU_03_2014_7'))
  } else {
    p <- ggplot(CDR3_fraction, aes(x = filename, 
                                   y = cloneFraction,
                                   fill = aaSeqCDR3,
                                   stratum = aaSeqCDR3,
                                   alluvium = aaSeqCDR3)) + 
      geom_alluvium(decreasing = FALSE) + 
      geom_stratum(decreasing = FALSE, stat = "alluvium") + 
      scale_fill_manual(breaks = names(mycolors[mycolors != "white"]),values = mycolors)
  }      
  
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
