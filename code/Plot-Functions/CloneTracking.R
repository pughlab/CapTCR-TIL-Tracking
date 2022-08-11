##################
# Clone Tracking #
##################

# Plots and tracks the redundant TIL clones
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param primary: desired sample to appear first (Baseline or TIL)

ClonetrackPlot <- function(patient, sampcohort, chain, clnefrc, primary){
  
  # Loading the data for the patient for a specific sampcohort, chain, and clone fraction
  Load_data(patient, sampcohort, chain, clnefrc)
  
  # Assigning the sample order to the predefined lists created in the data load section
  samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
  
  # Pulls relative gDNA TIL Infusion product if the sample does not have an infusion
  if(length(TIL_data$aaSeqCDR3) == 0){
    files <- list.files(paste(dir_clones ,"gDNA/CLONES_", chain, patient, "/", sep = ""))
    infusion <- files[grep('fusion', files)]
    len <- length(samporder)
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
    # Generating dataframe with added gDNA TIL Infusion product and updated y-axis order for plots
    CDR3_fraction <- rbind(Base_data, TIL_data, FW_data)
    samporder <- c(samporder[1],'TIL Infusion Product',samporder[2:len])
  }
  
  # Generates list of clones and assigns them to previously created distinct colour pallete
  CDR3_colors <- CDR3_colors[which(CDR3_colors$id==patient),]
  CDR3_presentcolors <- CDR3_colors[which(CDR3_colors$colored_clns %in% CDR3_fraction$aaSeqCDR3),]
  colored_clns <- CDR3_presentcolors$colored_clns
  non_clns <- CDR3_fraction$aaSeqCDR3[!CDR3_fraction$aaSeqCDR3 %in% colored_clns]
  message("Total number of coloured clonotypes:")
  print(length(colored_clns))
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
  if(patient == "TLML_1_" & sampcohort == "gDNA"){
    p <- ggplot(data=CDR3_fraction, aes(x=factor(filename, level=c('apheresis_2013_9', 'infusion_2013_8', 'FU_01_2014_1', 'FU_02_2014_4', 'FU_03_2014_7')), y=cloneFraction, fill=aaSeqCDR3, stratum=aaSeqCDR3, alluvium=aaSeqCDR3)) +
      geom_alluvium(decreasing=FALSE) +
      geom_stratum(decreasing = FALSE, stat = "alluvium") + 
      scale_fill_manual(breaks = names(mycolors[mycolors != "white"]), values=mycolors) +
      scale_x_discrete(limits=c('apheresis_2013_9', 'infusion_2013_8', '4 week sample', 'FU_01_2014_1', 'FU_02_2014_4', 'FU_03_2014_7'), labels=c('apheresis_2013_9', 'infusion_2013_8', '4 week sample', 'FU_01_2014_1', 'FU_02_2014_4', 'FU_03_2014_7'))
  }
  else {
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
