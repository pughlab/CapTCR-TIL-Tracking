############################
# cfDNA Clone Correlation #
###########################

# Plots the correlation between the cfDNA and gDNA repertoires 
# @param patients: list of patients for plot alignment, could be high, med, or low 
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param clones: Desired clones to analyze (TIL or Background)

# Creating R^2 function
rsq <- function (x, y) cor(x, y) ^ 2

cfDNA_correl <- function(patient, chain, clnefrc, clones){
  
  # Loading both gDNA and cfDNA data and assigning the 4 week and baseline data to differing data sets
  Load_data(patient, "gDNA", chain, clnefrc)
  gDNA_FW <- FW_data
  gDNA_Base <- Base_data
  gDNA_Total <- CDR3_fraction
  Load_data(patient, "cfDNA", chain, clnefrc)
  cfDNA_FW <- FW_data
  cfDNA_Base <- Base_data
  cfDNA_Total <- CDR3_fraction
  
  # Getting a list of redundant TIL clones and their colours
  TotalCDR3 <- unique(c(cfDNA_Total$aaSeqCDR3, gDNA_Total$aaSeqCDR3))
  
  # Creating the longitudinal dataframe containing all aaSeqCDR3s and their respective clonal fraction in both sample cohrots and timepoints
  LongitudinalComparison <- data.frame(matrix(nrow=length(TotalCDR3), ncol=5, data=NA))
  colnames(LongitudinalComparison) <- c("aaSeqCDR3", "gDNA_FW_cloneFraction", "cfDNA_FW_cloneFraction", 
                                        "gDNA_Base_cloneFraction", "cfDNA_Base_cloneFraction")
  LongitudinalComparison$aaSeqCDR3 <- TotalCDR3
  LongitudinalComparison$gDNA_FW_cloneFraction <- gDNA_FW$cloneFraction[match(LongitudinalComparison$aaSeqCDR3, gDNA_FW$aaSeqCDR3)]
  LongitudinalComparison$gDNA_Base_cloneFraction <- gDNA_Base$cloneFraction[match(LongitudinalComparison$aaSeqCDR3, gDNA_Base$aaSeqCDR3)]
  LongitudinalComparison$cfDNA_FW_cloneFraction <- cfDNA_FW$cloneFraction[match(LongitudinalComparison$aaSeqCDR3, cfDNA_FW$aaSeqCDR3)]
  LongitudinalComparison$cfDNA_Base_cloneFraction <- cfDNA_Base$cloneFraction[match(LongitudinalComparison$aaSeqCDR3, cfDNA_Base$aaSeqCDR3)]
  LongitudinalComparison[is.na(LongitudinalComparison)] <- 0
  
  # Creating colors dataframe containing the colors for the patient
  CDR3_colors_id <- CDR3_colors[which(CDR3_colors$id==patient),]
  col <- as.character(CDR3_colors_id$mycolors)
  names(col) <- as.character(CDR3_colors_id$colored_clns)
    
  # Creating Creating a dataframe for each timepoint
  Comparison_FWdata <- LongitudinalComparison[1:3]
  colnames(Comparison_FWdata) <- c("aaSeqCDR3", "gDNA_cloneFraction", "cfDNA_cloneFraction")
  Comparison_FWdata$Timepoint <- 2
  Comparison_FWdata$Color <- CDR3_colors_id$mycolors[match(Comparison_FWdata$aaSeqCDR3, CDR3_colors_id$colored_clns)]
  
  Comparison_Basedata <- LongitudinalComparison[-c(2:3)]
  colnames(Comparison_Basedata) <- c("aaSeqCDR3", "gDNA_cloneFraction", "cfDNA_cloneFraction")
  Comparison_Basedata$Timepoint <- 1
  Comparison_Basedata$Color <- CDR3_colors_id$mycolors[match(Comparison_Basedata$aaSeqCDR3, CDR3_colors_id$colored_clns)]
  
  # Seperating dataframe into only TILs or the background clones
  if(clones=="TIL"){
    Comparison_FWdata <- na.omit(Comparison_FWdata)
    Comparison_Basedata <- na.omit(Comparison_Basedata)
  }
  if(clones=="Background"){
    Comparison_FWdata <- Comparison_FWdata[is.na(Comparison_FWdata$Color),]
    Comparison_Basedata <- Comparison_Basedata[is.na(Comparison_Basedata$Color),]
  }
 
  # Calculating R^2 values for the baseline and 4 week samples
  FW_R2 <- rsq(Comparison_FWdata$gDNA_cloneFraction, Comparison_FWdata$cfDNA_cloneFraction)
  Base_R2 <- rsq(Comparison_Basedata$gDNA_cloneFraction, Comparison_Basedata$cfDNA_cloneFraction)
  
  # Creting scatter plot with R^2 value on it and assigning appropriate themes and log scale
  p <- ggplot(Comparison_FWdata, aes(x=log10(gDNA_cloneFraction), y=log10(cfDNA_cloneFraction))) + geom_point(aes(size=2,colour=aaSeqCDR3)) + scale_color_manual(values=col)
  myp_FW <- p + geom_text(size=7.5,x=-3, y=-1, label=paste0("~R^{2}: ", round(FW_R2, digits=3)), parse=TRUE) +
    coord_cartesian(ylim=c(-4, 0), xlim=c(-4, 0)) +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_blank(),
          axis.text.x = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "blank",
          plot.margin = unit(c(0.2,0,0,0),"cm"))
  p <- ggplot(Comparison_Basedata, aes(x=log10(gDNA_cloneFraction), y=log10(cfDNA_cloneFraction))) + geom_point(aes(size=0.25,colour=aaSeqCDR3)) + scale_color_manual(values=col)
  myp_Base <- p + geom_text(size=7.5, x=-3, y=-1, label=paste0("~R^{2}: ", round(Base_R2, digits=3)), parse=TRUE) +
    coord_cartesian(ylim=c(-4, 0), xlim=c(-4, 0)) +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_blank(),
          axis.text.x = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "blank",
          plot.margin = unit(c(0.2,0,0,0),"cm"))

  myp <<- plot_grid(myp_Base, myp_FW, ncol=2, align="h")
}
