##################
# Plot Alignment #
##################

# Aligns the clone tracking plots with a single x-axis
# @param patients: list of patients for plot alignment, could be high, med, or low 
# @param sampcohort: Desired sample cohort, could be DNA, RNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param figure: figure wanting to be created (clonetrack, diversity, relative, cfDNAcorrel, clonetrack_cfDNA or rel_div)
# @param primary: desired sample to appear first (Baseline or TIL)
# @param dir_output: directory to put the png output into
# @param file_output: desired filename to output 

alignment_fig <- function(patients, sampcohort, chain, figure, primary, dir_output, file_output, clones=NULL){
  
  # Creates a list of the number of samples for each patient and finds the patient with the maximum amount of samples
  lst <- list()
  for(i in 1:length(patients)){
    lst[[i]] <- length(eval(as.name(paste(patients[i], sampcohort, "_samporder", sep=""))))
  }
  max <- as.numeric(lst[(which.max(lst))])
  #x_labs <- order[1:max]
  
  # Creates a list of the required widths of each of the patients plots based on their number of samples compared to the patient with the maximum number of samples
  widths <- list()
  for(i in 1:length(patients)){
    if(patients[i] == 'TLML_1_' & sampcohort == 'DNA'){
      widths[i] <- (length(eval(as.name(paste(patients[i], sampcohort, "_samporder", sep=""))))+1)/max
    }
    else{
      widths[i] <- (length(eval(as.name(paste(patients[i], sampcohort, "_samporder", sep="")))))/(max)
    }
  }
  widths <- c(unlist(widths))
  
  width_df <- as.data.frame(widths)
  width_df$patients <- patients

  # Assigns patient name to the variable storing the x-axis deprived plot for each patient in the inputted list
  pltlst <- list()
  for(i in 1:length(patients)){
    if(figure=="clonetrack"){
      ClonetrackPlot(width_df$patients[i], sampcohort, chain, primary)
    }
    if(figure=="clonetrack_cfDNA"){
      cfDNA_clonetrack(width_df$patients[i], sampcohort, chain, primary)
    }
    if(figure=="cfDNAcorrel"){
      cfDNA_correl(width_df$patients[i], chain, clones)
      max <- 4.696
    }
    if(figure=="diversity"){
      DivPlot(width_df$patients[i], sampcohort, chain, primary, 500)
    }
    if(figure=="relative"){
      RelPlot(width_df$patients[i], sampcohort, chain, primary)
    }
    if(figure=="rel_div"){
      if(sampcohort == "DNA") {
        Overlay_RelDiv(width_df$patients[i], sampcohort, chain, primary, 500)
      }
      else {
        Overlay_RelDiv(width_df$patients[i], sampcohort, chain, primary, 50)
      }
    }
    
    myp1 <- myp + theme(axis.text.x=element_blank(), axis.title.y = element_blank(), axis.text.y=element_blank())
    
    assign(width_df$patients[i], myp1)
    pltlst[[i]] <- eval(as.name(width_df$patients[i]))
  }
  
  drawplot_command <- "draw_plot(eval(as.name(width_df$patients[X])), 0, (Z-X)/Z, width_df$widths[X], 1/Z)"
  drawplot_command <- chartr("Z", as.character(length(patients)), drawplot_command)
  commands <- "all_plots <- ggdraw()"
  for(i in 1:length(patients)){
    commands <- paste(commands, "+", chartr("X", as.character(i), drawplot_command))
  }
  eval(parse(text=commands))
  all_plots <<- all_plots
  # The final plot is printed as a png
  png(file = paste(dir_output, file_output, ".png", sep=""),
      width = 500*max,
      height = 800*length(patients),
      res=300)
  print(all_plots)
  dev.off()
}
