richness_boxplot <- function(GitHub_path) {
  # Get sample names ------------------- 
  all_variables <- ls(pattern = "samporder", envir=.GlobalEnv)
  all_variables <- c(all_variables[11:13], all_variables[14:16], all_variables[17:19], all_variables[1:2], all_variables[3:5], all_variables[23:24], all_variables[8:10], all_variables[20:22], all_variables[6:7])
  DNAdiversity <- c()
  RNAdiversity <- c()
  cfDNAdiversity <- c()
  for(x in all_variables){
    for(samp in eval(as.name(x))) {
      split_samp <- strsplit(x, "_")[[1]]
      patient <- paste0(split_samp[1], "_", split_samp[2], "_")
      cohort <- split_samp[3]
      timepoint <-
      if(cohort == "DNA") {
        value <- Richness(patient, "DNA", samp)
        DNAdiversity <- c(DNAdiversity, value)
      }
      else if(cohort == "RNA") {
        value <- Richness(patient, "RNA", samp)
        RNAdiversity <- c(RNAdiversity, value)
      }
      else if(cohort == "cfDNA") {
        value <- Richness(patient, "cfDNA", samp)
        cfDNAdiversity <- c(cfDNAdiversity, value)
      }
    }
  }
  Rich_df <- data.frame("Richness" = c(DNAdiversity, RNAdiversity, cfDNAdiversity), "Cohort" = c(rep("DNA", length(DNAdiversity)), rep("RNA", length(RNAdiversity)), rep("cfDNA", length(cfDNAdiversity))))
  
  png(file = paste(GitHub_path, "data/boxplot.png", sep=""),
      width = 564,
      height = 1154,
      res=300)
  print(ggplot(Rich_df, aes(x=Cohort, y=Richness, fill=Cohort)) + scale_y_log10() + 
    geom_boxplot(outlier.shape=NA) + theme_classic() + 
    geom_jitter(color="black", size=1.5, alpha=0.9, width=0.075) + 
    theme(legend.position = "none", plot.title = element_blank(), 
          axis.text = element_blank(), axis.title = element_blank()) + 
    scale_fill_manual(values=c("#FDC086", "#7FC97F", "#BEAED4")))
  dev.off()
}
