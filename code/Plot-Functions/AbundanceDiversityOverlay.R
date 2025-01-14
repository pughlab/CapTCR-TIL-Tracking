###################################
# Overlay Abundance and Diversity #
###################################

# Overlays the relative abundance and diversity plots
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be DNA, RNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param primary: desired sample to appear first (Baseline or TIL)

Overlay_RelDiv <- function(patient, sampcohort, chain, primary, max){

  # Loads diversity and relative abundance plots and creates dataframe combining the data from the 2
  DivPlot(patient, sampcohort, chain, primary, max)
  RelPlot(patient, sampcohort, chain, primary)
  #Rel_df$count <- Rel_df$count*100
  RelDiv_df <- cbind(Rel_df, Div_df)
  # Creates plot overalying the 2 with a double y-axis and colouring the second y-axis red 
  if(patient == "TLML_1_" & sampcohort == "DNA"){
    p <- ggplot(RelDiv_df, aes(x=factor(samporder, c("baseline_2013_9", "infusion_2013_8", "FU_01_2014_1", "FU_02_2014_4", "FU_03_2014_7")))) +
        scale_x_discrete(limits=c('baseline_2013_9', 'infusion_2013_8', '4 week sample', 'FU_01_2014_1', 'FU_02_2014_4', 'FU_03_2014_7'), labels=c('baseline_2013_9', 'infusion_2013_8', '4 week sample', 'FU_01_2014_1', 'FU_02_2014_4', 'FU_03_2014_7'))
  }
  else{
    p <- ggplot(RelDiv_df)
  }
  myp <<- p + geom_bar(aes(x=samporder, y=count, fill=type, alpha=0.5), stat="identity", width=1/3) + scale_fill_manual(values = c("steelblue","goldenrod")) +
          scale_y_continuous(labels=scales::percent, limits=c(0,1.0), sec.axis=sec_axis(~.*max)) + ylab(patient) +
          theme(axis.title.y = element_text(size = 13, angle=0, vjust=0.5),
                axis.title.x = element_blank(),
                axis.line = element_line(color = "black"),
                axis.text = element_text(size = 13),
                axis.text.x = element_text(angle = 50, hjust = 1)) +
          theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            legend.position = "blank",
            plot.margin = unit(c(0.2,0,0,0),"cm")) + 
          geom_line(aes(x=Filename, y=Diversity/max, group=1), stat="identity", color="#BF5549", size=1.5) +
          geom_point(aes(x=Filename, y=Diversity/max, shape=factor(Shape)), colour="#BF5549", size = 4.0)
}
