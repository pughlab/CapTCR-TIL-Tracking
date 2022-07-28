###################################
# Overlay Abundance and Diversity #
###################################

# Plots and tracks the redundant TIL clones
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param dir_clones: parent directory where clone files are located
# @param dir_samplekeys: directory where the sample keys excel file are located
# @param file_samplekeys: file name of the sample keys 
# @param primary: desired sample to appear first (Baseline or TIL)

Overlay_RelDiv <- function(patient, sampcohort, chain, clnefrc, dir_clones, 
                           dir_samplekeys, file_samplekeys, primary, max){

# Loads diversity and relative abundance plots and creates dataframe combining the data from the 2
DivPlot(patient, sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys, primary, max)
RelPlot(patient, sampcohort, chain, clnefrc, dir_clones, dir_samplekeys, file_samplekeys, primary)
RelDiv_df <- cbind(Rel_df, Div_df)
# Creates plot overalying the 2 with a double y-axis and colouring the second y-axis red  
myp <<- ggplot(RelDiv_df) +
        geom_bar(aes(x=samporder, y=count, fill=type, alpha=0.5), stat="identity", width=1/3) + scale_fill_manual(values = c("steelblue","goldenrod")) +
        scale_y_continuous(limits=c(0,1.0), sec.axis=sec_axis(~.*500)) + ylab(patient) +
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
        geom_line(aes(x=Filename, y=Diversity/500, group=1), stat="identity", color="#BF5549", size=1.5) +
        geom_point(aes(x=Filename, y=Diversity/500, shape=factor(Shape)), colour="#BF5549", size = 4.0)
}