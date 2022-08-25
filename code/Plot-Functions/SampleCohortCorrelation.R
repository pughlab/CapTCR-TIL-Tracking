#############################
# Sample Cohort Correlation #
#############################

# Plots the correlative overlap between sample cohorts (DNA, RNA, cfDNA)
# @param sampcohort1: Desired first sample cohort, could be DNA, RNA, or cfDNA
# @param sampcohort2: Desired second sample cohort, could be DNA, RNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param dir_output: directory for the output png
# @param file_output: png file name

rsq <- function (x, y) cor(x, y) ^ 2

SampleCohort_Correl <- function(sampcohort1, sampcohort2, chain, clnefrc, dir_output, file_output){
    
    #List of patients required for analysis (requires manual change)
    patients <- c("TLML_4_", "TLML_16_", "TLML_20_", "TLML_22_",
                  "TLML_26_", "TLML_29_")

    # Loops through the patients and loads the data from the 2 desired sample cohorts in one dataframe
    for(patient in patients){
        # Load fist sample cohort data
        samporder_cohort1 <- eval(as.name(paste0(patient, sampcohort1, "_samporder")))
        Sample1_fraction <- eval(as.name(paste0(patient, sampcohort1, "_", samporder_cohort1[3])))[c("aaSeqCDR3", "cloneFraction")]
        Sample1_fraction$clone2_Fraction <- NA
        
        Sample2_fraction <- eval(as.name(paste0(patient, sampcohort2, "_", samporder_cohort1[3])))[c("aaSeqCDR3", "cloneFraction")]
        Sample1_fraction$clone2_Fraction <- Sample2_fraction$cloneFraction[match(Sample1_fraction$aaSeqCDR3, Sample2_fraction$aaSeqCDR3)]
        #Removing all rows with NA values
        Sample1_fraction <- na.omit(Sample1_fraction)
        
        assign(paste(patient, "merged", sep=""), Sample1_fraction)
    }

    # Binds all of the patient's data together
    Ov_data <- rbind(TLML_4_merged, TLML_16_merged, TLML_20_merged, TLML_22_merged, TLML_26_merged, TLML_29_merged)

    # Calculates and prints R^2 value
    print(rsq(Ov_data$cloneFraction, Ov_data$clone2_Fraction))
    
    # Creates plot assigned to 'myp' global variable
    p <- ggplot(Ov_data, aes(x=cloneFraction, y=clone2_Fraction)) + geom_point(size=2, color="steelblue")
    myp <<- p + 
        scale_y_continuous(n.breaks=4, trans="log10") + scale_x_continuous(n.breaks=4, trans="log10") + geom_smooth(method=lm) +
        theme(axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text = element_text(size = 13),
              axis.text.x = element_text(angle = 50, hjust = 1)) +
        theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "blank",
          plot.margin = unit(c(0.2,0,0,0),"cm"))
 
    # The final plot is printed as a png
    png(file = paste(dir_output, file_output, ".png", sep=""),
        width = 3000,
        height = 3000,
        res=300)
   print(myp)
   dev.off() 
} 
