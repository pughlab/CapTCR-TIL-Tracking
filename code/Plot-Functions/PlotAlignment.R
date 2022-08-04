###############
# Plot Alignment #
###############

# Aligns the clone tracking plots with a single x-axis
# @param patients: list of patients for plot alignment, could be high, med, or low 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param dir_clones: parent directory where clone files are located
# @param dir_samplekeys: directory where the sample keys excel file are located
# @param file_samplekeys: file name of the sample keys 
# @param figure: figure wanting to be created (clonetrack, diversity, relative, or rel_div)
# @param primary: desired sample to appear first (Baseline or TIL)
# @param dir_output: directory to put the png output into
# @param file_output: desired filename to output 

alignment_fig <- function(patients, sampcohort, chain, clnefrc, dir_clones, dir_samplekeys,
                          file_samplekeys, figure, primary, dir_output, file_output){
    
    # Creates a list of the number of samples for each patient and finds the patient with the maximum amount of samples
    lst <- list()
    for(i in 1:length(patients)){
        lst[[i]] <- length(eval(as.name(paste(patients[i], sampcohort, sep=""))))
    }
    max <- as.numeric(lst[(which.max(lst))])
    #x_labs <- order[1:max]
    
    # Creates a list of the required widths of each of the patients plots based on their number of samples compared to the patient with the maximum number of samples
    widths <- list()
    for(i in 1:length(patients)){
        if(patients[i] == 'TLML_1_' & sampcohort == 'gDNA'){
          widths[i] <- (length(eval(as.name(paste(patients[i], sampcohort, sep=""))))+1)/max
        }
        else{
          widths[i] <- (length(eval(as.name(paste(patients[i], sampcohort, sep="")))))/(max)
        }
    }
    widths <- c(unlist(widths))
    
    width_df <- as.data.frame(widths)
    width_df$patients <- patients
    width_df <- width_df[order(width_df$widths),]
            
    # Assigns patient name to the variable storing the x-axis deprived plot for each patient in the inputted list
    pltlst <- list()
    for(i in 1:length(patients)){
        if(figure=="clonetrack"){
            ClonetrackPlot(width_df$patients[i], sampcohort, chain, clnefrc, dir_clones, 
                           dir_samplekeys, file_samplekeys, primary)
          }
        if(figure=="diversity"){
            DivPlot(width_df$patients[i], sampcohort, chain, clnefrc, dir_clones, 
                           dir_samplekeys, file_samplekeys, primary, 500)
          }
        if(figure=="relative"){
            RelPlot(width_df$patients[i], sampcohort, chain, clnefrc, dir_clones, 
                           dir_samplekeys, file_samplekeys, primary)
        }
        if(figure=="rel_div"){
            Overlay_RelDiv(width_df$patients[i], sampcohort, chain, clnefrc, dir_clones, 
                           dir_samplekeys, file_samplekeys, primary, 500)
        }

        myp1 <- myp + theme(axis.text.x=element_blank(), axis.title.y = element_blank(), axis.text.y=element_blank())
      
        assign(width_df$patients[i], myp1)
        pltlst[[i]] <- eval(as.name(width_df$patients[i]))
    }
    
    # For each possible length of list of patients, the plots are aligned on the same x-axis plot  
    if(length(patients)==5){
        all_plots <- ggdraw() + 
            draw_plot(eval(as.name(width_df$patients[1])), 0, 4/5, width_df$widths[1], 1/5) + 
            draw_plot(eval(as.name(width_df$patients[2])), 0, 3/5, width_df$widths[2], 1/5) +
            draw_plot(eval(as.name(width_df$patients[3])), 0, 2/5, width_df$widths[3], 1/5) +
            draw_plot(eval(as.name(width_df$patients[4])), 0, 1/5, width_df$widths[4], 1/5) + 
            draw_plot(eval(as.name(width_df$patients[5])), 0, 0, width_df$widths[5], 1/5) 
    }
    if(length(patients)==4){
        all_plots <- ggdraw() + 
            draw_plot(eval(as.name(width_df$patients[1])), 0, 0.75, width_df$widths[1], 0.25) + 
            draw_plot(eval(as.name(width_df$patients[2])), 0, 0.5, width_df$widths[2], 0.25) +
            draw_plot(eval(as.name(width_df$patients[3])), 0, 0.25, width_df$widths[3], 0.25) +
            draw_plot(eval(as.name(width_df$patients[4])), 0, 0, width_df$widths[4], 0.25)
    }
    if(length(patients)==3){
        all_plots <- ggdraw() + 
            draw_plot(eval(as.name(width_df$patients[1])), 0, 2/3, width_df$widths[1], 1/3, scale=1) + 
            draw_plot(eval(as.name(width_df$patients[2])), 0, 1/3, width_df$widths[2], 1/3, scale=1) +
            draw_plot(eval(as.name(width_df$patients[3])), 0, 0, width_df$widths[3], 1/3, scale=1)
    }
    if(length(patients)==2){
        all_plots <- ggdraw() + 
            draw_plot(eval(as.name(width_df$patients[1])), 0, 0.5, width_df$widths[1], 0.5) + 
            draw_plot(eval(as.name(width_df$patients[2])), 0, 0, width_df$widths[2], 0.5) 
    }
    all_plots <<- all_plots
    # The final plot is printed as a png
    png(file = paste(dir_output, file_output, ".png", sep=""),
        width = 400*max,
        height = 800*length(patients),
        res=300)
   print(all_plots)
   dev.off()
}