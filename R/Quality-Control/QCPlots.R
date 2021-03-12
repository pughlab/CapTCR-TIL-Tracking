#########################
# Quality Control Plots #
#########################

# Uses tabular log dataset from MiXCR to create and bind multiple quality control plots together showing different data points
# @param alignstatsfile: Name of the tabular file containing the align stats 
# @param assemblestatsfile: Name of the tabular file containing the assemble stats file
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param plotname: Name of the png to be printed with the aligned quality control plots
# @param inpath: The path to the align and assemble stats files
# @param plotpath: The path the png is to be printed into

mixcrQC.fx <- function(alignstatsfile, assemblestatsfile, sampcohort, 
                       plotname, inpath, plotpath){
    align_plots1 <- function (...) {
    pl <- list(...)
    stopifnot(do.call(all, lapply(pl, inherits, "gg")))
    gl <- lapply(pl, ggplotGrob)
    bind2 <- function(x, y) gtable:::rbind_gtable(x, y, "first")
    combined <- Reduce(bind2, gl[-1], gl[[1]])
    wl <- lapply(gl, "[[", "widths")
    combined$widths <- do.call(grid::unit.pmax, wl)
    grid::grid.newpage()
    grid::grid.draw(combined)
  }
  
    ## Plot Theme ##

    # Creates the theme 'mytheme' which features no panel or coloured background. #

    mytheme <- theme(axis.title.y = element_text(size = 28),
                       axis.title.x = element_blank(),
                       axis.line = element_line(color = "black"),
                       axis.text = element_text(size = 28),
                       axis.text.x = element_text(angle = 75, hjust = 1, size = 20)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "transparent",colour = NA),  
              legend.key = element_rect(fill = "white", colour = "white"),
              plot.margin = unit(c(0,0,0,0),"cm"))    


     ## Reading Log Files from MiXCR ##

    # Reads the outputed align and assemble stats from the MiXCR software. #

      alignstats <- read.csv(file = paste(inpath,alignstatsfile, sep =""),
                             sep = ",", header = T)    
      assemblestats <- read.csv(file = paste(inpath, assemblestatsfile, sep =""),
                                sep = ",", header = T)  
      
    # Creates a list of samples for the desired sample cohort
    myfilenames <- alignstats$SampleID
    samplelist <- myfilenames[grepl(sampcohort, myfilenames)]

    ## Sample List Matching. #

    # Matches the sample names to the previously generated SampleID column in the align and assemb
    # files. NOTE: The SampleID column was added manually in order to ensure proper x axis ordering
    # and size of the labels. 

      for(i in samplelist){
        alignstats$samplename[grepl(i, alignstats$SampleID)] <- i
      }    
      alignstats <- alignstats[!is.na(alignstats$samplename),]

      for(j in samplelist){
        assemblestats$samplename[grepl(j, assemblestats$SampleID)] <- j
      }
      assemblestats <- assemblestats[!is.na(assemblestats$samplename),]

    # Ensures that the samplename is in the same order as the samplelist. 
      alignstats$samplename <- factor(alignstats$samplename, levels = samplelist)
      assemblestats$samplename <- factor(assemblestats$samplename, levels = samplelist)

    ## Plot Creation from Quality Control Data ##
    
    # 'myplot_totalseq_align_seq' creates the plot that uses both the 'Total.sequencing.reads' and the 'Successfully.aligned.reads'
    # from the align stats file.              
    
    myplot_totalseq_align_seq <- ggplot(aes(x = samplename), 
                                data = alignstats) +
        geom_point(aes(y = Total.sequencing.reads),size = 7) +
        geom_point(aes(y = Successfully.aligned.reads),size = 7,color="#FF0000") +
        mytheme + theme(axis.text.x = element_blank()) +
        scale_y_continuous("Total Reads",sec.axis = dup_axis(name="Total Aligned Reads\n"))
    myplot_totalseq_align_seq <- myplot_totalseq_align_seq + theme(axis.title.y.right = element_text(color="#FF0000"))
    
    # 'myplot_percentalignedreads' creates the plot which uses the data 
    
    myplot_percentalignedreads <- ggplot(data=alignstats, aes(x=samplename, y=Percent.Aligned.Reads)) + 
        geom_bar(stat="identity") + mytheme + theme(axis.text.x=element_blank()) + ylim(0,100)
    
    

    #'myplot_averagereads' creates the plot that uses the data 'Average.number.of.reads.per.clonotype' 
    # from the assemble stats file.                                                                    

      myplot_averagereads <- ggplot(aes(x = samplename, y = Average.number.of.reads.per.clonotype), 
                                    data = assemblestats) + 
        geom_point(size = 7) + 
        mytheme +
        theme(axis.text.x = element_blank())

    #'myplot_readsused' creates the plot that uses the data 'Reads.used.in.clonotypes.before.clustering 
    # ..percent.of.total' from the assemble stats file.                                                                                                                  #

      myplot_readsused <- ggplot(aes(x = samplename, y = Reads.used.in.clonotypes.before.clustering..percent.of.total),       
                                 data = assemblestats) + 
        geom_point(size = 7) + 
        mytheme + 
        theme(axis.text.x = element_blank())

    #'myplot_totalclon' creates the plot that uses the data 'Final.clonotype.count' from the assemble 
    # stats file                                                                                      

      myplot_totalclon <- ggplot(aes(x = samplename, y = Final.clonotype.count),        
                                 data = assemblestats) + 
        geom_point(size = 7) + 
        mytheme

    # Exports all plots as pdfs combined using the 'align_plots1' function created earlier 
      png(file = paste(plotpath, plotname, sep = ""),
          width = 13000, 
          height = 8000,
          res = 300)
      align_plots1(myplot_totalseq_align_seq,
                   myplot_percentalignedreads + ylab("Percentage \nAligned Reads"),
                   myplot_readsused + ylab("Reads used \nin clonotypes"),
                   myplot_averagereads + ylab("Average \nreads/clonotype"),
                   myplot_totalclon + ylab("Total clonotypes"))    
      dev.off()  
    }                       