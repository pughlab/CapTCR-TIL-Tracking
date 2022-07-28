####################
# VJ Usage Step 1  #
####################

# Code modified from: github.com/pughlab/capTCR_seqAnalysis/blob/main/R/VJCassetteUsage/Step1_aaCDR3perVJ_input_generation.R 

#Required libraries: 
library(tidyverse)

VJUsage_Step1 <- function(file_path, file_path_output){
    # Create file list of TRB Clones MiXCR output while removing the TLDC 1 patient's samples
    file_list <- as.vector(list.files(path = file_path , recursive=TRUE, 
                            pattern = "^CLONES_TRB(.*)txt$"))
    file_list <- file_list[grepl("TLDC", file_list, fixed=TRUE) == FALSE]

    #Code description:
    ###Preparation of Mixcr output .txt files for VJ-cassette usage calculations. These prepared files will be fed into JS_aaCDR3perVJ_v2.5.py. 
    ###JS_aaCDR3perVJ_v2.5.py parses the amino acid CDR3 motifs encoded by each VJ cassette, and 
    ###outputs the number, entropy (H), and normalized entropy (Hnorm) of the aaCDR3 repertoire encoded by each VJ cassette.

    #You go to the 'clones' directory of your mixcr output and run the following commands to open R:
    #Module load R/3.5.0
    #R

    #Set the paths. Don't include '/' at the end of your path:

    #Depending on whether you want to work with TRB or TRA chains, the pattern you pass to this command, will differ. 
    #TRB analysis is always more reliable. Because MIXCR still lacks enough power to annotate TRA chains captured by capTCR-seq. So, TRA VJ-cassette 
    #analysis is not reliable enough to interpret and we'll stick to TRB chains for the rest of this analysis.

    for (fle in file_list) {

            input_file <-  read.delim(file.path (file_path , 
                                                 fle ) , 
                                      quote = "" , 
                                      row.names = NULL, 
                                      stringsAsFactors = FALSE) %>% 
                    #We just need a few columns from each .txt MIXCR output
                    select(nSeqCDR3 , aaSeqCDR3 ,
                           cloneCount , 
                           allVHitsWithScore , allJHitsWithScore)%>% 
                    #We only work with TRBs with amino-acid length of less than 30 to make sure 
                    #we're not including any erroneous clonotype. 
                    filter (str_length(aaSeqCDR3) <= 30) 
            #Since capTCR-seq only captures CDR3 and most of the times leaves germinal V and J sequences,
            #MIXCR has a hard time assigning V and J genes with high confidence. As a result, for some of the 
            #clonotypes it assigns more than 1 V or J gene to a CDR3 and gives them confidence scores.
            #We only pick the first V and J genes that have the highest scores and leave the rest.
            input_file$Vcassette <-  unlist(
                                      sapply(
                                        strsplit(input_file$allVHitsWithScore, split = "*", fixed = TRUE), 
                                        function(x) x[[1]][1], simplify=FALSE)
                                             ) 

            input_file$Jcassette <-  unlist(
                                      sapply(
                                        strsplit(input_file$allJHitsWithScore, split = "*", fixed = TRUE), 
                                        function(x) x[[1]][1], simplify=FALSE)
                                              ) 
            input_file$Counts <- input_file$cloneCount
            input_file$aaCDR3_filtered <- input_file$aaSeqCDR3
            input_file$ntCDR3 <- input_file$nSeqCDR3
            #We have V and J genes separately. We need to combine them:
            input_file$VJcombo <- stringr::str_c(input_file$Vcassette , 
                                                 "." , 
                                                 input_file$Jcassette)

            #We select the files in the same order JS_aaCDR3perVJ_v2.5.py expects:
            input_file <- input_file %>% select(VJcombo , Counts ,
                                                Vcassette , Jcassette,
                                                aaCDR3_filtered , ntCDR3 )


            #Writing out the output to a .tsv file required for JS_aaCDR3perVJ_v2.5.py:
            write_tsv(input_file , 
                      paste( file_path_output , "/" ,
                             #You need to modify the line below based on the name of your files.
                             #For example, if your .txt file name includes DNA in it (as it is my case),
                             #This line splits your file name by "-DNA" and only keeps the first part.
                             #If you don't have a preference for the name of your output file, 
                             #you can replace the line below with:
                             #strsplit(fle , split = ".txt$")[[1]][1]
                             #which only drops .txt from the ending of your input file name
                             #Suggestion: The best scenario is to only keep
                             #patient_id , cycle and sample type 
                             strsplit(fle , split = "-DNA", fixed = TRUE)[[1]][1] ,
                             "_only.productive.tsv" , 
                             sep = ""))

    }

    ###Concatenating the outputs for analysis:

    ###go to the directore where you've stored '_aaCDR3perVJ.tsv' files that you've generated via JS_aaCDR3perVJ_v2.5.py.
    file_path <- file_path_output
    file_list <- list.files(path = file_path , 
                            pattern = "aaCDR3perVJ.tsv")

    #In the following code, we want to store patient_id and cycle in 2 different columns.
    #If you have clean file names, where patient_id and cycle start at specific coordination,
    #you can easily pass the start and end coordinations to the following code, as I've done.
    #but if your file names are not neat, you gotta find a way to extract patient_id
    #and cycle and store them in 2 different columns:
    ##the 4 variables below are integers. The integers are the coordination of where 
    ##p_id or cycle start and end in the file name.
    ##for example, in 'CLONES_TRBLIB-14-0018-T3-M-DNA_only.productive.tsv' 
    ##p_id_start is 11, p_id_end is 21, cycle_start is 23 and cycle_end is 24.

    for (i in c (1 : length(file_list))) {
            p_id_start  <- as.numeric(unlist(gregexpr("TLML", file_list[i])))
            p_id_end <- p_id_start + 7
            #cycle_end <- nchar(file_list_TRB[j])
            #Extracting the patient_id from the file name:
            patient_id <- substr(x = file_list[i] , 
                                 start = p_id_start , 
                                 stop = p_id_end)
            #Extracting the sample cycle from the file name:
            # Shirin's Code
            #cycle <- substr(x = file_list_TRB[j] , 
            #                start = cycle_start , 
            #                stop = cycle_end)
            # Cam's Code
            file_norep <- file_list[i]
            if(grepl("TLML_", file_norep, fixed=TRUE)==TRUE){
                cycle1 <- unlist(strsplit(file_norep, "_"))[5]
                cohort <- unlist(strsplit(tail(unlist(strsplit(file_norep, "TCR")), n=1), "_"))[2]
            }
            else if((grepl("TLML-", file_norep, fixed=TRUE)==TRUE) & (grepl("17_02", file_norep, fixed=TRUE)==FALSE)){
                cycle1 <- unlist(strsplit(file_norep, "_"))[3]
                cohort <- unlist(strsplit(file_norep, "_"))[5]
            }
            else if((grepl("TLML-", file_norep, fixed=TRUE)==TRUE) & (grepl("17_02", file_norep, fixed=TRUE)==TRUE)){
                cycle1 <- unlist(strsplit(unlist(strsplit(file_norep, "-"))[4], "_"))[1]
                cohort <- unlist(strsplit(tail(unlist(strsplit(file_norep, "TCR")), n=1), "_"))[2]
            }
            if ( i == 1) {
                    comprehensive <- read_tsv(file_list[i])
                    comprehensive$Patient_id <- patient_id
                    comprehensive$Cycle <- cycle1
                    comprehensive$Cohort <- cohort

            }
            else {
                    transient <- read_tsv(file_list[i])    
                    transient$Patient_id <- patient_id
                    transient$Cycle <- cycle1
                    transient$Cohort <- cohort
                    comprehensive <- rbind(comprehensive , transient)
                    remove(transient)
            }

    }                                           


    ###write down the final file:
    ###file name should end with .tsv or .csv depending on the format you like to save it:
    file_name <- "VJUsage.csv"

    write.csv(comprehensive ,
              file.path(
                      file_path ,
                      file_name )
             )
}