###################
# VJ Usage Step 2 #
###################

# Code modified from: github.com/pughlab/capTCR_seqAnalysis/blob/main/R/VJCassetteUsage/Step2_TRBclonotypes_txtFile_concatenation.R     
 
# Creates excel file containing the data for VJ usage on each clonotype
# @param data_path: Path towards clonotype data
# @param file_path_output: Path towards desired output folder
# @param file_name: desired output file name

VJUsage_Step2 <- function(data_path, file_path_output, file_name){
    
    file_list_TRB <- as.vector(list.files(path = data_path , recursive=TRUE, 
                            pattern = "^CLONES_TRB(.*)txt$"))
    file_list_TRB <- file_list_TRB[grepl("TLDC", file_list_TRB, fixed=TRUE) == FALSE]

    i <- 1
    for (j in c(1:length(file_list_TRB))) {
            #Extracting the patient_id from the file name:
            patient_id <- unlist(strsplit(unlist(strsplit(file_list_TRB[j], "/"))[2], "TRB"))[2]
            #Extracting the sample cycle from the file name:
            # Shirin's Code
            #cycle <- substr(x = file_list_TRB[j] , 
            #                start = cycle_start , 
            #                stop = cycle_end)
            # Cam's Code
            file_norep <- unlist(strsplit(file_list_TRB[j], "/"))[3]
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
            #Reading in the TRB clones:
            mixcr_TRB <- read.delim(
                    file.path(data_path, 
                              file_list_TRB[j] ),  
                    header = TRUE, 
                    sep = "\t",
                    stringsAsFactors = FALSE,
                    na.strings = c("", "NA"))%>%
            #Selecting the required columns:
                    dplyr::select(cloneCount , cloneFraction , 
                                  aaSeqCDR3 , 
                                  nSeqCDR3 ,
                                  allVHitsWithScore , allJHitsWithScore)%>% 
            #We only work with TRBs with amino-acid length of less than 30 to make sure 
            #we're not including any erroneous clonotypes.
                    filter (str_length(aaSeqCDR3) <= 30) 

            #Since capTCR-seq only captures CDR3 and most of the times leaves germinal V and J sequences,
            #MIXCR has a hard time assigning V and J genes with high confidence. As a result, for some of the 
            #clonotypes it assigns more than 1 V or J gene to a CDR3 and gives them confidence scores.
            #We only pick the first V and J genes that have the highest scores and leave the rest.
            mixcr_TRB$Vcassette  <-  unlist(
                    sapply(
                            strsplit(mixcr_TRB$allVHitsWithScore, split = "*", fixed = TRUE), 
                            function(x) x[[1]][1], simplify=FALSE)
                    )

            mixcr_TRB$Jcassette  <-  unlist(
                    sapply(
                            strsplit(mixcr_TRB$allJHitsWithScore, split = "*", fixed = TRUE), 
                            function(x) x[[1]][1], simplify=FALSE)
                    ) 
            #We have V and J genes separately. We need to combine them and put them in 1 column:
            mixcr_TRB$VJcombo    <- stringr::str_c(mixcr_TRB$Vcassette , "." , mixcr_TRB$Jcassette)

            #Adding the p_is and cycle columns:
            mixcr_TRB$Patient_id <- patient_id
            mixcr_TRB$Cycle      <- cycle1
            mixcr_TRB$Cohort <- cohort

            if(i == 1){
                    comprehensive_TRB <- mixcr_TRB %>% 
                            select(Patient_id , Cycle , Cohort,
                                   cloneCount , cloneFraction , 
                                   aaSeqCDR3 , nSeqCDR3 ,
                                   VJcombo)
                    i <- i + 1   
            }

            else{
                    sub_file_TRB <- mixcr_TRB %>% 
                            select(Patient_id , Cycle , Cohort,
                                   cloneCount , cloneFraction , 
                                   aaSeqCDR3 , nSeqCDR3 ,
                                   VJcombo)
                    comprehensive_TRB <- rbind(comprehensive_TRB, sub_file_TRB)
                    rm(sub_file_TRB , mixcr_TRB)
            }

    }

    ###write down the final file:
    write.csv(comprehensive_TRB ,
              file.path(
                      file_path_output ,
                      file_name ))
}
