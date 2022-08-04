################################
# Unproductive Clones Analysis #
################################

# Analyzes the unproductive data consumption of the patient's clones
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param datapath: Path to desired files
# @param filename: Desired name of the output unproductive clones csv file

clonstatfx <- function(chain, datapath, filename, dir_output){
    
flelst <- list.files(datapath,
                     recursive = TRUE,
                     pattern = paste("CLONES", chain, sep = "_"))

# Creates a list of files and presents them to the user
flelst <- flelst[grepl(chain, flelst)]
    message("list of files:")
    print(flelst)
  
  # Creates the dataframe to store the unproductive clones data
  clones <- cbind.data.frame(flelst,
                    NA, NA, NA, NA)
  colnames(clones) <- c("filename", "total",
                        "out_of_frame", "stopcodon",
                       "productive")

  # For each file in the clones dataframe, the data is loaded and the amount of unproductive clones are recorded and data is stored in the clones dataframe
  for(i in 1:nrow(clones)){
    f <- clones$filename[i]
    mixcrfle <- read.table(paste0(datapath, f),
                           header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE,
                           na.strings = c("", "NA"))

    clones[i,colnames(clones) == "total"] <- length(mixcrfle$aaSeqCDR3)
    clones[i,colnames(clones) == "out_of_frame"] <- length(mixcrfle$aaSeqCDR3[grepl("_", mixcrfle$aaSeqCDR3)])
    clones[i,colnames(clones) == "stopcodon"] <- length(mixcrfle$aaSeqCDR3[grepl("[*]", mixcrfle$aaSeqCDR3) &
                                                                             !grepl("_", mixcrfle$aaSeqCDR3)])
    clones[i,colnames(clones) == "productive"] <-  length(mixcrfle$aaSeqCDR3[!grepl("[*]", mixcrfle$aaSeqCDR3) &
                                                                             !grepl("_", mixcrfle$aaSeqCDR3)])

  }
  
  # A CSV file is written to store the clones and exported
  write.csv(clones,
         file = paste(dir_output, "clonstats_", chain,
                       filename, ".csv", sep = ""),
         row.names = F)
}