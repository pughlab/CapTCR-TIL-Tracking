# Loading scripts
GitHub_path <- paste0(getwd(), "/")
dir_main <- paste(GitHub_path, "data/", sep="")

# Generating list of relevant R scripts from GitHub repo and then loading in functions and libraries
R_scripts <- list.files(path=paste(GitHub_path, "code/", sep=""), pattern="\\.R$", recursive=TRUE)
R_scripts <- R_scripts[!grepl("Quality-Control|VJUsage_Step1|VJUsage_Step2|Reproduce", R_scripts)]
for(script in R_scripts){
  source(paste(GitHub_path, "code/", script, sep=""))
}

# Loading MiXCR, sample keys, and CDR3 colors data
Initialize(dir_main, c("baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05"))
# Reproduces all data from the supplementary data 
ReproduceSupplementaryData <- function(GitHub_path){
  
  # Generating Table S6 of the supplementary material
  message("Generating Table S6")
  TableS6_df <- data.frame(matrix(NA, nrow=9, ncol=5))
  colnames(TableS6_df) <- c("Patient", "DNA baseline", "DNA first post-infusion", "cfDNA baseline", "cfDNA first post-infusion")
  TableS6_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 9,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS6_df$Patient)){
    DNAsamporder <- eval(as.name(paste0(TableS6_df$Patient[row], "DNA_samporder")))
    DNAbase <- persistence(TableS6_df$Patient[row],"DNA",DNAsamporder[1])
    DNApost <- persistence(TableS6_df$Patient[row],"DNA",DNAsamporder[3])
    if(TableS6_df$Patient[row] %in% c("TLML_4_", "TLML_20_", "TLML_16_", "TLML_22_", "TLML_29_", "TLML_26_")) {
      cfDNAsamporder <- eval(as.name(paste0(TableS6_df$Patient[row], "cfDNA_samporder")))
      cfDNAbase <- persistence(TableS6_df$Patient[row],"cfDNA",cfDNAsamporder[1])
      cfDNApost <- persistence(TableS6_df$Patient[row],"cfDNA",cfDNAsamporder[2])
    }
    else {
      cfDNAbase <- NA
      cfDNApost <- NA
    }
    TableS6_df[row,] <- c(TableS6_df$Patient[row], DNAbase, DNApost, cfDNAbase, cfDNApost)
    i <- i + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Generating Table S7, S8, S9 of the supplementary material
  message("Generating Table S7, S8, S9")
  TableS7_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  TableS8_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  TableS9_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS7_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS7_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS8_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS8_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS9_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS9_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 17,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS7_df$Patient)) {
    samporder <- eval(as.name(paste0(TableS7_df$Patient[row], "DNA_samporder")))
    for(column in 2:ncol(TableS7_df)) {
      if((column-1) <= length(samporder)) {
        count <- threshold_stats(TableS7_df$Patient[row], "DNA", samporder[column-1], "numberabove")
        percentage <- threshold_stats(TableS7_df$Patient[row], "DNA", samporder[column-1], "percentage")
        abundance <- threshold_stats(TableS7_df$Patient[row], "DNA", samporder[column-1], "abundance")
        TableS7_df[row, column] <- percentage
        TableS8_df[row, column] <- count
        TableS9_df[row, column] <- abundance
        i <- i + 1
        setTxtProgressBar(pb, i)
      }
    }
    if(TableS7_df$Patient[row] == "TLML_1_"){
      TableS7_df[row,] <- c("TLML_1_", TableS7_df[row, 2], TableS7_df[row, 3], TableS7_df[row,4], TableS7_df[row,4], TableS7_df[row, 5], TableS7_df[row, 6], NA, NA)
      TableS8_df[row,] <- c("TLML_1_", TableS8_df[row, 2], TableS8_df[row, 3], TableS8_df[row,4], TableS8_df[row,4], TableS8_df[row, 5], TableS8_df[row, 6], NA, NA)
      TableS9_df[row,] <- c("TLML_1_", TableS9_df[row, 2], TableS9_df[row, 3], TableS9_df[row,4], TableS9_df[row,4], TableS9_df[row, 5], TableS9_df[row, 6], NA, NA)
    }
  }
  close(pb)
  
  # Generating Table S10 of the supplementary material
  message("Creating TableS10")
  TableS10_df <- data.frame(matrix(NA, nrow=9, ncol=3))
  colnames(TableS10_df) <- c("Patient", "In baseline", "In infusion")
  TableS10_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 9,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS10_df$Patient)){
    samporder <- eval(as.name(paste0(TableS10_df$Patient[row], "DNA_samporder")))
    baseline <- expanded_reference(TableS10_df$Patient[row], "DNA", samporder[3], 
                                   TableS10_df$Patient[row], "DNA", samporder[1], "total")
    infusion <- expanded_reference(TableS10_df$Patient[row], "DNA", samporder[3], 
                                   TableS10_df$Patient[row], "DNA", samporder[2], "total")
    TableS10_df[row,] <- c(TableS10_df$Patient[row], baseline, infusion)
    i <- i + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Generating Table S11 of the supplementary material
  message("Creating TableS11")
  TableS11_df <- data.frame(matrix(NA, nrow=9, ncol=3))
  colnames(TableS11_df) <- c("Patient", "Baseline", "Infusion")
  TableS11_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 9,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS11_df$Patient)){
    samporder <- eval(as.name(paste0(TableS11_df$Patient[row], "DNA_samporder")))
    baseline <- relativerisk(TableS11_df$Patient[row], "DNA", samporder[1])$RelativeRisk
    infusion <- relativerisk(TableS11_df$Patient[row], "DNA", samporder[2])$RelativeRisk
    TableS11_df[row,] <- c(TableS11_df$Patient[row], baseline, infusion)
    i <- i + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Generating Table S12, S13, S14, S15, S16, S17
  message("Creating TableS12, S13, S14, S15, S16, S17")
  TableS12_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  TableS13_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  TableS14_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  TableS15_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  TableS16_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  TableS17_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS12_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS12_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS13_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS13_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS14_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS14_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS15_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS15_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS16_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS16_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS17_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS17_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 17,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS12_df$Patient)) {
    samporder <- eval(as.name(paste0(TableS12_df$Patient[row], "DNA_samporder")))
    for(column in 2:ncol(TableS12_df)) {
      if((column - 1) <= length(samporder)) {
        DNAsimpson <- 1/SimpsonIndex(TableS12_df$Patient[row], "DNA", samporder[column-1])
        DNAshannon <- ShannonIndex(TableS12_df$Patient[row], "DNA", samporder[column-1])
        DNArichness <- Richness(TableS12_df$Patient[row], "DNA", samporder[column-1])
        DNAclonality <- 1 - Eveness(TableS12_df$Patient[row], "DNA", samporder[column-1])
        DNAexpanded1 <- HighExpandedClone(TableS12_df$Patient[row], "DNA", samporder[column-1], 0.01)
        DNAexpanded5 <- HighExpandedClone(TableS12_df$Patient[row], "DNA", samporder[column-1], 0.05)
        TableS12_df[row, column] <- DNAsimpson
        TableS13_df[row, column] <- DNAshannon
        TableS14_df[row, column] <- DNArichness
        TableS15_df[row, column] <- DNAclonality
        TableS16_df[row, column] <- DNAexpanded1
        TableS17_df[row, column] <- DNAexpanded5
        i <- i + 1
        setTxtProgressBar(pb, i)
      }
    }
    if(TableS12_df$Patient[row] == "TLML_1_"){
      TableS12_df[row,] <- c("TLML_1_", TableS12_df[row, 2], TableS12_df[row, 3], TableS12_df[row,4], TableS12_df[row,4], TableS12_df[row, 5], TableS12_df[row, 6], NA, NA)
      TableS13_df[row,] <- c("TLML_1_", TableS13_df[row, 2], TableS13_df[row, 3], TableS13_df[row,4], TableS13_df[row,4], TableS13_df[row, 5], TableS13_df[row, 6], NA, NA)
      TableS14_df[row,] <- c("TLML_1_", TableS14_df[row, 2], TableS14_df[row, 3], TableS14_df[row,4], TableS14_df[row,4], TableS14_df[row, 5], TableS14_df[row, 6], NA, NA)
      TableS15_df[row,] <- c("TLML_1_", TableS15_df[row, 2], TableS15_df[row, 3], TableS15_df[row,4], TableS15_df[row,4], TableS15_df[row, 5], TableS15_df[row, 6], NA, NA)
      TableS16_df[row,] <- c("TLML_1_", TableS16_df[row, 2], TableS16_df[row, 3], TableS16_df[row,4], TableS16_df[row,4], TableS16_df[row, 5], TableS16_df[row, 6], NA, NA)
      TableS17_df[row,] <- c("TLML_1_", TableS17_df[row, 2], TableS17_df[row, 3], TableS17_df[row,4], TableS17_df[row,4], TableS17_df[row, 5], TableS17_df[row, 6], NA, NA)
    }
  }
  close(pb)
  
  # Generating Table S18, S19, S20, S21, S22
  message("Creating TableS18, S19, S20, S21, S22")
  TableS18_df <- data.frame(matrix(NA, nrow=6, ncol=3))
  TableS19_df <- data.frame(matrix(NA, nrow=6, ncol=3))
  TableS20_df <- data.frame(matrix(NA, nrow=6, ncol=3))
  TableS21_df <- data.frame(matrix(NA, nrow=6, ncol=3))
  TableS22_df <- data.frame(matrix(NA, nrow=6, ncol=3))
  colnames(TableS18_df) <- c("Patient", "baseline", "4W")
  TableS18_df$Patient <- c("TLML_4_", "TLML_16_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS19_df) <- c("Patient", "baseline", "4W")
  TableS19_df$Patient <- c("TLML_4_", "TLML_16_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS20_df) <- c("Patient", "baseline", "4W")
  TableS20_df$Patient <- c("TLML_4_", "TLML_16_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS21_df) <- c("Patient", "baseline", "4W")
  TableS21_df$Patient <- c("TLML_4_", "TLML_16_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  colnames(TableS22_df) <- c("Patient", "baseline", "4W")
  TableS22_df$Patient <- c("TLML_4_", "TLML_16_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 12,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS18_df$Patient)) {
    samporder <- eval(as.name(paste0(TableS18_df$Patient[row], "cfDNA_samporder")))
    for(column in 2:ncol(TableS18_df)) {
      cfDNAsimpson <- 1/SimpsonIndex(TableS18_df$Patient[row], "cfDNA", samporder[column-1])
      cfDNAshannon <- ShannonIndex(TableS18_df$Patient[row], "cfDNA", samporder[column-1])
      cfDNArichness <- Richness(TableS18_df$Patient[row], "cfDNA", samporder[column-1])
      cfDNAclonality <- 1 - Eveness(TableS18_df$Patient[row], "cfDNA", samporder[column-1])
      cfDNAexpanded5 <- HighExpandedClone(TableS18_df$Patient[row], "cfDNA", samporder[column-1], 0.05)
      TableS18_df[row, column] <- cfDNAsimpson
      TableS19_df[row, column] <- cfDNAshannon
      TableS20_df[row, column] <- cfDNArichness
      TableS21_df[row, column] <- cfDNAclonality
      TableS22_df[row, column] <- cfDNAexpanded5
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    if(TableS18_df$Patient[row] == "TLML_1_"){
      TableS18_df[row] <- c("TLML_1_", TableS18_df[row, 2], TableS18_df[row, 3], TableS18_df[row,4], TableS18_df[row,4], TableS18_df[row, 5], TableS18_df[row, 5], NA, NA)
      TableS19_df[row] <- c("TLML_1_", TableS19_df[row, 2], TableS19_df[row, 3], TableS19_df[row,4], TableS19_df[row,4], TableS19_df[row, 5], TableS19_df[row, 5], NA, NA)
      TableS20_df[row] <- c("TLML_1_", TableS20_df[row, 2], TableS20_df[row, 3], TableS20_df[row,4], TableS20_df[row,4], TableS20_df[row, 5], TableS20_df[row, 5], NA, NA)
      TableS21_df[row] <- c("TLML_1_", TableS21_df[row, 2], TableS21_df[row, 3], TableS21_df[row,4], TableS21_df[row,4], TableS21_df[row, 5], TableS21_df[row, 5], NA, NA)
      TableS22_df[row] <- c("TLML_1_", TableS22_df[row, 2], TableS22_df[row, 3], TableS22_df[row,4], TableS22_df[row,4], TableS22_df[row, 5], TableS22_df[row, 5], NA, NA)
    }
  }
  close(pb)
  
  # Exporting supplementary data excel file
  message("Exporting supplementary excel file")
  pb <- txtProgressBar(min = 0,max = 14,style = 3,   width = 50,char = "=")
  library(openxlsx)
  SupplementaryData <- createWorkbook()
  addWorksheet(SupplementaryData, "Table S6")
  addWorksheet(SupplementaryData, "Table S7")
  addWorksheet(SupplementaryData, "Table S8")
  addWorksheet(SupplementaryData, "Table S9")
  addWorksheet(SupplementaryData, "Table S10")
  addWorksheet(SupplementaryData, "Table S11")
  addWorksheet(SupplementaryData, "Table S12")
  addWorksheet(SupplementaryData, "Table S13")
  addWorksheet(SupplementaryData, "Table S14")
  addWorksheet(SupplementaryData, "Table S15")
  addWorksheet(SupplementaryData, "Table S16")
  addWorksheet(SupplementaryData, "Table S17")
  addWorksheet(SupplementaryData, "Table S18")
  addWorksheet(SupplementaryData, "Table S19")
  addWorksheet(SupplementaryData, "Table S20")
  addWorksheet(SupplementaryData, "Table S21")
  addWorksheet(SupplementaryData, "Table S22")
  
  writeData(SupplementaryData, sheet="Table S6", x=TableS6_df)
  setTxtProgressBar(pb, 4)
  writeData(SupplementaryData, sheet="Table S7", x=TableS7_df)
  setTxtProgressBar(pb, 5)
  writeData(SupplementaryData, sheet="Table S8", x=TableS8_df)
  setTxtProgressBar(pb, 6)
  writeData(SupplementaryData, sheet="Table S9", x=TableS9_df)
  setTxtProgressBar(pb, 7)
  writeData(SupplementaryData, sheet="Table S10", x=TableS10_df)
  setTxtProgressBar(pb, 8)
  writeData(SupplementaryData, sheet="Table S11", x=TableS11_df)
  setTxtProgressBar(pb, 9)
  writeData(SupplementaryData, sheet="Table S12", x=TableS12_df)
  setTxtProgressBar(pb, 10)
  writeData(SupplementaryData, sheet="Table S13", x=TableS13_df)
  setTxtProgressBar(pb, 11)
  writeData(SupplementaryData, sheet="Table S14", x=TableS14_df)
  setTxtProgressBar(pb, 12)
  writeData(SupplementaryData, sheet="Table S15", x=TableS15_df)
  setTxtProgressBar(pb, 13)
  writeData(SupplementaryData, sheet="Table S16", x=TableS16_df)
  setTxtProgressBar(pb, 14)
  writeData(SupplementaryData, sheet="Table S17", x=TableS17_df)
  setTxtProgressBar(pb, 14)
  writeData(SupplementaryData, sheet="Table S18", x=TableS18_df)
  setTxtProgressBar(pb, 14)
  writeData(SupplementaryData, sheet="Table S19", x=TableS19_df)
  setTxtProgressBar(pb, 14)
  writeData(SupplementaryData, sheet="Table S20", x=TableS20_df)
  setTxtProgressBar(pb, 14)
  writeData(SupplementaryData, sheet="Table S21", x=TableS21_df)
  setTxtProgressBar(pb, 14)
  writeData(SupplementaryData, sheet="Table S22", x=TableS22_df)
  setTxtProgressBar(pb, 14)
  saveWorkbook(SupplementaryData, paste(GitHub_path, "results/SupplementaryData_reproduced.xlsx", sep=""))
  close(pb)
}
ReproduceSupplementaryData(GitHub_path)

ReproduceClonetracking <- function(GitHub_path) { 
  
  message("Creating clone tracking plot")
  pb <- txtProgressBar(min = 0, max=10, style=3, width=50, char="=")
  # Fetching clone tracking plots
  DNApatients <- c("TLML_22_", "TLML_26_", "TLML_29_", "TLML_1_", "TLML_16_",
                   "TLML_7_", "TLML_20_", "TLML_4_", "TLML_18_")
  cfDNApatients <- c("TLML_22_", "TLML_26_", "TLML_29_", "TLML_16_", "TLML_20_",
                     "TLML_4_")
  setTxtProgressBar(pb, 1)
  alignment_fig(DNApatients, "DNA", "TRB", "clonetrack", "Baseline", 
                paste0(GitHub_path, "data/",sep=""), "clonetrack_DNA")
  setTxtProgressBar(pb, 2)
  for(i in 1:length(cfDNApatients)){
    alignment_fig(cfDNApatients[i], "cfDNA", "TRB", "clonetrack_cfDNA", 
                  "Baseline", paste0(GitHub_path, "data/", sep=""), 
                  paste0("clonetrack_", cfDNApatients[i]))
    setTxtProgressBar(pb, i + 2)
  }
  
  # Reading exported clone tracking plots and exporting into figure format
  DNAclonetrack <- image_read(paste0(GitHub_path, "data/clonetrack_DNA.png"))
  TLML22_clonetrack <- image_read(paste0(GitHub_path, "data/clonetrack_TLML_22_.png"))
  TLML26_clonetrack <- image_read(paste0(GitHub_path, "data/clonetrack_TLML_26_.png"))
  TLML29_clonetrack <- image_read(paste0(GitHub_path, "data/clonetrack_TLML_29_.png"))
  TLML16_clonetrack <-image_read(paste0(GitHub_path, "data/clonetrack_TLML_16_.png"))
  TLML20_clonetrack <-image_read(paste0(GitHub_path, "data/clonetrack_TLML_20_.png"))
  TLML4_clonetrack <- image_read(paste0(GitHub_path, "data/clonetrack_TLML_4_.png"))
  FigOverlay <- image_read(paste(GitHub_path, "data/Overlay/clonetracking_overlay.png", sep=""))
  
  setTxtProgressBar(pb, 9)
  
  png(paste(GitHub_path, "results/clonetracking_reproduced.png", sep=""), width = 6700, height = 7600, units = "px")
  plot.new()
  plot.window(xlim=c(0, 6700), ylim=c(0, 7600))
  rasterImage(DNAclonetrack, xleft=1464, ybottom=851, xright=4270  , ytop=7161)
  rasterImage(TLML22_clonetrack, xleft=4936 , xright=5989 , ybottom=6459 , ytop=7161)
  rasterImage(TLML26_clonetrack, xleft=4936 , xright=5989 , ybottom=5757 , ytop=6459)
  rasterImage(TLML29_clonetrack, xleft=4936 , xright=5989 , ybottom=5065 , ytop=5767)
  rasterImage(TLML16_clonetrack, xleft=4936 , xright=5989 , ybottom=3663 , ytop=4365)
  rasterImage(TLML20_clonetrack, xleft=4936 , xright=5989 , ybottom=2248 , ytop=2950)
  rasterImage(TLML4_clonetrack, xleft=4936 , xright=5989 , ybottom=1553 , ytop=2255)
  rasterImage(FigOverlay, xleft=0 , ybottom=0 , xright=6700 , ytop=7600)
  dev.off()
  
  unlink(paste(GitHub_path, "data/clonetrack_DNA.png", sep=""))
  unlink(paste(GitHub_path, "data/clonetrack_TLML_22_.png", sep=""))
  unlink(paste(GitHub_path, "data/clonetrack_TLML_26_.png", sep=""))
  unlink(paste(GitHub_path, "data/clonetrack_TLML_29_.png", sep=""))
  unlink(paste(GitHub_path, "data/clonetrack_TLML_16_.png", sep=""))
  unlink(paste(GitHub_path, "data/clonetrack_TLML_20_.png", sep=""))
  unlink(paste(GitHub_path, "data/clonetrack_TLML_4_.png", sep=""))
  
  setTxtProgressBar(pb, 10)
  close(pb)
  
}
ReproduceDiversity <- function(GitHub_path) {
  
  message("Creating diversity plot")
  pb <- txtProgressBar(min = 0, max=10, style=3, width=50, char="=")
  # Fetching clone tracking plots
  DNApatients <- c("TLML_22_", "TLML_26_", "TLML_29_", "TLML_1_", "TLML_16_",
                   "TLML_7_", "TLML_20_", "TLML_4_", "TLML_18_")
  cfDNApatients <- c("TLML_22_", "TLML_26_", "TLML_29_", "TLML_16_", "TLML_20_",
                     "TLML_4_")
  setTxtProgressBar(pb, 1)
  alignment_fig(DNApatients, "DNA", "TRB", "rel_div", "Baseline", 
                paste0(GitHub_path, "data/",sep=""), "diversity_DNA")
  setTxtProgressBar(pb, 2)
  for(i in 1:length(cfDNApatients)){
    alignment_fig(cfDNApatients[i], "cfDNA", "TRB", "rel_div", 
                  "Baseline", paste0(GitHub_path, "data/", sep=""), 
                  paste0("diversity_", cfDNApatients[i]))
    setTxtProgressBar(pb, i + 2)
  }
  
  # Reading exported clone tracking plots and exporting into figure format
  DNAclonetrack <- image_read(paste0(GitHub_path, "data/diversity_DNA.png"))
  TLML22_clonetrack <- image_read(paste0(GitHub_path, "data/diversity_TLML_22_.png"))
  TLML26_clonetrack <- image_read(paste0(GitHub_path, "data/diversity_TLML_26_.png"))
  TLML29_clonetrack <- image_read(paste0(GitHub_path, "data/diversity_TLML_29_.png"))
  TLML16_clonetrack <-image_read(paste0(GitHub_path, "data/diversity_TLML_16_.png"))
  TLML20_clonetrack <-image_read(paste0(GitHub_path, "data/diversity_TLML_20_.png"))
  TLML4_clonetrack <- image_read(paste0(GitHub_path, "data/diversity_TLML_4_.png"))
  FigOverlay <- image_read(paste(GitHub_path, "data/Overlay/diversity_overlay.png", sep=""))
  
  setTxtProgressBar(pb, 9)
  
  png(paste(GitHub_path, "results/diversity_reproduced.png", sep=""), width = 6500, height = 7500, units = "px")
  plot.new()
  plot.window(xlim=c(0, 6500), ylim=c(0, 7500))
  rasterImage(DNAclonetrack, xleft=1256, ybottom=736, xright=4051  , ytop=7026)
  rasterImage(TLML22_clonetrack, xleft=4888 , xright=5764 , ybottom=6325 , ytop=7025)
  rasterImage(TLML26_clonetrack, xleft=4888 , xright=5764 , ybottom=5623 , ytop=6323)
  rasterImage(TLML29_clonetrack, xleft=4888 , xright=5764 , ybottom=4917 , ytop=5617)
  rasterImage(TLML16_clonetrack, xleft=4888 , xright=5764 , ybottom=3531 , ytop=4231)
  rasterImage(TLML20_clonetrack, xleft=4888 , xright=5764 , ybottom=2119 , ytop=2819)
  rasterImage(TLML4_clonetrack, xleft=4888 , xright=5764 , ybottom=1419 , ytop=2119)
  rasterImage(FigOverlay, xleft=0 , ybottom=0 , xright=6500 , ytop=7500)
  dev.off()
  
  unlink(paste(GitHub_path, "data/diversity_DNA.png", sep=""))
  unlink(paste(GitHub_path, "data/diversity_TLML_22_.png", sep=""))
  unlink(paste(GitHub_path, "data/diversity_TLML_26_.png", sep=""))
  unlink(paste(GitHub_path, "data/diversity_TLML_29_.png", sep=""))
  unlink(paste(GitHub_path, "data/diversity_TLML_16_.png", sep=""))
  unlink(paste(GitHub_path, "data/diversity_TLML_20_.png", sep=""))
  unlink(paste(GitHub_path, "data/diversity_TLML_4_.png", sep=""))
  
  setTxtProgressBar(pb, 10)
  close(pb)
}
ReproduceOverlap <- function(outputmatrix, GitHub_path) {
  
  message("Creating overlap plot")
  pb <- txtProgressBar(min = 0, max=6, style=3, width=50, char="=")
  
  # Fetching clone tracking plots
  cfDNApatients <- c("TLML_22_", "TLML_26_", "TLML_29_", "TLML_16_", "TLML_20_",
                     "TLML_4_")
  setTxtProgressBar(pb, 1)
  overlap(outputmatrix, GitHub_path)
  setTxtProgressBar(pb, 2)
  alignment_fig(cfDNApatients, "cfDNA", "TRB", "cfDNAcorrel", "Baseline", 
                paste0(GitHub_path, "data/",sep=""), "cfDNAcorrel_TIL", clones = "TIL")
  setTxtProgressBar(pb, 3)
  alignment_fig(cfDNApatients, "cfDNA", "TRB", "cfDNAcorrel", "Baseline",
                paste0(GitHub_path, "data/", sep=""), "cfDNAcorrel_background", clones = "Background")
  setTxtProgressBar(pb, 4)
  richness_boxplot(GitHub_path)
  
  heatmap <- image_read(paste0(GitHub_path, "data/heatmap.png"))
  tils <- image_read(paste0(GitHub_path, "data/cfDNAcorrel_TIL.png"))
  background <- image_read(paste0(GitHub_path, "data/cfDNAcorrel_background.png"))
  richness <- image_read(paste0(GitHub_path, "data/boxplot.png"))
  FigOverlay <- image_read(paste0(GitHub_path, "data/Overlay/overlap_overlay.png"))
  
  setTxtProgressBar(pb, 5)
  
  png(paste(GitHub_path, "results/overlap_reproduced.png", sep=""), width = 12500, height = 4500, units = "px")
  plot.new()
  plot.window(xlim=c(0, 12500), ylim=c(0, 4500))
  rasterImage(heatmap, xleft=50, ybottom=-247, xright=5591  , ytop=4253)
  rasterImage(background, xleft=5991 , xright=7726 , ybottom=468 , ytop=3940)
  rasterImage(tils, xleft=8066 , xright=9801 , ybottom=468 , ytop=3938)
  rasterImage(richness, xleft=10387 , ybottom=424 , xright=12078 , ytop=3882)
  rasterImage(FigOverlay, xleft=0 , ybottom=0 , xright=12500 , ytop=4500)
  dev.off()
  
  unlink(paste(GitHub_path, "data/heatmap.png", sep=""))
  unlink(paste(GitHub_path, "data/boxplot.png", sep=""))
  unlink(paste(GitHub_path, "data/cfDNAcorrel_TIL.png", sep=""))
  unlink(paste(GitHub_path, "data/cfDNAcorrel_background.png", sep=""))
  
  setTxtProgressBar(pb, 6)
  close(pb)
}
ReproduceOverlap(overlapmatrix, GitHub_path)
ReproduceDiversity(GitHub_path)
ReproduceClonetracking(GitHub_path)