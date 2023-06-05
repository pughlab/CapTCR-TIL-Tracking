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
  
  # Generating Table S3 of the supplementary material
  message("Creating Table S3")
  TableS3_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS3_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS3_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 72,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS3_df$Patient)){
    for(column in 2:length(colnames(TableS3_df))){
      value <- CDR3perVJ(TableS3_df$Patient[row], "DNA", "TRB", colnames(TableS3_df)[column])
      TableS3_df[row, column] <- value
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    if(TableS3_df$Patient[row] == "TLML_1_"){
      value <- CDR3perVJ(TableS3_df$Patient[row], "DNA", "TRB", colnames(TableS3_df)[5])
      TableS3_df[row, 4] <- value
    }
  }
  close(pb)
  
  # Generating Table S4 of the supplementary material
  message("Creating Table S4")
  TableS4_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS4_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS4_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 72 ,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS4_df$Patient)){
    for(column in 2:length(colnames(TableS4_df))){
      if((column-1) <= length(eval(as.name(paste0(TableS4_df$Patient[row], "DNA_samporder"))))){
        value <- Convergence(TableS4_df$Patient[row], "DNA", "TRB", colnames(TableS4_df)[column])
        TableS4_df[row, column] <- value
      }
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    if(TableS4_df$Patient[row] == "TLML_1_"){
      value <- Convergence(TableS4_df$Patient[row], "DNA", "TRB", colnames(TableS4_df)[5])
      TableS4_df[row, 4] <- value
      value <- Convergence(TableS4_df$Patient[row], "DNA", "TRB", colnames(TableS4_df[7]))
      TableS4_df[row, 7] <- value
    }
  }
  close(pb)
  
  # Generating Table S5 of the supplementary material
  message("Creating Table S5")
  TableS5_df <- data.frame(matrix(NA, nrow=9, ncol=3))
  colnames(TableS5_df) <- c("Patient","baseline", "4W")
  TableS5_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 18,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS5_df$Patient)){
    samples <- eval(as.name(paste0(TableS5_df$Patient[row], "DNA_samporder")))
    baseline <- Persistence(TableS5_df$Patient[row], "DNA", "TRB", samples[1])
    i <- i + 1
    setTxtProgressBar(pb, i) 
    fweek <- Persistence(TableS5_df$Patient[row], "DNA", "TRB", samples[3])
    i <- i + 1
    setTxtProgressBar(pb, i)  
    TableS5_df[row,2:3] <- c(baseline, fweek)
  }
  close(pb)
  
  # Generating Table S6 of the supplementry material
  message("Creating Table S6")
  TableS6_df <- data.frame(matrix(NA, nrow=9, ncol=5))
  colnames(TableS6_df) <- c("Patient","Amount of >5% clones baseline", "Amount of >0.5% clones baseline", "Amount of >5% clones infusion", "Amount of >5% clones 4W")
  TableS6_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 36,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS6_df$Patient)){
    samples <- eval(as.name(paste0(TableS5_df$Patient[row], "DNA_samporder")))
    base_hyper <- HighExpandedClone(TableS6_df$Patient[row], "DNA", "TRB", samples[1], 0.05)
    i <- i + 1
    setTxtProgressBar(pb, i) 
    base_large <- HighExpandedClone(TableS6_df$Patient[row], "DNA", "TRB", samples[1], 0.005)
    i <- i + 1
    setTxtProgressBar(pb, i) 
    infusion_hyper <- HighExpandedClone(TableS6_df$Patient[row], "DNA", "TRB", samples[2], 0.05)
    i <- i + 1
    setTxtProgressBar(pb, i) 
    fw_hyper <- HighExpandedClone(TableS6_df$Patient[row], "DNA", "TRB", samples[3], 0.05)
    i <- i + 1
    setTxtProgressBar(pb, i) 
    TableS6_df[row, 2:5] <- c(base_hyper, base_large, infusion_hyper, fw_hyper)
  }
  close(pb)

  # Generating Table S7 of the supplementary material
  message("Creating Table S7")
  TableS7_df <- data.frame(matrix(NA, nrow=6, ncol=7))
  colnames(TableS7_df) <- c("Patient", "TIL % in baseline", "TIL % in 4 week", "cfDNA % in gDNA baseline",
                            "cfDNA % in gDNA 4 week", "Baseline cfDNA clone count", "4 week cfDNA clone count")
  TableS7_df$Patient <- c("TLML_4_", "TLML_16_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = length(TableS7_df$Patient),style = 3,   width = 50,char = "=")
  for(row in 1:length(TableS7_df$Patient)){
    First3Time_cfDNA <- eval(as.name(paste(TableS7_df$Patient[row], "cfDNA", "_samporder", sep="")))
    cfDNA_calc(TableS7_df$Patient[row], "TRB", "TotalTIL&Background")
    TableS7_df[row,2] <- result$Base_TIL
    TableS7_df[row,3] <- result$FW_TIL
    cfDNA_calc(TableS7_df$Patient[row], "TRB", "gDNAOverlap")
    TableS7_df[row, 4] <- result$gDNA_Base
    TableS7_df[row, 5] <- result$gDNA_FW
    Timepoint_char(TableS7_df$Patient[row], "cfDNA", "TRB", 0, "TCRCount", First3Time_cfDNA[1])
    TableS7_df[row, 6] <- TCRCount
    Timepoint_char(TableS7_df$Patient[row], "cfDNA", "TRB", 0, "TCRCount", First3Time_cfDNA[2])
    TableS7_df[row, 7] <- TCRCount
    setTxtProgressBar(pb, row)
  }
  close(pb)

  # Generating Table S8 of the supplementary material
  message("Creating Table S8")
  TableS8_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS8_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS8_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 9,style = 3,   width = 50,char = "=")
  for(row in 1:length(TableS8_df$Patient)){
    samporder <- eval(as.name(paste0(TableS8_df$Patient[row], "DNA_samporder")))
    for(column in 2:ncol(TableS8_df)){
      if((column - 1) <= length(samporder)){
        value <- BergerParkerIndex(TableS8_df$Patient[row], "DNA", "TRB", samporder[column-1])
        TableS8_df[row, column] <- value
      }
      setTxtProgressBar(pb, row)
    }
    if(TableS8_df$Patient[row] == "TLML_1_"){
      TableS8_df[5:7] <- TableS8_df[4:6]
    }
  }
  close(pb)

  # Generating Table S9 of the supplementary material
  message("Creating Table S9")
  TableS9_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS9_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS9_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 72,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS9_df$Patient)){
    samporder <- eval(as.name(paste0(TableS9_df$Patient[row], "DNA_samporder")))
    for(column in 2:ncol(TableS9_df)){
      if((column - 1) <= length(samporder)){
        value <- FractionExpandedClones(TableS9_df$Patient[row], "DNA", "TRB", samporder[column-1])
        TableS9_df[row, column] <- value
      }
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    if(TableS9_df$Patient[row] == "TLML_1_"){
      TableS9_df[5:7] <- TableS9_df[4:6]
    }
  }
  close(pb)

  # Generating TableS10 of the supplementary material
  message("Creating Table S10")
  TableS10_df <- data.frame(matrix(NA, nrow=9, ncol=7))
  colnames(TableS10_df) <- c("Patient", "Expanded baseline", "Background baseline", "Not in baseline", "Expanded infusion", "Background infusion", "Not in infusion")
  TableS10_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 18,style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS10_df$Patient)){
    samporder <- eval(as.name(paste0(TableS9_df$Patient[row], "DNA_samporder")))  
    expanded_baseline <- Timepointreference_changepoint(TableS10_df$Patient[row], "DNA", "TRB", samporder[1], "Expanded")
    background_baseline <- Timepointreference_changepoint(TableS10_df$Patient[row], "DNA", "TRB", samporder[1], "Background")
    notin_baseline <- 1 - expanded_baseline - background_baseline
    i <- i + 1
    setTxtProgressBar(pb, i)
    expanded_infusion <- Timepointreference_changepoint(TableS10_df$Patient[row], "DNA", "TRB", samporder[2], "Expanded")
    background_infusion <- Timepointreference_changepoint(TableS10_df$Patient[row], "DNA", "TRB", samporder[2], "Background")
    notin_infusion <- 1 - expanded_infusion - background_infusion
    i <- i + 1
    setTxtProgressBar(pb, i)
    TableS10_df[row,2:7] <- c(expanded_baseline, background_baseline, notin_baseline, expanded_infusion, background_infusion, notin_infusion)
  }
  close(pb)
  
  # Generating Table S11 of the supplementary material
  message("Creating Table S11")
  TableS11_df <- data.frame(matrix(NA, nrow=9, ncol=4))
  colnames(TableS11_df) <- c("Patient", "Baseline/Infusion", "Baseline/4 week", "4 week/Infusion")
  TableS11_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = length(TableS9_df$Patient),style = 3,   width = 50,char = "=")
  for(row in 1:length(TableS11_df$Patient)){
    samporder <- eval(as.name(paste0(TableS11_df$Patient[row], "DNA_samporder")))
    baseline_infusion <- MorisitaIndex(TableS11_df$Patient[row], "DNA", "TRB", samporder[1], samporder[2])
    baseline_fw <- MorisitaIndex(TableS11_df$Patient[row], "DNA", "TRB", samporder[1], samporder[3])
    fw_infusion <- MorisitaIndex(TableS11_df$Patient[row], "DNA", "TRB", samporder[2], samporder[3])
    TableS11_df[row, 2:4] <- c(baseline_infusion, baseline_fw, fw_infusion)
    setTxtProgressBar(pb, row)
  }
  close(pb)
  
  # Generating Table S12 of the supplementary material
  message("Creating Table S12")
  TableS12_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS12_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS12_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 72, style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS12_df$Patient)){
    samporder <- eval(as.name(paste0(TableS12_df$Patient[row], "DNA_samporder")))
    for(column in 2:length(colnames(TableS12_df))){
      if((column - 1) <= length(samporder)){
        value <- 1/SimpsonIndex(TableS12_df$Patient[row], "DNA", "TRB", samporder[column-1])
        TableS12_df[row, column] <- value
      }
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    if(TableS12_df$Patient[row] == "TLML_1_"){
      TableS12_df[row, 5:7] <- TableS12_df[row, 4:6]
    }
  }
  close(pb)
  
  # Generating Table S13 of the supplementary material
  message("Creating Table S13")
  TableS13_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS13_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS13_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 72, style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS13_df$Patient)){
    samporder <- eval(as.name(paste0(TableS13_df$Patient[row], "DNA_samporder")))
    for(column in 2:length(colnames(TableS13_df))){
      if((column - 1) <= length(samporder)){
        value <- ShannonIndex(TableS13_df$Patient[row], "DNA", "TRB", samporder[column-1])
        TableS13_df[row, column] <- value
      }
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    if(TableS13_df$Patient[row] == "TLML_1_"){
      TableS13_df[row, 5:7] <- TableS13_df[row, 4:6]
    }
  }
  close(pb)
  
  # Generating Table S14 of the supplementary material
  message("Creating Table S14")
  TableS14_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS14_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS14_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 72, style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS14_df$Patient)){
    samporder <- eval(as.name(paste0(TableS14_df$Patient[row], "DNA_samporder")))
    for(column in 2:length(colnames(TableS14_df))){
      if((column - 1) <= length(samporder)){
        value <- Richness(TableS14_df$Patient[row], "DNA", "TRB", samporder[column-1])
        TableS14_df[row, column] <- value
      }
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    if(TableS14_df$Patient[row] == "TLML_1_"){
      TableS14_df[row, 5:7] <- TableS14_df[row, 4:6]
    }
  }
  close(pb)

  # Generating Table S15 of the supplementary material
  message("Creating Table S15")
  TableS15_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS15_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS15_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 72, style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS15_df$Patient)){
    samporder <- eval(as.name(paste0(TableS15_df$Patient[row], "DNA_samporder")))
    for(column in 2:length(colnames(TableS15_df))){
      if((column - 1) <= length(samporder)){
        value <- Eveness(TableS15_df$Patient[row], "DNA", "TRB", samporder[column-1])
        TableS15_df[row, column] <- value
      }
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    if(TableS15_df$Patient[row] == "TLML_1_"){
      TableS15_df[row, 5:7] <- TableS15_df[row, 4:6]
    }
  }
  close(pb)
  
  # Generating Table S16 of the supplementary material
  message("Creating Table S16")
  TableS16_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS16_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS16_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = 72, style = 3,   width = 50,char = "=")
  i <- 0
  for(row in 1:length(TableS16_df$Patient)){
    samporder <- eval(as.name(paste0(TableS16_df$Patient[row], "DNA_samporder")))
    for(column in 2:length(colnames(TableS16_df))){
      if((column - 1) <= length(samporder)){
        value <- BergerParkerIndex(TableS16_df$Patient[row], "DNA", "TRB", samporder[column-1])
        TableS16_df[row, column] <- value
      }
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    if(TableS16_df$Patient[row] == "TLML_1_"){
      TableS16_df[row, 5:7] <- TableS16_df[row, 4:6]
    }
  }
  close(pb)
  
  # Exporting supplementary data excel file
  message("Exporting supplementary excel file")
  pb <- txtProgressBar(min = 0,max = 14,style = 3,   width = 50,char = "=")
  library(openxlsx)
  SupplementaryData <- createWorkbook()
  addWorksheet(SupplementaryData, "Table S3")
  addWorksheet(SupplementaryData, "Table S4")
  addWorksheet(SupplementaryData, "Table S5")
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
  writeData(SupplementaryData, sheet="Table S3", x=TableS3_df)
  setTxtProgressBar(pb, 1)
  writeData(SupplementaryData, sheet="Table S4", x=TableS4_df)
  setTxtProgressBar(pb, 2)
  writeData(SupplementaryData, sheet="Table S5", x=TableS5_df)
  setTxtProgressBar(pb, 3)
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
  saveWorkbook(SupplementaryData, paste(GitHub_path, "results/SupplementaryData_reproduced.xlsx", sep=""))
  close(pb)
}
ReproduceSupplementaryData(GitHub_path)

# Reproduces all figures from the manuscript
ReproduceFigures <- function(GitHub_path){
  message("Creating clone tracking and clonality plots")
  pb <- txtProgressBar(min = 0,max = 6,style = 3,   width = 50,char = "=")
  # Fetching clone tracking and relative abundance/diversity plots
  low_exp <- c("TLML_18_", "TLML_4_", "TLML_7_", "TLML_20_")
  high_exp <- c("TLML_29_", "TLML_1_", "TLML_26_", "TLML_22_", "TLML_16_")
  alignment_fig(low_exp, "DNA", "TRB", "clonetrack", "Baseline", paste(GitHub_path, "data/",sep=""), "Clonetrack_low")
  setTxtProgressBar(pb, 1)
  alignment_fig(high_exp, "DNA", "TRB", "clonetrack", "Baseline", paste(GitHub_path, "data/",sep=""), "Clonetrack_high")
  setTxtProgressBar(pb, 2)
  alignment_fig(low_exp, "DNA", "TRB", "rel_div", "Baseline", paste(GitHub_path, "data/",sep=""), "RelDiv_low")
  setTxtProgressBar(pb, 3)
  alignment_fig(high_exp, "DNA", "TRB", "rel_div", "Baseline", paste(GitHub_path, "data/",sep=""), "RelDiv_high")
  setTxtProgressBar(pb, 4)
  
  # Reading exported clone tracking and relative abundance/diversity plots and exporting into figure format
  ClonetrackHigh <- image_read(paste(GitHub_path, "data/Clonetrack_high.png", sep=""))
  ClonetrackLow <- image_read(paste(GitHub_path, "data/Clonetrack_low.png", sep=""))
  RelDivHigh <- image_read(paste(GitHub_path, "data/RelDiv_high.png", sep=""))
  RelDivLow <- image_read(paste(GitHub_path, "data/RelDiv_low.png", sep=""))
  FigOverlay <- image_read(paste(GitHub_path, "data/Clone_RelDivFigOverlay.png", sep=""))
  setTxtProgressBar(pb, 5)
  
  png(paste(GitHub_path, "results/CloneTrack_RelDiv_reproduced.png", sep=""), width = 8000, height = 8000, units = "px")
  plot.new()
  plot.window(xlim=c(0, 8000), ylim=c(0, 8000))
  rasterImage(ClonetrackHigh, xleft=1361, ybottom=870, xright=4273  , ytop=4510)
  rasterImage(ClonetrackLow, xleft=1343 , ybottom=4754, xright=3162 , ytop=7664)
  rasterImage(RelDivLow, xleft=4726 , ybottom=4777 , xright=6540 , ytop=7679)
  rasterImage(RelDivHigh, xleft=4719 , ybottom=892 , xright=7586 , ytop=4453)
  rasterImage(FigOverlay, xleft=0 , ybottom=0 , xright=8000 , ytop=8000)
  dev.off()
  
  unlink(paste(GitHub_path, "data/Clonetrack_high.png", sep=""))
  unlink(paste(GitHub_path, "data/Clonetrack_low.png", sep=""))
  unlink(paste(GitHub_path, "data/RelDiv_high.png", sep=""))
  unlink(paste(GitHub_path, "data/RelDiv_low.png", sep=""))
  
  setTxtProgressBar(pb, 6)
  close(pb)
  
  message("Creating VJ Usage plots")
  pb <- txtProgressBar(min = 0,max = 11,style = 3,   width = 50,char = "=")
  # Fetching VJ usage plots
  timepoint_order <- c("baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  VJUsage_Step3("TLML_1_", "DNA", timepoint_order[1:6], 1, 6, dir_main)
  setTxtProgressBar(pb, 1)
  VJUsage_Step3("TLML_4_", "DNA", timepoint_order[1:4], 1, 4, dir_main)
  setTxtProgressBar(pb, 2)
  VJUsage_Step3("TLML_7_", "DNA", timepoint_order[1:5], 1, 5, dir_main)
  setTxtProgressBar(pb, 3)
  VJUsage_Step3("TLML_16_", "DNA", timepoint_order[1:8], 1, 8, dir_main)
  setTxtProgressBar(pb, 4)
  VJUsage_Step3("TLML_18_", "DNA", timepoint_order[1:3], 1, 3, dir_main)
  setTxtProgressBar(pb, 5)
  VJUsage_Step3("TLML_20_", "DNA", timepoint_order[1:5], 1, 5, dir_main)
  setTxtProgressBar(pb, 6)
  VJUsage_Step3("TLML_22_", "DNA", timepoint_order[1:6], 1, 6, dir_main)
  setTxtProgressBar(pb, 7)
  VJUsage_Step3("TLML_26_", "DNA", timepoint_order[1:5], 1, 5, dir_main)
  setTxtProgressBar(pb, 8)
  VJUsage_Step3("TLML_29_", "DNA", timepoint_order[1:3], 1, 3, dir_main)
  setTxtProgressBar(pb, 9)
  
  # Reading exported clone tracking and relative abundance/diversity plots and exporting into figure format
  VJTLML_1_ <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_1__.png", sep=""))
  VJTLML_4_ <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_4__.png", sep=""))
  VJTLML_7_ <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_7__.png", sep=""))
  VJTLML_16 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_16__.png", sep=""))
  VJTLML_18 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_18__.png", sep=""))
  VJTLML_20 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_20__.png", sep=""))
  VJTLML_22 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_22__.png", sep=""))
  VJTLML_26 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_26__.png", sep=""))
  VJTLML_29 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_29__.png", sep=""))
  FigOverlay <- image_read(paste(GitHub_path, "data/VJFigOverlay.png", sep=""))
  setTxtProgressBar(pb, 10)
  
  png(paste(GitHub_path, "results/VJUsage_reproduced.png", sep=""), width = 5500, height = 6000, units = "px")
  plot.new()
  plot.window(xlim=c(0, 5500), ylim=c(0, 6000))
  rasterImage(VJTLML_1_, xleft=920, ybottom=1830, xright=4205  , ytop=2378)
  rasterImage(VJTLML_4_, xleft=920, ybottom=4730, xright=3112  , ytop=5278)
  rasterImage(VJTLML_7_, xleft=920, ybottom=4182, xright=3658  , ytop=4730)
  rasterImage(VJTLML_16, xleft=920, ybottom=726, xright=5304  , ytop=1274)
  rasterImage(VJTLML_18, xleft=920, ybottom=5270, xright=2563  , ytop=5818)
  rasterImage(VJTLML_20, xleft=920, ybottom=3634, xright=3658  , ytop=4182)
  rasterImage(VJTLML_22, xleft=920, ybottom=1257, xright=4205  , ytop=1805)
  rasterImage(VJTLML_26, xleft=920, ybottom=2366, xright=3658  , ytop=2914)
  rasterImage(VJTLML_29, xleft=920, ybottom=2909, xright=2563  , ytop=3457)
  rasterImage(FigOverlay, xleft=0 , ybottom=0 , xright=5500 , ytop=6000)
  dev.off()
  
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_1__.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_4__.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_7__.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_16__.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_18__.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_20__.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_22__.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_26__.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_29__.png", sep=""))
  setTxtProgressBar(pb, 11)
  close(pb)
  
  message("Creating cfDNA plots")
  pb <- txtProgressBar(min = 0,max = 10,style = 3,   width = 50,char = "=")
  # Fetching cfDNA plots
  low_exp <- c("TLML_4_", "TLML_20_")
  high_exp <- c("TLML_29_", "TLML_26_", "TLML_22_", "TLML_16_")
  alignment_fig(low_exp, "cfDNA", "TRB", "clonetrack_cfDNA", "Baseline", paste(GitHub_path, "data/",sep=""), "clonetrack_cfDNA_lowexp")
  setTxtProgressBar(pb, 1)
  alignment_fig(high_exp, "cfDNA", "TRB", "clonetrack_cfDNA", "Baseline", paste(GitHub_path, "data/",sep=""), "clonetrack_cfDNA_highexp")
  setTxtProgressBar(pb, 2)
  alignment_fig(low_exp, "DNA", "TRB", "clonetrack_cfDNA", "Baseline", paste(GitHub_path, "data/",sep=""), "clonetrack_gDNA_lowexp")
  setTxtProgressBar(pb, 3)
  alignment_fig(high_exp, "DNA", "TRB", "clonetrack_cfDNA", "Baseline", paste(GitHub_path, "data/",sep=""), "clonetrack_gDNA_highexp")
  setTxtProgressBar(pb, 4)
  alignment_fig(low_exp, "cfDNA", "TRB", "cfDNAcorrel", "Baseline", paste(GitHub_path, "data/",sep=""), "scatter_TIL_lowexp", clones ="TIL")
  setTxtProgressBar(pb, 5)
  alignment_fig(high_exp, "cfDNA", "TRB", "cfDNAcorrel", "Baseline", paste(GitHub_path, "data/",sep=""), "scatter_TIL_highexp", clones = "TIL")
  setTxtProgressBar(pb, 6)
  alignment_fig(low_exp, "cfDNA", "TRB", "cfDNAcorrel", "Baseline", paste(GitHub_path, "data/",sep=""), "scatter_back_lowexp", clones = "Background")
  setTxtProgressBar(pb, 7)
  alignment_fig(high_exp, "cfDNA", "TRB", "cfDNAcorrel", "Baseline", paste(GitHub_path, "data/",sep=""), "scatter_back_highexp", clones = "Background")
  setTxtProgressBar(pb, 8)
  
  # Reading and exporting cfDNA plots
  clonetrack_cfDNA_lowexp <- image_read(paste0(GitHub_path, "data/clonetrack_cfDNA_lowexp.png"))
  clonetrack_cfDNA_highexp <- image_read(paste0(GitHub_path, "data/clonetrack_cfDNA_highexp.png"))
  clonetrack_gDNA_lowexp <- image_read(paste0(GitHub_path, "data/clonetrack_gDNA_lowexp.png"))
  clonetrack_gDNA_highexp <- image_read(paste0(GitHub_path, "data/clonetrack_gDNA_highexp.png"))
  scatter_TIL_lowexp <- image_read(paste0(GitHub_path, "data/scatter_TIL_lowexp.png"))
  scatter_TIL_highexp <- image_read(paste0(GitHub_path, "data/scatter_TIL_highexp.png"))
  scatter_back_lowexp <- image_read(paste0(GitHub_path, "data/scatter_back_lowexp.png"))
  scatter_back_highexp <- image_read(paste0(GitHub_path, "data/scatter_back_highexp.png"))
  FigOverlay <- image_read(paste0(GitHub_path, "data/cfDNAFigOverlay.png"))
  setTxtProgressBar(pb, 9)
  
  png(paste(GitHub_path, "results/cfDNA_reproduced.png", sep=""), width = 11000, height = 6000, units = "px")
  plot.new()
  plot.window(xlim=c(0, 11000), ylim=c(0, 6000))
  rasterImage(clonetrack_cfDNA_lowexp, xleft=1321, ybottom=4224, xright=2043  , ytop=5668)
  rasterImage(clonetrack_cfDNA_highexp, xleft=1321, ybottom=945, xright=2050  , ytop=3861)
  rasterImage(clonetrack_gDNA_lowexp, xleft=2614, ybottom=4224, xright=4436  , ytop=5682)
  rasterImage(clonetrack_gDNA_highexp, xleft=2625, ybottom=941, xright=5532  , ytop=3847)
  rasterImage(clonetrack_cfDNA_lowexp, xleft=1321, ybottom=4224, xright=2043  , ytop=5668)
  rasterImage(scatter_TIL_lowexp, xleft=5988, ybottom=4246, xright=8132  , ytop=5676)
  rasterImage(scatter_TIL_highexp, xleft=5988, ybottom=945, xright=8126 , ytop=3796)
  rasterImage(scatter_back_lowexp, xleft=8643, ybottom=4274, xright=10784 , ytop=5701)
  rasterImage(scatter_back_highexp, xleft=8636, ybottom=924, xright=10784 , ytop=3788)
  rasterImage(FigOverlay, xleft=0 , ybottom=0 , xright=11000 , ytop=6000)
  dev.off()
  
  unlink(paste(GitHub_path, "data/clonetrack_cfDNA_lowexp.png", sep=""))
  unlink(paste(GitHub_path, "data/clonetrack_cfDNA_highexp.png", sep=""))
  unlink(paste(GitHub_path, "data/clonetrack_gDNA_lowexp.png", sep=""))
  unlink(paste(GitHub_path, "data/clonetrack_gDNA_highexp.png", sep=""))
  unlink(paste(GitHub_path, "data/scatter_TIL_lowexp.png", sep=""))
  unlink(paste(GitHub_path, "data/scatter_TIL_highexp.png", sep=""))
  unlink(paste(GitHub_path, "data/scatter_back_lowexp.png", sep=""))
  unlink(paste(GitHub_path, "data/scatter_back_highexp.png", sep=""))
  setTxtProgressBar(pb, 10)
  close(pb)
}
ReproduceFigures(GitHub_path)
