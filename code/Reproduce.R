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
  pb <- txtProgressBar(min = 0,max = length(TableS3_df$Patient),style = 3,   width = 50,char = "=")
  for(row in 1:length(TableS3_df$Patient)){
    CDR3perVJ(TableS3_df$Patient[row], "DNA", colnames(TableS3_df)[2:ncol(TableS3_df)])
    sorted_df <- total_df %>% arrange(factor(Cycle, levels = colnames(TableS3_df)[2:ncol(TableS3_df)]))
    if(TableS3_df$Patient[row] == "TLML_1_"){
      TableS3_df[row, 2:(length(sorted_df$Cycle)+2)] <- c(sorted_df$AverageCDR3[1:2], NA, sorted_df$AverageCDR3[3:(length(sorted_df$AverageCDR3))])
    } else {
      TableS3_df[row, 2:(length(sorted_df$Cycle)+1)] <- sorted_df$AverageCDR3
    }
    setTxtProgressBar(pb, row)
  }
  close(pb)
  
  # Generating Table S4 of the supplementary material
  message("Creating Table S4")
  TableS4_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS4_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS4_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = length(TableS4_df$Patient),style = 3,   width = 50,char = "=")
  for(row in 1:length(TableS4_df$Patient)){
    TCRConvergence(TableS4_df$Patient[row], "DNA", colnames(TableS4_df)[2:ncol(TableS4_df)])
    if(TableS4_df$Patient[row] != "TLML_1_"){
      sorted_df <- TCRConverg[colnames(TableS4_df)[2:(length(TCRConverg)+1)]]
      TableS4_df[row, 2:(length(sorted_df)+1)] <- sorted_df
    }
    if(TableS4_df$Patient[row] == "TLML_1_"){
      sorted_df <- TCRConverg[c("baseline", "infusion", colnames(TableS4_df)[5:(length(TCRConverg)+2)])]
      TableS4_df[row, 2:(length(sorted_df)+2)] <- c(sorted_df[1], sorted_df[2], NA, sorted_df[3], sorted_df[4], sorted_df[5])
    }
    setTxtProgressBar(pb, row)
  }
  close(pb)
  
  # Generating Table S5 of the supplementary material
  message("Creating Table S5")
  TableS5_df <- data.frame(matrix(NA, nrow=9, ncol=4))
  colnames(TableS5_df) <- c("Patient","baseline", "4W", "4W_expanded")
  TableS5_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = length(TableS5_df$Patient),style = 3,   width = 50,char = "=")
  for(row in 1:length(TableS5_df$Patient)){
    First3Time <- eval(as.name(paste(TableS5_df$Patient[row], "DNA", "_samporder", sep="")))[1:3]
    TIL_calc(TableS5_df$Patient[row], "DNA", "TRB", First3Time[1], expanded=FALSE)
    TableS5_df[row, 2] <- response
    TIL_calc(TableS5_df$Patient[row], "DNA", "TRB", First3Time[3], expanded=FALSE)
    TableS5_df[row, 3] <- response
    TIL_calc(TableS5_df$Patient[row], "DNA", "TRB", First3Time[3], expanded=TRUE)
    TableS5_df[row, 4] <- response
    setTxtProgressBar(pb, row)
  }
  close(pb)
  
  # Generating Table S6 of the supplementry material
  message("Creating Table S6")
  TableS6_df <- data.frame(matrix(NA, nrow=9, ncol=5))
  colnames(TableS6_df) <- c("Patient","Amount of >5% clones baseline", "Amount of >0.5% clones baseline", "Amount of >5% clones infusion", "Amount of >5% clones 4W")
  TableS6_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = length(TableS6_df$Patient),style = 3,   width = 50,char = "=")
  for(row in 1:length(TableS6_df$Patient)){
    First3Time <- eval(as.name(paste(TableS6_df$Patient[row], "DNA", "_samporder", sep="")))[1:3]
    Timepoint_char(TableS6_df$Patient[row], "DNA", "TRB", 0.05, "Oligo", First3Time[1])
    TableS6_df[row, 2] <- Sum_Oligo
    Timepoint_char(TableS6_df$Patient[row], "DNA", "TRB", 0.005, "Oligo", First3Time[1])
    TableS6_df[row, 3] <- Sum_Oligo
    Timepoint_char(TableS6_df$Patient[row], "DNA", "TRB", 0.05, "Oligo", First3Time[2])
    TableS6_df[row, 4] <- Sum_Oligo
    Timepoint_char(TableS6_df$Patient[row], "DNA", "TRB", 0.05, "Oligo", First3Time[3])
    TableS6_df[row, 5] <- Sum_Oligo
    setTxtProgressBar(pb, row)
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
  TableS8_df <- data.frame(matrix(NA, nrow=9, ncol=8))
  colnames(TableS8_df) <- c("Patient", "Top 50% of TIL to 4W", "Top 50% of 4W to TIL", "Top 50% of BL to 4W",
                            "Top 50% of 4W to BL", "Average basline clone fraction of top 10 4W clones",
                            "Average infusion clone fraction of top 10 4W clones", "Average 4W clone fraction of top 10 clones")
  TableS8_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = length(TableS8_df$Patient),style = 3,   width = 50,char = "=")
  for(row in 1:length(TableS8_df$Patient)){
    First3Time <- eval(as.name(paste(TableS8_df$Patient[row], "DNA", "_samporder", sep="")))[1:3]
    Top50.fx(TableS8_df$Patient[row], "DNA", "TRB", First3Time[2], First3Time[3])
    TableS8_df[row, 2] <- ReferenceFrc
    Top50.fx(TableS8_df$Patient[row], "DNA", "TRB", First3Time[3], First3Time[2])
    TableS8_df[row, 3] <- ReferenceFrc
    Top50.fx(TableS8_df$Patient[row], "DNA", "TRB", First3Time[1], First3Time[3])
    TableS8_df[row, 4] <- ReferenceFrc
    Top50.fx(TableS8_df$Patient[row], "DNA", "TRB", First3Time[3], First3Time[1])
    TableS8_df[row, 5] <- ReferenceFrc
    Top104WComp(TableS8_df$Patient[row], "DNA", "TRB")
    TableS8_df[row, 6] <- Base_cloneFraction
    TableS8_df[row, 7] <- TIL_cloneFraction
    TableS8_df[row, 8] <- Top104W_cloneFraction
    setTxtProgressBar(pb, row)
  }
  close(pb)
  
  # Generating Table S9 of the supplementary material
  message("Creating Table S9")
  TableS9_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS9_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS9_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = length(TableS9_df$Patient),style = 3,   width = 50,char = "=")
  for(row in 1:length(TableS9_df$Patient)){
    DivCalc(TableS9_df$Patient[row], "DNA", "TRB")
    if(TableS9_df$Patient[row] != "TLML_1_"){
      TableS9_df[row,2:(length(Div_df$Diversity)+1)] <- Div_df$Diversity
    }
    if(TableS9_df$Patient[row] == "TLML_1_"){
      TableS9_df[row,2:(length(Div_df$Diversity)+2)] <- c(Div_df$Diversity[1], Div_df$Diversity[2], NA, Div_df$Diversity[3], Div_df$Diversity[4], Div_df$Diversity[5])
    }
    setTxtProgressBar(pb, row)
  }
  close(pb)
  
  # Generating Table S10, S11, S12, and S13 of the supplementary material
  message("Creating Table S10-S13")
  TableS10_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS10_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS10_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  TableS11_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS11_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS11_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  TableS12_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS12_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS12_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  TableS13_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS13_df) <- c("Patient", "baseline", "infusion", "4W", "FU_01", "FU_02", "FU_03", "FU_04", "FU_05")
  TableS13_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16_", "TLML_18_", "TLML_20_", "TLML_22_", "TLML_26_", "TLML_29_")
  pb <- txtProgressBar(min = 0,max = length(TableS10_df$Patient),style = 3,   width = 50,char = "=")
  for(row in 1:length(TableS10_df$Patient)){
    RichEvCalc(TableS10_df$Patient[row], "DNA", "TRB")
    Div_df$Diversity <- abs(Div_df$Diversity)
    if(TableS9_df$Patient[row] != "TLML_1_"){
      TableS10_df[row, 2:(length(Div_df$Diversity)+1)] <- Div_df$Diversity
      TableS11_df[row, 2:(length(Div_df$Richness)+1)] <- Div_df$Richness
      TableS12_df[row, 2:(length(Div_df$Evenness)+1)] <- Div_df$Evenness
      TableS13_df[row, 2:(length(Div_df$Clonality)+1)] <- Div_df$Clonality
    }
    if(TableS10_df$Patient[row] == "TLML_1_"){
      TableS10_df[row, 2:(length(Div_df$Diversity)+2)] <- c(Div_df$Diversity[1], Div_df$Diversity[2], NA, Div_df$Diversity[3], Div_df$Diversity[4], Div_df$Diversity[5])
      TableS11_df[row, 2:(length(Div_df$Richness)+2)] <-c(Div_df$Richness[1], Div_df$Richness[2], NA, Div_df$Richness[3], Div_df$Richness[4], Div_df$Richness[5])
      TableS12_df[row, 2:(length(Div_df$Evenness)+2)] <- c(Div_df$Evenness[1], Div_df$Evenness[2], NA, Div_df$Evenness[3], Div_df$Evenness[4], Div_df$Evenness[5])
      TableS13_df[row, 2:(length(Div_df$Clonality)+2)] <- c(Div_df$Clonality[1], Div_df$Clonality[2], NA, Div_df$Clonality[3], Div_df$Clonality[4], Div_df$Clonality[5])
    }
    setTxtProgressBar(pb, row)
  }
  close(pb)
  
  # Exporting supplementary data excel file
  message("Exporting supplementary excel file")
  pb <- txtProgressBar(min = 0,max = 11,style = 3,   width = 50,char = "=")
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
