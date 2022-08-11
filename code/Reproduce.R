# Produces serperate files from mixcr output file

# Loading packages and scripts
GitHub_path <- "/Users/cameronkerr/CapTCRRepo/CapTCR-TIL-Tracking/"
dir_main <- paste(GitHub_path, "data/", sep="")
file_samplekeys <- "TLML_samples_keys.xlsx"

# Required packages:
install.packages("ggplot2")
install.packages("ggalluvial")
install.packages("readxl")
install.packages("openxlsx")
install.packages("cowplot")
install.packages("tidyverse")
install.packages("ggraph")
install.packages("igraph")
install.packages("gridExtra")
install.packages("bioseq")
install.packages("magick")

# Generating list of relevant R scripts from GitHub repo and then loading in functions and libraries
R_scripts <- list.files(path=paste(GitHub_path, "code/", sep=""), pattern="\\.R$", recursive=TRUE)
R_scripts <- R_scripts[!grepl("Quality-Control|VJUsage_Step1|VJUsage_Step2|Reproduce", R_scripts)]
for(script in R_scripts){
  source(paste(GitHub_path, "code/", script, sep=""))
}

# Reproduces all data from the supplementary data 
ReproduceSupplementaryData <- function(GitHub_path){
  
  # Generating Table S3 of the supplementary material
  TableS3_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS3_df) <- c("Patient","baseline_apheresis", "infusion", "4W", "FU1", "FU2", "FU3", "FU4", "FU5")
  TableS3_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16", "TLML_18", "TLML_20", "TLML_22", "TLML_26", "TLML_29")
  for(row in 1:length(TableS3_df$Patient)){
    CDR3perVJ(TableS3_df$Patient[row], "DNA", colnames(TableS3_df)[2:ncol(TableS3_df)], "VJTreemaps.csv", dir_main)
    sorted_df <- total_df %>% arrange(factor(Cycle, levels = colnames(TableS3_df)[2:ncol(TableS3_df)]))
    TableS3_df[row, 2:length(TableS3_df)] <- sorted_df$AverageCDR3
    if(TableS3_df$Patient[row] == "TLML_1_"){
      TableS3_df[row, 2:length(TableS3_df)] <- c(sorted_df$AverageCDR3[1:2], NA, sorted_df$AverageCDR3[3:(length(sorted_df$AverageCDR3)-1)])
    }
  }
  # Generating Table S4 of the supplementary material
  TableS4_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS4_df) <- c("Patient","baseline_apheresis", "infusion", "4W", "FU1", "FU2", "FU3", "FU4", "FU5")
  TableS4_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16", "TLML_18", "TLML_20", "TLML_22", "TLML_26", "TLML_29")
  for(row in 1:length(TableS4_df$Patient)){
    TCRConvergence(TableS4_df$Patient[row], "DNA", colnames(TableS4_df)[2:ncol(TableS4_df)], "VJTreemaps.csv", dir_main)
    if(TableS4_df$Patient[row] != "TLML_1_"){
      sorted_df <- TCRConverg[colnames(TableS4_df)[2:(length(TCRConverg)+1)]]
      TableS4_df[row, 2:(length(sorted_df)+1)] <- sorted_df
    }
    if(TableS4_df$Patient[row] == "TLML_1_"){
      sorted_df <- TCRConverg[c("baseline_apheresis", "infusion", colnames(TableS4_df)[5:(length(TCRConverg)+2)])]
      TableS4_df[row, 2:(length(sorted_df)+2)] <- c(sorted_df[1], sorted_df[2], NA, sorted_df[3], sorted_df[4], sorted_df[5])
    }
  }
  
  # Generating Table S5 of the supplementary material
  TableS5_df <- data.frame(matrix(NA, nrow=9, ncol=4))
  colnames(TableS5_df) <- c("Patient","baseline_apheresis", "4W", "4W_expanded")
  TableS5_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16", "TLML_18", "TLML_20", "TLML_22", "TLML_26", "TLML_29")
  for(row in 1:length(TableS5_df$Patient)){
    First3Time <- eval(as.name(paste(TableS5_df$Patient[row], "gDNA", sep="")))[1:3]
    TIL_calc(TableS5_df$Patient[row], "gDNA", "TRB", 0, First3Time[1], expanded=FALSE, dir_main, dir_main, file_samplekeys)
    TableS5_df[row, 2] <- response
    print(response)
    TIL_calc(TableS5_df$Patient[row], "gDNA", "TRB", 0, First3Time[3], expanded=FALSE, dir_main, dir_main, file_samplekeys)
    TableS5_df[row, 3] <- response
    print(response)
    TIL_calc(TableS5_df$Patient[row], "gDNA", "TRB", 0, First3Time[3], expanded=TRUE, dir_main, dir_main, file_samplekeys)
    TableS5_df[row, 4] <- response
    print(response)
    }
  
  # Generating Table S6 of the supplementry material
  TableS6_df <- data.frame(matrix(NA, nrow=9, ncol=5))
  colnames(TableS6_df) <- c("Patient","Amount of >5% clones baseline", "Amount of >0.5% clones baseline", "Amount of >5% clones infusion", "Amount of >5% clones 4W")
  TableS6_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16", "TLML_18", "TLML_20", "TLML_22", "TLML_26", "TLML_29")
  for(row in 1:length(TableS6_df$Patient)){
    First3Time <- eval(as.name(paste(TableS6_df$Patient[row], "gDNA", sep="")))[1:3]
    Timepoint_char(TableS6_df$Patient[row], "gDNA", "TRB", 0.05, "Oligo", First3Time[1], dir_main, dir_main, file_samplekeys)
    TableS6_df[row, 2] <- Sum_Oligo
    Timepoint_char(TableS6_df$Patient[row], "gDNA", "TRB", 0.005, "Oligo", First3Time[1], dir_main, dir_main, file_samplekeys)
    TableS6_df[row, 3] <- Sum_Oligo
    Timepoint_char(TableS6_df$Patient[row], "gDNA", "TRB", 0.05, "Oligo", First3Time[2], dir_main, dir_main, file_samplekeys)
    TableS6_df[row, 4] <- Sum_Oligo
    Timepoint_char(TableS6_df$Patient[row], "gDNA", "TRB", 0.05, "Oligo", First3Time[3], dir_main, dir_main, file_samplekeys)
    TableS6_df[row, 5] <- Sum_Oligo
  }
  
  # Generating Table S7 of the supplementary material
  TableS7_df <- data.frame(matrix(NA, nrow=9, ncol=8))
  colnames(TableS7_df) <- c("Patient", "Top 50% of TIL to 4W", "Top 50% of 4W to TIL", "Top 50% of BL to 4W",
                            "Top 50% of 4W to BL", "Average basline clone fraction of top 10 4W clones",
                            "Average infusion clone fraction of top 10 4W clones", "Average 4W clone fraction of top 10 clones")
  TableS7_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16", "TLML_18", "TLML_20", "TLML_22", "TLML_26", "TLML_29")
  for(row in 1:length(TableS7_df$Patient)){
    First3Time <- eval(as.name(paste(TableS7_df$Patient[row], "gDNA", sep="")))[1:3]
    Top50.fx(TableS7_df$Patient[row], "gDNA", "TRB", 0, First3Time[2], First3Time[3], dir_main, dir_main, file_samplekeys)
    TableS7_df[row, 2] <- ReferenceFrc
    Top50.fx(TableS7_df$Patient[row], "gDNA", "TRB", 0, First3Time[3], First3Time[2], dir_main, dir_main, file_samplekeys)
    TableS7_df[row, 3] <- ReferenceFrc
    Top50.fx(TableS7_df$Patient[row], "gDNA", "TRB", 0, First3Time[1], First3Time[3], dir_main, dir_main, file_samplekeys)
    TableS7_df[row, 4] <- ReferenceFrc
    Top50.fx(TableS7_df$Patient[row], "gDNA", "TRB", 0, First3Time[3], First3Time[1], dir_main, dir_main, file_samplekeys)
    TableS7_df[row, 5] <- ReferenceFrc
    Top104WComp(TableS7_df$Patient[row], "gDNA", "TRB", 0, dir_main, dir_main, file_samplekeys)
    TableS7_df[row, 6] <- Base_cloneFraction
    TableS7_df[row, 7] <- TIL_cloneFraction
    TableS7_df[row, 8] <- Top104W_cloneFraction
  }
  
  # Generating Table S8 of the supplementary material
  TableS8_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS8_df) <- c("Patient", "Apheresis / Baseline", "TIL Infusion Product", "4 Week sample", "FU1", "FU2", "FU3", "FU4", "FU5")
  TableS8_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16", "TLML_18", "TLML_20", "TLML_22", "TLML_26", "TLML_29")
  for(row in 1:length(TableS8_df$Patient)){
    DivCalc(TableS8_df$Patient[row], "gDNA", "TRB", 0, dir_main, dir_main, file_samplekeys)
    if(TableS8_df$Patient[row] != "TLML_1_"){
      TableS8_df[row,2:(length(Div_df$Diversity)+1)] <- Div_df$Diversity
    }
    if(TableS8_df$Patient[row] == "TLML_1_"){
      TableS8_df[row,2:(length(Div_df$Diversity)+2)] <- c(Div_df$Diversity[1], Div_df$Diversity[2], NA, Div_df$Diversity[3], Div_df$Diversity[4], Div_df$Diversity[5])
    }
    
  }
  
  # Generating Table S9, S10, S11, and S12 of the supplementary material
  TableS9_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS9_df) <- c("Patient", "Apheresis / Baseline", "TIL Infusion Product", "4 Week sample", "FU1", "FU2", "FU3", "FU4", "FU5")
  TableS9_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16", "TLML_18", "TLML_20", "TLML_22", "TLML_26", "TLML_29")
  TableS10_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS10_df) <- c("Patient", "Apheresis / Baseline", "TIL Infusion Product", "4 Week sample", "FU1", "FU2", "FU3", "FU4", "FU5")
  TableS10_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16", "TLML_18", "TLML_20", "TLML_22", "TLML_26", "TLML_29")
  TableS11_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS11_df) <- c("Patient", "Apheresis / Baseline", "TIL Infusion Product", "4 Week sample", "FU1", "FU2", "FU3", "FU4", "FU5")
  TableS11_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16", "TLML_18", "TLML_20", "TLML_22", "TLML_26", "TLML_29")
  TableS12_df <- data.frame(matrix(NA, nrow=9, ncol=9))
  colnames(TableS12_df) <- c("Patient", "Apheresis / Baseline", "TIL Infusion Product", "4 Week sample", "FU1", "FU2", "FU3", "FU4", "FU5")
  TableS12_df$Patient <- c("TLML_1_", "TLML_4_", "TLML_7_", "TLML_16", "TLML_18", "TLML_20", "TLML_22", "TLML_26", "TLML_29")
  for(row in 1:length(TableS9_df$Patient)){
    RichEvCalc(TableS9_df$Patient[row], "gDNA", "TRB", 0, dir_main, dir_main, file_samplekeys)
    Div_df$Diversity <- abs(Div_df$Diversity)
    if(TableS9_df$Patient[row] != "TLML_1_"){
      TableS9_df[row, 2:(length(Div_df$Diversity)+1)] <- Div_df$Diversity
      TableS10_df[row, 2:(length(Div_df$Richness)+1)] <- Div_df$Richness
      TableS11_df[row, 2:(length(Div_df$Evenness)+1)] <- Div_df$Evenness
      TableS12_df[row, 2:(length(Div_df$Clonality)+1)] <- Div_df$Clonality
    }
    if(TableS9_df$Patient[row] == "TLML_1_"){
      TableS9_df[row, 2:(length(Div_df$Diversity)+2)] <- c(Div_df$Diversity[1], Div_df$Diversity[2], NA, Div_df$Diversity[3], Div_df$Diversity[4], Div_df$Diversity[5])
      TableS10_df[row, 2:(length(Div_df$Richness)+2)] <-c(Div_df$Richness[1], Div_df$Richness[2], NA, Div_df$Richness[3], Div_df$Richness[4], Div_df$Richness[5])
      TableS11_df[row, 2:(length(Div_df$Evenness)+2)] <- c(Div_df$Evenness[1], Div_df$Evenness[2], NA, Div_df$Evenness[3], Div_df$Evenness[4], Div_df$Evenness[5])
      TableS12_df[row, 2:(length(Div_df$Clonality)+2)] <- c(Div_df$Clonality[1], Div_df$Clonality[2], NA, Div_df$Clonality[3], Div_df$Clonality[4], Div_df$Clonality[5])
    }
  }
  
  # Exporting supplementary data excel file
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
  writeData(SupplementaryData, sheet="Table S3", x=TableS3_df)
  writeData(SupplementaryData, sheet="Table S4", x=TableS4_df)
  writeData(SupplementaryData, sheet="Table S5", x=TableS5_df)
  writeData(SupplementaryData, sheet="Table S6", x=TableS6_df)
  writeData(SupplementaryData, sheet="Table S7", x=TableS7_df)
  writeData(SupplementaryData, sheet="Table S8", x=TableS8_df)
  writeData(SupplementaryData, sheet="Table S9", x=TableS9_df)
  writeData(SupplementaryData, sheet="Table S10", x=TableS10_df)
  writeData(SupplementaryData, sheet="Table S11", x=TableS11_df)
  writeData(SupplementaryData, sheet="Table S12", x=TableS12_df)
  saveWorkbook(SupplementaryData, paste(GitHub_path, "results/SupplementaryData_reproduced.xlsx", sep=""))
}
ReproduceSupplementaryData(GitHub_path)

# Reproduces all figures from the manuscript
ReproduceFigures <- function(GitHub_path){
  library(magick)
  # Fetching clone tracking and relative abundance/diversity plots
  low_exp <- c("TLML_18", "TLML_4_", "TLML_7_", "TLML_20")
  high_exp <- c("TLML_29", "TLML_1_", "TLML_26", "TLML_22", "TLML_16")
  alignment_fig(low_exp, "gDNA", "TRB", 0, dir_main, dir_main, file_samplekeys,
                "clonetrack", "Baseline", paste(GitHub_path, "data/",sep=""), "Clonetrack_low")
  alignment_fig(high_exp, "gDNA", "TRB", 0, dir_main, dir_main, file_samplekeys,
                "clonetrack", "Baseline", paste(GitHub_path, "data/",sep=""), "Clonetrack_high")
  alignment_fig(low_exp, "gDNA", "TRB", 0, dir_main, dir_main, file_samplekeys,
                "rel_div", "Baseline", paste(GitHub_path, "data/",sep=""), "RelDiv_low")
  alignment_fig(high_exp, "gDNA", "TRB", 0, dir_main, dir_main, file_samplekeys,
                "rel_div", "Baseline", paste(GitHub_path, "data/",sep=""), "RelDiv_high")
  
  # Reading exported clone tracking and relative abundance/diversity plots and exporting into figure format
  ClonetrackHigh <- image_read(paste(GitHub_path, "data/Clonetrack_high.png", sep=""))
  ClonetrackLow <- image_read(paste(GitHub_path, "data/Clonetrack_low.png", sep=""))
  RelDivHigh <- image_read(paste(GitHub_path, "data/RelDiv_high.png", sep=""))
  RelDivLow <- image_read(paste(GitHub_path, "data/RelDiv_low.png", sep=""))
  FigOverlay <- image_read(paste(GitHub_path, "data/Clone_RelDivFigOverlay.png", sep=""))
  
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
  
  # Fetching VJ usage plots
  timepoint_order <- c("baseline_apheresis", "infusion", "4W", "FU1", "FU2", "FU3", "FU4", "FU5")
  VJUsage_Step3("TLML_1_", "DNA", timepoint_order[1:6], 1, 6, dir_main, "VJUsage.csv", 
                dir_main, "VJTreemaps.csv", dir_main, "CDR3_colors.csv", dir_main)
  VJUsage_Step3("TLML_4_", "DNA", timepoint_order[1:4], 1, 4, dir_main, "VJUsage.csv", 
                dir_main, "VJTreemaps.csv", dir_main, "CDR3_colors.csv", dir_main)
  VJUsage_Step3("TLML_7_", "DNA", timepoint_order[1:5], 1, 5, dir_main, "VJUsage.csv", 
                dir_main, "VJTreemaps.csv", dir_main, "CDR3_colors.csv", dir_main)
  VJUsage_Step3("TLML_16", "DNA", timepoint_order[1:8], 1, 8, dir_main, "VJUsage.csv", 
                dir_main, "VJTreemaps.csv", dir_main, "CDR3_colors.csv", dir_main)
  VJUsage_Step3("TLML_18", "DNA", timepoint_order[1:3], 1, 3, dir_main, "VJUsage.csv", 
                dir_main, "VJTreemaps.csv", dir_main, "CDR3_colors.csv", dir_main)
  VJUsage_Step3("TLML_20", "DNA", timepoint_order[1:5], 1, 5, dir_main, "VJUsage.csv", 
                dir_main, "VJTreemaps.csv", dir_main, "CDR3_colors.csv", dir_main)
  VJUsage_Step3("TLML_22", "DNA", timepoint_order[1:6], 1, 6, dir_main, "VJUsage.csv", 
                dir_main, "VJTreemaps.csv", dir_main, "CDR3_colors.csv", dir_main)
  VJUsage_Step3("TLML_26", "DNA", timepoint_order[1:5], 1, 5, dir_main, "VJUsage.csv", 
                dir_main, "VJTreemaps.csv", dir_main, "CDR3_colors.csv", dir_main)
  VJUsage_Step3("TLML_29", "DNA", timepoint_order[1:3], 1, 3, dir_main, "VJUsage.csv", 
                dir_main, "VJTreemaps.csv", dir_main, "CDR3_colors.csv", dir_main)
  
  # Reading exported clone tracking and relative abundance/diversity plots and exporting into figure format
  VJTLML_1_ <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_1__.png", sep=""))
  VJTLML_4_ <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_4__.png", sep=""))
  VJTLML_7_ <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_7__.png", sep=""))
  VJTLML_16 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_16_.png", sep=""))
  VJTLML_18 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_18_.png", sep=""))
  VJTLML_20 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_20_.png", sep=""))
  VJTLML_22 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_22_.png", sep=""))
  VJTLML_26 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_26_.png", sep=""))
  VJTLML_29 <- image_read(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_29_.png", sep=""))
  FigOverlay <- image_read(paste(GitHub_path, "data/VJFigOverlay.png", sep=""))
  
  
  
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
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_16_.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_18_.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_20_.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_22_.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_26_.png", sep=""))
  unlink(paste(GitHub_path, "data/TLML_VJUsage_DNA_TLML_29_.png", sep=""))

  }
ReproduceFigures(GitHub_path)
