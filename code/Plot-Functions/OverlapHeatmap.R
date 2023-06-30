###################
# Overlap heatmap #
###################

overlap <- function(outputmatrix, GitHub_path) {
  annotation <- as.data.frame(matrix(NA, nrow = nrow(outputmatrix)))
  rownames(annotation) <- rownames(outputmatrix)
  annotation$Patient <- NA
  annotation$Cohort <- NA
  annotation$Timepoint <- NA
  for(i in 1:length(rownames(annotation))){
    split_r <- strsplit(rownames(annotation)[i], "_")[[1]]
    annotation$Patient[i] <- paste0("TLML_", split_r[2], "_")
    annotation$Cohort[i] <- split_r[3]
    annotation$Timepoint[i] <- split_r[4]
    if(annotation$Timepoint[i] == "4W"){
      annotation$Timepoint[i] <- "FW"
    }
    else if(annotation$Timepoint[i] == "FU"){
      annotation$Timepoint[i] <- paste0("FU_", split_r[5])
    }
    else if(annotation$Timepoint[i]=="12W"){
      annotation$Timepoint[i] <- "FU_01"
    }
    else if(annotation$Timepoint[i]=="17W"){
      annotation$Timepoint[i] <- "FU_02"
    }
  }
  annotation <- annotation[c(4:2)]
  ann_colors <- list(
    Patient = c(TLML_1_="#FFFF00", TLML_16_="#1CE6FF", TLML_18_="#FF34FF", TLML_20_="#FF4A46", 
                TLML_22_="#008941", TLML_26_="#006FA6", TLML_29_="#A30059", TLML_4_="#FFDBE5", TLML_7_="#7A4900"),
    Cohort = c(DNA="#7FC97F", RNA="#BEAED4", cfDNA="#FDC086"),
    Timepoint = c(baseline="yellow", infusion="blue", FW="darkred", FU_01="#B30000", FU_02="#E34A33", FU_03="#FC8D59", FU_04="#FDCC8A", FU_05="#FEF0D9", FU_07="white", FU_10="white", FU_12="white", FU_13="white", FU_16="white")
  )
  plot <- pheatmap(outputmatrix, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
           gaps_row=head(as.numeric(cumsum(table(annotation$Patient)[c(5,6,7,1,2,9,4,8,3)])), -1), gaps_col=head(as.numeric(cumsum(table(annotation$Patient)[c(5,6,7,1,2,9,4,8,3)])), -1),
           annotation_col=annotation, annotation_row=annotation, annotation_colors = ann_colors, cellheight=6, cellwidth = 6)

  png(file = paste0(GitHub_path, "data/heatmapraw.png"), width = 4000, height = 3244, units = "px", res=300)
  print(plot)
  dev.off() 
  
  heatmap <- image_read(paste0(GitHub_path, "data/heatmapraw.png"))
  FigOverlay <- image_read(paste(GitHub_path, "data/Overlay/heatmap_overlay.png", sep=""))
  whitespace <- image_read(paste(GitHub_path, "data/Overlay/whitespace.png", sep=""))
  
  png(paste(GitHub_path, "data/heatmap.png", sep=""), width = 4000, height = 3244, units = "px")
  plot.new()
  plot.window(xlim=c(0, 4000), ylim=c(0, 3244))
  rasterImage(heatmap, xleft=-247 , ybottom=-64 , xright=3753 , ytop=3180)
  rasterImage(FigOverlay, xleft=0 , ybottom=0 , xright=4000 , ytop=3244)
  rasterImage(whitespace, xleft=2948, ybottom=-64, xright=3537, ytop=300)
  dev.off()
  
  unlink(paste(GitHub_path, "data/heatmapraw.png", sep=""))
}
