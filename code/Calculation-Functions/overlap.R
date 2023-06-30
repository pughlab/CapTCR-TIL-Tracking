# Functions relating to calculating overlap and sample correlation
#
morisitaindex <- function(patient1, sampcohort1, timepoint1, 
                          patient2, sampcohort2, timepoint2) {
  # Computes the Morisita's overlap index between two samples
  #
  # Args:
  #   patient1: patient ID of sample 1
  #   sampcohort1: sample cohort (DNA, RNA, cfDNA) of sample 1
  #   timepoint1: full name of sample 1
  #   patient2: patient ID of sample 2
  #   sampcohort2: sample cohort (DNA, RNA, cfDNA) of sample 2
  #   timepoint2: full name of sample 2
  # Returns: 
  #   Morisita's overlap index  
  
  # Load data -------------------
  CDR3_df1 <- eval(as.name(paste0(patient1, sampcohort1, "_", timepoint1)))
  CDR3_df2 <- eval(as.name(paste0(patient2, sampcohort2, "_", timepoint2)))
  simpson1 <- SimpsonIndex(patient1, sampcohort1, timepoint1)
  simpson2 <- SimpsonIndex(patient2, sampcohort2, timepoint2)
  # Calculate index -------------------
  total_unique <- unique(c(CDR3_df1$aaSeqCDR3, CDR3_df2$aaSeqCDR3))
  sum_numerator <- 0
  for(CDR3 in total_unique){
    frac1 <- CDR3_df1[which(CDR3_df1$aaSeqCDR3==CDR3),]$cloneFraction
    frac2 <- CDR3_df2[which(CDR3_df2$aaSeqCDR3==CDR3),]$cloneFraction
    if(length(frac1)!=0 & length(frac2)!=0){
      sum_numerator <- sum_numerator + frac1*frac2
    }
  }
  2*sum_numerator/(simpson1+simpson2)
} 
get_overlapmatrix <- function(stat) {
  # Generates matrix of overlap between all samples
  #
  # Args:
  #   stat: "morisita" or "correlation"
  # Returns: 
  #   Overlap matrix 
  
  # Get sample names ------------------- 
  all_variables <- ls(pattern = "samporder", envir=.GlobalEnv)
  all_variables <- c(all_variables[11:13], all_variables[14:16], all_variables[17:19], all_variables[1:2], all_variables[3:5], all_variables[23:24], all_variables[8:10], all_variables[20:22], all_variables[6:7])
  samples <- c()
  for(x in all_variables){
    split_samp <- strsplit(x, "_")[[1]]
    patient <- paste0(split_samp[1], "_", split_samp[2], "_")
    cohort <- split_samp[3]
    samples <- c(samples, paste0(patient, cohort, "_", eval(as.name(x))))
  }
  # Generate overlap matrix -------------------
  n_samp <- length(samples)
  distancematrix <- as.data.frame(matrix(data = NA, nrow = n_samp, ncol=n_samp))
  colnames(distancematrix) <- samples
  rownames(distancematrix) <- samples
  pb <- txtProgressBar(min = 0,max = n_samp^2,style = 3,   width = 100,char = "=")
  i <- 0
  for(n in 1:n_samp){
    for(m in 1:n_samp){
      i <- i + 1
      setTxtProgressBar(pb, i)
      if(n < m){
        samp_n_split <- strsplit(samples[n], "_")[[1]]
        patient_n <- paste0(samp_n_split[1], "_", samp_n_split[2], "_")
        cohort_n <- samp_n_split[3]
        timepoint_n <- paste(samp_n_split[4:length(samp_n_split)], collapse="_")
        
        samp_m_split <- strsplit(samples[m], "_")[[1]]
        patient_m <- paste0(samp_m_split[1], "_", samp_m_split[2], "_")
        cohort_m <- samp_m_split[3]
        timepoint_m <- paste(samp_m_split[4:length(samp_m_split)], collapse="_")
        if(stat == "morisita") {
          value <- morisitaindex(patient_n, cohort_n, timepoint_n, patient_m, cohort_m, timepoint_m)
        }
        else if(stat == "correlation"){
          value <- samplecorrelation(patient_n, cohort_n, timepoint_n, patient_m, cohort_m, timepoint_m)
        }
        distancematrix[n,m] <- value
        distancematrix[m,n] <- value
      }
    }
  }
  # Return matrix -------------------
  close(pb)
  outputmatrix <- as.matrix(distancematrix)
  diag(outputmatrix) <- 1
  
  outputmatrix <<- outputmatrix
}
meanoverlap <- function(overlapmatrix, patient, sampcohorts, timepoints) {
  # Computes the mean Morisita's overlap over a category
  #
  # Args: 
  #   overlapmatrix: output from get_overlapmatrix
  #   patients: patient for comparison (empty string if any patient)
  #   sampcohorts: list of 1 or 2 sample cohorts (empty list if any cohort)
  #   timepoints: list of 1 or 2 timepoints (empty list if any timepoint)
  # Returns: 
  #   Mean overlap
  
  # Loop over all values -------------------
  values <- c()
  for(n in 1:nrow(overlapmatrix)){
    for(m in 1:nrow(overlapmatrix)){
      if(n > m){
        # Get patient, sampcohort, timepoint -------------------
        name_n <- rownames(outputmatrix)[n]
        name_m <- colnames(outputmatrix)[m]
        n_split <- strsplit(name_n, "_")[[1]]
        m_split <- strsplit(name_m, "_")[[1]]
        patient_n <- n_split[2]
        patient_m <- m_split[2]
        cohort_n <- n_split[3]
        cohort_m <- m_split[3]
        timepoint_n <- paste(n_split[4:length(n_split)], collapse="_")
        timepoint_m <- paste(m_split[4:length(m_split)], collapse="_")
        # Determine if value meets requiremnts -------------------
        # Patient requirements
        if(nchar(patient) > 0) {
          is_patient <- patient_n == patient_m & patient_n == patient
        }
        else {
          is_patient <- patient_n == patient_m
        }
        # Sample cohort requirements
        if(length(sampcohorts) == 2) {
          is_cohort <- cohort_n %in% sampcohorts & cohort_m %in% sampcohorts & cohort_n != cohort_m
        }
        else if(length(sampcohorts) == 1) {
          is_cohort <- cohort_n == cohort_m & cohort_n == sampcohorts[1]
        }
        else {
          is_cohort <- cohort_n == cohort_m
        }
        # Timepoint requirements
        if(length(timepoints) == 2) {
          is_timepoint <- timepoint_n %in% timepoints & timepoint_m %in% timepoints & timepoint_n != timepoint_m
        }
        else if(length(timepoints) == 1 & (timepoints[1] != "same" & timepoints[1] != "different")) {
          is_timepoint <- timepoint_n == timepoint_m & timepoint_n == timepoints[1]
        }
        else if(timepoints[1] == "same") {
          is_timepoint <- timepoint_n == timepoint_m
        }
        else if(timepoints[1] == "different"){
          is_timepoint <- TRUE
        }
        if(is_patient & is_cohort & is_timepoint){
          values <- c(values, overlapmatrix[n,m])
        }
      }
    }
  }
  mean(values)
}
samplecorrelation <- function(patient1, sampcohort1, timepoint1, patient2,
                              sampcohort2, timepoint2, TIL = FALSE, 
                              Background = FALSE) {
  # Computes the clone percentage correlation between two samples
  #
  # Args:
  #   patient1: patient ID of sample 1
  #   sampcohort1: sample cohort (DNA, RNA, cfDNA) of sample 1
  #   timepoint1: full name of sample 1
  #   patient2: patient ID of sample 2
  #   sampcohort2: sample cohort (DNA, RNA, cfDNA) of sample 2
  #   timepoint2: full name of sample 2
  #   TIL: TRUE if you want to analyze only redundant TILs
  #   Background: TRUE if you want to analyze clones that are not redundant TILs
  # Returns: 
  #   Clone percentage correlation
  
  # Load data -------------------
  CDR3_df1 <- eval(as.name(paste0(patient1, sampcohort1, "_", timepoint1)))
  CDR3_df2 <- eval(as.name(paste0(patient2, sampcohort2, "_", timepoint2)))
  total_unique <- unique(c(CDR3_df1$aaSeqCDR3, CDR3_df2$aaSeqCDR3))
  if(TIL) {
    total_unique <- CDR3_colors$colored_clns[which(CDR3_colors$id == patient1)]
  }
  else if(Background) {
    redTIL <- CDR3_colors$colored_clns[which(CDR3_colors$id == patient1)]
    total_unique <- total_unique[!(total_unique %in% redTIL)]
  }
  fraction_df <- data.frame("aaSeqCDR3" = total_unique, "cloneFrac1" = NA, "cloneFrac2" = NA)
  fraction_df$cloneFrac1 <- CDR3_df1$cloneFraction[match(fraction_df$aaSeqCDR3, CDR3_df1$aaSeqCDR3)]
  fraction_df$cloneFrac2 <- CDR3_df2$cloneFraction[match(fraction_df$aaSeqCDR3, CDR3_df2$aaSeqCDR3)]
  fraction_df[is.na(fraction_df)] <- 0
  # Calculate correlation -------------------
  cor(fraction_df$cloneFrac1, fraction_df$cloneFrac2) ^ 2
}
