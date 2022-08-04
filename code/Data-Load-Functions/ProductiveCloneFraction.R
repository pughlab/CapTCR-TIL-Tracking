###########################################################################
# Creates new cloneFraction column based on a subset of the original data #
###########################################################################

# Readjusts clone fraction data based on spliced data frame
# @param input: Clone data loaded from the 'Load_data' function

ProdcloneFrac <- function(input){ 
    
    # Creating list of files used in input data
    filenames <- unique(input$filename)
    
    # Adjusting clone fraction data of a splied dataframe
    output_df <- data.frame()
    for(i in 1:length(filenames)){
        subset_fraction <- dplyr::filter(input, grepl(filenames[i], filename)) 
        totalClone <- sum(subset_fraction$cloneCount)
        subset_fraction$cloneFraction <- as.numeric(as.character(subset_fraction$cloneCount))/totalClone
        output_df <-rbind(output_df, subset_fraction)
    }
    
    # Creating global variable for further analysis
    output_df <<- output_df
}