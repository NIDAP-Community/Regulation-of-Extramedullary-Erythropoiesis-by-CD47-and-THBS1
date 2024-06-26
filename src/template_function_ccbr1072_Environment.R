# [R] Get Package Info  (bdc9cdf9-c835-42e1-a084-bcf72b756ec1): v7
ccbr1072_Environment <- function() {

    # Function to Extract Session Information to a Dataframe
extract_session_info <- function() {
  # Get session information
  si <- sessionInfo()
  
  # Initialize the dataframe
  df <- data.frame(
    Package = character(),
    Version = character(),
    Attached = logical(),
    stringsAsFactors = FALSE
  )
  
  # Extract information for each package loaded
  for (package_info in si$otherPkgs) {
    row <- data.frame(
      Package = package_info$Package,
      Version = package_info$Version,
      Attached = TRUE,
      stringsAsFactors = FALSE
    )
    df <- rbind(df, row)
  }
  
  # Extract information for the base packages
  for (package in si$basePkgs) {
    row <- data.frame(
      Package = package,
      Version = NA,  # Version information is not available for base packages
      Attached = TRUE,
      stringsAsFactors = FALSE
    )
    df <- rbind(df, row)
  }
  
  # Extract information for the packages loaded via a namespace
  for (package in names(si$loadedOnly)) {
    row <- data.frame(
      Package = package,
      Version = si$loadedOnly[[package]]$Version,
      Attached = FALSE,
      stringsAsFactors = FALSE
    )
    df <- rbind(df, row)
  }
  
  # Return the dataframe
  return(df)
}

# Extract session information to dataframe
session_info_df <- extract_session_info()

si <- sessionInfo()
print(si)

# View the dataframe
return(session_info_df)

}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_ccbr1072_Environment.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
invisible(graphics.off())
var_ccbr1072_Environment<-ccbr1072_Environment()
invisible(graphics.off())
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
saveRDS(var_ccbr1072_Environment, paste0(rds_output,"/var_ccbr1072_Environment.rds"))
