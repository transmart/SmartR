#########################################################################################
## R-script to produce summary statistics (mean, median, etc.,) and static boxplot image
## for the data loading and data preprocessing tab
#########################################################################################
######
# Expected input data:
# Only data from high dimensional nodes. Either:
# * variable "loaded_variables" - a list of data.frames, containing one or more data.frames. Present as a variable in 
#                                 the global environment.
#   Data for multiple high dimensional nodes can be provided. Per high dimensional node 1 or 2 dataframes are to be
#   passed on, depending on whether 1 or 2 patient subsets are created.
#   Descriptive names (labels) are given to the data.frames so that it can be recognized which data.frame was derived 
#   from which data node and subset. 
#   Format of the labels: an unique identifier for the node (numerical identifier appended behind the letter "n"),
#   followed by _s1 or _s2 depending on subset (actually, subset number can be anything, as long as it is numerical), 
#   format: n<numerical_id_node>_s<numerical_id_subset>, e.g. n0_s1, n0_s2, n1_s1, n1_s2. 
#   
# or
# * variable "preprocessed" -  containing one single data.frame with the preprocessed data. Present as a variable in 
#   the global environment. It contains all the preprocessed data from the different nodes merged into a single 
#   data.frame.
# 
# The data.frames (coming from high dimensional nodes) in 'loaded_variables' and 'preprocessed' have columns: 
# 'Row.Label' and 'Bio.marker' (Bio.marker column isoptional), followed by the sample columns with the measurements 
# (for example columns named ASSAY_0001, ASSAY_0002 ...)
#
#
# Expected arguments to main() function:
# * phase parameter. Expected argument: "fetch" or "preprocess".  This parameter specifies whether the script is run for
#   the 'fetch data' or 'preprocess data' tab, and it is used to fetch the corresponding data from the appropriate  
#   variable (loaded_variables/preprocessed) and in the names of the output files of this script
#
# * projection of the data: "default_real_projection" (= intensity values/ counts) or "log_intensity" 
#   (log2 of intensity values/counts).
######
# Output:
# * 1 boxplot image per data node, png format. Name: <phase>_box_plot_node_<node identifier>.png.
# * 1 textfile per node with the summary statistics per subset, in json format. 
#     Name: <phase>_summary_stats_node_<node identifier>.json.
# for the "fetch" data output files the node identifier is extracted from the labels of loaded_variables, format:
#   n<numerical node identifier> (e.g. "n0"), for "preprocess data" the data from all nodes has been merged and
#   the word "all" is used for the node identifier in the output file names
#
# Note: If the data node is not high dimensional or the dataset is empty, no boxplot will be returned - only an image 
#       with the text "No data points to plot", also no mean, median etc. will be returned in the summary statistics: 
#       only variableLabel, node name and subset name are returned, and totalNumberOfValues = 1 and 
#       numberOfMissingValues = 1.
######
# If this script gets extended for clinical data:
#     ** right now this summary stats script is only adapted for high dimensional data nodes, later the functionality 
#         might be extended for clinical data. In that case it might be possible to recognize if it is high or low dim 
#         data based on the column names of the data.frame (assuming low dim data will also be passed on in the form of 
#         data.frames and will not contain a "Row.Label" column)
#
# NOTE:  the script is now assuming that the first column of high dim nodes is called "Row.Label", and that data coming
#         from low dim nodes NEVER have a first column named "Row.Label"".
#
###########


### SE: Currently it only generates summary statistics for high dimensional data!!!

library(jsonlite)
library(gplots)

## SE: For debug
#setwd("/Users/serge/GitHub/SmartR")

if (!exists("remoteScriptDir")) {  #  Needed for unit-tests
  remoteScriptDir <- "web-app/HeimScripts"
}

## Loading functions ##
summaryStatistics <- paste(remoteScriptDir, "/_shared_functions/SummaryStatistics/summaryStatistics.R", sep="")
source(summaryStatistics)



main <- function(phase = "", projection = "")
{
  input_data <- get_input_data(phase)
  
  
  ## Filtering out the high dimensional nodes on which this script will work
  if(phase == "fetch"){
    hd_idx.vec = grep("highDimensional", names(input_data))
    input_data = input_data[hd_idx.vec]
  }
  
  
  check_input(input_data, projection)
  
  data_measurements <- extract_measurements(input_data)
  
  
  # for preprocessed data the data from the two different subsets is merged
  # the summary stats and boxplot functions expect one data.frame per subset per node.
  # preprocessed data.frame therefore needs to be split
  if (phase == "preprocess")
  {
    data_measurements <- split_on_subsets(data_measurements)
  }
  
  summary_stats_json <- produce_summary_stats(data_measurements, phase)
  write_summary_stats(summary_stats_json)
  produce_boxplot(data_measurements, phase, projection)
  
  return(list(messages = "Finished successfully"))
}




## SE: for debug
# b = main(phase = "preprocess", projection = "log_intensity")
# a = main(phase = "fetch", projection = "log_intensity")





