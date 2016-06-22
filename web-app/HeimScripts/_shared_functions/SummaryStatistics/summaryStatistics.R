###############################################################
### Functions related to generating the summary statistics  ###
###############################################################

## Retrieve high dimensional input data according to phase
## parameter provided as input:
## - fetch -> loaded_variables
## - preprocess -> preprocessed  
get_input_data <- function(phase)
{
  #check if the phase parameter
  if (is.na(phase))
  {
    stop(
      "Supply phase parameter to function \'main()\'. Expected input: \'fetch\' or \'preprocess\'"
    )
  }
  
  if (phase != "fetch" & phase != "preprocess" & !is.na(phase))
  {
    error_msg <-
      paste(
        "Incorrect value for phase parameter. Given input: \'", phase, "\'.",
        "Expected input: \'fetch\' or \'preprocess\'", sep = ""
      )
    stop(error_msg)
  }
  
  # for fetch data tab the data from 'loaded_variables' should be used
  if (phase == "fetch")
  {
    if (!exists("loaded_variables"))
    {
      stop(
        "Summary stats is run for phase \'fetch\', but variable \'loaded_variables\' does not exist in the environment."
      )
    }
    input_data <- loaded_variables
  }
  
  # for preprocess data tab the data from 'preprocessed' should be used
  if (phase == "preprocess")
  {
    if (!exists("preprocessed"))
    {
      stop(
        "Summary stats is run for phase \'preprocess\', but variable \'preprocessed\' does not exist in the environment."
      )
    }
    input_data <- list("preprocessed" = preprocessed$HD)
  }
  
  return(input_data)
}


#check if provided variables and phase info are in line with expected input as described at top of this script
check_input <- function(datasets, projection)
{
  #expected input: list of data.frames
  items_list <- sapply(datasets, class)
  
  if (class(datasets) != "list" | !all(items_list == "data.frame"))
    #for a data.frame is.list() also returns TRUE. Class returns "data.frame" in that case
  {
    stop("Unexpected input. Expected input: a list, containing one or more data.frames")
  }
  
  # all items in the list are expected to have some unique identifier for the node (numerical identifier appended behind
  # the letter "n"), followed by an underscore and a subset identifier - s<numerical id> - depending on subset, 
  # e.g. n0_s1, n0_s2, n1_s1, n1_s2.
  dataset_names <- names(datasets)
  #expected_format_names <- "^n[[:digit:]]+_s[[:digit:]]+$"
  expected_format_names <- "^.+_n[0-9]+_s[1-2]{1}$"
  names_in_correct_format <-
    grepl(expected_format_names, dataset_names)
  if (any(!names_in_correct_format &
          dataset_names !=  "preprocessed"))
  {
    stop(
      paste(
        "One or more labels of the datasets do not have the expected format.",
        "\nExpected format: either the label should be \'preprocessed\' or it should be an unique numerical identifier",
        "for the node appended behind the letter \n\',",
        "\nfollowed by an underscore and an unique numerical identifier for the subset appended behind an \'s\',",
        "\ne.g. n0_s1, n0_s2, n1_s1, n1_s2. "
      )
    )
  }
  
  #all labels of the data.frames should be unique
  dataset_names <- names(datasets)
  if (any(duplicated(dataset_names)))
  {
    stop(
      "Not all data.frame labels are unique; one or more labels in the \'loaded_variables\' or \'preprocessed\' variable are duplicated"
    )
  }
  
  # check if data.frame "preprocessed" has the right column names
  # expected format of column names for sample columns is <sample_name>_n<numerical node id>_s<numerical subset id>, 
  # eg. sample234_n1_s1
  if (dataset_names[1] == "preprocessed")
    #for preprocessed data there is only 1 data.frame
  {
    column_names <- colnames(datasets$preprocessed)
    expected_format_names <- "^.+_n[0-9]+_s[1-2]{1}$"
    names_in_correct_format <-
      grepl(expected_format_names, column_names)
    if (!all(
      names_in_correct_format |
      column_names == "Row.Label" | column_names == "Bio.marker"
    ))
    {
      stop(
        paste(
          "The column names of the sample columns of the data.frame \'preprocessed\' do not have the expected format.",
          "\nExpected column names for the feature columns (probes/genes/proteins/etc): \'Row.Label\'", 
          "and (optional) \'Bio.marker\' ",
          "\nSample columns should adhere to the format: <sample_name>_n<numerical node id>_s<numerical subset id>,",
          "eg. sample234_n1_s1"
        )
      )
    }
  }
  
  
  if (is.na(projection))
  {
    stop(
      "Supply projection parameter to function \'main()\'. Expected input:  \'default_real_projection\' or \'log_intensity\'"
    )
  }
  
  if (!is.na(projection) &
      projection != "default_real_projection" &
      projection != "log_intensity")
  {
    stop(
      "Incorrect value for projection parameter - expected input:  \'default_real_projection\' or \'log_intensity\'"
    )
  }
  
  
  
  
  
}


# preprocessed data contains one big data.frame. Summary stats and boxplot should be made per subset, so the dataset 
# gets split on subsets
split_on_subsets <- function(preprocessed_measurements)
{
  preprocessed_measurements <- preprocessed_measurements$preprocessed
  
  split_measurements <- list()
  used_columns <- c()
  
  
  # get possible subset identifiers. Format of column names for sample columns is 
  #   <sample_name>_n<numerical node id>_s<numerical subset id>, eg. sample234_n1_s1
  column_names <- colnames(preprocessed_measurements)
  
  subsets_match_indices <- regexpr("s[1-2]{1}$", column_names)
  subsets_matches <-
    regmatches(x = column_names, m = subsets_match_indices)
  subset_names <- unique (subsets_matches)
  
  if (length(subset_names) == 0)
  {
    stop(
      "Something went wrong when splitting the preprocessed dataset into separate datasets per cohort (subset).
      Check the format of the column names of variable \'preprocessed\'.
      Expected format for column labels (of sample columns) of data.frame \'preprocessed\' is:
      <sample_name>_n<numerical identifier for the node>_s<numerical identifier for the subset>
      for example sample23341_n1_s1"
    )
    
  }
  
  #extract for each subset the corresponding data
  for (subset_name in subset_names)
  {
    subset_columns <-
      grep(paste(subset_name, "$", sep = ""), column_names) #strict matching to avoid accidental matching 
    #(e.g. if sample name contains "s1" as well)
    
    item_name <- paste("preprocessed_", subset_name, sep = "")
    
    split_measurements[[item_name]] <-
      preprocessed_measurements[, subset_columns]
    used_columns <- c(used_columns, subset_columns)
  }
  
  # each column should only be present in one of the two split datasets, 
  # and all columns should be present in one or the other split dataset
  if (any(duplicated(used_columns)) |
      !all(1:ncol(preprocessed_measurements) %in% used_columns))
  {
    stop(
      "Something went wrong when splitting the preprocessed dataset into separate datasets per cohort (subset).
      Maybe the column names of variable \'preprocessed\' are not of the correct format?
      Expected format for column labels of data.frame \'preprocessed\' is:",
      "<sample_name>_n<numerical identifier for the node>_s<numerical identifier for the subset>, e.g. sample1_n1_s1"
    )
  }
  return(split_measurements)
}



# Extract the measurement values from the data.frames.
#  The data.frames (coming from high dimensional nodes) have columns: 
#   Row.Label, Bio.marker (optional), ASSAY_0001, ASSAY_0002 ...
#  All columns except Row.Label and Bio.marker contain the measurement values.
extract_measurements <- function(datasets)
{
  for (i in 1:length(datasets))
  {
    dataset <- datasets[[i]]
    dataset_id <- names(datasets)[[i]]
    colNames <- colnames(dataset)
    
    # test if the data.frame contains data from a high dimensional data node
    #   the first column of the data from high dimensional nodes should be "Row.Label", 
    #   a second column called "Bio.marker" is optional.
    #   For low dimensional data these columns will contain different names
    is_highDim <- (colNames[1] == "Row.Label")
    
    # if it is not high dim, then it is clinical data. Right now this script is only written to handle high dimensional 
    #   data.
    if (!is_highDim)
    {
      datasets[[i]] <- NA
    }
    
    # only keep the measurement values, by removing the row label and biomarker columns
    if (is_highDim)
    {
      non_measurement_columns <-
        which(colNames %in% c("Row.Label","Bio.marker"))
      datasets[[i]] <-
        dataset[,-non_measurement_columns, drop = F]
      
      #drop empty data columns for check if measurements are numeric.
      tmp_measurements <- dataset[,-non_measurement_columns, drop = F]
      empty_columns <- apply(tmp_measurements, allNA, MARGIN = 2)
      tmp_measurements <- tmp_measurements[,!empty_columns, drop = F]
      ##TEST IF GOES WELL IF ALL EMPTY
      if (ncol(tmp_measurements) > 0 & !all(sapply(tmp_measurements, FUN = is.numeric)))
      {
        stop(
          paste(
            "Correct extraction of data columns was not possible for dataset ",dataset_id,
            ". It seems that, aside from the Row.Label and Bio.marker column,", 
            " there are one or more non numeric data columns in the data.frame.", sep = ""
          )
        )
      }
      
    }
  }
  
  return(datasets)
}

allNA <- function(v1)
{
  return(all(is.na(v1)))
}

convertNodeNames <- function(nodeNames) {
  nodes <- sub("_s[1-2]{1}", "", nodeNames)
  names <- sapply(nodes, function(el) fetch_params$ontologyTerms[[el]]$name)
  variableLabel <- sapply(1:length(names), function(i) sub(".*_", paste(names[i], "_", sep=""), nodeNames[i]))
  variableLabel
}

# Function to produce one JSON file per data node containing the summary stats per subset for that node. 
# Summary stats include: the number of values, number of missing values quartiles, min, max, mean, std deviation, median
# NOTE: this function converts the data.frame to a vector containing all datapoints from the dataframe and calculates 
#       the statistics over this entire vector. 
#       The statistics are thus not calculated per data column or row, but over the whole data.frame
# NOTE 2: quartiles, mean and other statistics can only be calculated on numeric data. 
#         If this method is extended for clinical data:
#         * Are the desired summary statistics the same for clinical data (ie. mean, sd, median, etc.)  
#           or should different stats be calculated?
#         * extend this method to recognize low dim data and split up the data.frame or apply this method to each data 
#           column separately to calculate the statistics per clinical variable separately 
#           (ie. supply the clinical data as a list containing a separate data.frame/vector for each clinical variable)
#         * build in  a test to determine if a variable is numeric or categorical and only calculate the statistics if 
#           numeric.
#         * If a variable is categorical: are missing values in that case NA or  "" (empty string?)
produce_summary_stats <- function(measurement_tables, phase)
{
  # construct data.frame to store the results from the summary statistics in
  table_columns <-
    c(
      "variableLabel", "node", "subset", "totalNumberOfValuesIncludingMissing", "numberOfMissingValues",
      "numberOfSamples", "min", "max", "mean", "standardDeviation", "q1", "median", "q3"
    )
  
  result_table <-
    as.data.frame(matrix(
      NA, nrow = length(measurement_tables), ncol = length(table_columns),
      dimnames = list(names(measurement_tables), table_columns)
    ))
  
  # add information about node and subset identifiers
  result_table$subset <- gsub(".*_","", rownames(result_table)) #take everything after _
  result_table$node <- discardSubset(rownames(result_table))

  nodes <-  result_table$node
  result_table$node[nodes == "preprocessed"] <-
    "preprocessed_allNodesMerged"

  # calculate summary stats per data.frame
  for (i in 1:length(measurement_tables))
  {
    # get the name of the data.frame, identifying the node and subset
    identifier <- names(measurement_tables)[i]
    result_table[identifier, "variableLabel"] <- ifelse(grepl("preprocessed", identifier), identifier, convertNodeNames(identifier))

    if (!all(is.na(measurement_tables[[i]])))
    {
      cols <- ncol(measurement_tables[[i]])
      result_table[identifier, "numberOfSamples"] <- ifelse(is.null(cols), 1, cols)
    }
    
    # convert data.frame to a vector containing all values of that data.frame, a vector remains a vector
    measurements <- unlist(measurement_tables[[i]])
    
    # determine total number of values and number of missing values
    is.missing <- is.na(measurements)
    result_table[identifier, "totalNumberOfValuesIncludingMissing"] <-
      length(measurements)
    result_table[identifier, "numberOfMissingValues"] <-
      length(which(is.missing))
    
    # calculate descriptive statistics, only for numerical data.
    # the 50% quantile is the median. 0 and 100% quartiles are min and max respectively
    
    result_table[identifier, "mean"] <- mean(measurements, na.rm = T)
    result_table[identifier, "standardDeviation"] <-
      sd(measurements, na.rm = T)
    
    quartiles <- quantile(measurements, na.rm = T)
    result_table[identifier, "min"] <- quartiles["0%"]
    result_table[identifier, "max"] <- quartiles["100%"]
    result_table[identifier, "q1"] <- quartiles["25%"]
    result_table[identifier, "median"] <- quartiles["50%"]
    result_table[identifier, "q3"] <- quartiles["75%"]
  }
  
  # remove rownames so that rownames are not included in the JSON output.
  rownames(result_table) <- 1:nrow(result_table)
  
  # write the summary statistics for each node to a separate file, in JSON format
  summary_stats_all_nodes <- list()
  
  unique_nodes <- unique(result_table$node)
  for (node in unique_nodes)
  {
    partial_table <-
      result_table[which(result_table$node == node), ,drop = F]
    rownames(partial_table) <-
      1:nrow(partial_table) #does not influence json result, however is needed for unit testing (matching rownumbers).
    if (node != "preprocessed_allNodesMerged") {
      fileName <-
        paste(phase,"_summary_stats_node_", node, ".json", sep = "")
    }
    if (node == "preprocessed_allNodesMerged") {
      fileName <- paste(phase,"_summary_stats_node_all.json", sep = "")
    }
    summary_stats_all_nodes[[fileName]] <- partial_table
  }

  return(summary_stats_all_nodes)
}



write_summary_stats <- function(summary_stats)
{
  for (i in 1:length(summary_stats))
  {
    summary_stats_JSON <- toJSON(summary_stats[[i]], dataframe = "rows", pretty = T, digits = I(17))
    fileName <- names(summary_stats)[[i]]
    write(summary_stats_JSON, fileName)
  }
}


# Function that outputs one box plot image per data node
produce_boxplot <- function(measurement_tables, phase, projection)
{
  #get node identifiers
  nodes <- unique( discardSubset(names(measurement_tables)) )
  
  if (projection == "default_real_projection") {
    projection <- "intensity"
  }
  if (projection == "log_intensity") {
    projection <- "log2(intensity)"
  }
  
  # convert the tables to vectors for use with the boxplot function
  # this converts a data.frame to a vector containing all values from the data.frame, a vector remains a vector
  measurement_vectors <- measurement_tables
  for (i in 1:length(measurement_vectors))
  {
    measurement_vectors[[i]] <- unlist(measurement_vectors[[i]])
  }
  
  boxplot_results_all_nodes <- list()
  
  
  #make a separate boxplot for each node and write to a PNG file.
  for (node in nodes)
  {
    # grab the data.frames corresponding to the selected node
    identifiers_single_node <-
      grep(paste("^",node, sep = ""), names(measurement_vectors), value = T)
    single_node_data <- measurement_vectors[identifiers_single_node]
    
    #remove node prefix from the names (labels) of the data.frame
    names(single_node_data) <-
      gsub(".*_","",names(single_node_data))
    
    #make sure the subsets are always ordered s1, s2, ...
    single_node_data <-
      single_node_data[order(names(single_node_data))]
    
    ## create box plot, output to PNG file
    if (phase != "preprocess") {
      fileName <-
        paste(phase, "_box_plot_node_", node, ".png", sep = "")
    }
    if (phase == "preprocess") {
      fileName <- paste(phase, "_box_plot_node_all.png", sep = "")
    }
    png(filename = fileName)
    
    # in case there is data present: create box plot
    if (!all(is.na(single_node_data)))
    {
      if (phase != "preprocess") {
        plot_title <- paste("Box plot node:", fetch_params$ontologyTerms[[node]]$name)
      }
      if (phase == "preprocess") {
        plot_title <- "Box plot node: preprocessed - all nodes merged"
      }
      
      boxplot_results_all_nodes[[fileName]] <-
        boxplot(
          single_node_data, col = "grey", show.names = T, ylab = projection,
          main = plot_title, outline = F, pch = 20, cex =
            0.2
        )
    }
    
    # if there are no data values: create image with text "No data points to plot"
    if (all(is.na(single_node_data)))
    {
      sinkplot("start")
      write("No data points\n\   to plot","")
      sinkplot("plot")
      box("outer", lwd = 2)
      boxplot_results_all_nodes[[fileName]] <-
        "No data points to plot"
    }
    
    dev.off()
  }
  return(boxplot_results_all_nodes)
}

discardSubset <- function(labels) {
  gsub("_s[1-2]{1}","",labels)
}

