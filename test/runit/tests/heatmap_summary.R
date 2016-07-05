source("web-app/HeimScripts/heatmap/summary.R")


### main test data ###
test_set <-
  data.frame(
    Row.Label = c("1007_s_at", "1053_at", "117_at", "121_at", "1255_g_at"), 
    Bio.marker = letters[1:5], 
    GSM210004 = c(179.26, NA, 0, 0.00734,-23.0234),
    GSM210005 = c(185.841, 42.2334, 16.7999, 64.5686, 3.82495), 
    GSM210006 = c(179.15, 35.5867, 33.9371, 105.912, 1.15958), 
    GSM210007 = c(367.436, 60.1003, 31.3962, 140.325, 5.08933), 
    stringsAsFactors = F
  )

test_set_measurements <-
  test_set[, c("GSM210004", "GSM210005", "GSM210006", "GSM210007")]

test_set2 <- test_set
test_set2[,3:6] <- test_set2[,3:6] + 1
test_set2_measurements <-
  test_set2[, c("GSM210004", "GSM210005", "GSM210006", "GSM210007")]

test_data <- list("highDimensional_n0_s1" = test_set)
test_data_measurements <- list("highDimensional_n0_s1" = test_set_measurements)

test_set_preprocessed  <- test_set
colnames(test_set_preprocessed) <-
  c("Row.Label" , "Bio.marker","GSM210004_n0_s1", "GSM210005_n0_s1", "GSM210006_n1_s1", "GSM210007_n1_s1")
test_set_preprocessed <-
  cbind(test_set_preprocessed, test_set2_measurements, stringsAsFactors = F)
colnames(test_set_preprocessed) <-
  c("Row.Label" , "Bio.marker","GSM210004_n0_s1", "GSM210005_n0_s1", "GSM210006_n1_s1", "GSM210007_n1_s1",
    "GSM210004_n0_s2", "GSM210005_n0_s2", "GSM210006_n1_s2", "GSM210007_n1_s2")


## Creating fetch_params test data ##
test_data_ontologyTerms_highDimensional = list(
  key = "\\\\Demo Data\\Demo Data\\Sorlie(2003) GSE4382\\Biomarker Data\\Gene Expression\\Custom Array Stanford Microarray Database\\Breast\\",
  level = 5,
  fullName = "\\Demo Data\\Sorlie(2003) GSE4382\\Biomarker Data\\Gene Expression\\Custom Array Stanford Microarray Database\\Breast\\",
  name = "Breast",
  tooltip = "\\Demo Data\\Sorlie(2003) GSE4382\\Biomarker Data\\Gene Expression\\Custom Array Stanford Microarray Database\\Breast\\",
  visualAttributes = c("LEAF", "ACTIVE", "HIGH_DIMENSIONAL"),
  metadata = NULL,
  dimensionCode = "\\Demo Data\\Sorlie(2003) GSE4382\\Biomarker Data\\Gene Expression\\Custom Array Stanford Microarray Database\\Breast\\",
  dimensionTableName = "CONCEPT_DIMENSION"
)

test_set_fetch_params = list(
  ontologyTerms = list(highDimensional_n0 = test_data_ontologyTerms_highDimensional),
  resultInstanceIds = c(31249, 31250),
  assayConstraints = NULL,
  dataConstraints = NULL,
  projection = "log_intensity"
)

.setUp <- function() {
  # A .png file is created by produce_boxplot. Should be written to temporary directory
  assign("origDirectory", getwd(), envir = .GlobalEnv)
  dir <- tempdir()
  dir.create(dir)
  setwd(dir)
  
  #loaded_variables and preprocessed are sometimes removed by a test function
  assign("loaded_variables", test_data, envir = .GlobalEnv)
  assign("preprocessed", list(preprocessed=test_set_preprocessed), envir = .GlobalEnv)
  assign("fetch_params", test_set_fetch_params, envir = .GlobalEnv)
  
}

.tearDown <- function() {
  dir <- getwd()
  setwd(origDirectory)
  unlink(dir, recursive = T, force = T)
}



### unit tests for function get_input_data ###

## get_input_data should return value of loaded_variables if phase = "fetch"
test.get_input_data.fetch <- function() {
  fetched_data <- get_input_data("fetch")
  checkEquals(loaded_variables, fetched_data)
}

## get_input_data should take variable preprocessed and return its contents as a named list if phase = "preprocess"
test.get_input_data.preprocess <- function() {
  checkEquals( preprocessed, get_input_data("preprocess"))
}

## get_input_data should throw error if phase is NA or anything other than "fetch" and "preprocess
test.get_input_data.phaseNA <- function() {
  checkException(get_input_data(NA))
}

test.get_input_data.phaseIncorrect <- function() {
  checkException(get_input_data("Incorrect_phase"))
}

## get_input_data should throw error if phase is fetch and loaded_variables does not exist
test.get_input_data.no_loaded_variables <- function() {
  rm(loaded_variables, envir = .GlobalEnv)
  checkException(get_input_data("fetch"))
}

## get_input_data should throw error if phase is preprocess and preprocessed does not exist
test.get_input_data.no_preprocessed <- function() {
  rm(preprocessed, envir = .GlobalEnv)
  checkException(get_input_data("preprocess"))
}


### unit tests for function extract_measurements ###

#test if function works if data has a "Row.Label" and a "Bio.marker" column
test.extract_measurements.simplecase1 <- function() {
  checkEquals(test_data_measurements, extract_measurements(test_data))
}

#test if it works if data has "Row.Label" column but no "Bio.marker" column
test.extract_measurements.simplecase2 <- function() {
  test_data_tmp <-
    list("highDimensional_n0_s1" = test_set[, c("Row.Label", "GSM210004", "GSM210005", "GSM210006", "GSM210007")])
  checkEquals(test_data_measurements, extract_measurements(test_data_tmp))
}

#test if it works for preprocessed data
test.extract_measurements.preprocesseddata <- function() {
  checkEquals(
    list(preprocessed = test_set_preprocessed[,-c(1:2)]), 
    extract_measurements(list(preprocessed = test_set_preprocessed))
    )
}


#should return a list with a missing value if the data has no "Row.Label" column
test.extract_measurements.noRowLabel <- function() {
  tmp <- test_data
  colnames(tmp[[1]])[1] <- "No.Row.Label"
  checkEquals(list("highDimensional_n0_s1" = NA), extract_measurements(tmp))
}


# function should return error if the dataframe has non-numeric columns in addition to "Row.Label" and "Bio.marker"
test.extract_measurements.nonNumeric <- function() {
  tmp <-
    list("highDimensional_n0_s1" = cbind(
      test_set, nonNumericCol = letters[1:5], stringsAsFactors = F
    ))
  checkException(extract_measurements(tmp))
}

#should work if data for multiple nodes is provided
test.extract_measurements.multiplenodes <- function() {
  test_data_multiple_nodes <-
    list("highDimensional_n0_s1" = test_set, "highDimensional_n1_s1" = test_set2, "highDimensional_n2_s1" = test_set)
  test_data_measurements_multiple_nodes <-
    list("highDimensional_n0_s1" = test_set_measurements,"highDimensional_n1_s1" = test_set2_measurements, "highDimensional_n2_s1" = test_set_measurements)
  checkEquals(
    test_data_measurements_multiple_nodes, extract_measurements(test_data_multiple_nodes)
  )
}

#should work if data for multiple subsets is provided
test.extract_measurements.multiplesubsets <- function() {
  test_data_multiple_subsets <-
    list("highDimensional_n0_s1" = test_set2, "highDimensional_n0_s2" = test_set, "highDimensional_n1_s1" = test_set)
  test_data_measurements_multiple_subsets <-
    list("highDimensional_n0_s1" = test_set2_measurements,"highDimensional_n0_s2" = test_set_measurements, "highDimensional_n1_s1" = test_set_measurements)
  checkEquals(
    test_data_measurements_multiple_subsets, extract_measurements(test_data_multiple_subsets)
  )
}


#what if number of samples = 1 or number of genes = 1?
test.extract_measurements.1sample <- function() {
  test_data_one_sample <-
    list("highDimensional_n0_s1" = test_set[,c("Row.Label", "Bio.marker", "GSM210004")])
  checkEquals(list("highDimensional_n0_s1" = test_set[,c("GSM210004"), drop = F]), extract_measurements(test_data_one_sample))
}

test.extract_measurements.1probe <- function() {
  test_data_one_probe <- list("highDimensional_n0_s1" = test_set[1,])
  checkEquals(
    list("highDimensional_n0_s1" = test_set[1,c("GSM210004", "GSM210005", "GSM210006", "GSM210007")]), 
    extract_measurements(test_data_one_probe)
    )
}

#1 probe, 1 sample
test.extract_measurements.1sample1probe <- function() {
  test_data_one_probe_sample <-
    list("highDimensional_n0_s1" = test_set[1,c("Row.Label", "Bio.marker", "GSM210004")])
  checkEquals(list("highDimensional_n0_s1" = test_set[1,c("GSM210004"), drop = F]), extract_measurements(test_data_one_probe_sample))
}

### unit tests for function split_on_subsets ###

#should return one data.frame if there is only one subset
test.split_on_subsets.onesubset <- function() {
  one_subset_data <-
    test_set_preprocessed[,c("GSM210004_n0_s1", "GSM210005_n0_s1", "GSM210006_n1_s1", "GSM210007_n1_s1")]
  preprocessed_measurements <- list(preprocessed = one_subset_data)
  split_preprocessed_measurements <-
    list(preprocessed_s1 = one_subset_data)
  
  checkEquals(split_preprocessed_measurements, split_on_subsets(preprocessed_measurements))
}

#should return two data.frames if there are two subsets
test.split_on_subsets.twosubsets <- function() {
  preprocessed_measurements <-
    list(preprocessed = 
           test_set_preprocessed[,c("GSM210004_n0_s1", "GSM210005_n0_s2", "GSM210006_n1_s1", "GSM210007_n1_s2")])
  split_preprocessed_measurements <-
    list(preprocessed_s1 = test_set_preprocessed[,c("GSM210004_n0_s1", "GSM210006_n1_s1")],
         preprocessed_s2 = test_set_preprocessed[,c("GSM210005_n0_s2", "GSM210007_n1_s2")])
  checkEquals(split_preprocessed_measurements, split_on_subsets(preprocessed_measurements))
}

# should also work in case sample name contains '_'
test.split_on_subsets.preprocessed_samplename_underscore <-
  function() {
    preprocessed_measurements <-
      list(preprocessed = 
             test_set_preprocessed[,c("GSM210004_n0_s1", "GSM210005_n0_s2", "GSM210006_n1_s1", "GSM210007_n1_s2")])
    colnames(preprocessed_measurements[[1]]) <-
      c("sample _ 1_n1_s1", "sample _ 2_n1_s1","sample _ 3_n1_s2","sample _ 4_n1_s2")
    split_preprocessed_measurements <-
      list(preprocessed_s1 = preprocessed_measurements[[1]][,c("sample _ 1_n1_s1", "sample _ 2_n1_s1")],
           preprocessed_s2 = preprocessed_measurements[[1]][,c("sample _ 3_n1_s2","sample _ 4_n1_s2")])
    checkEquals(split_preprocessed_measurements, split_on_subsets(preprocessed_measurements))
  }


### unit tests for function produce_summary_stats ###
# input is numeric, as output extract_measurements is numeric

phase <- "fetch"

#summary stats corresponding to test_set
summary_stats_table <- data.frame(
  "variableLabel" = "highDimensional_n0_s1",
  "node" = "highDimensional_n0",
  "subset" = "s1",
  "totalNumberOfValuesIncludingMissing" = 20,
  "numberOfMissingValues" = 1,
  "numberOfSamples" = 4,
  "min" = -23.0234,
  "max" = 367.436,
  "mean" = 75.242316,
  "standardDeviation" = 97.450718,
  "q1" = 4.45714,
  "median" = 35.5867,
  "q3" = 123.1185,
  stringsAsFactors = F
)

#summary stats corresponding to test_set2
summary_stats_table2 <- summary_stats_table
summary_stats_table2[,c("min", "max", "mean", "q1", "median", "q3")] <-
  summary_stats_table2 [,c("min", "max", "mean", "q1", "median", "q3")] +
  1


# 1 node, 1 subset
test.produce_summary_stats.simplecase <- function() {
  checkEquals(
    list("fetch_summary_stats_node_highDimensional_n0.json" = summary_stats_table), produce_summary_stats(test_data_measurements, phase)
  )
}

#multiple nodes, multiple subsets
test.produce_summary_stats.multiplenodesandsubsets <- function() {
  test_data_measurements_multiple_nodes_subsets <-
    list(
      "highDimensional_n0_s1" = test_set_measurements,"highDimensional_n0_s2" = test_set2_measurements,
      "highDimensional_n1_s1" = test_set2_measurements, "highDimensional_n1_s2" = test_set_measurements
    )
  
  expected_result_n0 <-
    rbind(summary_stats_table,summary_stats_table2)
  expected_result_n0[2 , c("variableLabel","subset")] <-
    c("highDimensional_n0_s2", "s2")
  
  expected_result_n1 <-
    rbind(summary_stats_table2,summary_stats_table)
  expected_result_n1[1 , c("variableLabel", "node","subset")] <-
    c("highDimensional_n1_s1", "highDimensional_n1", "s1")
  expected_result_n1[2 , c("variableLabel", "node","subset")] <-
    c("highDimensional_n1_s2", "highDimensional_n1", "s2")
  
  checkEquals(
    list(
      "fetch_summary_stats_node_highDimensional_n0.json" = expected_result_n0, "fetch_summary_stats_node_highDimensional_n1.json" = expected_result_n1
    ), produce_summary_stats(test_data_measurements_multiple_nodes_subsets, phase)
  )
}

# preprocessed data
test.produce_summary_stats.preprocessed <- function() {
  split_preprocessed_measurements <-
    list(preprocessed_s1 = 
           test_set_preprocessed[,c("GSM210004_n0_s1", "GSM210005_n0_s1", "GSM210006_n1_s1", "GSM210007_n1_s1")],
         preprocessed_s2 = 
           test_set_preprocessed[,c("GSM210004_n0_s2", "GSM210005_n0_s2", "GSM210006_n1_s2", "GSM210007_n1_s2")])
  expected_result <-
    rbind(summary_stats_table,summary_stats_table2)
  expected_result[1 , c("variableLabel","subset", "node")] <-
    c("preprocessed_s1", "s1", "preprocessed_allNodesMerged")
  expected_result[2 , c("variableLabel","subset", "node")] <-
    c("preprocessed_s2", "s2", "preprocessed_allNodesMerged")
  
  checkEquals(
    list("preprocess_summary_stats_node_all.json" = expected_result), 
    produce_summary_stats(split_preprocessed_measurements, "preprocess")
  )
}




# what if list item = NA
test.produce_summary_stats.itemNA <- function() {
  test_NA_set <- list("highDimensional_n0_s1" = NA)
  
  expected_result <- summary_stats_table
  expected_result[,c("totalNumberOfValuesIncludingMissing", "numberOfMissingValues")] <-
    c(1,1)
  expected_result[,c("min", "max", "mean", "standardDeviation", "q1", "median", "q3")] <-
    as.numeric(NA) #numeric
  expected_result[,"numberOfSamples"] <- NA #logical
  
  checkEquals(
    list("fetch_summary_stats_node_highDimensional_n0.json" = expected_result), produce_summary_stats(test_NA_set, phase)
  )
}

# what if number of samples = 1 , or number of genes = 1?
test.produce_summary_stats.1sample <- function() {
  test_data_one_sample <- list("highDimensional_n0_s1" = data.frame(sample = 1:5))
  expected_result <- summary_stats_table
  expected_result[,c(
    "totalNumberOfValuesIncludingMissing", "numberOfMissingValues","numberOfSamples"
  )] <- c(5, 0, 1)
  expected_result[,c("min", "max", "mean", "standardDeviation", "q1", "median", "q3")] <-
    c(1, 5, 3, 1.58113883, 2, 3, 4)
  
  checkEquals(
    list("fetch_summary_stats_node_highDimensional_n0.json" = expected_result), produce_summary_stats(test_data_one_sample, phase)
  )
}

test.produce_summary_stats.1probe <- function() {
  test_data_one_probe <-
    data.frame(
      sample1 = 1, sample2 = 2, sample3 = 3, sample4 = 4, sample5 = 5
    )
  test_data_one_probe <- list("highDimensional_n0_s1" = test_data_one_probe)
  
  expected_result <- summary_stats_table
  expected_result[,c(
    "totalNumberOfValuesIncludingMissing", "numberOfMissingValues","numberOfSamples"
  )] <- c(5, 0, 5)
  expected_result[,c("min", "max", "mean", "standardDeviation", "q1", "median", "q3")] <-
    c(1, 5, 3, 1.58113883, 2, 3, 4)
  
  checkEquals(
    list("fetch_summary_stats_node_highDimensional_n0.json" = expected_result), produce_summary_stats(test_data_one_probe, phase)
  )
}

#1 probe, 1 sample
test.produce_summary_stats.1probe1sample <- function() {
  test_data_one_probe_sample <- list("highDimensional_n0_s1" = data.frame(sample = 1))
  expected_result <- summary_stats_table
  expected_result[,c(
    "totalNumberOfValuesIncludingMissing", "numberOfMissingValues","numberOfSamples"
  )] <- c(1, 0, 1)
  expected_result[,c("min", "max", "mean", "standardDeviation", "q1", "median", "q3")] <-
    c(1, 1, 1, as.numeric(NA), 1, 1, 1)
  
  checkEquals(
    list("fetch_summary_stats_node_highDimensional_n0.json" = expected_result), produce_summary_stats(test_data_one_probe_sample, phase)
  )
}


### unit tests for function produce_boxplot ###
# input = numeric, as output extract_measurements is numeric (tested in one of the earlier unit tests)

# boxplot output corresponding to test_set
boxplot_table <- list(
  stats = matrix(c(
    -23.02340, 4.45714, 35.58670, 123.11850, 185.84100
  ),5,1),
  n = 19,
  conf = matrix(c(-7.42529712,78.59869712), 2, 1),
  out = c(GSM2100071 = 367.436),
  group = 1,
  names = "s1"
)

# boxplot output corresponding to test_set2
boxplot_table_2 <- list(
  stats = matrix(c(
    -22.02340, 5.45714, 36.58670, 124.11850, 186.84100
  ),5,1),
  n = 19,
  conf = matrix(c(-6.42529712,79.59869712), 2, 1),
  out = c(GSM2100071 = 368.436),
  group = 1,
  names = "s1"
)

# boxplot output corresponding to 1 node with two subsets, each subset containing the data from test_set
boxplot_table_2subsets <- list(
  stats = matrix(
    c(-23.02340, 4.45714, 35.58670, 123.11850, 185.84100),5,2, byrow = F
  ),
  n = c(19, 19),
  conf = matrix(c(-7.42529712,78.59869712), 2, 2, byrow = F),
  out = c(GSM2100071 = 367.436, GSM2100071 = 367.436),
  group = c(1, 2),
  names = c("s1", "s2")
)

# boxplot output corresponding to 1 node with two subsets, the first subset containing the data from test_set
#  and the second from test_set2
boxplot_table_2subsets2 <- list(
  stats = matrix(
    c(
      -23.02340, 4.45714, 35.58670, 123.11850, 185.84100,-22.02340, 5.45714, 36.58670, 124.11850, 186.84100
    ),5,2, byrow = F
  ),
  n = c(19, 19),
  conf = matrix(
    c(-7.42529712,78.59869712,-6.42529712,79.59869712), 2, 2, byrow = F
  ),
  out = c(GSM2100071 = 367.436, GSM2100071 = 368.436),
  group = c(1, 2),
  names = c("s1", "s2")
)

# 1 node, 1 subset
test.produce_boxplot.simplecase <- function() {
  checkEquals(
    list("fetch_box_plot_node_highDimensional_n0.png" = boxplot_table), produce_boxplot(test_data_measurements, phase, "log_intensity")
  )
}



#multiple nodes, multiple subsets
test.produce_boxplot.multiplenodesandsubsets <- function() {
  test_data_measurements_multiple_nodes_subsets <-
    list(
      "highDimensional_n0_s1" = test_set_measurements, "highDimensional_n0_s2" = test_set2_measurements,
      "highDimensional_n1_s1" = test_set_measurements, "highDimensional_n1_s2" = test_set_measurements
    )
  checkEquals(
    list(
      "fetch_box_plot_node_highDimensional_n0.png" = boxplot_table_2subsets2, "fetch_box_plot_node_highDimensional_n1.png" = boxplot_table_2subsets
    ), produce_boxplot(
      test_data_measurements_multiple_nodes_subsets, phase, "log_intensity"
    )
  )
}

#preprocessed data
test.produce_boxplot.preprocessed <- function() {
  split_preprocessed_measurements <-
    list(preprocessed_s1 = 
           test_set_preprocessed[, c("GSM210004_n0_s1", "GSM210005_n0_s1", "GSM210006_n1_s1", "GSM210007_n1_s1")],
         preprocessed_s2 = 
           test_set_preprocessed[, c("GSM210004_n0_s2", "GSM210005_n0_s2", "GSM210006_n1_s2", "GSM210007_n1_s2")])
  expected_result <-
    list("preprocess_box_plot_node_all.png" = boxplot_table_2subsets2)
  names(expected_result$preprocess_box_plot_node_all.png$out) <-
    c("GSM210007_n1_s11", "GSM210007_n1_s21")
  checkEquals(
    expected_result, produce_boxplot(
      split_preprocessed_measurements, "preprocess", "log_intensity"
    )
  )
}

# what if list item = NA
test.produce_boxplot.itemNA <- function() {
  test_NA_set <- list("highDimensional_n0_s1" = NA)
  expected_result <- "No data points to plot"
  checkEquals(
    list("fetch_box_plot_node_highDimensional_n0.png" = expected_result), produce_boxplot(test_NA_set, phase,"log_intensity")
  )
}

# what if number of samples = 1 , or number of genes = 1?
test.produce_boxplot.1sample <- function() {
  # make a data.frame with only one column that contains all the measurements of the test_set, but now in one column. 
  # the statistics then remain the same
  test_data_one_sample <-
    list("highDimensional_n0_s1" = data.frame(sample = unlist(test_data_measurements))) 
  expected_result <- boxplot_table
  names(expected_result$out) <- "sample16"
  checkEquals(
    list("fetch_box_plot_node_highDimensional_n0.png" = expected_result), produce_boxplot(test_data_one_sample, phase, "log_intensity")
  )
}

test.produce_boxplot.1probe <- function() {
  # make a data.frame with only one row that contains all the measurements of the test_set, but now in one row. 
  # the statistics then remain the same
  test_data_one_probe <- list("highDimensional_n0_s1" = as.data.frame(matrix(
    unlist(test_data_measurements), 1, 20, dimnames = list(1, letters[1:20])
  ))) 
  expected_result <- boxplot_table
  names(expected_result$out) <- "p"
  checkEquals(
    list("fetch_box_plot_node_highDimensional_n0.png" = expected_result), produce_boxplot(test_data_one_probe, phase, "log_intensity")
  )
}

#1 probe, 1 sample
test.produce_boxplot.1probe1sample <- function() {
  test_data_one_probe_sample <- list("highDimensional_n0_s1" = data.frame(sample = 1))
  expected_result <- list(
    stats = matrix(rep(1, 5), 5, 1),
    n = 1,
    conf = matrix(c(1,1), 2, 1),
    out = as.numeric(NULL),
    group = as.numeric(NULL),
    names = "s1"
  )
  checkEquals(
    list("fetch_box_plot_node_highDimensional_n0.png" = expected_result), 
    produce_boxplot(test_data_one_probe_sample, phase, "log_intensity")
  )
}
