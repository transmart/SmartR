##########################################################################
## This file contains more generic functionalities when working with    ##
## "df-type"-data frame containing Microarray GEX data  as well as      ##
## GEX matrix variables containing measurements and sample ids as       ##
## colnames.                                                            ##
## "df-type": a data frame with GEX matrix data as well as Limmy diff   ##
## expr statistics or other statistics used to order and filter the df  ##
## GEX data frame (variance, coefficient of variation, ... )            ##
##########################################################################

if (!exists("remoteScriptDir")) {  #  Needed for unit-tests
  remoteScriptDir <- "web-app/HeimScripts"
}

## Loading functions ##
utils <- paste(remoteScriptDir, "/_shared_functions/Generic/utils.R", sep="")
source(utils)



## Ranking/ordering the data frame according to
## selected statistics and retaining
## the accordingly highest ranked probes
applyRanking <- function (df, ranking, max_rows) {
  nrows = min(max_rows, nrow(df))

  if (ranking %in% c("ttest", "logfold")) {
    df["SIGNIFICANCE_ABS"] <- abs(df["SIGNIFICANCE"])
    df <- df[with(df, order(-SIGNIFICANCE_ABS)), ]
    df["SIGNIFICANCE_ABS"] <- NULL
    df <- df[1:nrows, ]
    df <- df[with(df, order(SIGNIFICANCE, decreasing=TRUE)), ]
  } else if(ranking %in% c("pval", "adjpval")) {
    df <- df[with(df, order(SIGNIFICANCE)), ]
    df <- df[1:nrows, ]
  } else {
    df <- df[with(df, order(-SIGNIFICANCE)), ]
    df <- df[1:nrows, ]
  }

  df
}


## Sorting of measurement columns
## in input data frame
applySorting <- function(df,sorting) {
  
  measurements <- getMeasurements(df)
  colNames <- names(measurements)
  subsets <- getSubset(colNames)
  nodes <- getNode(colNames)
  subjects <- getSubject(colNames)

  
  timelineValues <- getTimelineValues(nodes, fetch_params$ontologyTerms)

  
   ## Converting subjects vector to integer if required for ordering
   subjects_are_integer.logical = all(grepl("^\\d+$", subjects, perl = TRUE))
 
   if(subjects_are_integer.logical)
     subjects = as.integer(subjects)

   
  ## Performing the sorting of columns
    ## - according to nodes
  if (sorting == "nodes") {
    inds <- order(subsets, timelineValues, nodes, subjects)
  
    ## - or subjects
  } else if (sorting == "subjects") {
    inds <- order(subsets, subjects, timelineValues, nodes)
    ## - else just stop the script as sth is wrong in the input params for this function
  } else{
    stop(paste("applySorting: Sorting can only be performed according to nodes/subjects. Please check your input:", sorting))
  }

  measurements <- measurements[, inds]

  ## Returning the data frame with reordered columns
  cbind(df[, c(1,2)], measurements)
}



## Generates a filtered data frame containing the GEX matrix data
## as well as the statistics used to order and filter the df GEX data frame.
## Ordering and filtering can be performed based on Limma-based
## statistics (B value, P-value, Adjusted P-value) but also other statistical
## functions like variance, mean, median
addStats <- function(df, sorting, ranking, max_rows) {
  
  # In order to prevent displaying the table from previous run.
  cleanUpLimmaOutput(markerTableJson = markerTableJson)
  
  
  measurements  <- getMeasurements(df)

  
  rankingMethod <- getRankingMethod(ranking)
  twoSubsets    <- hasTwoSubsets(measurements)

  logfold.values <- data.frame(LOGFOLD=numeric())
  ttest.values <- data.frame(TTEST=numeric())
  pval.values <- data.frame(PVAL=numeric())
  adjpval.values <- data.frame(ADJPVAL=numeric())
  bval.values <- data.frame(BVAL=numeric())
  rankingScore <- data.frame(SIGNIFICANCE=numeric())
  means = data.frame(MEAN=numeric())
  sdses = data.frame(SD=numeric())
  
  
  # rankingMethod is either a function or character string "limma"
  useLimma <- !is.function(rankingMethod) 
  
  # Is the GEX matrix valid for differential expression analysis this means containing enough not NA values?
  validLimmaMeasurements <- isValidLimmaMeasurements(measurements)


  #this is the case for more than 1 sample available in the data frame
  if (ncol(df) > 3) {


    df <- applySorting(df,sorting)  # no need for sorting in one sample case, we do it here only

    

    validLimmaMeasurements <- isValidLimmaMeasurements(measurements)


    if (twoSubsets && validLimmaMeasurements) {
      markerTable  <- getDEgenes(df)
      
      logfold.values <- markerTable["logfold"]
      ttest.values <- markerTable["ttest"]
      pval.values <- markerTable["pval"]
      adjpval.values <- markerTable["adjpval"]
      bval.values <- markerTable["bval"]
    }                                 
    
    

    ## Cannot use limma diff expr analysis for a single subset only.
    if (useLimma  && !twoSubsets) { 
      stop( paste("Illegal ranking method: ", ranking, " two subsets needed.") )
    }

    ## rankingScore provided based on Limma stat output
    if (useLimma && validLimmaMeasurements)  {

      if (!ranking %in% colnames(markerTable) )
        stop(paste("Illegal ranking method selected: ", ranking) )
      
      # Obtain the rankingScore based on selected diff expr statistics
      rankingScore <- markerTable[ranking]
      
    } 
    
 
    ## In case that the user does not want limma
    if(!useLimma){
        rankingScore <- apply(measurements, 1, rankingMethod, na.rm = TRUE )  # Calculating ranking per probe (per row)
    }
    
  
    
    ## Copy of markerTable that will be used for file dump
    if(exists("markerTable")){
      markerTable_forFileDump = markerTable
      markerTable_forFileDump["SIGNIFICANCE"] <- rankingScore
      markerTable_forFileDump                 <- applyRanking(markerTable_forFileDump, ranking, max_rows)
      markerTable_forFileDump["SIGNIFICANCE"] <- NULL
      writeMarkerTable(markerTable_forFileDump, markerTableJson = markerTableJson)
    }
    

    # this is just an auxiliary column - it will not be used for JSON.
    sdses <- apply(measurements,1 ,sd , na.rm = T)  
  
  } else{
    
    ## In case there is only one sample some ranking methods won't work
    if(ranking %in% c("mean", "median")){
      rankingScore <- apply(measurements, 1, rankingMethod, na.rm = TRUE )  # Calculating ranking per probe (per row)
    } else{
      stop(paste("This dataset contains only one sample. Illegal ranking method selected: ", ranking))
    }
 
    
  }
  
  # this is just an auxiliary column - it will not be used for JSON.
  means <- rowMeans(measurements, na.rm = T)
  
  
  df["MEAN"]         <- means
  df["SD"]           <- sdses
  df["SIGNIFICANCE"] <- rankingScore
  df["LOGFOLD"] <- logfold.values
  df["TTEST"] <- ttest.values
  df["PVAL"] <- pval.values
  df["ADJPVAL"] <- adjpval.values
  df["BVAL"] <- bval.values

  if (useLimma & !validLimmaMeasurements) {
    # don't apply ranking but throw an error message
    stop(paste("The GEX matrix does not contain enough valid measurements to",
                "perform differential expression analysis.",
                "Please check the high dimensional input data"))
  } else {
    # else apply ranking
    df <- applyRanking(df, ranking, max_rows)
  }
  
  return(df)
}


## This function generates a data frame containing the complete set
## of  the following statistics:
## -- coefficient of variation
## -- variance
## -- mean
## -- median
## -- log Fold change
## -- P-value
## -- Adjusted P-value
## -- B value
## Statistics are only computed if data allows it. Else the corresponding values
## for this  column in the results data frame remains set to NA.
getAllStatForExtDataFrame = function(df){
  
  ## We expect in this case that Limma has been performed and the
  ## data frame contains valid statistical data for Limma
  hasLimma =  ifelse(isValidLimmaMeasurements(getMeasurements(df)), TRUE, FALSE)
  

  ## Checking that the stat colnames exist for Limma analysis
  colnamesStat = c("LOGFOLD", "TTEST", "PVAL", "ADJPVAL", "BVAL")
  containsAllColnamesStat = ifelse(length(which(colnames(df) %in% colnamesStat))==length(colnamesStat), TRUE, FALSE)
  
  
  ## SE: HERE WE CAN CHECK IF limma-specific columns contain numeric data!!!!!
  
  
  ## Checking if coefficient of variation and variance can be calculated
  performCoef = isValidCoefMeasurements(getMeasurements(df))
  performVariance = isValidVarianceMeasurements(getMeasurements(df))
  performRange = isValidRangeMeasurements(getMeasurements(df))
  
  
  ## Statistics that are provided in the output ranking data frame
  ## coef, variance, range, means, median, logfold, ttest, pval, adjpval and bval
  ## so altogether 11 columns including uid
  
  Stat.df = data.frame(ROWNAME = df$ROWNAME, COEF = NA, VARIANCE = NA,
                              RANGE = NA, MEAN = NA,
                              MEDIAN = NA, TTEST = NA,
                              LOGFOLD = NA, PVAL = NA,
                              ADJPVAL = NA, BVAL = NA)
  
  
  
  ### STAT INFO ... ###
  
  ## - for limma
  if(hasLimma && containsAllColnamesStat){
    
#     ## Vectors storing the ordered item position
#     ## in input vector in drecreasing order ...
#     idx_logfold_rank.vec = order(df$LOGFOLD, decreasing = TRUE)
#     idx_ttest_rank.vec = order(df$TTEST, decreasing = TRUE)
#     idx_pval_rank.vec = order(df$PVAL, decreasing = TRUE)
#     idx_adjpval_rank.vec = order(df$ADJPVAL, decreasing = TRUE)
#     idx_bval_rank.vec = order(df$BVAL, decreasing = TRUE)
    
#     ## ... and applying this data to obtain the corresponding
#     ## ranks according to UID
#     for(i in 1:length(idx_logfold_rank.vec)){
#       rankingStat.df$LOGFOLD[idx_logfold_rank.vec[i]] = i
#       rankingStat.df$TTEST[idx_ttest_rank.vec[i]] = i
#       
#       rankingStat.df$PVAL[idx_pval_rank.vec[i]] = i
#       rankingStat.df$ADJPVAL[idx_adjpval_rank.vec[i]] = i
#       rankingStat.df$BVAL[idx_bval_rank.vec[i]] = i
#     }
    
    
    Stat.df$LOGFOLD = df$LOGFOLD
    Stat.df$TTEST = df$TTEST
    Stat.df$PVAL = df$PVAL
    Stat.df$ADJPVAL = df$ADJPVAL
    Stat.df$BVAL = df$BVAL
  }
    
  ## - for other statistics (coefficient of variance, variance,
  ##   range, mean, median)
    
  stat_tests = logical( length = 5)
  stat_tests[stat_tests==FALSE] = TRUE
    
  names(stat_tests) = c("coeffVar", "var", "normRange", "mean", "median")
    
    
  stat_tests2col = c("COEF", "VARIANCE", "RANGE", "MEAN", "MEDIAN")
  names(stat_tests2col) = names(stat_tests)
    
  ## Checking if Coefficient of variation and Variance can
  ## be computed
  if(!performCoef)
    stat_tests["coeffVar"] = FALSE
    
  if(!performVariance)
    stat_tests["var"] = FALSE
    

  if(!performRange)
    stat_tests["normRange"] = FALSE
  
      
  ## Getting the measurements data frame
  ## to perform ranking on
  measurements.df = getMeasurements(df)
    
  ## Performing all statistics tests and store the values in the output data frame
  ## (besides Limma specific ones)
  for(i in 1:length(stat_tests)){
    if(stat_tests[i]){
      stat_test = names(stat_tests)[i]
        
      ## Calculating values for the different stats
      stat_test_value.vec = apply(measurements.df, 1, stat_test)
      
      ## Storing the different stats as col in results data frame
      Stat.df[[stat_tests2col[ stat_test ]]] = stat_test_value.vec
    }
  }
  
  return(Stat.df)
}




## Returns data frame containing
## the intensity/expression measures as matrix
## and the header which corresponds to the sample
## names. In case the extended data frame generated
## by function addStats gets provided the corresponding
## stat columns get removed.
getMeasurements <- function(df) {
  
  
  ## This is the first col containing measurements
  idx_first_data_col = 3
  
  ## Removing the statistics columns in case the
  ## extended df containing stat info has been provided as input
  idx = -which(colnames(df) %in% c("MEAN", "SD",
                                   "SIGNIFICANCE", "LOGFOLD",
                                   "TTEST", "PVAL", "ADJPVAL", "BVAL"))
  
  
  ## In case it is an extended data frame containing stat info
  if(length(idx)>0){
    ## Removing the stat columns
    df = df[,idx]
    ## first col containing measurements is getting
    ## adapted to this data frame type
    idx_first_data_col = 2
  }
     
  ## Data columns get returned
  if (ncol(df) > idx_first_data_col){
    return(df[, idx_first_data_col:ncol(df)])
  } else {
    return( df[idx_first_data_col] )
  }
}



## This function removes rows with duplicate probe id (). Only the data row with first
## occurence of the probe id is kept. Probe id and gene symbol get merged to uid variable.
mergeDuplicates <- function(df) {

  ## Which probes have multiple entries in df GEX matrix?
  ## This situation can happen when non unique probe ids get
  ## used like for example gene symbols
  dupl.where <- duplicated(df$Row.Label)
  dupl.rows <- df[dupl.where, ]
  df <- df[! dupl.where, ]

  uids <- paste(df$Row.Label, df$Bio.marker, sep="--")
  
  uids[df$Row.Label == dupl.rows$Row.Label] <- paste(uids[df$Row.Label == dupl.rows$Row.Label],
                                                     dupl.rows$Bio.marker[df$Row.Label == dupl.rows$Row.Label],
                                                     sep="--")
  
  df <- cbind(ROWNAME=uids, df[, -c(1,2)])

  df
}


# nodeID has usually this format: 'X123_highDimensional_n0_s1)
# this method pretifies it with the actual node label like this: '123_BreastCancer'
idToNodeLabel <- function(ids, ontologyTerms) {
  # extract patientID (123)
  patientIDs <- sub("_.+_n[0-9]+_s[1-2]{1}$", "", ids, perl=TRUE) # remove the _highDimensional_n0_s1
  patientIDs <- sub("^X", "", patientIDs, perl=TRUE) # remove the X
  # extract subset (s1)
  subsets <- substring(ids, first=nchar(ids)-1, last=nchar(ids))
  # extract node label (Breast)
  nodes <- sub("^.+?_", "", ids, perl=TRUE) # remove the X123_
  nodes <- as.vector(substring(nodes, first=1, last=nchar(nodes)-3)) # remove the _s1
  nodeLabels <- as.vector(sapply(nodes, function(node) return(ontologyTerms[[node]]$name)))
  # put everything together (123, Breast, s1)
  paste(patientIDs, nodeLabels, subsets, sep="_")
}


## Checking if a variable called preprocessed exists in R
## workspace, else loaded_variables is used to create data frame df.
## Column names in data frame get modified by replacing matrix id
## (e.g.: n0, n1, ...) by corresponding name in fetch_params list var
parseInput<- function() {
  
  ## Retrieving the input data frame
  if (exists("preprocessed")) {
    
    ## Retrieving low and high dim data into separate vars
    df = preprocessed$HD
    ld = preprocessed$LD
    
  } else {
    
    ld_var.idx = grep("^(categoric)|(numeric)", names(loaded_variables), perl = TRUE)
    hd_var.idx = grep("^highDimensional", names(loaded_variables), perl = TRUE)
    
    ## Merging high dim data from different nodes
    df <- mergeFetchedData(loaded_variables[hd_var.idx])
    
    ## Either there is low dim data available ...
    if(length(ld_var.idx)>0){
      ld = loaded_variables[ld_var.idx]
    ## ... or not
    } else{
      ld = NULL
    }
  }
  
  ## Renaming the column names in the data frame:
  ## - removing "X" as prefix
  ## - replacing node id by node name
  ## (e.g. 'X144_n0_s1' -> '144_Breast_s2')
  colnames(df)[c(-1,-2)] = idToNodeLabel(colnames(df)[c(-1,-2)], fetch_params$ontologyTerms)
  
  return(list(HD = df, LD = ld))
}


# ## Checking if a variable called preprocessed exists in R
# ## workspace, else loaded_variables is used to create data frame df.
# ## Column names in data frame get modified by replacing matrix id
# ## (e.g.: n0, n1, ...) by corresponding name in fetch_params list var
# parseInput <- function() {
#   
#   ## Retrieving the input data frame
#   if (exists("preprocessed")) {
#     df <- preprocessed
#   } else {
#     df <- mergeFetchedData(loaded_variables)
#   }
#   
#   ## Renaming the column names in the data frame:
#   ## - removing "X" as prefix
#   ## - replacing node id by node name
#   ## (e.g. 'X144_n0_s1' -> '144_Breast_s2')
#   colnames(df)[c(-1,-2)] = idToNodeLabel(colnames(df)[c(-1,-2)], fetch_params$ontologyTerms)
#   
#   return(df)
# }



mergeFetchedData <- function(listOfHdd){
  df <- listOfHdd[[1]]

  #test if the different data.frames all contain the exact same set of probe IDs/metabolites/etc, independent of order.
  row.Labels<- df$Row.Label

  for(i in 1:length(listOfHdd)){
    if(!all(listOfHdd[[i]]$Row.Label %in% row.Labels) | !all(row.Labels %in% listOfHdd[[i]]$Row.Label) ){
      assign("errors", "Mismatched probe_ids - different platform used?", envir = .GlobalEnv)
    }
  }

  #merge data.frames
  expected.rowlen <- nrow(df)
  labels <- names(listOfHdd)
  df <- add.subset.label(df,labels[1])

  if(length(listOfHdd) > 1){
    for(i in 2:length(listOfHdd)){
      df2 <- listOfHdd[[i]]
      label <- labels[i]
      df2 <- add.subset.label(df2,label)
      df <- merge(df, df2 ,by = c("Row.Label","Bio.marker"), all = T)
      if(nrow(df) != expected.rowlen){
        assign("errors", "Mismatched probe_ids - different platform used?", envir = .GlobalEnv)
      }
    }
  }
  return(df)
}



add.subset.label <- function(df,label) {
  sample.names <- c("")
  if (ncol(df) == 3) {
    sample.names <-
      colnames(df)[3] # R returns NA instead of column name
    # for colnames(df[,3:ncol(df)])
  }else{
    measurements <- getMeasurements(df)
    sample.names <- colnames(measurements)
  }
  for (sample.name in sample.names) {
    new.name <- paste(sample.name,label,sep = "_")
    colnames(df)[colnames(df) == sample.name] <- new.name
  }
  return(df)
}



## Removing existing statistic data from
## from df object
cleanUp <- function(df) {
  df["MEAN"] <- NULL
  df["SD"] <- NULL
  df["SIGNIFICANCE"] <- NULL
  df["LOGFOLD"] <- NULL
  df["TTEST"] <- NULL
  df["PVAL"] <- NULL
  df["ADJPVAL"] <- NULL
  df["BVAL"] <- NULL
  df
}



## Ranking/ordering the data frame according to
## selected statistics and retaining
## the accordingly highest ranked probes
applyRanking <- function (df, ranking, max_rows) {
  nrows = min(max_rows, nrow(df))

  if (ranking %in% c("ttest", "logfold")) {
    df["SIGNIFICANCE_ABS"] <- abs(df["SIGNIFICANCE"])
    df <- df[with(df, order(-SIGNIFICANCE_ABS)), ]
    df["SIGNIFICANCE_ABS"] <- NULL
    df <- df[1:nrows, ]
    df <- df[with(df, order(SIGNIFICANCE, decreasing=TRUE)), ]
  } else if(ranking %in% c("pval", "adjpval")) {
    df <- df[with(df, order(SIGNIFICANCE)), ]
    df <- df[1:nrows, ]
  } else {
    df <- df[with(df, order(-SIGNIFICANCE)), ]
    df <- df[1:nrows, ]
  }

  df
}


## Creating an unpivoted-type data frame based on high dim analysis data frame
## providing essentially as output the value, and zscore for given colname and rowname
buildFields <- function(df) {
  

  
  df <- melt(df, na.rm=T, id=c("ROWNAME", "MEAN", "SD", "SIGNIFICANCE", "LOGFOLD", "TTEST", "PVAL", "ADJPVAL", "BVAL"))
  # melt implicitly casts
  # characters to factors to make your life more exciting,
  # in order to encourage more adventures it does not
  # have characters.as.factors=F param.
  df$variable <- as.character(df$variable)
  
  ## Caclulating Z-score for sample data
  ## - If there is more then one sample
  num_samples = length(unique(df$variable))
  
  if(num_samples>1){
    ZSCORE      <- (df$value - df$MEAN) / df$SD
  } else if(num_samples==1) {
    ## - If there is only one sample we use the sample data itself to calculate the z-score
    ## for the values
    ZSCORE      <- (df$value - mean(df$value)) / sd(df$value)
  } else {
    
    stop(paste("The GEX matrix does not contain enough samples to",
               "perform z-score calculation.",
               "Please check the high dimensional input data"))
    
    
  }
  
  
  idx_exclude.vec = which(colnames(df) %in% c("MEAN", "SD", "SIGNIFICANCE", "LOGFOLD", "TTEST", "PVAL", "ADJPVAL", "BVAL"))
  df = df[, -idx_exclude.vec]

  names(df)   <- c("ROWNAME", "COLNAME", "VALUE")
  df["PATIENTID"] <- getSubject(df$COLNAME)
  df["ZSCORE"] <- ZSCORE
  df["SUBSET"] <- getSubset(df$COLNAME)

  return(df[, c("PATIENTID", "ROWNAME", "COLNAME", "VALUE", "ZSCORE", "SUBSET")])
}



## Creates a data frame containing the restructured high dim data.
## providing the cohort/subset information for each colname and rowname as value.
## Rowname is set here to "Cohort" instead of probeid to help trace down
## what sample id/colname corresponds to which cohort
buildExtraFieldsHighDim <- function(df) {
  
  ROWNAME <- rep("Cohort", nrow(df))
  COLNAME <- as.character(df$COLNAME)
  PATIENTID = as.integer(getSubject(COLNAME))
  SUBSET = as.integer(getSubset(COLNAME))
  TYPE <- rep("subset",nrow(df))
  VALUE <- df$SUBSET
  ZSCORE <- rep(NA, length(VALUE))
  
  extraFields <- unique(data.frame(PATIENTID, COLNAME, ROWNAME, VALUE, ZSCORE, TYPE, SUBSET, stringsAsFactors = FALSE))

}




## Creates a data frame containing the restructured low dim data. This relies essentially an unpivot operation 
## for the input but using as input a list of data tables/data frames.
## The output column names and content are as follows:
## - FEATURE: The low dim node name (e.g. "\\Demo Data\\Sorlie(2003) GSE4382\\Subjects\\Demographics\\Age\\")
## - PATIENTID: The patient id (e.g. "112")
## - TYPE: Value type ("numeric", "categoric")
## - VALUE: The value  (e.g.: 55, "Cause of Death, Other", NA )
##
## If a value is not existing, the corresponding cell in the returned data frame
## is set to NA.
buildExtraFieldsLowDim <- function(ld.list) {
  
  if(is.null(ld.list)){
    return(NULL)
  }
  
  varNames_with_subset.vec = unlist(names(ld.list))

  varNames_without_subset.vec = paste(strsplit2(varNames_with_subset.vec, "_")[,1],
                                      strsplit2(varNames_with_subset.vec, "_")[,2], sep="_")
  
  
  type.vec = strsplit2(varNames_with_subset.vec, "_")[,1]
  
  subset.vec = strsplit2(varNames_with_subset.vec, "_")[,3]

    
  fullnames = character(length=length(varNames_without_subset.vec))
  
  ## Declaring the variables to store the extra field data
  ROWNAME.vec = character(length = 0)
  PATIENTID.vec = character(length = 0)  
  VALUE.vec = character(length = 0)
  COLNAME.vec = character(length = 0)
  TYPE.vec = character(length = 0)
  SUBSET.vec = character(length = 0)
  ZSCORE.vec = character(length = 0)


  ## Iterating over full low dim variable names 
  for(i in 1:length(varNames_without_subset.vec)){
    
    ROWNAME = get("fullName", get(varNames_without_subset.vec[i], fetch_params$ontologyTerms))
    NAME = get("name", get(varNames_without_subset.vec[i], fetch_params$ontologyTerms))
    
    ## Accessing data frame for given full low dim variable name
    ## This data frame provides for each subject the corresponding
    ## value for given variable
    ld_var.df = ld.list[[varNames_with_subset.vec[i]]]

    for(j in 1:dim(ld_var.df)[1]){
      if (ld_var.df[j, 2] == "" || is.na(ld_var.df[j, 2])) next

      ROWNAME.vec = c(ROWNAME.vec, ROWNAME)
      PATIENTID.vec = c(PATIENTID.vec, ld_var.df[j,1])
      TYPE.vec = c(TYPE.vec, type.vec[i])
      SUBSET.vec = c(SUBSET.vec, subset.vec[i])
      VALUE.vec = c(VALUE.vec, ld_var.df[j,2])
      COLNAME.vec = c(COLNAME.vec, paste(ld_var.df[j,1], NAME , subset.vec[i], sep="_"))
    }
  }

  res.df = data.frame( PATIENTID = as.integer(PATIENTID.vec),
                       COLNAME = COLNAME.vec,
                       ROWNAME = ROWNAME.vec,
                       VALUE = VALUE.vec,
                       ZSCORE = rep(NA, length(PATIENTID.vec)),
                       TYPE = TYPE.vec,
                       SUBSET = as.integer(gsub("^s", "", SUBSET.vec)), stringsAsFactors=FALSE)

  # z-score computation must be executed on both cohorts, hence it happens after all the data are in res.df
  rownames <- unique(res.df$ROWNAME)
  for (rowname in rownames) {
      sub.res.df <- res.df[res.df$ROWNAME == rowname, ]
      if (sub.res.df[1,]$TYPE == "numeric") {
          values <- as.numeric(sub.res.df$VALUE)
          ZSCORE.values <- (values - mean(values)) / sd(values)
          res.df[res.df$ROWNAME == rowname, ]$ZSCORE <- ZSCORE.values
      }
  }

  return(res.df)
}




## Filtering out the rows in df
## where no Gene Symbol/Biomarker is
## defined
dropEmptyGene <- function(d){
  d[!(d$Bio.marker == ""|
        is.null(d$Bio.marker) |
        is.na(d$Bio.marker) |
        is.nan(d$Bio.marker)),]
}



getSubset1Length <- function(measurements) {
  sum(grepl(pattern = SUBSET1REGEX, x = colnames(measurements))) # returns number of column names satisfying regexp.
}



hasTwoSubsets <- function(measurements) {
  getSubset1Length(measurements) < ncol(measurements)
}


## Checking if the GEX matrix contains
## at least two samples
hasMinTwoSamples <- function(measurements) {
  if(dim(measurements)[2]>1){
    return(TRUE)
  } else{
    return(FALSE)
  }
}


## Wrapper function to check if GEX matrix is valid
## to calculate coefficient of variation
isValidCoefMeasurements <- function(measurements){
  return(hasMinTwoSamples(measurements))
}

## Wrapper function to check if GEX matrix is valid
## to calculate variance
isValidVarianceMeasurements <- function(measurements){
  return(hasMinTwoSamples(measurements))
}

## Wrapper function to check if GEX matrix is valid
## to calculate range
isValidRangeMeasurements <- function(measurements){
  return(hasMinTwoSamples(measurements))
}


getSubset <- function(patientIDs) {
  splittedIds <- strsplit(patientIDs,"_s") # During merge,
  # which is always run we append subset id, either
  # _s1 or _s2 to PATIENTID.
  subsets <- sapply(splittedIds, FUN = tail_elem) # In proper patienid

}


## Returns the nth last element from a vector
## Integers are casted from string to integer if necessary
## This function is helpful after performing strsplit 
tail_elem <- function(vect, n = 1) {
  
 value = vect[length(vect) - n + 1]
 
 if(length(grep("^\\d+$", value, value = FALSE, perl = TRUE))>0){
  as.integer(value)
 } else{
   return(value)
 }
 
 
 
}


# to check if a subset contains at least one non missing value
subsetHasNonNA <- function (subset, row) {
  # select the  measurements of a subset by matching the subset label with column names
  # each measurements column has subset information as suffix
  # eg:??X1000314002_n0_s2, X1000314002_n0_s1
  subsetMeasurement <- row[grep(paste(c(subset, '$'), collapse=''), names(row))]
  # check if there's non missing values
  sum(!is.na(subsetMeasurement)) > 0
}



# to check if the row contains 4 non missing values and if each subset contains at least one
# non missing value
validMeasurementsRow <- function (row) {
  
  ## SE: Minimum 3 valid measures are needed for Limma
  ##sum(!is.na(row)) > 3 & subsetHasNonNA('s1', row) &  subsetHasNonNA('s2', row) 
  sum(!is.na(row)) > 2 & subsetHasNonNA('s1', row) &  subsetHasNonNA('s2', row)
  
}



## Convert GEX matrix to Z-scores
toZscores <- function(measurements) {
  
  ## Row-wise Z-score if more then one sample 
  if(ncol(measurements)>1){
    measurements = scale(t(measurements))
    measurements = t(measurements)
    
  } else{
    ## Column-wise z-score when only one sample
    measurements <- scale(measurements)
  }
  
  return(measurements)
}



# Coefficient of variation
coeffVar     <- function(x, na.rm = TRUE) {
  c_sd <- sd(x, na.rm = na.rm)
  c_mean <- mean(x, na.rm = na.rm)
  if (c_mean == 0) {c_mean <- 0.0001} # override mean with 0.0001 when it's zero
  c_sd/c_mean
}



# Specific implementation of range.
normRange <- function(x, na.rm = TRUE) {
  # NAs are already removed but we keep this parameter so that apply does not have to be called differently
  # for normRange
  x <- removeOutliers(x)
  max(x) - min(x)
}



# We define outliers as measurements falling outside .25 or .75
# quantiles by more than 1.5 of interquantile range.
removeOutliers <- function(x) {
  qnt    <- quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H      <- 1.5 * IQR(x, na.rm = TRUE)
  result <- x[!(x < (qnt[1] - H))]  # Below .25 by more than 1.5 IQR
  result <- result[!(result > (qnt[2] + H))]  # Above .75 by more than 1.5 IQR
  na.omit(result)
}



## Returns the statistical ranking function based
## on the provided input string except for Diff expr analysis
## related input strings where string "Limma" is simply returned.
getRankingMethod <- function(rankingMethodName) {
  if (rankingMethodName == "variance") {
    return(var)
  } else if (rankingMethodName == "coef") {
    return(coeffVar)
  } else if (rankingMethodName == "range") {
    return(normRange)
  } else if (rankingMethodName == "mean") {
    return(mean)
  } else if (rankingMethodName == "median") {
    return(median)
  } else {
    return("Limma")
  }
}


writeDataForZip <- function(df, zScores, pidCols) {
  df      <- df[ , -which(names(df) %in% pidCols)]  # Drop patient columns
  df      <- cbind(df,zScores)                      # Replace with zScores
  df      <- df[ , -which(colnames(df) == "SIGNIFICANCE")]
  write.table(
    df,
    "heatmap_data.tsv",
    sep = "\t",
    na = "",
    row.names = FALSE,
    col.names = TRUE
  )
}



## Probe signal aggregation based on maxMean for initial data frame containing probe id, biomarker
## and sample measurements. Probes are merged according to maxMean, this means the row with highest mean
## intensity for same biomarker will be retained.
aggregate.probes <- function(df) {
  
  if(nrow(df)<2){
    stop("Cannot aggregate probes: there only is data for a single probe (ie. only one row of data) or 
         there is insufficient bio.marker information for the selected probes to be able to match the probes to 
         biomarkers for aggregation (e.g. in case of micro-array data to match probe ID to gene symbol). 
         Suggestion: skip probe aggregation.")
  }
  
  measurements <- df[,3:ncol(df)]
  
  
  row.names(measurements) <- df[,1]
  collapsed <- collapseRows(measurements, df[,2], df[,1], "MaxMean",
                            connectivityBasedCollapsing = FALSE, #in Rmodules = TRUE. In our spec, not required
                            methodFunction = NULL, # It only needs to be specified if method="function"
                            #connectivityPower = 1, # ignored when connectivityBasedCollapsing = FALSE
                            selectFewestMissing = FALSE)
  
  #collapsed returns 0 if collapseRows does not work, otherwise it returns a list
  if(is.numeric(collapsed))
  {
    stop("Probe aggregation is not possible: there are too many missing data points in the data for succesfull probe
         aggregation. Skip probe aggregation or select more genes/proteins/metabolites/... 
         
         Note: for aggregation only the data from probes that have accompanying bio.marker information is used
         (e.g. in case of micro-array data only the probes are used that have gene symbol information coupled to it).
         This could mean the dataset used for aggregation is smaller than the initially selected dataset.")
  }
  collapsedMeasurements <- collapsed$datETcollapsed
  Bio.marker <- collapsed$group2row[,1] # first column of this matrix always contains gene
  Row.Label <- collapsed$group2row[,2]  # second column of this matrix always contains probe_id
  lastColIndex <- ncol(df)
  lastbutOne <- lastColIndex -1
  df <- data.frame(collapsedMeasurements)
  df["Bio.marker"] <- Bio.marker
  df["Row.Label"] <- Row.Label
  row.names(df) <- NULL # WGCNA adds row.names. We do not need them to be set
  return(df[,c(lastColIndex, lastbutOne , 1:(lastbutOne-1))])
  
  }

