#!/usr/bin/env Rscript

library(reshape2)
library(limma)
library(jsonlite)

## SE: Just to get things working for dev purposes
# rm(list = ls())
# load("/Users/serge/Documents/Projects/SmartR/Development_env_Input_workspace/R_workspace_objects/Heatmap/data.Rda")
# load("/Users/serge/Documents/Projects/SmartR/Development_env_Input_workspace/R_workspace_objects/Heatmap/fetchParams.Rda")
# load("/Users/serge/Documents/Projects/SmartR/Development_env_Input_workspace/R_workspace_objects/Heatmap/loaded_variables_withLDD.Rda")
# load("/Users/serge/Documents/Projects/SmartR/Development_env_Input_workspace/R_workspace_objects/Heatmap/fetch_params_withLDD.Rda")
# setwd("/Users/serge/GitHub/SmartR")

if (!exists("remoteScriptDir")) {  #  Needed for unit-tests
  remoteScriptDir <- "web-app/HeimScripts"
}

## Loading functions ##
utils <- paste(remoteScriptDir, "/_shared_functions/Generic/utils.R", sep="")
limmaUtils <- paste(remoteScriptDir, "/_shared_functions/GEX/limmaUtils.R", sep="")
dataFrameUtils <- paste(remoteScriptDir, "/_shared_functions/GEX/DataFrameAndGEXmatrixUtils.R", sep="")


source(utils)
source(limmaUtils)
source(dataFrameUtils)


SUBSET1REGEX <- "_s1$"  # Regex identifying columns of subset 1.
markerTableJson <- "markerSelectionTable.json" # Name of the json file with limma outputs




main <- function() {


  ## Get Gene Expression Matrix as data frame
  data.list = parseInput()
  df = data.list$HD
  
  
  ## Defining the content for these variables to
  ## perform differential expression analysis:
  max_rows = dim(df)[1]
  sorting = "nodes"
  ranking = "pval"


  ## File containing the original
  ## GEX values that can be downloaded
  ## by the user
  write.table(
    df,
    "volcanoplot_orig_values.tsv",
    sep = "\t",
    na = "",
    row.names = FALSE,
    col.names = TRUE
  )

  df          <- addStats(df, sorting, ranking, max_rows)
  df          <- mergeDuplicates(df)
  df          <- df[1:min(max_rows, nrow(df)), ]  #  apply max_rows

  fields      <- buildFields(df)
  uids        <- df[, 1]
  patientIDs  <- unique(fields["PATIENTID"])[,1]

  negativeLog10PvalValues = -log10(df["PVAL"][,1])

  ## Output json object containing results
  jsn <- list(
    "uids"                    = uids,
    "logfoldValues"           = df["LOGFOLD"][,1],
    "pvalValues"              = df["PVAL"][,1],
    "negativeLog10PvalValues" = negativeLog10PvalValues,
    "patientIDs"              = patientIDs,
    "warnings"                = c() # initiate empty vector
  )
  jsn <- toJSON(jsn)

  writeRunParams()

  measurements <- cleanUp(df)  # temporary stats like SD and MEAN need to be removed for clustering to work

  write(jsn, file = "volcanoplot.json")
  # json file be served the same way
  # like any other file would - get name via
  # /status call and then /download

  msgs <- c("Finished successfuly")
  list(messages = msgs)
  
  #return(df)
}


#####################
#####################


## #SE: For dev purposes we call the function here
##out = main()


