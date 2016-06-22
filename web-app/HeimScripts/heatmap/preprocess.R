library(WGCNA)

# Pre-processing data saved in preprocessed
# Pre-processing parameters saved in global preprocessing_params
# They can be considered fresh if the variable preprocessed exists


# SE: Just to get things working for dev purposes
#rm(list = ls())
#load("/Users/serge/Documents/Projects/SmartR/Development_env_Input_workspace/R_workspace_objects/Heatmap/data.Rda")
#load("/Users/serge/Documents/Projects/SmartR/Development_env_Input_workspace/R_workspace_objects/Heatmap/fetchParams.Rda")
#load("/Users/serge/Documents/Projects/SmartR/Development_env_Input_workspace/R_workspace_objects/Heatmap/loaded_variables_withLDD.Rda")
#load("/Users/serge/Documents/Projects/SmartR/Development_env_Input_workspace/R_workspace_objects/Heatmap/fetch_params_withLDD.Rda")
#setwd("/Users/serge/GitHub/SmartR")
#######


if (!exists("remoteScriptDir")) {  #  Needed for unit-tests
  remoteScriptDir <- "web-app/HeimScripts"
}

## Loading functions ##
utils <- paste(remoteScriptDir, "/_shared_functions/Generic/utils.R", sep="")
dataFrameUtils <- paste(remoteScriptDir, "/_shared_functions/GEX/DataFrameAndGEXmatrixUtils.R", sep="")

source(utils)
source(dataFrameUtils)


## SE: Handles Low and high dimensional data. Creates preprocessed variable as a list storing HD and LD data
## In case no LD data is available this list item is set to NULL
main <- function(aggregate=FALSE) {
  msgs = c("")
  
  df <- loaded_variables[grep("highDimensional", names(loaded_variables))]

  if (aggregate && is.null(dim(df))) {
      stop("Cannot aggregate probes: there only is data for a single probe (ie. only one row of data) or there is insufficient bio.marker information for the selected probes to be able to match the probes to biomarkers for aggregation (e.g. in case of micro-array data to match probe ID to gene symbol). Suggestion: skip probe aggregation.")
  } 

  ## High dimensional data
  df <- mergeFetchedData(loaded_variables[grep("highDimensional", names(loaded_variables))])
  ## Low dimensional data
  ld.idx = grep( "^(categoric)|(numeric)", names(loaded_variables))
  
  if(length(ld.idx)==0)
    ld = c()
  else
    ld = loaded_variables[ld.idx]
  
  
  ## SE: Modified cutoff for samples to 2 instead of 3
  good.input <- ncol(df) > 2 && nrow(df) > 0 && sum(df$Bio.marker != "") > 0 # min one sample, contains any rows, non empty Bio.marker column.
  
  if(aggregate && good.input){
    df <- dropEmptyGene(df)
    
    aggr  <- aggregate.probes(df)
    
   
    
    ## SE: modified assign value to list type to handle also low dim data
    assign("preprocessed", list(HD = aggr, LD=ld), envir = .GlobalEnv)
    
    discarded.rows <- nrow(df) - nrow(aggr)
    msgs <- paste("Total discarded rows:",discarded.rows)
  }
  else if(aggregate && !good.input){
    stop("Incorrect subset - in order to perform probe aggregation more than one sample is needed.")
  }else{
    msgs <- c("No preprocessing applied.")
    
    ## SE: modified assign value to list type to handle also low dim data
      assign("preprocessed", list(HD = df, LD = ld), envir = .GlobalEnv)
  }

  assign("preprocessing_params", list(aggregate=aggregate), envir = .GlobalEnv)

  
  list(finished=T,messages=msgs)
 
  ## SE: For debug 
  #return(df)
}




## SE: For debug
#main(aggregate=T)
# print(dim(df))
