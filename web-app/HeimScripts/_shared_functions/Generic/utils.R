######################################################################
## This file contains  generic functionalities when working with    ##
## SmartR data objects                                              ##                           
######################################################################



## Returns a vector providing the data tree node names 
## as a named vector where a node name identifier like "numeric_n0" is mapped to
## its real name in the tree (e.g. "Alive")
node2name = function(){
  
  names = character(length = length(names(fetch_params$ontologyTerms)))
  node_names = character(length = length(names(fetch_params$ontologyTerms)))
  
  for(i in 1:length(names(fetch_params$ontologyTerms))){
    node_name = names(fetch_params$ontologyTerms)[i]
    name = get(names(fetch_params$ontologyTerms)[i], fetch_params$ontologyTerms )$name
    
    names[i] = name
    node_names[i] = node_name
  }
  names(names) = node_names
  
  return(names)
}



## Retrieve High dim node for vector of patientIDs
getNode <- function(patientIDs) {

    splittedIds <- strsplit(patientIDs,"_") # During merge, which is always
    # run we append subset id, either
    # _s1 or _s2 to PATIENTID.
    
    sapply(splittedIds, FUN = tail_elem,n = 2) # In proper patienid subset will
    # always be  at the end.
    # This select last but one elemnt
    # - the node
}



getTimelineValues <- function(nodes, ontologyTerms) {
    sapply(nodes, function(n) {
        metaValue <- ontologyTerms[[n]]$metadata$seriesMeta$value
        t <- !is.na(as.numeric(metaValue))
        if (length(t) && t) {
            as.numeric(metaValue)
        } else {
            Inf
        }
    }, USE.NAMES = FALSE)
}

getSubject <- function(str) {
    sub("_.*", "", str)
}

## This function dumps the run parameters into a json format file
## The function can be called either
## ... with three unnamed params providing the following variables/info in this order:
## -- max_rows => integer()
## -- sorting => character("nodes", "subjects")
## -- ranking => charcter("logfold", "ttest", "pval", "adjpval", "mean", "median", "coef", "range")
## e.g.: "writeRunParams(max_rows, sorting, ranking)"
## ... or with a group of named params:
## e.g.: "writeRunParams(some_parameter=some_parameter, another_param=another_param, ... )"
writeRunParams <- function(...) {
  
  ## Get the data objects passed
  params <- list(...)
  
  ## The default way of calling the function, this
  ## means providing 3 args and without named params
  if(length(names(params))==0 & length(params)==3) 
    names(params) = c("max_rows", "sorting", "ranking")  
  
  ## Get name of this function
  calledFun = as.list(sys.call())[[1]]
  
  ## Checking if params is a named list
  ## this means named parameters have been
  ## provided as input
  if(length(names(params))==0 & length(params)!=0)
    print(paste("Please provided named function parameters as input for ", calledFun, sep=": "))
  
  params <- c(params, fetch_params)
  
  if (exists("preprocessed") && exists("preprocessing_params")) {
    params <- c(params, preprocessing_params)
  }
  write(toJSON(params, pretty = TRUE), 'params.json')
}




