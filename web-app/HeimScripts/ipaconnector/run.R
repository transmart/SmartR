library(reshape2)
library(limma)
library(jsonlite)

if (!exists("remoteScriptDir")) {  #  Needed for unit-tests
    remoteScriptDir <- "web-app/HeimScripts"
}


## Loading functions ##
utils <- paste(remoteScriptDir, "/_shared_functions/Generic/utils.R", sep="")
limmaUtils <- paste(remoteScriptDir, "/_shared_functions/GEX/limmaUtils.R", sep="")
dataFrameUtils <- paste(remoteScriptDir, "/_shared_functions/GEX/DataFrameAndGEXmatrixUtils.R", sep="")
heatmapUtils <- paste(remoteScriptDir, "/_shared_functions/Clustering/heatmapUtils.R", sep="")

source(utils)
source(limmaUtils)
source(dataFrameUtils)
source(heatmapUtils)
#######################


# Regex identifying columns of subset 1
SUBSET1REGEX <- "_s1$"

# Name of the json output file with DEA limma outputs
markerTableJson <- "ipaconnector.json"


## Main function to generate json output file for limma-based DEA
## Input params:
## - significanceMeasure: the measure to define statistical significance. Possible values are 'pval' and 'adjpval'.
##                        'pval' corresponds to the raw p-value. 'adjpval' corresponds to BH FDR in limma.
## - significanceCutoff: statistical cutoff for filtering the DEG list
## - fcCutoff: log2 fold change cutoff. The absolute value of this will be used for filtering the list of DEG
##
## Output: a json format object file storing the following info as vector:
## - probeID
## - symbol
## - logFC
## - pval
## - adjPval


main <- function(significanceMeasure = 'pval', significanceCutoff = 0.05, fcCutoff = 1.5) {

    ### Input arg check for numeric data##

    BadSignificanceMeasure = FALSE
    BadSignificanceCutoff  = FALSE
    BadFcCutoff            = FALSE

    if( ! significanceMeasure %in% c('pval', 'adjpval')){
      BadSignificanceMeasure = TRUE
      significanceMeasure = 'pval'
    }

    if( ! is.numeric(significanceCutoff) )
      BadSignificanceCutoff = TRUE

    if( significanceCutoff > 1 )
      BadSignificanceCutoff = TRUE

    if( significanceCutoff < 0 )
      BadSignificanceCutoff = TRUE

    if( ! is.numeric(fcCutoff) ){
      BadFcCutoff = TRUE
      fcCutoff    = 0
    }


    ## Returns a list containing two variables named HD and LD
    data.list <- parseInput()

    ## Splitting up input into low dim and high dim vars
    hd.df = data.list$HD

    # Just the number of probes tested when performing DEA
    max_rows = dim(hd.df)[1]

    ## Creating the extended diff expr analysis data frame containing besides the input data,
    ## a set of statistics. The returned data frame is ranked according to provided ranking statistic
    hd.df          <- addStats(hd.df, significanceMeasure, max_rows)

    hd.df          <- mergeDuplicates(hd.df)

    ## A df containing the computed values for
    ## all possible statistical methods
    statistics_hd.df = getAllStatForExtDataFrame(hd.df)


    probe_annot = matrix(unlist(strsplit(as.character(statistics_hd.df[,"ROWNAME"]), "--")), ncol = 2, byrow = TRUE)
    probeID = as.character(probe_annot[,1])
    symbol = as.character(probe_annot[,2])

    statistics_hd.df = cbind(SYMBOL = symbol, statistics_hd.df)
    statistics_hd.df = cbind(PROBEID = probeID, statistics_hd.df)

    statistics_hd.df =  statistics_hd.df[, c("PROBEID", "SYMBOL", "LOGFOLD", "PVAL", "ADJPVAL")]


    ## Filtering significant probes ##
    fc_cutoff_sig.idx = which(abs(statistics_hd.df[, "LOGFOLD"]) >= abs(fcCutoff))

    if(significanceMeasure == "adjpval"){
      pval_cutoff_sig.idx = which(abs(statistics_hd.df[, "ADJPVAL"]) < significanceCutoff)
    } else if(significanceMeasure == "pval"){
      pval_cutoff_sig.idx = which(abs(statistics_hd.df[, "PVAL"]) < significanceCutoff)
    }

    stat_and_biol_sigProbes.idx = intersect(pval_cutoff_sig.idx, fc_cutoff_sig.idx)
    statistics_hd.df = statistics_hd.df[stat_and_biol_sigProbes.idx, ]


    ## Ordering the data frame according to p-value/statistical significance
    statistics_hd.df = statistics_hd.df[order(statistics_hd.df[, "PVAL"]), ]


    ## The returned jsn object that will be dumped to file
    jsn <- list(
        'header'    = c('probeID', 'symbol', 'logFC', 'pval', 'adjPval'),
        'columnTypes' = c('STRING', 'STRING', 'NUMERIC', 'NUMERIC', 'NUMERIC'),
        'probeID'   = statistics_hd.df[, "PROBEID"],
        'symbol'    = statistics_hd.df[, "SYMBOL"],
        'logFC'     = statistics_hd.df[, "LOGFOLD"],
        'pval'      = statistics_hd.df[, "PVAL"],
        'adjPval'   = statistics_hd.df[, "ADJPVAL"],
        'warnings'  = c() # initiate empty vector
    )


    ## Checking if there are any DEGs available
    if(dim(statistics_hd.df)[1]==0)
      noSignificanceValues = TRUE
    else
      noSignificanceValues = FALSE


    ## If no significanceValues are available throw a warning:
    if (noSignificanceValues)
        jsn$warnings <- append(jsn$warnings, c("No differentially expressed genes detected."))


    ## If statistical significance measure input argument is not correct throw a warning:
    if(BadSignificanceMeasure)
      jsn$warnings <- append(jsn$warnings, c("Statistical significance measure parameter is not correct. Using default significance measure (raw p-value)."))


    ## If statistical cutoff input argument is not correct throw a warning:
    if(BadSignificanceCutoff)
      jsn$warnings <- append(jsn$warnings, c("Statistical significance cutoff parameter is not correct. DEG list could not be filtered based on this parameter."))


    ## If fold change cutoff input argument is not correct throw a warning:
    if(BadFcCutoff)
      jsn$warnings <- append(jsn$warnings, c("Fold change cutoff parameter is not correct. DEG list could not be filtered based on this parameter."))


    ## Transforming the output list to json format
    jsn <- toJSON(jsn, pretty = TRUE, digits = I(17))

    write(jsn, file = markerTableJson)

    msgs <- c("Finished successfuly")
    list(messages = msgs)
}
