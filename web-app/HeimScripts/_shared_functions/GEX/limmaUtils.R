################################################################################
## This file contains more specific limma and limma-related functionalities   ##
## to check the input data frame content and perform diff expr analysis for   ##
## Microarray GEX data                                                        ##
################################################################################




## Loading functions ##
dataFrameUtils <- paste(remoteScriptDir, "/_shared_functions/GEX/DataFrameAndGEXmatrixUtils.R", sep="")

source(dataFrameUtils)


## Differential expression analysis comparing Subset 1 to subset 2 
## using Limma package
getDEgenes <- function(df) {
  measurements    <- getMeasurements(df)
  design          <- getDesign(measurements)
  contrast.matrix <- makeContrasts( S2-S1, levels = design )
  fit             <- lmFit(measurements, design)
  fit             <- contrasts.fit(fit, contrast.matrix)
  fit             <- eBayes(fit)
  contr           <- 1  #  We need a vector, not a df, so we'll do [, contr] on all stats
  top.fit         <- data.frame (
    logfold = fit$coefficients[, contr],
    ttest   = fit$t[, contr],
    pval    = fit$p.value[, contr],
    adjpval = p.adjust(
    p       = fit$p.value[, contr],
    method  ='fdr'),
    bval    = fit$lods[, contr]
  )
  
  cbind(df[,1:2], top.fit)
}



## Generate the design matrix
getDesign <- function(measurements) {
  
  subsets <- getSubset(colnames(measurements)) #s1 = 1, s2 = 2
  classVectorS1 <- subsets             #s1 = 1, s2 = 2
  classVectorS2 <- - subsets + 3       #s1 = 2, s2 = 1
  
  DesignMatrix = cbind(S1=subsets, S2=subsets)
  DesignMatrix[which(DesignMatrix[, "S1"]==2), "S1"]=0
  DesignMatrix[which(DesignMatrix[, "S2"]==1), "S2"]=0
  DesignMatrix[which(DesignMatrix[, "S2"]==2), "S2"]=1

  return(DesignMatrix)

}



# checking valid measurements on a matrix level
isValidLimmaMeasurements <- function (measurements) {
    sum(apply(measurements, 1, validMeasurementsRow)) > 0
}



## Writing jsn file containing the Limma diff expr analysis results
writeMarkerTable <- function(markerTable, markerTableJson = "markerSelectionTable.json"){
  colnames(markerTable) <- c("rowLabel", "biomarker",
                             "log2FoldChange", "t", "pValue", "adjustedPValue", "B")
  jsn                   <- toJSON(markerTable, pretty = TRUE, digits = I(17))
  write(jsn, file = markerTableJson)
}



## Deleting the json-format topTable
## from last Limma diff expr analysis run
cleanUpLimmaOutput <- function(markerTableJson = "markerSelectionTable.json") {
  if (file.exists(markerTableJson)) {
    file.remove(markerTableJson)
  }
}






