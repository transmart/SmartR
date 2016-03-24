library(reshape2)
library(limma)

if (!exists("remoteScriptDir")) {  #  Needed for unit-tests
    remoteScriptDir <- "web-app/HeimScripts"
}

utils <- paste(remoteScriptDir, "/heatmap/utils.R", sep="")
limmaUtils <- paste(remoteScriptDir, "/heatmap/limmaUtils.R", sep="")
source(utils)
source(limmaUtils)


SUBSET1REGEX <- "_s1$"  # Regex identifying columns of subset 1.
markerTableJson <- "markerSelectionTable.json" # Name of the json file with limma outputs

main <- function(max_rows = 100,
                 sorting = "nodes",
                 ranking = "coef",
                 pdMapLogin = "",
                 pdMapPassword = "",
                 pdMapLinkEnabled = false) {
    max_rows <- as.numeric(max_rows)
    verifyInput(max_rows, sorting)
    df <- parseInput()
    write.table(
        df,
        "heatmap_orig_values.tsv",
        sep = "\t",
        na = "",
        row.names = FALSE,
        col.names = TRUE
    )
    df          <- addStats(df, sorting, ranking, max_rows)
    df          <- mergeDuplicates(df)
    df          <- df[1:min(max_rows, nrow(df)), ]  #  apply max_rows
    fields      <- buildFields(df)
    extraFields <- buildExtraFields(fields)
    uids        <- df[, 1]
    patientIDs  <- unique(fields["PATIENTID"])[,1]

    significanceValues <- df["SIGNIFICANCE"][,1]
    features <- unique(extraFields["FEATURE"])[,1]
    jsn <- list(
        "fields"             = fields,
        "patientIDs"         = patientIDs,
        "uids"               = uids,
        "significanceValues" = significanceValues,
        "logfoldValues"      = df["LOGFOLD"][,1],
        "ttestValues"        = df["TTEST"][,1],
        "pvalValues"         = df["PVAL"][,1],
        "adjpvalValues"      = df["ADJPVAL"][,1],
        "bvalValues"         = df["BVAL"][,1],
        "ranking"            = ranking,
        "features"           = features,
        "extraFields"        = extraFields,
        "maxRows"            = max_rows,
        "pdMapLinkEnabled"   = pdMapLinkEnabled,
        "pdMapLogin"         = pdMapLogin,
        "pdMapPassword"      = pdMapPassword,
        "warnings"           = c() # initiate empty vector
    )
    writeRunParams(max_rows, sorting, ranking)
    measurements <- cleanUp(df)  # temporary stats like SD and MEAN need
    # to be removed for clustering to work

    # discard UID column
    if (ncol(df) > 2){
        measurements <- measurements[, 2:ncol(measurements)]
    } else {
        measurements <- measurements[2]
    }

    measurements <- toZscores(measurements)
    if (is.na(significanceValues)) {
        jsn$warnings <- append(jsn$warnings, c("Significance sorting could not be done due to insufficient data"))
    }
    jsn <- addClusteringOutput(jsn, measurements) #
    jsn <- toJSON(jsn, pretty = TRUE, digits = I(17))
    writeDataForZip(df, measurements, patientIDs)  # for later zip generation
    write(jsn, file = "heatmap.json")
    # json file be served the same way
    # like any other file would - get name via
    # /status call and then /download

    msgs <- c("Finished successfuly")
    list(messages = msgs)
}

writeMarkerTable <- function(markerTable){
    colnames(markerTable) <- c("rowLabel", "biomarker",
    "log2FoldChange", "t", "pValue", "adjustedPValue", "B")
    jsn                   <- toJSON(markerTable, pretty = TRUE, digits = I(17))
    write(jsn, file = markerTableJson)
}

cleanUpLimmaOutput <- function() {
    if (file.exists(markerTableJson)) {
        file.remove(markerTableJson)
    }
}

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

addStats <- function(df, sorting, ranking, max_rows) {
    measurements  <- getMeasurements(df)
    rankingMethod <- getRankingMethod(ranking)
    twoSubsets    <- hasTwoSubsets(measurements)

    logfold.values <- data.frame(LOGFOLD=numeric())
    ttest.values <- data.frame(TTEST=numeric())
    pval.values <- data.frame(PVAL=numeric())
    adjpval.values <- data.frame(ADJPVAL=numeric())
    bval.values <- data.frame(BVAL=numeric())
    rankingScore <- data.frame(SIGNIFICANCE=numeric())

    if (ncol(df) > 3) {
        #this is the case for more than 1 sample
        df <- applySorting(df,sorting)  # no need for sorting in one sample case, we do it here only

        # is it valid measurements ?
        validLimmaMeasurements <- isValidLimmaMeasurements(measurements)

        if (twoSubsets && validLimmaMeasurements) {
            markerTable  <- getDEgenes(df)  # relies on columns being sorted,
        }                                           # with subset1 first, this is because of the way
        # design matrix is being constructed


        useLimma <- !is.function(rankingMethod)  # rankingMethod is either a function or a character "limma"
        if (useLimma  && !twoSubsets) {  # cannot use any of the limma methods for single subset.
            stop( paste("Illegal ranking method: ", ranking, " two subsets needed.") )
        }

        if (useLimma && validLimmaMeasurements) {
            if (!ranking %in% colnames(markerTable) ) {
                stop(paste("Illegal ranking method selected: ", ranking) )
            }
            logfold.values <- markerTable["logfold"]
            ttest.values <- markerTable["ttest"]
            pval.values <- markerTable["pval"]
            adjpval.values <- markerTable["adjpval"]
            bval.values <- markerTable["bval"]

            rankingScore <- markerTable[ranking]
        } else if (useLimma && !validLimmaMeasurements)  {
            # when differential expression rank criteria
            # is selected while it's not valid limma measurements,
            # provide empty data frame to the significance columns.
        } else {
            rankingScore <- apply(measurements, 1, rankingMethod, na.rm = TRUE )  # Calculating
        }
        # ranking per probe (per row)
        means <- rowMeans(measurements, na.rm = T)  # this is just an
                                                    # auxiliary column - it will not be
                                                    # used for JSON.
        sdses <- apply(measurements,1 ,sd , na.rm = T)  # this is just
                                                        # an auxiliary column -
                                                        # it will not be used for JSON.

        cleanUpLimmaOutput()  # In order to prevent displaying the table from previous run.

        if (twoSubsets && validLimmaMeasurements) {
            markerTable["SIGNIFICANCE"] <- rankingScore
            markerTable                 <- applyRanking(markerTable, ranking, max_rows)
            markerTable["SIGNIFICANCE"] <- NULL
            writeMarkerTable(markerTable)
        }
    } else {
        #one sample
        variance <- 1.0  # We cannot get any significance measure for just
        # one sample
        if ( nrow(measurements) > 1) {
            meanS <- mean(df[,3], na.rm = T)  # For just one sample we take mean
            # over all genes (whole column), not per gene
            sds <- sd(df[,3],na.rm = T)  # For just one sample we take mean
            # over all genes (whole column), not per gene
        } else {
            meanS <- 1
            sds  <- 1
        }

        means <- rep(meanS, nrow(df))
        sdses <- rep(sds, nrow(df))
        rankingScore <- rep(variance, nrow(df))
    }

    df["MEAN"]         <- means
    df["SD"]           <- sdses
    df["SIGNIFICANCE"] <- rankingScore
    df["LOGFOLD"] <- logfold.values
    df["TTEST"] <- ttest.values
    df["PVAL"] <- pval.values
    df["ADJPVAL"] <- adjpval.values
    df["BVAL"] <- bval.values

    if (useLimma & !validLimmaMeasurements) {
        # dont apply ranking
    } else {
        # else apply
        df <- applyRanking(df, ranking, max_rows)
    }

    df
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


toZscores <- function(measurements) {
    measurements <- scale(t(measurements))
    t(measurements)
}

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

verifyInput <- function(max_rows, sorting) {
    if (max_rows <= 0) {
        stop("Max rows argument needs to be higher than zero.")
    }
    if (!(sorting == "nodes" || sorting == "subjects")) {
        stop("Unsupported sorting type. Only nodes and subjects allowed")
    }
}

applySorting <- function(df,sorting) {
    measurements <- getMeasurements(df)
    colNames <- names(measurements)
    subsets <- getSubset(colNames)
    nodes <- getNode(colNames)
    subjects <- getSubject(colNames)

    timelineValues <- getTimelineValues(nodes, fetch_params$ontologyTerms)

    if (sorting == "nodes") {
        inds <- order(subsets, timelineValues, nodes, subjects)
    } else {
        inds <- order(subsets, subjects, timelineValues, nodes)
    }

    measurements <- measurements[, inds]
    cbind(df[, c(1,2)], measurements)
}

buildFields <- function(df) {
    df <- melt(df, na.rm=T, id=c("UID", "MEAN", "SD", "SIGNIFICANCE", "LOGFOLD", "TTEST", "PVAL", "ADJPVAL", "BVAL"))
    # melt implicitly casts
    # characters to factors to make your life more exciting,
    # in order to encourage more adventures it does not
    # have characters.as.factors=F param.
    df$variable <- as.character(df$variable)
    ZSCORE      <- (df$value - df$MEAN) / df$SD
    df["MEAN"]  <- NULL
    df["SD"]    <- NULL

    names(df)   <- c("UID", "SIGNIFICANCE", "LOGFOLD", "TTEST", "PVAL", "ADJPVAL", "BVAL", "PATIENTID", "VALUE")
    df["ZSCORE"] <- ZSCORE
    df["SUBSET"] <- getSubset(df$PATIENTID)

    df$PATIENTID <- replaceNodeIDNodeLabel(df$PATIENTID, fetch_params$ontologyTerms)

    return(df)
}

writeDataForZip <- function(df, zScores, pidCols) {
    df      <- df[ , -which(names(df) %in% pidCols)]  # Drop patient columns
    df      <- cbind(df,zScores)                      # Replace with zScores
    write.table(
        df,
        "heatmap_data.tsv",
        sep = "\t",
        na = "",
        row.names = FALSE,
        col.names = TRUE
    )
}

writeRunParams <- function(max_rows, sorting, ranking) {
    params <- list(max_rows = max_rows,
    sorting  = sorting,
    ranking  = ranking)
    params <- c(params, fetch_params)
    if (exists("preprocessed") && exists("preprocessing_params")) {
        params <- c(params, preprocessing_params)
    }
    write(toJSON(params, pretty = TRUE), 'params.json')
}


buildExtraFields <- function(df) {
    FEATURE <- rep("Cohort", nrow(df))
    PATIENTID <- as.character(df$PATIENTID)
    TYPE <- rep("subset",nrow(df))
    VALUE <- df$SUBSET
    extraFields <- data.frame(FEATURE, PATIENTID, TYPE, VALUE, stringsAsFactors = FALSE)
}

computeDendrogram <- function(distances, linkageMethod) {
    as.dendrogram(hclust(distances, method = linkageMethod))
}

dendrogramToJSON <- function(d) {
    totalMembers <- attributes(d)$members
    add_json <- function(x, start) {
        members <- attributes(x)$members
        height <- attributes(x)$height
        index <- (start - 1):(start + members - 2)
        index <- paste(index, collapse = ' ')
        jsonString <<- paste(
            jsonString,
            sprintf('{"height":"%s", "index":"%s", "children":[', height, index)
        )
        if (is.leaf(x)) {
            jsonString <<- paste(jsonString, ']}')
        } else {
            add_json(x[[1]], start)
            jsonString <<- paste(jsonString, ",")
            leftMembers <- attributes(x[[1]])$members
            add_json(x[[2]], start + leftMembers)
            jsonString <<- paste(jsonString, "]}")
        }
    }
    jsonString <- ""
    add_json(d, 1)
    return(jsonString)
}


addClusteringOutput <- function(jsn, measurements_arg) {
    # we need to discard rows and columns without at least two values
    # determine rows without two non-NAs
    thresholdRows <- floor(ncol(measurements_arg) / 2 ) + 1  #  Half of the lentgh +1 has to be filled otherwise
    thresholdCols <- floor(nrow(measurements_arg) / 2 ) + 1  #  there might not be enough overlap for clustering
    logicalSelection <- apply(measurements_arg, 1, function(row) length(which(!is.na(row))) >= thresholdRows)
    measurements_rows <- measurements_arg[logicalSelection, ]
    # and for columns
    logicalSelection <- apply(measurements_arg, 2, function(col) length(which(!is.na(col))) >= thresholdCols)
    measurements_cols <- t(measurements_arg[ , logicalSelection])
    jsn$numberOfClusteredRows <- nrow(measurements_rows)
    jsn$numberOfClusteredColumns <- nrow(measurements_cols)  # still nrow (transposed)
    if (is.null(jsn$numberOfClusteredRows)) jsn$numberOfClusteredRows <- 0
    if (is.null(jsn$numberOfClusteredColumns)) jsn$numberOfClusteredColumns <- 0
    if (jsn$numberOfClusteredRows < 2 | jsn$numberOfClusteredColumns < 2 ) {  # Cannot cluster less than 2x2 matrix
        jsn$warnings <- append(jsn$warnings, c("Clustering could not be done due to insufficient data"))
        return(jsn)
    }

    euclideanDistancesRow <- dist(measurements_rows, method = "euclidean")
    manhattanDistancesRow <- dist(measurements_rows, method = "manhattan")
    euclideanDistancesCol <- dist(measurements_cols, method = "euclidean")
    manhattanDistancesCol <- dist(measurements_cols, method = "manhattan")

    colDendrogramEuclideanComplete <- computeDendrogram( euclideanDistancesCol, 'complete')
    colDendrogramEuclideanSingle <- computeDendrogram( euclideanDistancesCol, 'single')
    colDendrogramEuclideanAverage <- computeDendrogram( euclideanDistancesCol, 'average')
    rowDendrogramEuclideanComplete <- computeDendrogram( euclideanDistancesRow, 'complete')
    rowDendrogramEuclideanSingle <- computeDendrogram( euclideanDistancesRow, 'single')
    rowDendrogramEuclideanAverage <- computeDendrogram( euclideanDistancesRow, 'average')

    colDendrogramManhattanComplete <- computeDendrogram( manhattanDistancesCol, 'complete')
    colDendrogramManhattanSingle <- computeDendrogram( manhattanDistancesCol, 'single')
    colDendrogramManhattanAverage <- computeDendrogram( manhattanDistancesCol, 'average')
    rowDendrogramManhattanComplete <- computeDendrogram( manhattanDistancesRow, 'complete')
    rowDendrogramManhattanSingle <- computeDendrogram( manhattanDistancesRow, 'single')
    rowDendrogramManhattanAverage <- computeDendrogram( manhattanDistancesRow, 'average')

    calculateOrderInTermsOfIndexes <- function(dendrogram, originalOrderedLabels) {
        clusterOrderedLabels  <- labels(dendrogram)

        allIndexes <- 1 : length(originalOrderedLabels)
        orderOfClustered <- match(clusterOrderedLabels, originalOrderedLabels)

        # put the elements that were not clustered at the end
        notIncluded <- allIndexes[!is.element(allIndexes, orderOfClustered)]
        c(orderOfClustered, notIncluded) - 1  # start at 0
    }

    columnOrder <- function(dendrogram) {
        calculateOrderInTermsOfIndexes(dendrogram, colnames(measurements_arg))
    }
    rowOrder <- function(dendrogram) {
        calculateOrderInTermsOfIndexes(dendrogram, rownames(measurements_arg))
    }


    jsn$hclustEuclideanComplete <- list(
        columnOrder(colDendrogramEuclideanComplete),
        rowOrder(rowDendrogramEuclideanComplete),
        dendrogramToJSON(colDendrogramEuclideanComplete),
        dendrogramToJSON(rowDendrogramEuclideanComplete)
    )

    jsn$hclustEuclideanSingle <- list(
        columnOrder(colDendrogramEuclideanSingle),
        rowOrder(rowDendrogramEuclideanSingle),
        dendrogramToJSON(colDendrogramEuclideanSingle),
        dendrogramToJSON(rowDendrogramEuclideanSingle)
    )

    jsn$hclustEuclideanAverage <- list(
        columnOrder(colDendrogramEuclideanAverage),
        rowOrder(rowDendrogramEuclideanAverage),
        dendrogramToJSON(colDendrogramEuclideanAverage),
        dendrogramToJSON(rowDendrogramEuclideanAverage)
    )

    jsn$hclustManhattanComplete <- list(
        columnOrder(colDendrogramManhattanComplete),
        rowOrder(rowDendrogramManhattanComplete),
        dendrogramToJSON(colDendrogramManhattanComplete),
        dendrogramToJSON(rowDendrogramManhattanComplete)
    )

    jsn$hclustManhattanSingle <- list(
        columnOrder(colDendrogramManhattanSingle),
        rowOrder(rowDendrogramManhattanSingle),
        dendrogramToJSON(colDendrogramManhattanSingle),
        dendrogramToJSON(rowDendrogramManhattanSingle)
    )

    jsn$hclustManhattanAverage <- list(
        columnOrder(colDendrogramManhattanAverage),
        rowOrder(rowDendrogramManhattanAverage),
        dendrogramToJSON(colDendrogramManhattanAverage),
        dendrogramToJSON(rowDendrogramManhattanAverage)
    )
    return(jsn)
}
