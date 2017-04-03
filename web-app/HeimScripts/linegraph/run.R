library(jsonlite)
library(reshape2)

main <- function() {
    df <- buildCrossfilterCompatibleDf(loaded_variables, fetch_params)
    checkTimeNameSanity(df)

    numeric.stats.df <- getStatsForNumericType(df)
    df <- merge(df, numeric.stats.df, by=c("bioMarker", "timeInteger", "subset"), all=TRUE)

    output <- list()
    output$data_matrix <- df

    # make available to other scripts
    runResults <<- output

    json <- toJSON(output, pretty=TRUE, digits=I(17))
    write(json, file="linegraph.json")
    list(messages="Finished successfully")
}

# returns character vector (e.g. c("\\Demo Study\\Vital Status\\Alive\\Week1", "\\Demo Study\\Vital Status\\Alive\\Week2", ...))
getFullNames <- function(loaded_variables, fetch_params) {
    names.without.subset <- sub("_s[1-2]{1}$", "", names(loaded_variables))
    fullNames <- sapply(names.without.subset, function(el) fetch_params$ontologyTerms[[el]]$fullName)
    as.character(as.vector(fullNames))
}

# returns integer vector (e.g. c(1,2,2,2,1,1,2))
getSubsets <- function(loaded_variables) {
    subsets <- sub("^.*_s", "", names(loaded_variables))
    as.integer(subsets)
}

# returns character vector (e.g. c("numeric", "numeric", "categoric", ...))
getTypes <- function(loaded_variables) {
    types <- sub("_.*$", "", names(loaded_variables))
    types[types == "highData"] <- "highDimensional"
    types[types == "numData"] <- "numeric"
    types[types == "catData"] <- "categoric"
    as.character(types)
}

# returns integer of extracted time if possible, NULL otherwise
extractTime <- function(string) {
    match <- regmatches(string, regexpr("((?>- )|-)?\\d+(?>\\.\\d+)?", string, perl=TRUE))
    timeInteger <- as.numeric(match[1])
    if (!is.na(timeInteger) && timeInteger%%1 == 0) {
        return(timeInteger)
    }
    NULL
}

# returns df that is compatible with crossfilter.js
buildCrossfilterCompatibleDf <- function(loaded_variables, fetch_params) {
    # gather information
    allSubsets <- getSubsets(loaded_variables)
    fullNames <- getFullNames(loaded_variables, fetch_params)
    types <- getTypes(loaded_variables)

    # initialize empty df
    df <- data.frame(patientID=character(),
                     value=c(), # can be string or integer
                     timeInteger=integer(),
                     timeString=character(),
                     bioMarker=character(),
                     type=character(),
                     subset=integer(),
                     ranking=integer(),
                     stringsAsFactors=FALSE)


    # maps name to numeric time (e.g. "Week 12" to 12)
    times <- list()
    # build big df step by step via binding row-wise every loaded variable
    for (i in 1:length(names(loaded_variables))) {
        variable <- loaded_variables[[i]]
        variable.df <- data.frame()

        if (types[i] == "highDimensional") {
            colnames(variable)[-(1:2)] <- sub("^X", "", colnames(variable[-(1:2)]))
            variable.df <- melt(variable, id.vars=c("Row.Label", "Bio.marker"), variable.name="patientID")

            variable.label.df <- data.frame()
            # create for each unique row label an own df
            for (rowLabel in unique(variable.df$Row.Label)) {
                variable.label.df <- variable.df[variable.df$Row.Label == rowLabel, ]
                split <- strsplit(fullNames[i], "\\\\")[[1]]
                subsets <- rep(allSubsets[i], nrow(variable.label.df))
                patientIDs <- as.numeric(as.vector(variable.label.df$patientID))
                patientIDs <- paste(patientIDs, subsets, sep="_")
                bioMarker <- paste(split[length(split) - 2], split[length(split) - 1], sep="//")
                bioMarker <- paste(bioMarker, rowLabel, variable.label.df$Bio.marker[1], sep="--")
                timeString <- split[length(split)]
                timeInteger <- times[timeString][[1]]
                # if timeString never occured before, assign it a new timeInteger
                if (is.null(timeInteger)) {
                    extractedTime <- extractTime(timeString)
                    if (is.null(extractedTime)) {
                        timeIntegers <- as.vector(unlist(times))
                        timeInteger <- 0
                        while (TRUE) {
                            if (timeInteger %in% timeIntegers) {
                                timeInteger <- timeInteger + 1
                            } else {
                                break;
                            }
                        }
                    } else {
                        timeInteger <- extractedTime
                    }
                    times[timeString] <- timeInteger
                }
                variable.label.df <- variable.label.df[, -c(1,2)]
                variable.label.df <- data.frame(patientID=patientIDs,
                                                value=as.numeric(as.vector(variable.label.df$value)),
                                                timeInteger=rep(timeInteger, nrow(variable.label.df)),
                                                timeString=rep(timeString, nrow(variable.label.df)),
                                                bioMarker=rep(bioMarker, nrow(variable.label.df)),
                                                type=rep('numeric', nrow(variable.label.df)),
                                                subset=subsets,
                                                ranking=rep(NA, nrow(variable.label.df)),
                                                stringsAsFactors=FALSE)
                # no value -> no interest
                variable.label.df <- variable.label.df[variable.label.df$value != "" & !is.na(variable.label.df$value), ]
                if (nrow(variable.label.df) == 0) next
                df <- rbind(df, variable.label.df)
            }
        } else {
            variable.df <- data.frame(patientID=as.character(variable[,1]), value=variable[,2])
            split <- strsplit(fullNames[i], "\\\\")[[1]]
            subsets <- rep(allSubsets[i], nrow(variable.df))
            patientIDs <- as.character(as.vector(variable.df$patientID))
            patientIDs <- paste(patientIDs, subsets, sep="_")
            bioMarker <- paste(split[length(split) - 2], split[length(split) - 1], sep="//")
            timeString <- split[length(split)]
            timeInteger <- times[timeString][[1]]
            values <- c()
            if (types[i] == "categoric") {
                values <- as.vector(variable.df$value)
                rankings <- rep(0, nrow(variable.df))
            } else {
                values <- as.numeric(as.vector(variable.df$value))
                rankings <- rep(NA, nrow(variable.df))
            }
            # if timeString never occured before, assign it a new timeInteger
            if (is.null(timeInteger)) {
                extractedTime <- extractTime(timeString)
                if (is.null(extractedTime)) {
                    timeIntegers <- as.vector(unlist(times))
                    timeInteger <- 0
                    while (TRUE) {
                        if (timeInteger %in% timeIntegers) {
                            timeInteger <- timeInteger + 1
                        } else {
                            break;
                        }
                    }
                } else {
                    timeInteger <- extractedTime
                }
                times[timeString] <- timeInteger
            }
            # attach additional information
            variable.df <- data.frame(patientID=patientIDs,
                                 value=values,
                                 timeInteger=rep(timeInteger, nrow(variable.df)),
                                 timeString=rep(timeString, nrow(variable.df)),
                                 bioMarker=rep(bioMarker, nrow(variable.df)),
                                 type=rep(types[i], nrow(variable.df)),
                                 subset=subsets,
                                 ranking=rankings,
                                 stringsAsFactors=FALSE)
            # no value -> no interest
            variable.df <- variable.df[variable.df$value != "" & !is.na(variable.df$value), ]
            if (nrow(variable.df) == 0) next
            df <- rbind(df, variable.df)
        }
    }

    # compute inverse freq of each event as a default ranking/order
    totalNumOfEvents <- nrow(df[df$type == "categoric", ])
    for (bioMarker in unique(df$bioMarker)) {
        occurences <- nrow(df[df$bioMarker == bioMarker, ])
        inverseFreq <- 1 / (occurences / totalNumOfEvents)
        df[df$bioMarker == bioMarker, ]$ranking <- inverseFreq
    }


    # before we are done we assign a unique id to every row to make it easier for the front-end
    df <- cbind(id=1:nrow(df), df)

    df
}

# compute several statistics for the given data frame and return a statistic df
getStatsForNumericType <- function(df) {
    numeric.df <- df[df$type == 'numeric', ]
    timeIntegers <- unique(numeric.df$timeInteger)
    bioMarkers <- unique(numeric.df$bioMarker)
    subsets <- unique(numeric.df$subset)

    stats.df <- data.frame()
    for (subset in subsets) {
        for (bioMarker in bioMarkers) {
            for (timeInteger in timeIntegers) {
                current.df <- numeric.df[numeric.df$timeInteger == timeInteger &
                                         numeric.df$bioMarker == bioMarker &
                                         numeric.df$subset == subset, ]
                                     if (nrow(current.df) == 0) next
                                     values <- as.numeric(current.df$value)
                                     mean <- mean(values)
                                     median <- median(values)
                                     sd <- sd(values)
                                     sem <- sd / sqrt(length(values))
                                     stats.df <- rbind(stats.df, data.frame(bioMarker=bioMarker,
                                                                            timeInteger=timeInteger,
                                                                            subset=subset,
                                                                            sd=sd,
                                                                            sem=sem,
                                                                            mean=mean,
                                                                            median=median))
            }
        }
    }
    stats.df
}

# time (e.g. 15) and name (e.g. Day 15) must have a 1:1 relationship
# It is not possible to have multiple times for one name or multiple names for one time
checkTimeNameSanity <- function(df) {
    df.without.duplicates <- unique(df[, c("timeInteger", "timeString")])
    timeSane = nrow(df.without.duplicates) == length(unique(df.without.duplicates$timeInteger))
    nameSane = nrow(df.without.duplicates) == length(unique(df.without.duplicates$timeString))
    if (! (timeSane && nameSane)) {
        stop("Node names and assigned time values must have a 1:1 relationship.
             E.g. two nodes Age/Week1, Bloodpressure/Week1, must have both the same assigned time (e.g. 1)")
    }
}
