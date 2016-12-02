library(jsonlite)
library(plyr)
library(reshape2)

main <- function(excludedPatientIDs = integer(), transformation="raw") {
    output <- list()
    output$transformation <- transformation

    df <- buildCrossfilterCompatibleDf(loaded_variables, fetch_params)
    if (transformation == "log2") {
        df$value <- log2(df$value)
    } else if (transformation == "log10") {
        df$value <- log10(df$value)
    }

    df <- df[!is.infinite(df$value), ]
    if (length(unique(df$bioMarker)) >= 2) {
        ano <- anova(lm(value ~ bioMarker, data=df))
        p <- ano[1,5]
        p <- ifelse(p < .Machine$double.eps, paste("<", as.character(.Machine$double.eps)), p)
        output$pValue <- p
    }
    output$dataMatrix <- df

    toJSON(output)
}

# returns df that is compatible with crossfilter.js
buildCrossfilterCompatibleDf <- function(loaded_variables, fetch_params) {
    # gather information
    subsets <- getSubsets(loaded_variables)
    names <- getNames(loaded_variables, fetch_params)
    types <- getTypes(loaded_variables)

    # initialize empty df
    df <- data.frame(patientID=integer(),
                     value=integer(),
                     bioMarker=character(),
                     type=character(),
                     subset=integer(),
                     stringsAsFactors=FALSE)

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
                bioMarker <- paste(names[i], rowLabel, variable.label.df$Bio.marker[1], sep="--")
                bioMarker <- paste(bioMarker, " s", subsets[i], sep="")
                variable.label.df <- variable.label.df[, -c(1,2)]
                variable.label.df <- data.frame(patientID=as.numeric(as.vector(variable.label.df$patientID)),
                                                value=as.numeric(as.vector(variable.label.df$value)),
                                                bioMarker=rep(bioMarker, nrow(variable.label.df)),
                                                type=rep('numeric', nrow(variable.label.df)),
                                                subset=rep(subsets[i], nrow(variable.label.df)),
                                                stringsAsFactors=FALSE)
                # no value -> no interest
                variable.label.df <- variable.label.df[variable.label.df$value != "" & !is.na(variable.label.df$value), ]
                if (nrow(variable.label.df) == 0) next
                df <- rbind(df, variable.label.df)
            }
        } else  if (types[i] == "numeric"){
            variable.df <- data.frame(patientID=as.integer(variable[,1]), value=variable[,2])
            bioMarker <- names[i]
            bioMarker <- paste(bioMarker, " s", subsets[i], sep="")
            values <- as.numeric(as.vector(variable.df$value))
            # attach additional information
            variable.df <- data.frame(patientID=as.numeric(as.vector(variable.df$patientID)),
                                 value=values,
                                 bioMarker=rep(bioMarker, nrow(variable.df)),
                                 type=rep(types[i], nrow(variable.df)),
                                 subset=rep(subsets[i], nrow(variable.df)),
                                 stringsAsFactors=FALSE)
            # no value -> no interest
            variable.df <- variable.df[variable.df$value != "" & !is.na(variable.df$value), ]
            if (nrow(variable.df) == 0) next
            df <- rbind(df, variable.df)
        }
    }

    # before we are done we assign a unique id to every row to make it easier for the front-end
    df <- cbind(id=1:nrow(df), df)

    groups <- loaded_variables[grep("^groups", names(loaded_variables))]
    if (length(groups) > 0 ) {
        subsets <- getSubsets(groups)
        for (i in 1:length(groups)) {
            group <- groups[[i]]
            subset <- subsets[i]
            na.omit(group)
            group <- group[group[,2] != "",]
            groupName <- group[1,2]
            patients <- group[, "Row.Label"]
            if (! any(patients %in% df$patientID)) next
            df[df$patientID %in% patients & df$subset == subset, ]$bioMarker <-
                paste(df[df$patientID %in% patients & df$subset == subset, ]$bioMarker, groupName, sep=" g:")
        }
    }
    df
}

# returns character vector (e.g. c("\\Demo Study\\Vital Status\\Alive\\Week1", "\\Demo Study\\Vital Status\\Alive\\Week2", ...))
getNames <- function(loaded_variables, fetch_params) {
    names.without.subset <- sub("_s[1-2]{1}$", "", names(loaded_variables))
    names <- sapply(names.without.subset, function(el) fetch_params$ontologyTerms[[el]]$name)
    as.character(as.vector(names))
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
