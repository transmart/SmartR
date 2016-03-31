
getMeasurements <- function(df) {
    if (ncol(df) > 3){
        return(df[, 3:ncol(df)])
    } else {
        return( df[3] )
    }
}

mergeDuplicates <- function(df) {
    dupl.where <- duplicated(df$Row.Label)
    dupl.rows <- df[dupl.where, ]
    df <- df[! dupl.where, ]
    uids <- paste(df$Row.Label, df$Bio.marker, sep="--")
    uids[df$Row.Label == dupl.rows$Row.Label] <- paste(uids[df$Row.Label == dupl.rows$Row.Label], dupl.rows$Bio.marker[df$Row.Label == dupl.rows$Row.Label], sep="--")
    df <- cbind(UID=uids, df[, -c(1,2)])
    df
}

parseInput <- function() {
    if (exists("preprocessed")) {
        df <- preprocessed
    } else {
        df <- mergeFetchedData(loaded_variables)
    }
    return(df)
}

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

# nodeID has usually this format: 'X123_highDimensional_n0_s1)
# this method pretifies it with the actual node label like this: '123_BreastCancer'
replaceNodeIDNodeLabel <- function(ids, ontologyTerms) {
    patientIDs <- sub("_.+_n[0-9]+_s[0-9]+", "", ids, perl=TRUE) # remove the _highDimensional_n0_s1
    patientIDs <- sub("^X", "", patientIDs, perl=TRUE) # remove the X
    nodes <- sub("_s[0-9]+", "", ids, perl=TRUE) # remove the _s1
    nodes <- sub(".+?_", "", nodes, perl=TRUE) # remove the X123_
    # replace highDimensional with Breast
    nodeLabels <- lapply(ontologyTerms[nodes], function(terms) return(terms$name))
    paste(patientIDs, nodeLabels, sep="_")
}

getSubject <- function(patientIDs) {
    splittedIds <- strsplit(patientIDs,"_")
    sapply(splittedIds, FUN = discardNodeAndSubject)
}

discardNodeAndSubject <- function(label) {
    label <- strsplit(label,"_")
    endOfSubject <-
    length(label) - 2  #last too elements are node and subset.
    label <- label[1:endOfSubject]
    paste(label, collapse = "_")
}

fixString <- function(str) {
    gsub("[^a-zA-Z0-9-]", "", str, perl = TRUE)
}

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

dropEmptyGene <- function(d){
    d[!(d$Bio.marker == ""|
        is.null(d$Bio.marker) |
        is.na(d$Bio.marker) |
        is.nan(d$Bio.marker)),]
}

hasTwoSubsets <- function(measurements) {
    getSubset1Length(measurements) < ncol(measurements)
}

getSubset <- function(patientIDs) {
    splittedIds <- strsplit(patientIDs,"_s") # During merge,
    # which is always run we append subset id, either
    # _s1 or _s2 to PATIENTID.
    subsets <- sapply(splittedIds, FUN = tail_elem) # In proper patienid

}

tail_elem <- function(vect, n = 1) {
    as.integer(vect[length(vect) - n + 1])
}
