main <- function() {
    output <- list()
    output$sourceIDs <- c()
    output$sourceSubsets <- c()
    output$targetIDs <- c()
    variables <- names(loaded_variables)
    source.vars <- names(loaded_variables)[grep("source", variables)]
    target.vars <- names(loaded_variables)[grep("target", variables)]
    target.vars <- sub("_s[1-2]{1}$", "", target.vars)
    target.vars <- unique(target.vars)
    for (var in source.vars) {
        data <- loaded_variables[var][[1]]
        data <- na.omit(data)
        data <- data[data[, 2] != "", ]
        if (nrow(data)) {
            source.id <- data[1,2]
            sourceSubset <- sub("^.*_s", "", var)
            output$sourceIDs <- c(output$sourceIDs, source.id)
            output$sourceSubsets <- c(output$sourceSubsets, sourceSubset)
        }
    }
    cohortNodes <- data.frame(key=c(),
                              level=c(),
                              fullName=c(),
                              name=c(),
                              tooltip=c(),
                              visualAttributes=c(),
                              subset=c())
    for (var in target.vars) {
        params <- fetch_params$ontologyTerms[var][[1]]
        subsets <- output$sourceSubsets[which(params$name == output$sourceIDs)]
        for (subset in subsets) {
            cohortNodes <- rbind(cohortNodes, data.frame(key=params$key,
                                                         level=params$level,
                                                         fullName=params$fullName,
                                                         name=params$name,
                                                         tooltip=params$tooltip,
                                                         visualAttributes=paste(params$visualAttributes, collapse=","),
                                                         subset=strtoi(subset)))
        }
    }
    output$cohortNodes <- cohortNodes
    toJSON(output)
}
