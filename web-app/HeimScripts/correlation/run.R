library(reshape2)

main <- function(method = "pearson", transformation = "raw", selectedPatientIDs = integer()) {

    df1 <- loaded_variables$datapoints_n0_s1
    df2 <- loaded_variables$datapoints_n1_s1

    if (nrow(df1) == 0) {
        stop(paste("Variable '", fetch_params$ontologyTerms$datapoints_n0$name, "' has no patients for subset 1"), sep="")
    }
    if (nrow(df2) == 0) {
        stop(paste("Variable '", fetch_params$ontologyTerms$datapoints_n1$name, "' has no patients for subset 1"), sep="")
    }
    num_data <- merge(df1, df2, by="Row.Label")
    colnames(num_data) <- c("patientID", "x", "y")

    if (transformation == "log2") {
        num_data$x <- log2(num_data$x)
        num_data$y <- log2(num_data$y)
    } else if (transformation == "log10") {
        num_data$x <- log10(num_data$x)
        num_data$y <- log10(num_data$y)
    }
    num_data <- na.omit(num_data)
    num_data <- num_data[!is.infinite(num_data$x), ]
    num_data <- num_data[!is.infinite(num_data$y), ]

    cat_data <- data.frame(patientID=integer(), annotation=character())
    filtered.loaded_variables <- get.loaded_variables.by.source("annotations", loaded_variables)
    if (length(filtered.loaded_variables) > 0) {
        merged.df <- Reduce(function(...) merge(..., by='Row.Label', all=T), filtered.loaded_variables)
        merged.df <- merged.df[, colSums(is.na(merged.df)) != nrow(merged.df)] # remove NA columns
        
        annotations <- apply(merged.df[,-1], 1, function(row) {
                row <- row[row != ""]
                paste(row, collapse="-AND-")
        })
        cat_data <- data.frame(
           patientID=as.integer(merged.df$Row.Label),
           annotation=as.character(annotations)
        )
    }

    df <- num_data
    if (nrow(cat_data) > 0) {
        df <- merge(df, cat_data, by="patientID")
    } else {
        df$annotation <- ''
    }

    colnames(df) <- c("patientID", "x", "y", "annotation")

    if (length(selectedPatientIDs) > 0) {
        df <- df[df$patientID %in% selectedPatientIDs, ]
    }

    corTest <- tryCatch({
        test <- cor.test(df$x, df$y, method=method)
        test$p.value <- ifelse(test$p.value == 0, paste("<", as.character(.Machine$double.eps)), as.character(test$p.value))
        test
    }, error = function(e) {
        ll <- list()
        ll$estimate <- as.numeric(NA)
        ll$p.value <- as.numeric(NA)
        ll
    })

    regLineSlope <- corTest$estimate * (sd(df$y) / sd(df$x))
    regLineYIntercept <- mean(df$y) - regLineSlope * mean(df$x)

    output <- list(
        correlation = corTest$estimate,
        pvalue = corTest$p.value,
        regLineSlope = regLineSlope,
        regLineYIntercept = regLineYIntercept,
        xArrLabel = fetch_params$ontologyTerms$datapoints_n0$fullName,
        yArrLabel = fetch_params$ontologyTerms$datapoints_n1$fullName,
        method = method,
        transformation = transformation,
        patientIDs = df$patientID,
        annotations = unique(df$annotation),
        points = df
    )
    toJSON(output)
}

conceptStrToFolderStr <- function(s) {
    splitString <- strsplit(s, "")[[1]]
    backslashs <- which(splitString == "\\")
    substr(s, 0, tail(backslashs, 2)[1])
}

