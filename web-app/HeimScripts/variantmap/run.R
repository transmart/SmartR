library(reshape2)

main <- function(variants) {
    output <- list()

    save(variants, file="/tmp/variants.Rda")
    save(loaded_variables, file="/tmp/loaded_variables.Rda")
    save(fetch_params, file="/tmp/fetch_params.Rda")

    if (exists("loaded_variables")) {
        hdd.idx <- grep("^highDimensional", names(loaded_variables))
        cat.idx <- grep("^categoric", names(loaded_variables))
        num.idx <- grep("^numeric", names(loaded_variables))

        hdd.data <- loaded_variables[hdd.idx]
        hdd.data <- lapply(hdd.data, function(x) {
             names(x) <- sub("^X", "", names(x))
             x <- melt(x, id.vars=c("Row.Label", "Bio.marker"))
             x
        })
        hdd.data <- do.call(rbind, hdd.data)
        names(hdd.data) <- c("Row.Label", "gene", "subject", "expr")
        data <- merge(variants, hdd.data[, -1], by=c("subject", "gene"), all.x=T)
        output$data <- data
    } else {
        output$data <- variants
    }

    json <- toJSON(output)
    return(json)
}
