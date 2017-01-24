main <- function(variants) {
    output <- list()
    output$variantsMatrix <- variants
    save(variants, file="/tmp/variants.Rda")
    save(loaded_variables, file="/tmp/loaded_variables.Rda")
    save(fetch_params, file="/tmp/fetch_params.Rda")

    json <- toJSON(output)
    return(json)
}
