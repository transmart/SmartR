
main <- function(lowRangeOperator = "&lt;", lowRangeValue = "", highRangeOperator = "&lt;", highRangeValue = "") {

    save(loaded_variables, file="/tmp/data/loaded_variables.Rda")
    save(fetch_params, file="/tmp/data/fetch_params.Rda")

    stop(paste("R Script: lowRangeOperator = '", lowRangeOperator, "'  lowRangeValue = '", lowRangeValue, "'  highRangeOperator = '", highRangeOperator, "'  highRangeValue = '", highRangeValue, "'"), sep="")

}

