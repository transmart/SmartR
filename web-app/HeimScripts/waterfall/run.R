
main <- function(lowRangeOperator = "&lt;", lowRangeValue = "", highRangeOperator = "&lt;", highRangeValue = "") {

#    save(loaded_variables, file="/tmp/data/loaded_variables.Rda")
#    save(fetch_params, file="/tmp/data/fetch_params.Rda")
#
#    stop(paste("R Script: lowRangeOperator = '", lowRangeOperator, "'  lowRangeValue = '", lowRangeValue, "'  highRangeOperator = '", highRangeOperator, "'  highRangeValue = '", highRangeValue, "'"), sep="")

    if (length(loaded_variables)==1) {
        numData <- loaded_variables$numData_n0_s1
    } else if (length(loaded_variables)==2) {
        numdata <- rbind(loaded_variables$numData_n0_s1, loaded_variables$numData_n0_s2)
    }

    numdata = numdata[order(numdata[2], numdata[1]),]

    colnames(numdata) <- c("patientID", fetch_params$ontologyTerms$numData_n0$fullName)

    output <- list(
    plot_data = numdata,
    xArrLabel = "Patient ID",
    plot_title = fetch_params$ontologyTerms$datapoints_n0$fullName,
    lowRangeOperator = lowRangeOperator,
    lowRangeValue = lowRangeValue,
    highRangeOperator = highRangeOperator,
    highRangeValue = highRangeValue
    )

    toJSON(output)
}

