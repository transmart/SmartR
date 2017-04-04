
main <- function(selLowRange = "&lt;", txtLowRange = "", selHighRange = "&lt;", txtHighRange = "") {

    print("BUILDING WATERFALL DATA")
    if (selLowRange == "&lt;" & txtLowRange == "" & selHighRange == "&lt;" & txtHighRange == "")
    {
        print("Dummy print")
    }

    num_data <- loaded_variables$numData_n0_s1

    #Read the input data.
    if (nrow(num_data) == 0) {
        stop(paste("Variable '", fetch_params$ontologyTerms$num_data_n0$name, "' has no patients for subset 1"), sep="")
    }

    require(MASS)
    #Write the final data file.
    # write.matrix(finalData,"outputfile.txt",sep = "\t")
    # Using write.table; write.matrix was leaving trailing white-space in the file - see JIRA issue TRANSREL-24.
    write.table(num_data,"~/num_data.dat", sep = "\t", quote = FALSE, row.names = FALSE)
}

