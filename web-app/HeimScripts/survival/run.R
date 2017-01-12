# Load libraries
library(survival)

# Constants for TimeConversions
dayMonthRelation <- 30
dayYearRelation <- 365
monthYearRelation <- 12

# Global variables
dummy_category <- "NO_CATEGORY_SELECTED"
censor_value <- "1"
censor_interpretation <- "negative"
selected_subsets <- c()
selected_subsets_count <- 0
selected_categories <- c()
selected_categories_count <- 0
subset_1 <- data.frame()
subset_2 <- data.frame()
data_sets <- list()
survival_data_sets <- list()
cox_regression_data_sets <- list()
cox_regression_result_1 <- data.frame()
cox_regression_result_2 <- data.frame()

# MAIN function which is called by SmartR automatically
main <- function(plotWidth = "800", plotHeight = "500", timeIn = "days", timeOut = "days", legendType = "inner", showRiskTable = TRUE, legendPosition = "bottom", mergeSubsets = FALSE, mergeCategories = FALSE, censorInterpretation = "negative") {
	
	# Fill global variables
	censor_interpretation <<- censorInterpretation
	if(censor_interpretation == "negative") {
		censor_value <<- "1"
	} else {
		censor_value <<- "0"
	}
	getSelectedSubsets()
	getSelectedCategories()
	
	# Make local settings
	x_label = fetch_params$ontologyTerms$time_n0$fullName
	if(mergeSubsets == 'TRUE') {
		mergeSubsets = TRUE
	} else {
		mergeSubsets = FALSE
	}
	if(mergeCategories == 'TRUE') {
		mergeCategories = TRUE
	} else {
		mergeCategories = FALSE
	}
	if(showRiskTable == 'TRUE') {
		showRiskTable = TRUE
	} else {
		showRiskTable = FALSE
	}
	
	# Subset 1
	getSubsetData(1)
	#log1 <- file("subset_1.log")
	#write.table(subset_1, log1, sep="\t", na="", row.names=FALSE)
	
	# Subset 2
	getSubsetData(2)
	#log2 <- file("subset_2.log")
	#write.table(subset_2, log2, sep="\t", na="", row.names=FALSE)
	
	# Convert TimeMeasurement of the previously defined subsets
	convertDataSetsTimeMeasurement(timeIn, timeOut)
	
	# Generate the different datasets that should be analyised according to settings
	generateDataSets(mergeSubsets, mergeCategories)
	
	# Make survival analysis for all the data_sets
	generateSurvivalData()
	
	# Make Cox Regression Data Sets
	generateCoxRegressionData(mergeSubsets, mergeCategories)
	
	# Make Cox Regression
	generateCoxRegressionResults(mergeSubsets, mergeCategories)
	
	# Generate output for visualization
	output <- 	list(
					survival_data = survival_data_sets,
					statistics1 = cox_regression_result_1,
					statistics2 = cox_regression_result_2,
					selected_subsets = selected_subsets,
					selected_categories = selected_categories,
					merge_subsets = mergeSubsets,
					merge_categories = mergeCategories,
					plot_width = plotWidth,
					plot_height = plotHeight,
					time_out = timeOut,
					legend_type = legendType,
					legend_position = legendPosition,
					show_risk_table = showRiskTable,
					x_label = x_label
				)
	
	toJSON(output)
	
}

# Function that extracts the name of the selected Subsets
# Currently not really working
# todo: get selected concepts for generation of subsets
getSelectedSubsets <- function() {
	if(!is.null(loaded_variables$time_n0_s1)) {
		selected_subsets <<- c("Subset 1")
		selected_subsets_count <<- 1
	}
	if(!is.null(loaded_variables$time_n0_s2)) {
		selected_subsets <<- c(selected_subsets, "Subset 2")
		selected_subsets_count <<- selected_subsets_count + 1
	}
}

# Function that extracts the name of the selected Categories
getSelectedCategories <- function() {
	
	if(!is.null(fetch_params$ontologyTerms$category_n0)) {
		selected_categories <<- fetch_params$ontologyTerms$category_n0$fullName
		selected_categories_count <<- 1
	}
	
	if(!is.null(fetch_params$ontologyTerms$category_n1)) {
		selected_categories <<- c(selected_categories, fetch_params$ontologyTerms$category_n1$fullName)
		selected_categories_count <<- selected_categories_count + 1
	}
	
	if(!is.null(fetch_params$ontologyTerms$category_n2)) {
		selected_categories <<- c(selected_categories, fetch_params$ontologyTerms$category_n2$fullName)
		selected_categories_count <<- selected_categories_count + 1
	}
	
	if(!is.null(fetch_params$ontologyTerms$category_n3)) {
		selected_categories <<- c(selected_categories, fetch_params$ontologyTerms$category_n3$fullName)
		selected_categories_count <<- selected_categories_count + 1
	}
	
}

# Function that gets the data for each subset and returns it as data.frame
getSubsetData <- function(subsetNumber) {
	
	# Subset
	subset <- data.frame()
	
	if(subsetNumber == 1 | subsetNumber == 2) {
	
		# Time (necessary) (max=1)
		if(subsetNumber == 1) {
			time <- loaded_variables$time_n0_s1
			if(nrow(time) == 0) {
				stop(paste("Variable '", fetch_params$ontologyTerms$time_n0$name, "' has no patients for subset ", subsetNumber), sep="")
			}
		} else {
			if(!is.null(loaded_variables$time_n0_s2)) {
				time <- loaded_variables$time_n0_s2
				if(nrow(time) == 0) {
					stop(paste("Variable '", fetch_params$ontologyTerms$time_n0$name, "' has no patients for subset ", subsetNumber), sep="")
				}
			} else {
				time <- data.frame()
			}
		}
		
		if(nrow(time)!=0) {
			
			subset <- time
			
			# Censoring (optional) (max=1)
			if(subsetNumber == 1) {
				if(!is.null(loaded_variables$censoring_n0_s1)) {
					censoring <- loaded_variables$censoring_n0_s1
					if(censor_interpretation == "negative") {
						censoring[censoring == ""] <- 0
					} else {
						censoring[censoring == ""] <- 1
					}
					if(nrow(censoring) == 0) {
						stop(paste("Variable '", fetch_params$ontologyTerms$censoring_n0$name, "' has no patients for subset ", subsetNumber), sep="")
					} else {
						subset <- merge(subset, censoring, by="Row.Label")
					}
				} else {
					#If no event was selected, we consider everyone to have had the event.
					subset <- cbind(subset, 1)
				}
				colnames(subset)[ncol(subset)] <- "CENSOR"
				if(censor_interpretation == "negative") {
					subset$CENSOR[!(subset$CENSOR=="0")] <- 1
				} else {
					subset$CENSOR[!(subset$CENSOR=="1")] <- 0
				}
				subset$CENSOR <- as.numeric(subset$CENSOR)
			} else {
				if(!is.null(loaded_variables$censoring_n0_s2)) {
					censoring <- loaded_variables$censoring_n0_s2
					if(censor_interpretation == "negative") {
						censoring[censoring == ""] <- 0
					} else {
						censoring[censoring == ""] <- 1
					}
					if(nrow(censoring) == 0) {
						stop(paste("Variable '", fetch_params$ontologyTerms$censoring_n0$name, "' has no patients for subset ", subsetNumber), sep="")
					} else {
						subset <- merge(subset, censoring, by="Row.Label")
					}
				} else {
					#If no event was selected, we consider everyone to have had the event.
					subset <- cbind(subset, 1)
				}
				colnames(subset)[3] <- "CENSOR"
				if(censor_interpretation == "negative") {
					subset$CENSOR[!(subset$CENSOR=="0")] <- 1
				} else {
					subset$CENSOR[!(subset$CENSOR=="1")] <- 0
				}
				subset$CENSOR <- as.numeric(subset$CENSOR)
			}
			
			colnames(subset)[2] <- "TIME"
			
		}
		
	}
	
	if(subsetNumber == 1) {
		subset_1 <<- subset
	}
	
	if(subsetNumber == 2) {
		subset_2 <<- subset
	}
	
}

generateDataSets <- function(mergeSubsets, mergeCategories) {
	
	this_data_set <- data.frame()
	
	if(selected_categories_count == 0) {
		# Subset 1
		if(nrow(subset_1) > 0) {
			this_data_set <- cbind(subset_1, 1)
			colnames(this_data_set)[4] <- "CATEGORY"
			this_data_set <- this_data_set[,c(1,2,4,3)]
			data_sets[[length(data_sets)+1]] <<- this_data_set
		}
		# Subset 2
		if(nrow(subset_2) > 0) {
			this_data_set <- cbind(subset_2, 1)
			colnames(this_data_set)[4] <- "CATEGORY"
			this_data_set <- this_data_set[,c(1,2,4,3)]
			data_sets[[length(data_sets)+1]] <<- this_data_set
		}
	}
	
	if(mergeCategories && selected_categories_count > 0) {
		# Subset 1
		if(nrow(subset_1) > 0) {
			subset_1_data_set <- subset_1
			# Category 1
			if(!is.null(loaded_variables$category_n0_s1)) {
				category <- loaded_variables$category_n0_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[4] <- "CATEGORY_1"
				subset_1_data_set <- subset_1_data_set[,c(1,2,4,3)]
			}
			# Category 2
			if(!is.null(loaded_variables$category_n1_s1)) {
				category <- loaded_variables$category_n1_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[5] <- "CATEGORY_2"
				subset_1_data_set <- subset_1_data_set[,c(1,2,3,5,4)]
			}
			# Category 3
			if(!is.null(loaded_variables$category_n2_s1)) {
				category <- loaded_variables$category_n2_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[6] <- "CATEGORY_3"
				subset_1_data_set <- subset_1_data_set[,c(1,2,3,4,6,5)]
			}
			# Category 4
			if(!is.null(loaded_variables$category_n3_s1)) {
				category <- loaded_variables$category_n3_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[7] <- "CATEGORY_4"
				subset_1_data_set <- subset_1_data_set[,c(1,2,3,4,5,7,6)]
			}
			if(selected_categories_count == 1) {
				subset_1_data_set <- subset_1_data_set[!(is.na(subset_1_data_set$CATEGORY_1) | subset_1_data_set$CATEGORY_1==""),]
			}
			if(selected_categories_count == 2) {
				subset_1_data_set <- subset_1_data_set[!(is.na(subset_1_data_set$CATEGORY_1) | subset_1_data_set$CATEGORY_1=="") | !(is.na(subset_1_data_set$CATEGORY_2) | subset_1_data_set$CATEGORY_2==""),]
			}
			if(selected_categories_count == 3) {
				subset_1_data_set <- subset_1_data_set[!(is.na(subset_1_data_set$CATEGORY_1) | subset_1_data_set$CATEGORY_1=="") | !(is.na(subset_1_data_set$CATEGORY_2) | subset_1_data_set$CATEGORY_2=="") | !(is.na(subset_1_data_set$CATEGORY_3) | subset_1_data_set$CATEGORY_3==""),]
			}
			if(selected_categories_count == 4) {
				subset_1_data_set <- subset_1_data_set[!(is.na(subset_1_data_set$CATEGORY_1) | subset_1_data_set$CATEGORY_1=="") | !(is.na(subset_1_data_set$CATEGORY_2) | subset_1_data_set$CATEGORY_2=="") | !(is.na(subset_1_data_set$CATEGORY_3) | subset_1_data_set$CATEGORY_3=="") | !(is.na(subset_1_data_set$CATEGORY_4) | subset_1_data_set$CATEGORY_4==""),]
			}
			subset_1_data_set <- data.frame(Row.Label=subset_1_data_set$Row.Label, TIME=subset_1_data_set$TIME, CENSOR=subset_1_data_set$CENSOR)
			subset_1_data_set <- cbind(subset_1_data_set, 1)
			colnames(subset_1_data_set)[4] <- "CATEGORY"
			subset_1_data_set <- subset_1_data_set[,c(1,2,4,3)]
			data_sets[[length(data_sets)+1]] <<- subset_1_data_set
		}
		# Subset 2
		if(nrow(subset_2) > 0) {
			subset_2_data_set <- subset_2
			# Category 1
			if(!is.null(loaded_variables$category_n0_s2)) {
				category <- loaded_variables$category_n0_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[4] <- "CATEGORY_1"
				subset_2_data_set <- subset_2_data_set[,c(1,2,4,3)]
			}
			# Category 2
			if(!is.null(loaded_variables$category_n1_s2)) {
				category <- loaded_variables$category_n1_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[5] <- "CATEGORY_2"
				subset_2_data_set <- subset_2_data_set[,c(1,2,3,5,4)]
			}
			# Category 3
			if(!is.null(loaded_variables$category_n2_s2)) {
				category <- loaded_variables$category_n2_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[6] <- "CATEGORY_3"
				subset_2_data_set <- subset_2_data_set[,c(1,2,3,4,6,5)]
			}
			# Category 4
			if(!is.null(loaded_variables$category_n3_s2)) {
				category <- loaded_variables$category_n3_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[7] <- "CATEGORY_4"
				subset_2_data_set <- subset_2_data_set[,c(1,2,3,4,5,7,6)]
			}
			if(selected_categories_count == 1) {
				subset_2_data_set <- subset_2_data_set[!(is.na(subset_2_data_set$CATEGORY_1) | subset_2_data_set$CATEGORY_1==""),]
			}
			if(selected_categories_count == 2) {
				subset_2_data_set <- subset_2_data_set[!(is.na(subset_2_data_set$CATEGORY_1) | subset_2_data_set$CATEGORY_1=="") | !(is.na(subset_2_data_set$CATEGORY_2) | subset_2_data_set$CATEGORY_2==""),]
			}
			if(selected_categories_count == 3) {
				subset_2_data_set <- subset_2_data_set[!(is.na(subset_2_data_set$CATEGORY_1) | subset_2_data_set$CATEGORY_1=="") | !(is.na(subset_2_data_set$CATEGORY_2) | subset_2_data_set$CATEGORY_2=="") | !(is.na(subset_2_data_set$CATEGORY_3) | subset_2_data_set$CATEGORY_3==""),]
			}
			if(selected_categories_count == 4) {
				subset_2_data_set <- subset_2_data_set[!(is.na(subset_2_data_set$CATEGORY_1) | subset_2_data_set$CATEGORY_1=="") | !(is.na(subset_2_data_set$CATEGORY_2) | subset_2_data_set$CATEGORY_2=="") | !(is.na(subset_2_data_set$CATEGORY_3) | subset_2_data_set$CATEGORY_3=="") | !(is.na(subset_2_data_set$CATEGORY_4) | subset_2_data_set$CATEGORY_4==""),]
			}
			subset_2_data_set <- data.frame(Row.Label=subset_2_data_set$Row.Label, TIME=subset_2_data_set$TIME, CENSOR=subset_2_data_set$CENSOR)
			subset_2_data_set <- cbind(subset_2_data_set, 1)
			colnames(subset_2_data_set)[4] <- "CATEGORY"
			subset_2_data_set <- subset_2_data_set[,c(1,2,4,3)]
			data_sets[[length(data_sets)+1]] <<- subset_2_data_set
		}
	}
	
	# make data_sets combining subset_data with selected_categories
	if(!mergeCategories && selected_categories_count > 0) {
		# Subset 1 && Category 1
		if(nrow(subset_1) > 0 && !is.null(loaded_variables$category_n0_s1)) {
			category <- loaded_variables$category_n0_s1
			this_data_set <- merge(subset_1, category, by="Row.Label")
			colnames(this_data_set)[4] <- "CATEGORY"
			this_data_set <- this_data_set[,c(1,2,4,3)]
			this_data_set <- this_data_set[!(is.na(this_data_set$CATEGORY) | this_data_set$CATEGORY==""),]
			data_sets[[length(data_sets)+1]] <<- this_data_set
		}
		# Subset 1 && Category 2
		if(nrow(subset_1) > 0 && !is.null(loaded_variables$category_n1_s1)) {
			category <- loaded_variables$category_n1_s1
			this_data_set <- merge(subset_1, category, by="Row.Label")
			colnames(this_data_set)[4] <- "CATEGORY"
			this_data_set <- this_data_set[,c(1,2,4,3)]
			this_data_set <- this_data_set[!(is.na(this_data_set$CATEGORY) | this_data_set$CATEGORY==""),]
			data_sets[[length(data_sets)+1]] <<- this_data_set
		}
		# Subset 1 && Category 3
		if(nrow(subset_1) > 0 && !is.null(loaded_variables$category_n2_s1)) {
			category <- loaded_variables$category_n2_s1
			this_data_set <- merge(subset_1, category, by="Row.Label")
			colnames(this_data_set)[4] <- "CATEGORY"
			this_data_set <- this_data_set[,c(1,2,4,3)]
			this_data_set <- this_data_set[!(is.na(this_data_set$CATEGORY) | this_data_set$CATEGORY==""),]
			data_sets[[length(data_sets)+1]] <<- this_data_set
		}
		# Subset 1 && Category 4
		if(nrow(subset_1) > 0 && !is.null(loaded_variables$category_n3_s1)) {
			category <- loaded_variables$category_n3_s1
			this_data_set <- merge(subset_1, category, by="Row.Label")
			colnames(this_data_set)[4] <- "CATEGORY"
			this_data_set <- this_data_set[,c(1,2,4,3)]
			this_data_set <- this_data_set[!(is.na(this_data_set$CATEGORY) | this_data_set$CATEGORY==""),]
			data_sets[[length(data_sets)+1]] <<- this_data_set
		}
		# Subset 2 && Category 1
		if(nrow(subset_2) > 0 && !is.null(loaded_variables$category_n0_s2)) {
			category <- loaded_variables$category_n0_s2
			this_data_set <- merge(subset_2, category, by="Row.Label")
			colnames(this_data_set)[4] <- "CATEGORY"
			this_data_set <- this_data_set[,c(1,2,4,3)]
			this_data_set <- this_data_set[!(is.na(this_data_set$CATEGORY) | this_data_set$CATEGORY==""),]
			data_sets[[length(data_sets)+1]] <<- this_data_set
		}
		# Subset 2 && Category 2
		if(nrow(subset_2) > 0 && !is.null(loaded_variables$category_n1_s2)) {
			category <- loaded_variables$category_n1_s2
			this_data_set <- merge(subset_2, category, by="Row.Label")
			colnames(this_data_set)[4] <- "CATEGORY"
			this_data_set <- this_data_set[,c(1,2,4,3)]
			this_data_set <- this_data_set[!(is.na(this_data_set$CATEGORY) | this_data_set$CATEGORY==""),]
			data_sets[[length(data_sets)+1]] <<- this_data_set
		}
		# Subset 2 && Category 3
		if(nrow(subset_2) > 0 && !is.null(loaded_variables$category_n2_s2)) {
			category <- loaded_variables$category_n2_s2
			this_data_set <- merge(subset_2, category, by="Row.Label")
			colnames(this_data_set)[4] <- "CATEGORY"
			this_data_set <- this_data_set[,c(1,2,4,3)]
			this_data_set <- this_data_set[!(is.na(this_data_set$CATEGORY) | this_data_set$CATEGORY==""),]
			data_sets[[length(data_sets)+1]] <<- this_data_set
		}
		# Subset 2 && Category 4
		if(nrow(subset_2) > 0 && !is.null(loaded_variables$category_n3_s2)) {
			category <- loaded_variables$category_n3_s2
			this_data_set <- merge(subset_2, category, by="Row.Label")
			colnames(this_data_set)[4] <- "CATEGORY"
			this_data_set <- this_data_set[,c(1,2,4,3)]
			this_data_set <- this_data_set[!(is.na(this_data_set$CATEGORY) | this_data_set$CATEGORY==""),]
			data_sets[[length(data_sets)+1]] <<- this_data_set
		}
	}
	
	if(mergeSubsets && selected_subsets_count == 2 && selected_categories_count == 0) {
		data_set_subset_1 <- as.data.frame(data_sets[1])
		data_set_subset_2 <- as.data.frame(data_sets[2])
		merge_data_set <- rbind(data_set_subset_1, data_set_subset_2)
		data_sets <<- list()
		data_sets[[length(data_sets)+1]] <<- merge_data_set
	}
	
	if(mergeSubsets && selected_subsets_count == 2 && selected_categories_count > 0) {
		merged_data_sets <- list()
		for(i in 1:(length(data_sets)/2)) {
			data_set_subset_1 <- as.data.frame(data_sets[i])
			data_set_subset_2 <- as.data.frame(data_sets[i + (length(data_sets) / 2)])
			merge_data_set <- rbind(data_set_subset_1, data_set_subset_2)
			merged_data_sets[[length(merged_data_sets)+1]] <- merge_data_set
		}
		data_sets <<- merged_data_sets
	}

}

generateSurvivalData <- function() {
	
	for(data_set in data_sets) {
		if(length(data_set$TIME)>0) {
			survival_fit_data = survfit(Surv(data_set$TIME, data_set$CENSOR)~1)
			survival_data <- data.frame(survival_fit_data$time, survival_fit_data$n.risk, survival_fit_data$n.event)
			colnames(survival_data) <- c("t", "n", "d")
		} else {
			survival_data <- data.frame(t=0, n=0, d=0)
		}
		survival_data_sets[[length(survival_data_sets)+1]] <<- survival_data
	}
	
}

generateCoxRegressionData <- function(mergeSubsets, mergeCategories) {
	
	this_data_set <- data.frame()
	
	if(selected_categories_count == 0) {
		# Subset 1
		if(nrow(subset_1) > 0) {
			this_data_set <- cbind(subset_1, 1)
			colnames(this_data_set)[4] <- "CATEGORY"
			cox_regression_data_sets[[length(cox_regression_data_sets)+1]] <<- this_data_set
		}
		# Subset 2
		if(nrow(subset_2) > 0) {
			this_data_set <- cbind(subset_2, 1)
			colnames(this_data_set)[4] <- "CATEGORY"
			cox_regression_data_sets[[length(cox_regression_data_sets)+1]] <<- this_data_set
		}
	}
	
	if(mergeCategories && selected_categories_count > 0) {
		# Subset 1
		if(nrow(subset_1) > 0) {
			subset_1_data_set <- subset_1
			# Category 1
			if(!is.null(loaded_variables$category_n0_s1)) {
				category <- loaded_variables$category_n0_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[4] <- "CATEGORY_1"
			}
			# Category 2
			if(!is.null(loaded_variables$category_n1_s1)) {
				category <- loaded_variables$category_n1_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[5] <- "CATEGORY_2"
			}
			# Category 3
			if(!is.null(loaded_variables$category_n2_s1)) {
				category <- loaded_variables$category_n2_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[6] <- "CATEGORY_3"
			}
			# Category 4
			if(!is.null(loaded_variables$category_n3_s1)) {
				category <- loaded_variables$category_n3_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[7] <- "CATEGORY_4"
			}
			subset_1_data_set <- data.frame(Row.Label=subset_1_data_set$Row.Label, TIME=subset_1_data_set$TIME, CENSOR=subset_1_data_set$CENSOR)
			subset_1_data_set <- cbind(subset_1_data_set, 1)
			colnames(subset_1_data_set)[4] <- "CATEGORY"
			cox_regression_data_sets[[length(cox_regression_data_sets)+1]] <<- subset_1_data_set
		}
		# Subset 2
		if(nrow(subset_2) > 0) {
			subset_2_data_set <- subset_2
			# Category 1
			if(!is.null(loaded_variables$category_n0_s2)) {
				category <- loaded_variables$category_n0_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[4] <- "CATEGORY_1"
			}
			# Category 2
			if(!is.null(loaded_variables$category_n1_s2)) {
				category <- loaded_variables$category_n1_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[5] <- "CATEGORY_2"
			}
			# Category 3
			if(!is.null(loaded_variables$category_n2_s2)) {
				category <- loaded_variables$category_n2_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[6] <- "CATEGORY_3"
			}
			# Category 4
			if(!is.null(loaded_variables$category_n3_s2)) {
				category <- loaded_variables$category_n3_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[7] <- "CATEGORY_4"
			}
			if(selected_categories_count == 1) {
				subset_2_data_set <- subset_2_data_set[!(is.na(subset_2_data_set$CATEGORY_1) | subset_2_data_set$CATEGORY_1==""),]
			}
			if(selected_categories_count == 2) {
				subset_2_data_set <- subset_2_data_set[!(is.na(subset_2_data_set$CATEGORY_1) | subset_2_data_set$CATEGORY_1=="") | !(is.na(subset_2_data_set$CATEGORY_2) | subset_2_data_set$CATEGORY_2==""),]
			}
			if(selected_categories_count == 3) {
				subset_2_data_set <- subset_2_data_set[!(is.na(subset_2_data_set$CATEGORY_1) | subset_2_data_set$CATEGORY_1=="") | !(is.na(subset_2_data_set$CATEGORY_2) | subset_2_data_set$CATEGORY_2=="") | !(is.na(subset_2_data_set$CATEGORY_3) | subset_2_data_set$CATEGORY_3==""),]
			}
			if(selected_categories_count == 4) {
				subset_2_data_set <- subset_2_data_set[!(is.na(subset_2_data_set$CATEGORY_1) | subset_2_data_set$CATEGORY_1=="") | !(is.na(subset_2_data_set$CATEGORY_2) | subset_2_data_set$CATEGORY_2=="") | !(is.na(subset_2_data_set$CATEGORY_3) | subset_2_data_set$CATEGORY_3=="") | !(is.na(subset_2_data_set$CATEGORY_4) | subset_2_data_set$CATEGORY_4==""),]
			}
			subset_2_data_set <- data.frame(Row.Label=subset_2_data_set$Row.Label, TIME=subset_2_data_set$TIME, CENSOR=subset_2_data_set$CENSOR)
			subset_2_data_set <- cbind(subset_2_data_set, 1)
			colnames(subset_2_data_set)[4] <- "CATEGORY"
			cox_regression_data_sets[[length(cox_regression_data_sets)+1]] <<- subset_2_data_set
		}
	}
	
	# make data_sets combining subset_data with selected_categories
	if(!mergeCategories && selected_categories_count > 0) {
		# Subset 1
		if(nrow(subset_1) > 0) {
			subset_1_data_set <- subset_1
			# Subset 1 && Category 1
			if(!is.null(loaded_variables$category_n0_s1)) {
				category <- loaded_variables$category_n0_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[4] <- "CATEGORY_1"
				subset_1_data_set$CATEGORY_1[is.na(subset_1_data_set$CATEGORY_1)] <- ""
			}
			# Subset 1 && Category 2
			if(!is.null(loaded_variables$category_n1_s1)) {
				category <- loaded_variables$category_n1_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[5] <- "CATEGORY_2"
				subset_1_data_set$CATEGORY_2[is.na(subset_1_data_set$CATEGORY_2)] <- ""
			}
			# Subset 1 && Category 3
			if(!is.null(loaded_variables$category_n2_s1)) {
				category <- loaded_variables$category_n2_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[6] <- "CATEGORY_3"
				subset_1_data_set$CATEGORY_3[is.na(subset_1_data_set$CATEGORY_3)] <- ""
			}
			# Subset 1 && Category 4
			if(!is.null(loaded_variables$category_n3_s1)) {
				category <- loaded_variables$category_n3_s1
				subset_1_data_set <- merge(subset_1_data_set, category, by="Row.Label")
				colnames(subset_1_data_set)[7] <- "CATEGORY_4"
				subset_1_data_set$CATEGORY_4[is.na(subset_1_data_set$CATEGORY_4)] <- ""
			}
			cox_regression_data_sets[[length(cox_regression_data_sets)+1]] <<- subset_1_data_set
		}
		
		# Subset 2
		if(nrow(subset_2) > 0) {
			subset_2_data_set <- subset_2
			# Subset 2 && Category 1
			if(!is.null(loaded_variables$category_n0_s2)) {
				category <- loaded_variables$category_n0_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[4] <- "CATEGORY_1"
				subset_2_data_set$CATEGORY_1[is.na(subset_2_data_set$CATEGORY_1)] <- ""
			}
			# Subset 2 && Category 2
			if(!is.null(loaded_variables$category_n1_s2)) {
				category <- loaded_variables$category_n1_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[5] <- "CATEGORY_2"
				subset_2_data_set$CATEGORY_2[is.na(subset_2_data_set$CATEGORY_2)] <- ""
			}
			# Subset 2 && Category 3
			if(!is.null(loaded_variables$category_n2_s2)) {
				category <- loaded_variables$category_n2_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[6] <- "CATEGORY_3"
				subset_2_data_set$CATEGORY_3[is.na(subset_2_data_set$CATEGORY_3)] <- ""
			}
			# Subset 2 && Category 4
			if(!is.null(loaded_variables$category_n3_s2)) {
				category <- loaded_variables$category_n3_s2
				subset_2_data_set <- merge(subset_2_data_set, category, by="Row.Label")
				colnames(subset_2_data_set)[7] <- "CATEGORY_4"
				subset_2_data_set$CATEGORY_4[is.na(subset_2_data_set$CATEGORY_4)] <- ""
			}
			cox_regression_data_sets[[length(cox_regression_data_sets)+1]] <<- subset_2_data_set
		}
		
	}
	
	# merge Subsets
	if(mergeSubsets && selected_subsets_count == 2 && selected_categories_count == 0) {
		data_set_subset_1 <- as.data.frame(cox_regression_data_sets[1])
		data_set_subset_2 <- as.data.frame(cox_regression_data_sets[2])
		merge_data_set <- rbind(data_set_subset_1, data_set_subset_2)
		cox_regression_data_sets <<- list()
		cox_regression_data_sets[[length(cox_regression_data_sets)+1]] <<- merge_data_set
	}
	
	if(mergeSubsets && selected_subsets_count == 2 && selected_categories_count > 0) {
		merged_data_sets <- list()
		for(i in 1:(length(cox_regression_data_sets)/2)) {
			data_set_subset_1 <- as.data.frame(cox_regression_data_sets[i])
			data_set_subset_2 <- as.data.frame(cox_regression_data_sets[i + (length(cox_regression_data_sets) / 2)])
			merge_data_set <- rbind(data_set_subset_1, data_set_subset_2)
			merged_data_sets[[length(merged_data_sets)+1]] <- merge_data_set
		}
		cox_regression_data_sets <<- merged_data_sets
	}
	
}

generateCoxRegressionResults <- function(mergeSubsets, mergeCategories) {
	
	rownames_result_1 = c("Number of Subjects", "Number of Events", "Concordance", "Rsquared", "Likelihood ratio test", "Wald test", "Score (logrank) test")
	columnnames_result_2 = c("1", "Category", "Cox Coefficient", "Hazards Ratio", "Standard Deviation", "z-Score", "Lower Range of Hazards Ratio", "Upper Range of Hazards Ratio")
	
	cox_regression_result_1 <- data.frame(ROWNAMES=rownames_result_1)
	cox_regression_result_2 <- data.frame(ROWNAMES=columnnames_result_2)
	
	current_subset = 1
	for(cox_regression_data_set in cox_regression_data_sets) {
		
		cox <- data.frame()
		
		if(mergeCategories) {
			cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY, data=cox_regression_data_set, method="efron", robust="F"))
		} else {
			if(selected_categories_count == 0) {
				cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY, data=cox_regression_data_set, method="efron", robust="F"))
			}
			if(selected_categories_count == 1) {
				classList <- as.vector(gsub("\\s","_",gsub("^\\s+|\\s+$", "",(cox_regression_data_set$CATEGORY_1))))
				#1
				if(length(unique(classList)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1, data=cox_regression_data_set, method="efron", robust="F"))
				}
			}
			if(selected_categories_count == 2) {
				classList1 <- as.vector(gsub("\\s","_",gsub("^\\s+|\\s+$", "",(cox_regression_data_set$CATEGORY_1))))
				classList2 <- as.vector(gsub("\\s","_",gsub("^\\s+|\\s+$", "",(cox_regression_data_set$CATEGORY_2))))
				#11
				if(length(unique(classList1)) > 1 && length(unique(classList2)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_2, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#10
				if(length(unique(classList1)) > 1 && length(unique(classList2)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#01
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_2, data=cox_regression_data_set, method="efron", robust="F"))
				}
			}
			if(selected_categories_count == 3) {
				classList1 <- as.vector(gsub("\\s","_",gsub("^\\s+|\\s+$", "",(cox_regression_data_set$CATEGORY_1))))
				classList2 <- as.vector(gsub("\\s","_",gsub("^\\s+|\\s+$", "",(cox_regression_data_set$CATEGORY_2))))
				classList3 <- as.vector(gsub("\\s","_",gsub("^\\s+|\\s+$", "",(cox_regression_data_set$CATEGORY_3))))
				#111
				if(length(unique(classList1)) > 1 && length(unique(classList2)) > 1 && length(unique(classList3)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_2 + CATEGORY_3, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#110
				if(length(unique(classList1)) > 1 && length(unique(classList2)) > 1 && length(unique(classList3)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_2, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#101
				if(length(unique(classList1)) > 1 && length(unique(classList2)) <= 1 && length(unique(classList3)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_3, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#100
				if(length(unique(classList1)) > 1 && length(unique(classList2)) <= 1 && length(unique(classList3)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#011
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) > 1 && length(unique(classList3)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_2 + CATEGORY_3, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#010
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) > 1 && length(unique(classList3)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_2, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#001
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) <= 1 && length(unique(classList3)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_3, data=cox_regression_data_set, method="efron", robust="F"))
				}
			}
			if(selected_categories_count == 4) {
				classList1 <- as.vector(gsub("\\s","_",gsub("^\\s+|\\s+$", "",(cox_regression_data_set$CATEGORY_1))))
				classList2 <- as.vector(gsub("\\s","_",gsub("^\\s+|\\s+$", "",(cox_regression_data_set$CATEGORY_2))))
				classList3 <- as.vector(gsub("\\s","_",gsub("^\\s+|\\s+$", "",(cox_regression_data_set$CATEGORY_3))))
				classList4 <- as.vector(gsub("\\s","_",gsub("^\\s+|\\s+$", "",(cox_regression_data_set$CATEGORY_4))))
				#1111
				if(length(unique(classList1)) > 1 && length(unique(classList2)) > 1 && length(unique(classList3)) > 1 && length(unique(classList4)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_2 + CATEGORY_3 + CATEGORY_4, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#1110
				if(length(unique(classList1)) > 1 && length(unique(classList2)) > 1 && length(unique(classList3)) > 1 && length(unique(classList4)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_2 + CATEGORY_3, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#1101
				if(length(unique(classList1)) > 1 && length(unique(classList2)) > 1 && length(unique(classList3)) <= 1 && length(unique(classList4)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_2 + CATEGORY_4, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#1100
				if(length(unique(classList1)) > 1 && length(unique(classList2)) > 1 && length(unique(classList3)) <= 1 && length(unique(classList4)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_2, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#1011
				if(length(unique(classList1)) > 1 && length(unique(classList2)) <= 1 && length(unique(classList3)) > 1 && length(unique(classList4)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_3 + CATEGORY_4, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#1010
				if(length(unique(classList1)) > 1 && length(unique(classList2)) <= 1 && length(unique(classList3)) > 1 && length(unique(classList4)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_3, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#1001
				if(length(unique(classList1)) > 1 && length(unique(classList2)) <= 1 && length(unique(classList3)) <= 1 && length(unique(classList4)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1 + CATEGORY_4, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#1000
				if(length(unique(classList1)) > 1 && length(unique(classList2)) <= 1 && length(unique(classList3)) <= 1 && length(unique(classList4)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_1, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#0111
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) > 1 && length(unique(classList3)) > 1 && length(unique(classList4)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_2 + CATEGORY_3 + CATEGORY_4, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#0110
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) > 1 && length(unique(classList3)) > 1 && length(unique(classList4)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_2 + CATEGORY_3, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#0101
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) > 1 && length(unique(classList3)) <= 1 && length(unique(classList4)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_2 + CATEGORY_4, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#0100
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) > 1 && length(unique(classList3)) <= 1 && length(unique(classList4)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_2, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#0011
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) <= 1 && length(unique(classList3)) > 1 && length(unique(classList4)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_3 + CATEGORY_4, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#0010
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) <= 1 && length(unique(classList3)) > 1 && length(unique(classList4)) <= 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~  CATEGORY_3, data=cox_regression_data_set, method="efron", robust="F"))
				}
				#0001
				if(length(unique(classList1)) <= 1 && length(unique(classList2)) <= 1 && length(unique(classList3)) <= 1 && length(unique(classList4)) > 1) {
					cox <- summary(coxph(Surv(TIME, CENSOR) ~ CATEGORY_4, data=cox_regression_data_set, method="efron", robust="F"))
				}
			}
		}
		
		# RESULTS FOR TABLE 1
		number_of_subjects = cox$n
		if(is.null(number_of_subjects)) {
			number_of_subjects = 0
		}
		number_of_events = cox$nevent
		if(is.null(number_of_events)) {
			number_of_events = 0
		}
		concordance = cox$concordance
		concordance = gsub(",", "", unlist(strsplit(toString(cox$concordance), "\\s+")))
		concordance = paste(paste(paste(toString(signif(as.numeric(concordance[1]), digits=3)), " (se = ", sep=""), toString(signif(as.numeric(concordance[2]), digits=3)), sep=""), ")", sep="")
		rsq = gsub(",", "", unlist(strsplit(toString(cox$rsq), "\\s+")))
		rsq = paste(paste(paste(toString(signif(as.numeric(rsq[1]), digits=3)), " (max possible = ", sep=""), toString(signif(as.numeric(concordance[2]), digits=3)), sep=""), ")", sep="")
		likelihood = gsub(",", "", unlist(strsplit(toString(cox$logtest), "\\s+")))
		likelihood = paste(paste(paste(paste(toString(signif(as.numeric(likelihood[1]), digits=4)), " on ", sep=""), toString(signif(as.numeric(likelihood[2]), digits=4)), sep=""), " df, p=", sep=""), toString(signif(as.numeric(likelihood[3]), digits=4)), sep="")
		waldtest = gsub(",", "", unlist(strsplit(toString(cox$waldtest), "\\s+")))
		waldtest = paste(paste(paste(paste(toString(signif(as.numeric(waldtest[1]), digits=4)), " on ", sep=""), toString(signif(as.numeric(waldtest[2]), digits=4)), sep=""), " df, p=", sep=""), toString(signif(as.numeric(waldtest[3]), digits=4)), sep="")
		score = gsub(",", "", unlist(strsplit(toString(cox$sctest), "\\s+")))
		score = paste(paste(paste(paste(toString(signif(as.numeric(score[1]), digits=4)), " on ", sep=""), toString(signif(as.numeric(score[2]), digits=4)), sep=""), " df, p=", sep=""), toString(signif(as.numeric(score[3]), digits=4)), sep="")
		
		sub_result_1 <- c(number_of_subjects, number_of_events, concordance, rsq, likelihood, waldtest, score)
		if(current_subset==1) {
			cox_regression_result_1[["SUB_1"]] = sub_result_1
		} else {
			cox_regression_result_1[["SUB_2"]] = sub_result_1
		}
		
		# RESULTS FOR TABLE 2
		subset = current_subset
		if(selected_categories_count > 0) {
			if(mergeCategories && selected_categories_count > 1) {
				category = paste(selected_categories, collapse=" & ")
				cox_coefficient = toString(signif(as.numeric(cox$coefficients[1]), digits=4))
				hazards_ratio = toString(signif(as.numeric(cox$coefficients[2]), digits=4))
				standard_deviation = toString(signif(as.numeric(cox$coefficients[3]), digits=4))
				z_score = toString(signif(as.numeric(cox$coefficients[4]), digits=4))
				lower_hazards_ratio = toString(signif(as.numeric(cox$conf.int[3]), digits=4))
				upper_hazards_ratio = toString(signif(as.numeric(cox$conf.int[4]), digits=4))
				sub_result_2 <- c(subset, category, cox_coefficient, hazards_ratio, standard_deviation, z_score, lower_hazards_ratio, upper_hazards_ratio)
				if(current_subset==1) {
					cox_regression_result_2[[1]] = sub_result_2
				} else {
					cox_regression_result_2[[2]] = sub_result_2
				}
			} else {
				index = 1
				for(selected_category in selected_categories) {
					category = selected_category
					cox_coefficient = toString(signif(as.numeric(cox$coefficients[1 + index - 1]), digits=4))
					hazards_ratio = toString(signif(as.numeric(cox$coefficients[2 + index - 1]), digits=4))
					standard_deviation = toString(signif(as.numeric(cox$coefficients[3 + index - 1]), digits=4))
					z_score = toString(signif(as.numeric(cox$coefficients[4 + index - 1]), digits=4))
					lower_hazards_ratio = toString(signif(as.numeric(cox$conf.int[3 + index - 1]), digits=4))
					upper_hazards_ratio = toString(signif(as.numeric(cox$conf.int[4 + index - 1]), digits=4))
					sub_result_2 <- c(subset, category, cox_coefficient, hazards_ratio, standard_deviation, z_score, lower_hazards_ratio, upper_hazards_ratio)
					if(current_subset==1) {
						cox_regression_result_2[[1+index-1]] = sub_result_2
					} else {
						cox_regression_result_2[[2+index-2+selected_categories_count]] = sub_result_2
					}
					index = index + 1
				}
			}
		} else {
			category = "NO_CATEGORY_SELECTED"
			cox_coefficient = "- - -"
			hazards_ratio = "- - -"
			standard_deviation = "- - -"
			z_score = "- - -"
			lower_hazards_ratio = "- - -"
			upper_hazards_ratio = "- - -"
			sub_result_2 <- c(subset, category, cox_coefficient, hazards_ratio, standard_deviation, z_score, lower_hazards_ratio, upper_hazards_ratio)
			cox_regression_result_2[[current_subset]] = sub_result_2
		}
		
		current_subset = current_subset + 1
		
	}
	
	cox_regression_result_1 <<- cox_regression_result_1
	cox_regression_result_2 <- t(cox_regression_result_2)
	cox_regression_result_2 <<- cox_regression_result_2
	
}

convertDataSetsTimeMeasurement <- function(timeIn, timeOut) {
	
	# check if conversion is necessary
	if(!(timeIn == timeOut)) {
		if(timeIn == "days" && timeOut == "months") {
			convertDataSetsFromDaysToMonths()
		}
		if(timeIn == "days" && timeOut == "years") {
			convertDataSetsFromDaysToYears()
		}
		if(timeIn == "months" && timeOut == "days") {
			convertDataSetsFromMonthsToDays()
		}
		if(timeIn == "months" && timeOut == "years") {
			convertDataSetsFromMonthsToYears()
		}
		if(timeIn == "years" && timeOut == "days") {
			convertDataSetsFromYearsToDays()
		}
		if(timeIn == "years" && timeOut == "months") {
			convertDataSetsFromYearsToMonths()
		}
	}
	
}

convertDataSetsFromDaysToMonths <- function() {
	# Subset 1
	if(nrow(subset_1) > 0) {
		subset_1$TIME <<- convertTimeFromDaysToMonths(subset_1$TIME)
	}
	# Subset 2
	if(nrow(subset_2) > 0) {
		subset_2$TIME <<- convertTimeFromDaysToMonths(subset_2$TIME)
	}
}

convertDataSetsFromDaysToYears <- function() {
	# Subset 1
	if(nrow(subset_1) > 0) {
		subset_1$TIME <<- convertTimeFromDaysToYears(subset_1$TIME)
	}
	# Subset 2
	if(nrow(subset_2) > 0) {
		subset_2$TIME <<- convertTimeFromDaysToYears(subset_2$TIME)
	}
}

convertDataSetsFromMonthsToDays <- function() {
	# Subset 1
	if(nrow(subset_1) > 0) {
		subset_1$TIME <<- convertTimeFromMonthsToDays(subset_1$TIME)
	}
	# Subset 2
	if(nrow(subset_2) > 0) {
		subset_2$TIME <<- convertTimeFromMonthsToDays(subset_2$TIME)
	}
}

convertDataSetsFromMonthsToYears <- function() {
	# Subset 1
	if(nrow(subset_1) > 0) {
		subset_1$TIME <<- convertTimeFromMonthsToYears(subset_1$TIME)
	}
	# Subset 2
	if(nrow(subset_2) > 0) {
		subset_2$TIME <<- convertTimeFromMonthsToYears(subset_2$TIME)
	}
}

convertDataSetsFromYearsToDays <- function() {
	# Subset 1
	if(nrow(subset_1) > 0) {
		subset_1$TIME <<- convertTimeFromYearsToDays(subset_1$TIME)
	}
	# Subset 2
	if(nrow(subset_2) > 0) {
		subset_2$TIME <<- convertTimeFromYearsToDays(subset_2$TIME)
	}
}

convertDataSetsFromYearsToMonths <- function() {
	# Subset 1
	if(nrow(subset_1) > 0) {
		subset_1$TIME <<- convertTimeFromYearsToMonths(subset_1$TIME)
	}
	# Subset 2
	if(nrow(subset_2) > 0) {
		subset_2$TIME <<- convertTimeFromYearsToMonths(subset_2$TIME)
	}
}

convertTimeFromDaysToMonths <- function(timeInDays) {
	return(timeInDays / dayMonthRelation)
}

convertTimeFromDaysToYears <- function(timeInDays) {
	return(timeInDays / dayYearRelation)
}

convertTimeFromMonthsToDays <- function(timeInMonths) {
	return(timeInMonths * dayMonthRelation)
}

convertTimeFromMonthsToYears <- function(timeInMonths) {
	return(timeInMonths / monthYearRelation)
}

convertTimeFromYearsToDays <- function(timeInYears) {
	return(timeInYears * dayYearRelation)
}

convertTimeFromYearsToMonths <- function(timeInYears) {
	return(timeInYears * monthYearRelation)
}