
main <- function(analysis_type = "unpaired_group") {

#    save(loaded_variables, file="/tmp/data/loaded_variables.Rda")
#    save(fetch_params, file="/tmp/data/fetch_params.Rda")

    if (is.null(fetch_params$assayConstraints)) {
        stop("No assays satisfy the provided criteria")
    }

    if (analysis_type=="unpaired_group") {
        print("Two group analysis (not paired):")

        ## conditionsFile	= read.delim(fileName,header=F,stringsAsFactors=F)
        ## countfiles 		= conditionsFile[,1]
        ## conditions 		= as.factor(conditionsFile[,2])
        ## sampleID		= conditionsFile[,3]
        ## print("Reading count files...")
        ## countTable 		= read.delim(countfiles[1],header=F,stringsAsFactors=F,row.names=1)
        ## for (b in 2:length(countfiles)) {
        ## 	countTable = cbind(countTable,read.delim(countfiles[b],header=F,stringsAsFactors=F,row.names=1))
        ## }

        countTable <- read.table(readcountfileName, header=TRUE, sep='\t', quote='"', as.is=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
        phenodata  <- read.table(phenodatafileName, header=TRUE, sep='\t', quote='"', strip.white=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

        # Make rownames equal to the regionname
        if ( 'regionname' %in% colnames(countTable) ) {
        rownames(countTable) <- countTable$regionname
        } else {
        stop("||FRIENDLY||Expecting readcountTable to at least have a column regionname. Please check your region variable selection and run again.")
        }

        # Filter phenodata for patients that have data in countTable
        phenodata <- phenodata[paste("readcount.",phenodata$PATIENT_NUM,sep="") %in% colnames(countTable), ]

        # Filter for HTSeq predefined counts:
        exclude_HTSeq = c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
        exclude_DEXSeq = c("_ambiguous","_empty","_lowaqual","_notaligned")
        exclude = match(c(exclude_HTSeq, exclude_DEXSeq),rownames(countTable))
        exclude = exclude[is.na(exclude)==0]
        if(length(exclude) != 0)  {
        countTable = countTable[-exclude,]
        }

        # Make sure that the order of the subjects/samples in the readcount columns is consistent with the order of the subjects/samples in the phenodata rows
        reionGeneSymbol <- countTable[c('regionname', 'genesymbol')]
        countTable = countTable[grep('^readcount.', colnames(countTable))]

        # If provided data set does not contain readcount data at all, stop further processing
        if (all(is.na(countTable))) stop("||FRIENDLY||R cannot perform edgeR analysis when NO readcount data is provided for any of the samples. Please check your variable or cohort selection and run again.")

        # If provided data set has missing readcount data, find out if there is some structure in the missing data values
        # try if removing samples (columns) or transcripts (rows) which do not contain any data (rows or columns with NA's only), leaves us a valid data subset
        if (any(is.na(countTable))) {
        # only keep rows which contain not only NA's (< ncol)
        ncol = ncol(countTable)
        countTable <- countTable[ rowSums(is.na(countTable)) != ncol , ]
        # only keep columns which contain not only NA's (< nrow)
        nrow = nrow(countTable)
        countTable <- countTable[ , colSums(is.na(countTable)) != nrow ]
        }
        # If data set does not contain readcount data for all transcripts and samples, stop further processing
        if (nrow(countTable)==0 | any(is.na(countTable))) stop("||FRIENDLY||R cannot perform edgeR analysis if not all readcount data is provided for all selected samples. Please check your variable or cohort selection and run again.")

        # Make row names equal to the sample id
        rownames(phenodata) <- phenodata[,"PATIENT_NUM"]

        # Extract sample list from RNASeq data column names for which readcounts have been observed
        samplelist <- sub("readcount.", "" , colnames(countTable))
        # Find out if readcounts are available for samples for which no group information is available. Ignore the readcounts for those samples.
        countTable <- countTable[ , samplelist %in% rownames(phenodata) ]
        # Update sample list
        samplelist <- sub("readcount.", "" , colnames(countTable))

        # Reorder phenodata rows to match the order in the RNASeq data columns
        phenodata <- phenodata[samplelist,,drop=FALSE]

        conditions <- phenodata$group

        dge 	= DGEList(counts=countTable,group=conditions,genes=rownames(countTable))
        print("Calculating normalization factors...")
        dge		= calcNormFactors(dge)
        print("Estimating common dispersion...")
        dge 	= estimateCommonDisp(dge)
        print("Estimating tagwise dispersion...")
        dge 	= estimateTagwiseDisp(dge)
    } else {
        print("Multi group analysis")
        ## conditionsFile 	= read.delim(fileName,header=F,stringsAsFactors=F)
        ## countfiles 		= conditionsFile[,1]
        ## conditions 		= as.factor(conditionsFile[,2])
        ## sampleID		= conditionsFile[,3]
        ## print("Reading count files...")
        ## countTable 		= read.delim(countfiles[1],header=F,stringsAsFactors=F,row.names=1)
        ## for (b in 2:length(countfiles)) {
        ## 	countTable = cbind(countTable,read.delim(countfiles[b],header=F,stringsAsFactors=F,row.names=1))
        ## }

        countTable <- read.table(readcountfileName, header=TRUE, sep='\t', quote='"', as.is=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
        phenodata  <- read.table(phenodatafileName, header=TRUE, sep='\t', quote='"', strip.white=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

        # Make rownames equal to the regionname
        if ( 'regionname' %in% colnames(countTable) ) {
            rownames(countTable) <- countTable$regionname
        } else {
            stop("||FRIENDLY||Expecting readcountTable to at least have a column regionname. Please check your region variable selection and run again.")
        }

        # Filter phenodata for patients that have data in countTable
        phenodata <- phenodata[paste("readcount.",phenodata$PATIENT_NUM,sep="") %in% colnames(countTable), ]

        # Make sure that the order of the subjects/samples in the readcount columns is consistent with the order of the subjects/samples in the phenodata rows
        # Extract sample list from RNASeq data column names for which readcounts have been observed
        reionGeneSymbol <- countTable[c('regionname', 'genesymbol')]
        countTable = countTable[grep('readcount.', colnames(countTable))]
        samplelist <- sub("readcount.", "" , colnames(countTable))

        # Make row names equal to the sample id
        rownames(phenodata) <- phenodata[,"PATIENT_NUM"]
        # Reorder phenodata rows to match the order in the RNASeq data columns
        phenodata <- phenodata[samplelist,,drop=FALSE]

        conditions <- sort(phenodata$group)
        ngrp <- length(unique(conditions))

        # Filter for HTSeq predifined counts:
        exclude_HTSeq = c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
        exclude_DEXSeq = c("_ambiguous","_empty","_lowaqual","_notaligned")
        exclude = match(c(exclude_HTSeq, exclude_DEXSeq),rownames(countTable))
        exclude = exclude[is.na(exclude)==0]
        if(length(exclude) != 0)  {
            countTable = countTable[-exclude,]
        }

        #colnames(countTable) 	= sampleID
        dge = DGEList(counts=countTable,genes=rownames(countTable))
        design = model.matrix(~conditions,data=dge$samples)
        rownames(design) = colnames(dge)
        print("Calculating normalization factors...")
        dge		= calcNormFactors(dge)
        print("Estimating common dispersion...")
        dge 	= estimateGLMCommonDisp(dge,design)
        print("Estimating trended dispersion...")
        dge 	= estimateGLMTrendedDisp(dge,design)
        print("Estimating tagwise dispersion...")
        dge 	= estimateGLMTagwiseDisp(dge,design)
    }

    output <- list(
        plot_title = fetch_params$ontologyTerms$datapoints_n0$fullName,
        dge = dge
    )
    print("Done!")

    toJSON(output)
}

