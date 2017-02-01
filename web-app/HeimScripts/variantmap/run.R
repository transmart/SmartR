library(reshape2)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(PolyPhen.Hsapiens.dbSNP131)
library(SIFT.Hsapiens.dbSNP132)


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
vcfFile <- "variants.vcf"

main <- function(variants) {
    output <- list()

    load("/tmp/json.Rda")
    return(json)

    save(variants, file="/tmp/variants.Rda")
    save(loaded_variables, file="/tmp/loaded_variables.Rda")
    save(fetch_params, file="/tmp/fetch_params.Rda")

    vcf <- variants
    names(vcf) <- c("subject", "gene", "POS", "REF", "ALT", "CHROM", "var", "ID", "frq", "subset")
    vcf <- cbind(vcf, QUAL=rep(".", nrow(vcf)),
                      FILTER=rep(".", nrow(vcf)),
                      INFO=rep(".", nrow(vcf)),
                      FORMAT=rep("GT", nrow(vcf)))
    vcf <- dcast(vcf, CHROM+POS+ID+REF+ALT+QUAL+FILTER+INFO+FORMAT~subject, fill="0|0", value.var="var")
    annotation.data <- vcf
    annotation.data <- cbind(QUERYID=rownames(annotation.data), annotation.data)
    NAMES <- names(vcf)
    NAMES[1] <- paste0("#", NAMES[1])
    cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n', file=vcfFile)
    write.table(vcf, file=vcfFile, row.names=F, col.names=NAMES, sep="\t", quote=F, append=T)
    vcf <- readVcf(vcfFile)
    seqlevels(vcf) = paste0("chr", seqlevels(vcf))
    rd <- rowRanges(vcf)
    loc <- locateVariants(rd, txdb, AllVariants())
    coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)
    nms <- names(coding)
    idx <- mcols(coding)$CONSEQUENCE != "synonymous"
    nonsyn <- coding[idx]
    names(nonsyn) <- nms[idx]
    rsids <- unique(names(nonsyn)[grep("rs", names(nonsyn), fixed=TRUE)])
    pp <- select(PolyPhen.Hsapiens.dbSNP131, keys=rsids, cols=c("TRAININGSET", "PREDICTION", "PPH2PROB"))
    prediction <- pp[!is.na(pp$PREDICTION), ][, c("RSID", "PREDICTION")]
    names(prediction) <- c("ID", "PREDICTION")
    locations <- aggregate(LOCATION ~ QUERYID, data=mcols(loc), FUN=function(x) toString(unique(x)))
    consequences <- aggregate(CONSEQUENCE ~ QUERYID, data=mcols(coding), FUN=function(x) toString(unique(x)))
    annotation.data <- merge(consequences, annotation.data, by="QUERYID", all.y=T)
    annotation.data <- merge(locations, annotation.data, by="QUERYID", all.y=T)
    annotation.data <- merge(prediction, annotation.data, by="ID", all.y=T)
    annotation.data <- melt(annotation.data, id.vars=c("QUERYID", "PREDICTION", "LOCATION", "CONSEQUENCE", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),
         variable.name="subject", value.name="var")
    annotation.data <- annotation.data[, names(annotation.data) %in% c("PREDICTION", "LOCATION", "CONSEQUENCE", "CHROM", "POS", "subject")]
    names(annotation.data) <- c("prediction", "location", "consequence", "chr", "pos", "subject")
    if (exists("loaded_variables")) {
        hdd.idx <- grep("^highDimensional", names(loaded_variables))
        cat.idx <- grep("^categoric", names(loaded_variables))
        num.idx <- grep("^numeric", names(loaded_variables))

        # HIGHDIMENSIONAL
        hdd.data <- loaded_variables[hdd.idx]
        hdd.data <- lapply(hdd.data, function(x) {
             names(x) <- sub("^X", "", names(x))
             x <- melt(x, id.vars=c("Row.Label", "Bio.marker"))
             x
        })
        hdd.data <- do.call(rbind, hdd.data)
        names(hdd.data) <- c("Row.Label", "gene", "subject", "expr")
        expression.data <- merge(variants, hdd.data[, -1], by=c("subject", "gene"), all.x=T)
        data <- merge(annotation.data, expression.data, by=c("subject", "chr", "pos"))

        zScore.data <- data[, c("subject", "gene", "expr")]
        zScore.data <- zScore.data[!duplicated(zScore.data[, 1:2]), ]
        zScore.data <- dcast(zScore.data, gene~subject)
        colNames <- colnames(zScore.data)
        genes <- zScore.data$gene
        zScore.data <- cbind(genes, t(apply(zScore.data[,-1], 1, scale)))
        colnames(zScore.data) <- colNames
        zScore.data <- as.data.frame(zScore.data)
        zScore.data <- melt(zScore.data, id.vars="gene", na.rm=T, variable.name="subject", value.name="zscore")
        data <- merge(data, zScore.data, by=c("subject", "gene"), all.x=T)

        # NUMERIC
        ldd.data <- loaded_variables[c(num.idx, cat.idx)]
        ldd.data <- lapply(seq_along(ldd.data), function(data, names, i) {
                   melt(data[i], id.vars="Row.Label")
        }, data=ldd.data, names=names(ldd.data))
        ldd.data <- do.call(rbind, ldd.data)
        ldd.data <- ldd.data[, -2]
        names(ldd.data) <- c("subject", "value", "node")
        fullNames <- do.call(rbind, lapply(fetch_params$ontologyTerms[sub("_s\\d$", "", ldd.data$node)], function(x) x$fullName))
        types <- sub("_.+", "", ldd.data$node)
        identifiers <- paste(types, "fullName:", fullNames)
        ldd.data$node <- identifiers
        names(ldd.data)[3] <- "identifier"
        ldd.data <- dcast(ldd.data, subject~...)
        data <- merge(data, ldd.data, by="subject", all.x=T)
    } else {
        data <- annotation.data
    }

    output$data <- data
    json <- toJSON(output)
    save(json, file="/tmp/json.Rda")
    return(json)
}
