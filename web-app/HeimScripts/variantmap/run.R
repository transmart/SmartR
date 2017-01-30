library(reshape2)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
vcfFile <- "variants.vcf"

main <- function(variants) {
    output <- list()

    # save(variants, file="/tmp/variants.Rda")
    # save(loaded_variables, file="/tmp/loaded_variables.Rda")
    # save(fetch_params, file="/tmp/fetch_params.Rda")

    vcf <- variants
    names(vcf) <- c("subject", "gene", "POS", "REF", "ALT", "CHROM", "var", "frq", "subset")
    vcf <- cbind(vcf, QUAL=rep(".", nrow(vcf)),
                      FILTER=rep(".", nrow(vcf)),
                      ID=rep(".", nrow(vcf)),
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
    locations <- aggregate(LOCATION ~ QUERYID, data=mcols(loc), FUN=function(x) toString(unique(x)))
    consequences <- aggregate(CONSEQUENCE ~ QUERYID, data=mcols(coding), FUN=function(x) toString(unique(x)))
    annotation.data <- merge(consequences, annotation.data, by="QUERYID", all.y=T)
    annotation.data <- merge(locations, annotation.data, by="QUERYID", all.y=T)
    annotation.data <- melt(annotation.data, id.vars=c("QUERYID", "LOCATION", "CONSEQUENCE", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),
         variable.name="subject", value.name="var")
    annotation.data <- annotation.data[, names(annotation.data) %in% c("LOCATION", "CONSEQUENCE", "CHROM", "POS", "subject")]
    names(annotation.data) <- c("location", "consequence", "chr", "pos", "subject")
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
        expression.data <- merge(variants, hdd.data[, -1], by=c("subject", "gene"), all.x=T)
        data <- merge(annotation.data, expression.data, by=c("subject", "chr", "pos"))
    } else {
        data <- annotation.data
    }

    output$data <- data
    json <- toJSON(output)
    return(json)
}
