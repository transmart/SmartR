library(reshape2)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
vcfFile <- "/tmp/variants.vcf"

main <- function(variants) {
    output <- list()

    save(variants, file="/tmp/variants.Rda")
    save(loaded_variables, file="/tmp/loaded_variables.Rda")
    save(fetch_params, file="/tmp/fetch_params.Rda")

    vcf <- variants
    names(vcf) <- c("subject", "gene", "POS", "REF", "ALT", "CHROM", "var", "frq", "subset")
    vcf <- cbind(vcf, QUAL=rep(".", nrow(vcf)),
                      FILTER=rep(".", nrow(vcf)),
                      ID=rep(".", nrow(vcf)),
                      INFO=rep(".", nrow(vcf)),
                      FORMAT=rep("GT", nrow(vcf)))
    vcf <- dcast(vcf, CHROM+POS+ID+REF+ALT+QUAL+FILTER+INFO+FORMAT~subject, fill="0|0", value.var="var")
    NAMES <- names(vcf)
    NAMES[1] <- paste0("#", NAMES[1])
    cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n', file=vcfFile)
    write.table(vcf, file=vcfFile, row.names=F, col.names=NAMES, sep="\t", quote=F, append=T)
    vcf <- readVcf(vcfFile)
    seqlevels(vcf) = paste0("chr", seqlevels(vcf))
    rd <- rowRanges(vcf)
    loc <- locateVariants(rd, txdb, CodingVariants())
    coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)

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
