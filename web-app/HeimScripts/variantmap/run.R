library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

main <- function(variants) {
    output <- list()
    output$variantsMatrix <- variants
    save(variants, file="/tmp/variants.Rda")
    save(loaded_variables, file="/tmp/loaded_variables.Rda")
    save(fetch_params, file="/tmp/fetch_params.Rda")

    vcf <- variants
    vcf <- cbind(vcf, QUAL=rep(".", nrow(vcf)),
                      FILTER=rep(".", nrow(vcf)),
                      ID=rep(".", nrow(vcf)),
                      INFO=rep(".", nrow(vcf)),
                      FORMAT=rep("GT", nrow(vcf)))
    vcf <- dcast(vcf, CHROM+POS+ID+REF+ALT+QUAL+FILTER+INFO+FORMAT~subject, fill="0|0", value.var="var")
    NAMES <- names(vcf)
    NAMES[1] <- paste0("#", NAMES[1])
    # write to first line: ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    write.table(vcf, file="/tmp/variants.vcf", row.names=F, col.names=NAMES, sep="\t", quote=F)
    vcf <- readVcf("/tmp/variants.vcf")
    seqlevels(vcf) = paste0("chr", seqlevels(vcf))
    rd <- rowRanges(vcf)
    loc <- locateVariants(rd, txdb, CodingVariants())
    coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)
}
