# Shared helpers for the test suite

ex_file <- function(fn)
    system.file("extdata", fn, package="GDSAnnotator", mustWork=TRUE)

favor_gds <- function() ex_file("favor_chr22_sub.gds")

# Build a tiny SeqArray genotype GDS over a given set of variants, so the
# write path of seqAnnotate(add_to_gds=) can be exercised.
make_geno_gds <- function(chr, pos, ref, alt)
{
    vcf <- tempfile(fileext=".vcf")
    hdr <- c("##fileformat=VCFv4.2",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "##contig=<ID=22>",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2")
    rows <- sprintf("%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t0/1\t0/0",
        chr, as.integer(pos), ref, alt)
    writeLines(c(hdr, rows), vcf)
    gds <- tempfile(fileext=".gds")
    SeqArray::seqVCF2GDS(vcf, gds, verbose=FALSE)
    unlink(vcf, force=TRUE)
    gds
}
