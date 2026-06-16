.n_info <- function(gds)
{
    f <- SeqArray::seqOpen(gds, allow.duplicate=TRUE)
    on.exit(SeqArray::seqClose(f))
    nm <- ls.gdsn(index.gdsn(f, "annotation/info"), recursive=TRUE,
        include.dirs=FALSE)
    sum(!grepl("^@|/@", nm))
}

.n_var <- function(gds)
{
    f <- SeqArray::seqOpen(gds, allow.duplicate=TRUE)
    on.exit(SeqArray::seqClose(f))
    objdesp.gdsn(index.gdsn(f, "variant.id"))$dim
}

test_that("seqToGDS_FAVOR converts CSV and attaches a digest", {
    out <- tempfile(fileext=".gds")
    on.exit(unlink(out, force=TRUE))
    seqToGDS_FAVOR(ex_file("favor_chr22_sub.csv.gz"), out, compress="ZIP",
        verbose=FALSE)
    expect_gt(.n_var(out), 0L)
    expect_gt(.n_info(out), 0L)
    f <- SeqArray::seqOpen(out, allow.duplicate=TRUE)
    on.exit(SeqArray::seqClose(f), add=TRUE)
    expect_true("md5" %in% names(get.attr.gdsn(index.gdsn(f, "position"))))
})

test_that("seqToGDS_gnomAD converts a VCF", {
    out <- tempfile(fileext=".gds")
    on.exit(unlink(out, force=TRUE))
    seqToGDS_gnomAD(ex_file("gnomad.genomes.v4.chr22_sub.vcf.gz"), out,
        compress="ZIP", verbose=FALSE)
    expect_gt(.n_var(out), 0L)
    expect_gt(.n_info(out), 0L)
})

test_that("seqToGDS_VEP splits CSQ and honours LZMA recompression", {
    out <- tempfile(fileext=".gds")
    on.exit(unlink(out, force=TRUE))
    seqToGDS_VEP(ex_file("example_wgs_sites_chr22_vep.vcf.gz"), out,
        compress="LZMA", verbose=FALSE)
    f <- SeqArray::seqOpen(out, allow.duplicate=TRUE)
    on.exit(SeqArray::seqClose(f), add=TRUE)
    csq <- ls.gdsn(index.gdsn(f, "annotation/info/CSQ.list"))
    expect_gt(length(csq), 0L)
    # recompression to LZMA must have actually happened
    nd <- index.gdsn(f, paste0("annotation/info/CSQ.list/", csq[1]))
    expect_match(objdesp.gdsn(nd)$compress, "^LZMA")
})

test_that("seqToGDS_SnpEff splits the ANN field", {
    out <- tempfile(fileext=".gds")
    on.exit(unlink(out, force=TRUE))
    seqToGDS_SnpEff(ex_file("example_wgs_sites_chr22_snpeff.vcf.gz"), out,
        compress="ZIP", verbose=FALSE)
    f <- SeqArray::seqOpen(out, allow.duplicate=TRUE)
    on.exit(SeqArray::seqClose(f), add=TRUE)
    expect_gt(length(ls.gdsn(index.gdsn(f, "annotation/info/ANN.list"))), 0L)
})
