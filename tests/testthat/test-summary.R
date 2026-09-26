test_that("seqValueCounts tabulates a categorical INFO field", {
    gds <- ex_file("gnomad.genomes.v4.chr22_sub.gds")
    f <- SeqArray::seqOpen(gds, allow.duplicate=TRUE)
    on.exit(SeqArray::seqClose(f))
    tab <- seqValueCounts(f, "variant_type", verbose=FALSE)
    expect_s3_class(tab, "table")
    expect_gt(sum(tab), 0L)
})

test_that("per_variant=TRUE counts each (variant, value) pair once", {
    gds <- ex_file("example_wgs_sites_chr22_snpeff.gds")
    a <- seqValueCounts(gds, "ANN.list/Annotation_Impact", verbose=FALSE)
    b <- seqValueCounts(gds, "ANN.list/Annotation_Impact", per_variant=TRUE,
        verbose=FALSE)
    expect_s3_class(b, "table")
    expect_setequal(names(b), names(a))
    expect_true(all(b[names(a)] <= a))
    expect_lt(sum(b), sum(a))
})

test_that("the counts of several files are merged", {
    gds <- ex_file("gnomad.genomes.v4.chr22_sub.gds")
    cp <- tempfile(fileext=".gds")
    file.copy(gds, cp)
    on.exit(unlink(cp, force=TRUE))
    a <- seqValueCounts(gds, "variant_type", verbose=FALSE)
    b <- seqValueCounts(c(gds, cp), "variant_type", verbose=FALSE)
    expect_equal(as.vector(b[names(a)]), 2L * as.vector(a))
})
