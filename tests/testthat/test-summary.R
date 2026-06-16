test_that("seqValueCounts tabulates a categorical INFO field", {
    gds <- ex_file("gnomad.genomes.v4.chr22_sub.gds")
    f <- SeqArray::seqOpen(gds, allow.duplicate=TRUE)
    on.exit(SeqArray::seqClose(f))
    tab <- seqValueCounts(f, "variant_type", verbose=FALSE)
    expect_s3_class(tab, "table")
    expect_gt(sum(tab), 0L)
})

test_that("seqUnitGroupAnnot works on a scalar field via a cond function", {
    unit <- seqUnitGroupAnnot(favor_gds(), varnm="cadd_phred",
        cond=list(function(x) ifelse(x >= 10, "damaging", "tolerated")),
        verbose=FALSE)
    expect_s3_class(unit, "SeqUnitListClass")
    expect_true(all(c("desp","index") %in% names(unit)))
    expect_true(all(c("chr","start","end") %in% names(unit$desp)))
    expect_setequal(unit$desp[[1]], c("damaging","tolerated"))
    # the group sizes must account for every annotated variant
    expect_gt(sum(lengths(unit$index)), 0L)
})

test_that("use_info=TRUE and use_info=FALSE give the same result", {
    a <- seqUnitGroupAnnot(favor_gds(), varnm="cadd_phred",
        use_info=TRUE, verbose=FALSE)
    b <- seqUnitGroupAnnot(favor_gds(), varnm="annotation/info/cadd_phred",
        use_info=FALSE, verbose=FALSE)
    expect_equal(lengths(a$index), lengths(b$index))
    expect_equal(nrow(a$desp), nrow(b$desp))
})

test_that("seqUnitGroupAnnot validates the 'by' argument", {
    expect_error(
        seqUnitGroupAnnot(favor_gds(), varnm="cadd_phred", by=5L,
            verbose=FALSE),
        "by")
})
