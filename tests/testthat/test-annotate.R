test_that("seqAnnotList returns a DataFrame describing the INFO fields", {
    ann <- seqAnnotList(favor_gds())
    expect_s4_class(ann, "DataFrame")
    expect_true(all(c("name","type","trait","description") %in% colnames(ann)))
    expect_gt(nrow(ann), 0L)
})

# Pull real variants from the example file to build matching queries.
.sample_variants <- function(n=5L)
{
    f <- SeqArray::seqOpen(favor_gds(), allow.duplicate=TRUE)
    on.exit(SeqArray::seqClose(f))
    chr <- seqGetData(f, "$chromosome")
    pos <- seqGetData(f, "position")
    ref <- seqGetData(f, "$ref")
    alt <- seqGetData(f, "$alt")
    val <- seqGetData(f, "annotation/info/cadd_phred")
    k <- head(which(is.finite(val)), n)
    list(chr=chr[k], pos=pos[k], ref=ref[k], alt=alt[k], cadd=val[k])
}

test_that("seqAnnotate works with a single annotation name (drop=FALSE)", {
    v <- .sample_variants(5L)
    snp <- paste(v$chr, v$pos, v$ref, v$alt, sep="-")
    # single varnm must not collapse the DataFrame to a vector
    res <- seqAnnotate(snp, favor_gds(), varnm="cadd_phred", verbose=FALSE)
    expect_s4_class(res, "DataFrame")
    expect_equal(ncol(res), 1L)
    expect_equal(as.numeric(res[["cadd_phred"]]), as.numeric(v$cadd),
        tolerance=1e-5)
})

test_that("seqAnnotate matches across character, data.frame and GRanges", {
    v <- .sample_variants(5L)
    snp <- paste(v$chr, v$pos, v$ref, v$alt, sep="-")
    a <- seqAnnotate(snp, favor_gds(), varnm=c("cadd_phred","linsight"),
        verbose=FALSE)
    df <- data.frame(chr=v$chr, pos=v$pos, ref=v$ref, alt=v$alt)
    b <- seqAnnotate(df, favor_gds(), varnm=c("cadd_phred","linsight"),
        verbose=FALSE)
    expect_equal(as.numeric(a[["cadd_phred"]]), as.numeric(b[["cadd_phred"]]),
        tolerance=1e-5)
})

test_that("a variant absent from the file yields NA", {
    res <- seqAnnotate("22-1-A-T", favor_gds(), varnm="cadd_phred",
        verbose=FALSE)
    expect_equal(nrow(res), 1L)
    expect_true(is.na(as.numeric(res[["cadd_phred"]])))
})

test_that("add_to_gds writes correctly aligned annotation values", {
    v <- .sample_variants(5L)
    # append a fake variant that is not in the annotation file
    geno <- make_geno_gds(c(v$chr, "22"), c(v$pos, 1L),
        c(v$ref, "A"), c(v$alt, "T"))
    on.exit(unlink(geno, force=TRUE))
    out <- tempfile(fileext=".gds")
    on.exit(unlink(out, force=TRUE), add=TRUE)

    seqAnnotateGDS(geno, favor_gds(), varnm="cadd_phred",
        add_to_gds=out, verbose=FALSE)

    g <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(g), add=TRUE)
    written <- seqGetData(g, "annotation/info/cadd_phred")
    expect_equal(as.numeric(written[seq_along(v$cadd)]),
        as.numeric(v$cadd), tolerance=1e-5)
    expect_true(is.na(as.numeric(written[length(written)])))  # fake variant
})
