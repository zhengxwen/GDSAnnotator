snpeff_gds <- function() ex_file("example_wgs_sites_chr22_snpeff.gds")
vep_gds <- function() ex_file("example_wgs_sites_chr22_vep.gds")

.n_variant <- function(gds)
{
    f <- SeqArray::seqOpen(gds, allow.duplicate=TRUE)
    on.exit(SeqArray::seqClose(f))
    objdesp.gdsn(index.gdsn(f, "variant.id"))$dim
}

test_that("seqAnnotStat summarises a SnpEff file, independently of bsize", {
    st <- seqAnnotStat(snpeff_gds(), verbose=FALSE)
    expect_s3_class(st, "SeqAnnotStat")
    expect_equal(st$source, "SnpEff")
    expect_equal(st$n_variant, .n_variant(snpeff_gds()))
    expect_equal(sum(st$counts$chrom), st$n_variant)
    expect_equal(sum(st$counts$var_class), sum(st$counts$n_allele) )
    expect_lte(sum(st$counts$cons_severe), st$n_variant)
    expect_equal(sum(st$counts$cons_severe), sum(st$counts$impact_severe))
    expect_true(all(vapply(st$counts, is.numeric, FALSE)))
    # the impact categories are ordered by severity
    expect_equal(names(st$counts$impact),
        intersect(c("HIGH", "MODERATE", "LOW", "MODIFIER"),
            names(st$counts$impact)))
    # a small block size gives exactly the same result
    st2 <- seqAnnotStat(snpeff_gds(), bsize=7L, verbose=FALSE)
    expect_identical(st2$counts, st$counts)
    expect_output(print(st), "SeqAnnotStat")
})

test_that("seqAnnotStat works on a VEP file and on a file without annotation", {
    st <- seqAnnotStat(vep_gds(), verbose=FALSE)
    expect_equal(st$source, "Ensembl VEP")
    expect_true(all(c("cons_all", "cons_severe", "impact", "impact_severe",
        "region", "biotype", "feature_type") %in% names(st$counts)))
    # the VEP version is taken from the VCF header
    expect_true("VEP" %in% st$header$id)
    expect_equal(names(st$counts)[1L], "chrom")
    st0 <- seqAnnotStat(ex_file("gnomad.genomes.v4.chr22_sub.gds"),
        verbose=FALSE)
    expect_equal(st0$source, "unknown")
    expect_null(st0$counts$cons_all)
    expect_equal(sum(st0$counts$chrom), st0$n_variant)
})

test_that("parallel processing gives the same counters", {
    st <- seqAnnotStat(vep_gds(), verbose=FALSE)
    st2 <- seqAnnotStat(vep_gds(), parallel=2L, verbose=FALSE)
    expect_identical(st2$counts, st$counts)
})

test_that("the statistics of several files are merged", {
    a <- seqAnnotStat(vep_gds(), verbose=FALSE)
    b <- seqAnnotStat(snpeff_gds(), verbose=FALSE)
    ab <- seqAnnotStat(c(vep_gds(), snpeff_gds()), verbose=FALSE)
    expect_equal(ab$n_variant, a$n_variant + b$n_variant)
    expect_length(ab$file, 2L)
    expect_equal(sum(ab$counts$chrom), ab$n_variant)
    expect_equal(sum(ab$counts$var_class),
        sum(a$counts$var_class) + sum(b$counts$var_class))
    expect_equal(sum(ab$counts$cons_all),
        sum(a$counts$cons_all) + sum(b$counts$cons_all))
})

test_that("seqAnnotGeneTable returns the per-gene counts", {
    st <- seqAnnotStat(snpeff_gds(), verbose=FALSE)
    tab <- seqAnnotGeneTable(st)
    expect_s3_class(tab, "data.frame")
    expect_equal(names(tab)[1:2], c("gene", "total"))
    expect_equal(sum(tab$total), sum(st$counts$gene_impact))
    expect_true(all(diff(tab$total) <= 0))
    expect_equal(unname(rowSums(tab[, -(1:2), drop=FALSE])), tab$total)
    # NULL when there is no gene annotation
    st0 <- seqAnnotStat(ex_file("gnomad.genomes.v4.chr22_sub.gds"),
        verbose=FALSE)
    expect_null(seqAnnotGeneTable(st0))
})

test_that("seqAnnotReport writes HTML, Markdown and R Markdown reports", {
    st <- seqAnnotStat(vep_gds(), verbose=FALSE)
    for (ext in c(".html", ".md", ".Rmd"))
    {
        fn <- tempfile(fileext=ext)
        out <- seqAnnotReport(st, fn, title="VEP test report", verbose=FALSE)
        expect_true(file.exists(fn))
        expect_equal(normalizePath(out), normalizePath(fn))
        txt <- readLines(fn, warn=FALSE)
        expect_true(any(grepl("VEP test report", txt, fixed=TRUE)))
        unlink(fn, force=TRUE)
    }
    # the HTML report is self-contained and the title is escaped
    fn <- tempfile(fileext=".html")
    seqAnnotReport(st, fn, title="a <b> & c", verbose=FALSE)
    txt <- paste(readLines(fn, warn=FALSE), collapse="\n")
    expect_match(txt, "a &lt;b&gt; &amp; c", fixed=TRUE)
    expect_match(txt, "<svg", fixed=TRUE)
    expect_false(grepl("<script", txt, fixed=TRUE))
    # the format can be forced, and a GDS file name is accepted directly
    fn2 <- tempfile(fileext=".txt")
    seqAnnotReport(st, fn2, format="markdown", verbose=FALSE)
    expect_match(readLines(fn2, n=1L, warn=FALSE), "^# ")
    fn3 <- tempfile(fileext=".md")
    seqAnnotReport(snpeff_gds(), fn3, verbose=FALSE)
    expect_true(any(grepl("SnpEff", readLines(fn3, warn=FALSE), fixed=TRUE)))
    unlink(c(fn, fn2, fn3), force=TRUE)
})

test_that("the R Markdown report can be rendered", {
    skip_if_not_installed("rmarkdown")
    skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
    st <- seqAnnotStat(snpeff_gds(), verbose=FALSE)
    rmd <- tempfile(fileext=".Rmd")
    seqAnnotReport(st, rmd, title="SnpEff test report", verbose=FALSE)
    html <- rmarkdown::render(rmd, output_file=tempfile(fileext=".html"),
        quiet=TRUE)
    on.exit(unlink(c(rmd, html), force=TRUE))
    expect_true(file.exists(html))
    expect_true(any(grepl("SnpEff test report", readLines(html, warn=FALSE),
        fixed=TRUE)))
})
