test_that("seqToGDS_Nirvana imports the example JSON", {
    skip_if_not_installed("jsonlite")
    json <- ex_file("example_nirvana.json.gz")
    out <- tempfile(fileext=".gds")
    on.exit(unlink(out, force=TRUE))
    seqToGDS_Nirvana(json, out, compress="ZIP", verbose=FALSE)

    f <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(f), add=TRUE)

    # multiallelic position 10510038 is split into two records => 4 variants
    expect_equal(objdesp.gdsn(index.gdsn(f, "variant.id"))$dim, 4L)
    expect_equal(as.character(seqGetData(f, "$chromosome")), rep("22", 4L))
    expect_equal(seqGetData(f, "position"),
        c(10510007L, 10510038L, 10510038L, 10510282L))
    expect_equal(seqGetData(f, "allele"),
        c("T,C", "T,C", "T,A", "G,C"))

    info <- ls.gdsn(index.gdsn(f, "annotation/info"))
    # a field that is entirely missing (oneKg) must be dropped
    expect_false("onekg_af" %in% info)
    expect_true(all(c("gnomad_af", "clinvar", "gene", "consequence") %in% info))

    # numeric typing & values
    expect_equal(seqGetData(f, "annotation/info/gnomad_an"),
        c(152312L, 150000L, NA, 151000L))
    expect_equal(seqGetData(f, "annotation/info/gnomad_af"),
        c(0.0012, 0.02, NaN, 0.5), tolerance=1e-6)
    # multiple values for one variant are joined by 'sep'
    expect_equal(seqGetData(f, "annotation/info/consequence")[3],
        "missense_variant;splice_region_variant")
    # the alt-specific ClinVar significance lands on the right record
    expect_equal(seqGetData(f, "annotation/info/clinvar"),
        c("benign", "", "pathogenic", ""))
})

test_that("Nirvana header provenance is stored on 'description'", {
    skip_if_not_installed("jsonlite")
    out <- tempfile(fileext=".gds")
    on.exit(unlink(out, force=TRUE))
    seqToGDS_Nirvana(ex_file("example_nirvana.json.gz"), out, compress="ZIP",
        verbose=FALSE)
    f <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(f), add=TRUE)
    a <- get.attr.gdsn(index.gdsn(f, "description"))
    expect_equal(a$genomeAssembly, "GRCh38")
    expect_match(a$annotator, "Nirvana")
    expect_match(a$dataSources, "ClinVar")
})

test_that("the Nirvana GDS can be queried with seqAnnotate", {
    skip_if_not_installed("jsonlite")
    out <- tempfile(fileext=".gds")
    on.exit(unlink(out, force=TRUE))
    seqToGDS_Nirvana(ex_file("example_nirvana.json.gz"), out, compress="ZIP",
        verbose=FALSE)
    res <- seqAnnotate(c("22-10510007-T-C", "22-10510038-T-A"), out,
        varnm=c("clinvar", "gene"), verbose=FALSE)
    expect_equal(res$clinvar, c("benign", "pathogenic"))
    expect_equal(res$gene, c("GENEA", "GENEB"))
})

test_that("fields can be selected by catalogue Name with catalogue types", {
    skip_if_not_installed("jsonlite")
    out <- tempfile(fileext=".gds")
    on.exit(unlink(out, force=TRUE))
    seqToGDS_Nirvana(ex_file("example_nirvana.json.gz"), out,
        fields=c("gene", "gnomad_af", "gnomad_an"), compress="ZIP",
        verbose=FALSE)
    f <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(f), add=TRUE)
    expect_setequal(ls.gdsn(index.gdsn(f, "annotation/info")),
        c("gene", "gnomad_af", "gnomad_an"))
    # gnomad_an is declared 'integer' in nirvana_output_format.csv
    expect_equal(as.character(
        objdesp.gdsn(index.gdsn(f, "annotation/info/gnomad_an"))$type),
        "Integer")
    # gnomad_af is declared 'float'
    expect_equal(as.character(
        objdesp.gdsn(index.gdsn(f, "annotation/info/gnomad_af"))$type),
        "Real")
})

test_that("the field catalogue is installed and well-formed", {
    fn <- system.file("extdata", "nirvana_output_format.csv",
        package="GDSAnnotator")
    expect_true(nzchar(fn))
    d <- read.csv(fn, stringsAsFactors=FALSE)
    expect_true(all(c("Name", "Path", "Type", "Default") %in% colnames(d)))
    expect_false(anyDuplicated(d$Name) > 0L)
})

test_that("a custom field selection works (user-supplied dot paths)", {
    skip_if_not_installed("jsonlite")
    out <- tempfile(fileext=".gds")
    on.exit(unlink(out, force=TRUE))
    seqToGDS_Nirvana(ex_file("example_nirvana.json.gz"), out,
        fields=c(af="gnomad.allAf", csq="transcripts.consequence"),
        compress="ZIP", verbose=FALSE)
    f <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(f), add=TRUE)
    expect_setequal(ls.gdsn(index.gdsn(f, "annotation/info")), c("af", "csq"))
})

test_that("a non-Nirvana file is rejected", {
    skip_if_not_installed("jsonlite")
    bad <- tempfile(fileext=".json.gz")
    con <- gzfile(bad, "wt"); writeLines('{"foo":1}', con); close(con)
    on.exit(unlink(bad, force=TRUE))
    expect_error(seqToGDS_Nirvana(bad, tempfile(fileext=".gds"),
        verbose=FALSE), "Nirvana")
})
