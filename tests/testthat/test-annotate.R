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

test_that("seqAnnotateVCF matches exact alleles in a sites-only VCF", {
    skip_if_not_installed("VariantAnnotation")
    v <- .sample_variants(1L)
    vcf <- tempfile(fileext=".vcf")
    on.exit(unlink(vcf, force=TRUE))
    hdr <- c("##fileformat=VCFv4.2", "##contig=<ID=22>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    row <- sprintf("%s\t%d\t.\t%s\t%s\t.\tPASS\t.",
        v$chr, v$pos, v$ref, v$alt)
    writeLines(c(hdr, row), vcf)

    res <- seqAnnotateVCF(vcf, favor_gds(), varnm="cadd_phred",
        verbose=FALSE)

    expect_equal(nrow(res), 1L)
    expect_equal(as.numeric(res[["cadd_phred"]]), as.numeric(v$cadd),
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

test_that("a small block size gives the same output as a single block", {
    v <- .sample_variants(20L)
    # append a fake variant that is not in the annotation file
    geno <- make_geno_gds(c(v$chr, "22"), c(v$pos, 1L),
        c(v$ref, "A"), c(v$alt, "T"))
    on.exit(unlink(geno, force=TRUE))
    # 'cadd_phred' is fixed-length, 'rsid' is a character variable
    nm <- c("cadd_phred", "rsid")

    get_out <- function(bsize)
    {
        out <- tempfile(fileext=".gds")
        on.exit(unlink(out, force=TRUE))
        seqAnnotateGDS(geno, favor_gds(), varnm=nm, add_to_gds=out,
            bsize=bsize, verbose=FALSE)
        g <- SeqArray::seqOpen(out)
        on.exit(SeqArray::seqClose(g), add=TRUE, after=FALSE)
        lapply(nm, function(s)
        {
            n <- index.gdsn(g, paste0("annotation/info/", s))
            # the storage mode should not depend on the block size either
            list(val=seqGetData(g, paste0("annotation/info/", s), .tolist=TRUE),
                storage=objdesp.gdsn(n)$storage)
        })
    }

    single <- get_out(1000000L)
    for (bs in c(1L, 3L, 20L, 21L))
        expect_equal(get_out(bs), single)
})

test_that("varnm can be missing (all annotations) or NULL (locate only)", {
    v <- .sample_variants(3L)
    snp <- paste(v$chr, v$pos, v$ref, v$alt, sep="-")
    ann <- seqAnnotList(favor_gds())
    df <- data.frame(chr=v$chr, pos=v$pos, ref=v$ref, alt=v$alt)
    # a missing 'varnm' means all the annotations of the file
    res <- seqAnnotate(df, favor_gds(), verbose=FALSE)
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), ann$name)
    expect_equal(nrow(res), 3L)
    res2 <- seqAnnotate(snp, favor_gds(), verbose=FALSE)
    expect_equal(colnames(res2), ann$name)
    # NULL: only the location of the variants
    loc <- seqAnnotate(snp, favor_gds(), varnm=NULL, verbose=FALSE)
    expect_setequal(colnames(loc), c("variant_idx", "file_idx"))
    expect_false(anyNA(loc$variant_idx))
    expect_equal(as.integer(loc$file_idx), rep(1L, 3L))
})

test_that("GRanges and IRanges queries return the variants in the ranges", {
    v <- .sample_variants(5L)
    rng <- range(v$pos)
    gr <- GenomicRanges::GRanges(seqnames="chr22",
        ranges=IRanges::IRanges(start=rng[1L], end=rng[2L]))
    a <- seqAnnotate(gr, favor_gds(), varnm="cadd_phred", verbose=FALSE)
    expect_s4_class(a, "DataFrame")
    expect_gte(nrow(a), 5L)
    expect_true(all(v$cadd %in% as.numeric(a$cadd_phred)))
    b <- seqAnnotate(IRanges::IRanges(start=rng[1L], end=rng[2L]),
        favor_gds(), varnm="cadd_phred", chr="22", verbose=FALSE)
    expect_equal(as.numeric(a$cadd_phred), as.numeric(b$cadd_phred))
    # missing 'varnm' with an IRanges input
    d <- seqAnnotate(IRanges::IRanges(start=rng[1L], end=rng[2L]),
        favor_gds(), chr="22", verbose=FALSE)
    expect_equal(nrow(d), nrow(a))
})

test_that("variants are dispatched to the file of their chromosome", {
    # a second annotation file on a fake chromosome '21', made from the CSV
    df <- read.csv(ex_file("favor_chr22_sub.csv.gz"))
    df$chromosome <- "21"
    df$variant_vcf <- sub("^22-", "21-", df$variant_vcf)
    csv21 <- tempfile(fileext=".csv")
    write.csv(df, csv21, row.names=FALSE)
    gds21 <- tempfile(fileext=".gds")
    on.exit(unlink(c(csv21, gds21), force=TRUE))
    seqToGDS_FAVOR(csv21, gds21, compress="ZIP", root="", verbose=FALSE)
    v <- .sample_variants(4L)
    # two variants on chr22, two on chr21 and one unknown variant
    snp <- c(paste("22", v$pos[1:2], v$ref[1:2], v$alt[1:2], sep="-"),
        paste("21", v$pos[3:4], v$ref[3:4], v$alt[3:4], sep="-"),
        "21-1-A-T")
    res <- seqAnnotate(snp, c(favor_gds(), gds21), varnm="cadd_phred",
        verbose=FALSE)
    expect_equal(nrow(res), 5L)
    expect_equal(as.numeric(res$cadd_phred[1:4]), as.numeric(v$cadd),
        tolerance=1e-5)
    expect_true(is.na(res$cadd_phred[5L]))
    # the same with the files in the opposite order
    res2 <- seqAnnotate(snp, c(gds21, favor_gds()), varnm="cadd_phred",
        verbose=FALSE)
    expect_equal(as.numeric(res2$cadd_phred), as.numeric(res$cadd_phred))
    # the location of the variants
    loc <- seqAnnotate(snp, c(favor_gds(), gds21), varnm=NULL, verbose=FALSE)
    expect_equal(as.integer(loc$file_idx)[1:4], c(1L, 1L, 2L, 2L))
    expect_true(is.na(loc$variant_idx[5L]))
})

test_that("the 'chr' prefix is ignored when matching the chromosome", {
    v <- .sample_variants(3L)
    a <- seqAnnotate(paste(v$chr, v$pos, v$ref, v$alt, sep="-"), favor_gds(),
        varnm="cadd_phred", verbose=FALSE)
    b <- seqAnnotate(paste0("chr", v$chr, ":", v$pos, ":", v$ref, ":", v$alt),
        favor_gds(), varnm="cadd_phred", verbose=FALSE)
    expect_equal(as.numeric(b$cadd_phred), as.numeric(a$cadd_phred))
    expect_false(anyNA(b$cadd_phred))
})

test_that("annot_gds can be a list of opened GDS files", {
    v <- .sample_variants(3L)
    snp <- paste(v$chr, v$pos, v$ref, v$alt, sep="-")
    f1 <- SeqArray::seqOpen(favor_gds(), allow.duplicate=TRUE)
    f2 <- SeqArray::seqOpen(ex_file("gnomad.genomes.v4.chr22_sub.gds"),
        allow.duplicate=TRUE)
    on.exit({ SeqArray::seqClose(f1); SeqArray::seqClose(f2) })
    a <- seqAnnotate(snp, list(f1), varnm="cadd_phred", verbose=FALSE)
    expect_equal(as.numeric(a$cadd_phred), as.numeric(v$cadd), tolerance=1e-5)
    # both files are on chromosome 22: a variant absent from the first file
    #   is looked up in the second one
    loc <- seqAnnotate(c(snp, "22-10510756-C-A"), list(f2, f1), varnm=NULL,
        verbose=FALSE)
    expect_false(anyNA(loc$variant_idx[1:3]))
    expect_true(all(as.integer(loc$file_idx[1:3]) %in% 1:2))
    # the files are still open afterwards
    expect_s4_class(f1, "SeqVarGDSClass")
    expect_equal(nrow(seqAnnotate(snp, list(f1), varnm="cadd_phred",
        verbose=FALSE)), 3L)
})

test_that("seqAnnotate works on CollapsedVCF and ExpandedVCF objects", {
    skip_if_not_installed("VariantAnnotation")
    v <- .sample_variants(3L)
    hdr <- c("##fileformat=VCFv4.2", "##contig=<ID=22>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    # a third allele, different from the reference and alternative alleles
    other <- setdiff(c("A", "C", "G", "T"), c(v$ref[2L], v$alt[2L]))[1L]
    rows <- c(
        sprintf("%s\t%d\t.\t%s\t%s\t.\tPASS\t.",
            v$chr[1L], v$pos[1L], v$ref[1L], v$alt[1L]),
        # a multi-allelic site gives one row per alternative allele
        sprintf("%s\t%d\t.\t%s\t%s,%s\t.\tPASS\t.",
            v$chr[2L], v$pos[2L], v$ref[2L], v$alt[2L], other))
    vcf <- tempfile(fileext=".vcf")
    on.exit(unlink(vcf, force=TRUE))
    writeLines(c(hdr, rows), vcf)
    obj <- VariantAnnotation::readVcf(vcf)
    expect_s4_class(obj, "CollapsedVCF")
    a <- seqAnnotate(obj, favor_gds(), varnm="cadd_phred", verbose=FALSE)
    expect_s4_class(a, "DataFrame")
    expect_equal(nrow(a), 3L)
    expect_equal(as.numeric(a$cadd_phred[1:2]), as.numeric(v$cadd[1:2]),
        tolerance=1e-5)
    # the same from the expanded VCF object and from the file name
    ex <- VariantAnnotation::expand(obj)
    expect_s4_class(ex, "ExpandedVCF")
    b <- seqAnnotate(ex, favor_gds(), varnm="cadd_phred", verbose=FALSE)
    expect_equal(as.numeric(b$cadd_phred), as.numeric(a$cadd_phred))
    d <- seqAnnotateVCF(vcf, favor_gds(), varnm="cadd_phred", verbose=FALSE)
    expect_equal(as.numeric(d$cadd_phred), as.numeric(a$cadd_phred))
    # 'varnm' can be missing
    expect_equal(ncol(seqAnnotate(obj, favor_gds(), verbose=FALSE)),
        nrow(seqAnnotList(favor_gds())))
})

test_that("a VCF row without alternative allele is skipped", {
    skip_if_not_installed("VariantAnnotation")
    v <- .sample_variants(2L)
    hdr <- c("##fileformat=VCFv4.2", "##contig=<ID=22>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    rows <- c(
        sprintf("%s\t%d\t.\t%s\t.\t.\tPASS\t.",
            v$chr[1L], v$pos[1L], v$ref[1L]),
        sprintf("%s\t%d\t.\t%s\t%s\t.\tPASS\t.",
            v$chr[2L], v$pos[2L], v$ref[2L], v$alt[2L]))
    vcf <- tempfile(fileext=".vcf")
    on.exit(unlink(vcf, force=TRUE))
    writeLines(c(hdr, rows), vcf)
    a <- seqAnnotateVCF(vcf, favor_gds(), varnm="cadd_phred", verbose=FALSE)
    expect_equal(nrow(a), 1L)
    expect_equal(as.numeric(a$cadd_phred), as.numeric(v$cadd[2L]),
        tolerance=1e-5)
})

test_that("an unsupported input class gives an informative error", {
    expect_error(seqAnnotate(1:3, favor_gds(), varnm="cadd_phred",
        verbose=FALSE), "No seqAnnotate\\(\\) method")
    expect_error(seqAnnotate(list(a=1), favor_gds(), verbose=FALSE),
        "class")
})
