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

# all the nodes under 'annotation/info' (including the hidden index nodes):
#   names in order, data, storage modes and attributes
.info_nodes <- function(fn)
{
    f <- openfn.gds(fn, allow.duplicate=TRUE)
    on.exit(closefn.gds(f))
    nd <- index.gdsn(f, "annotation/info")
    nm <- ls.gdsn(nd, recursive=TRUE, include.hidden=TRUE, include.dirs=FALSE)
    v <- lapply(nm, function(s)
    {
        n <- index.gdsn(nd, s)
        list(data=read.gdsn(n), storage=objdesp.gdsn(n)$storage,
            attr=get.attr.gdsn(n))
    })
    names(v) <- nm
    v
}

test_that("parallel writes the same nodes as serial processing", {
    snpeff <- ex_file("example_wgs_sites_chr22_snpeff.gds")
    f <- SeqArray::seqOpen(snpeff, allow.duplicate=TRUE)
    chr <- seqGetData(f, "chromosome"); pos <- seqGetData(f, "position")
    ref <- seqGetData(f, "$ref"); alt <- seqGetData(f, "$alt")
    SeqArray::seqClose(f)
    # a variant not in the annotation file, and one on another chromosome
    geno <- make_geno_gds(c(chr, "22", "X"), c(pos, 1L, 100L),
        c(ref, "A", "A"), c(alt, "T", "T"))
    out <- replicate(4L, tempfile(fileext=".gds"))
    on.exit(unlink(c(geno, out), force=TRUE))
    # a float, a flag, lists of strings, integers or floats with up to 25
    #   values per variant, a list with at most one value per variant, and
    #   a factor
    nm <- c("AF", "IMPUTED", "ANN.list/Gene_Name", "ANN.list/Rank",
        "NMD.list/Gene_Name", "LOF.list/Percent_of_transcripts_affected",
        ":filter")
    # to a new file under 'root', with several blocks
    cl <- parallel::makeCluster(2L)
    on.exit(parallel::stopCluster(cl), add=TRUE)
    seqAnnotateGDS(geno, snpeff, varnm=nm, add_to_gds=out[1L], root="ann",
        bsize=7L, verbose=FALSE)
    seqAnnotateGDS(geno, snpeff, varnm=nm, add_to_gds=out[2L], root="ann",
        bsize=7L, parallel=cl, verbose=FALSE)
    a <- .info_nodes(out[1L])
    expect_length(a, 2L * length(nm))
    expect_identical(.info_nodes(out[2L]), a)
    # in place, and again with the existing nodes replaced at their positions
    #   and the variables in another order
    file.copy(geno, out[3:4])
    for (v in list(nm, rev(nm)[-1L]))
    {
        seqAnnotateGDS(out[3L], snpeff, varnm=v, add_to_gds=TRUE, bsize=50L,
            verbose=FALSE)
        seqAnnotateGDS(out[4L], snpeff, varnm=v, add_to_gds=TRUE, bsize=50L,
            parallel=2L, verbose=FALSE)
    }
    a <- .info_nodes(out[3L])
    expect_length(a, 2L * length(nm))
    expect_identical(.info_nodes(out[4L]), a)
    # the processes open the files of a SeqAnnotGDSList object again, and
    #   the files stay open
    x <- seqAnnotOpen(snpeff, allow.duplicate=TRUE)
    on.exit(seqClose(x), add=TRUE)
    out5 <- tempfile(fileext=".gds")
    on.exit(unlink(out5, force=TRUE), add=TRUE)
    seqAnnotateGDS(geno, x, varnm=nm, add_to_gds=out5, root="ann",
        bsize=7L, parallel=cl, verbose=FALSE)
    expect_identical(.info_nodes(out5), .info_nodes(out[1L]))
    expect_s4_class(seqAnnotate("22-1-A-T", x, varnm=NULL, verbose=FALSE),
        "DataFrame")
    # no temporary file left
    expect_length(list.files(tempdir(), "\\.gds\\.(tmp[0-9]+_|map)"), 0L)
    expect_error(seqAnnotateGDS(geno, snpeff, varnm=nm, add_to_gds=out[1L],
        parallel="2", verbose=FALSE), "parallel")
})

test_that("the values are written as stored, and a factor stays a factor", {
    # read the variants and the variables 'nm' of a GDS file
    .read <- function(fn, nm=character())
    {
        f <- SeqArray::seqOpen(fn, allow.duplicate=TRUE)
        on.exit(SeqArray::seqClose(f))
        nm <- c(chr="chromosome", pos="position", ref="$ref", alt="$alt", nm)
        lapply(nm, function(s) seqGetData(f, s, .padNA=FALSE))
    }
    out <- replicate(3L, tempfile(fileext=".gds"))
    on.exit(unlink(out, force=TRUE))
    # SnpEff: all the variants and one variant not found
    snpeff <- ex_file("example_wgs_sites_chr22_snpeff.gds")
    s <- .read(snpeff, c(nmd="annotation/info/NMD.list/Gene_Name"))
    geno <- make_geno_gds(c(s$chr, "22"), c(s$pos, 1L), c(s$ref, "A"),
        c(s$alt, "T"))
    on.exit(unlink(geno, force=TRUE), add=TRUE)
    seqAnnotateGDS(geno, snpeff, varnm=c("AF", "NMD.list/Gene_Name"),
        add_to_gds=out[1L], verbose=FALSE)
    a <- .read(out[1L], c(nmd="annotation/info/NMD.list/Gene_Name",
        af="annotation/info/AF"))
    # a variant without value in the annotation file, or not found, has no
    #   value in a list variable, instead of a missing value
    expect_identical(a$nmd$length, c(s$nmd$length, 0L))
    expect_identical(a$nmd$data, s$nmd$data)
    # one value per variant in a vector variable, missing if not found
    expect_identical(a$af$length, rep(1L, length(s$pos) + 1L))
    expect_true(is.na(a$af$data[length(a$af$data)]))
    # the filters of two files with different levels, and a variant not found
    annot <- c(ex_file("gnomad.genomes.v4.chr22_sub.gds"), favor_gds())
    v <- lapply(annot, .read, nm=c(filter="annotation/filter"))
    i1 <- seq(1L, length(v[[1L]]$pos), by=10L)
    i2 <- seq(5L, length(v[[2L]]$pos), by=50L)
    geno2 <- make_geno_gds("22", c(v[[1L]]$pos[i1], v[[2L]]$pos[i2], 1L),
        c(v[[1L]]$ref[i1], v[[2L]]$ref[i2], "A"),
        c(v[[1L]]$alt[i1], v[[2L]]$alt[i2], "T"))
    on.exit(unlink(geno2, force=TRUE), add=TRUE)
    seqAnnotateGDS(geno2, annot, varnm=":filter", add_to_gds=out[2L],
        verbose=FALSE)
    seqAnnotateGDS(geno2, annot, varnm=":filter", add_to_gds=out[3L],
        bsize=3L, verbose=FALSE)
    expect_identical(.info_nodes(out[3L]), .info_nodes(out[2L]))
    flt <- .read(out[2L], c(filter="annotation/info/.filter"))$filter
    expect_s3_class(flt$data, "factor")
    expect_identical(levels(flt$data),
        union(levels(v[[1L]]$filter), levels(v[[2L]]$filter)))
    # the filter of each variant in the file where it is found
    f <- SeqArray::seqOpen(geno2, allow.duplicate=TRUE)
    loc <- seqAnnotate(f, annot, varnm=NULL, verbose=FALSE)
    SeqArray::seqClose(f)
    fi <- as.integer(loc$file_idx)
    expect_true(all(c(1L, 2L) %in% fi[!is.na(loc$variant_idx)]))
    expect_identical(as.character(flt$data), vapply(seq_along(fi),
        function(k) as.character(v[[fi[k]]]$filter[loc$variant_idx[k]]), ""))
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

test_that("seqAnnotOpen opens the files once for more than one call", {
    v <- .sample_variants(3L)
    snp <- paste(v$chr, v$pos, v$ref, v$alt, sep="-")
    gnomad <- ex_file("gnomad.genomes.v4.chr22_sub.gds")
    a <- seqAnnotate(snp, favor_gds(), varnm="cadd_phred", verbose=FALSE)
    x <- seqAnnotOpen(c(favor_gds(), gnomad))
    on.exit(seqClose(x))
    expect_s4_class(x, "SeqAnnotGDSList")
    expect_length(x, 2L)
    expect_s4_class(x[[1L]], "SeqVarGDSClass")
    m <- S4Vectors::mcols(x)
    expect_equal(m$filename, normalizePath(c(favor_gds(), gnomad)))
    expect_equal(as.list(m$chromosome), list("22", "22"))
    expect_equal(m$n_variant[1L],
        seqSummary(x[[1L]], "genotype", verbose=FALSE)$dim[3L])
    expect_output(show(x), "2 annotation GDS file")
    # the same as with the file names, and the files stay open
    expect_equal(seqAnnotate(snp, x, varnm="cadd_phred", verbose=FALSE), a)
    expect_equal(seqAnnotate(snp, x, varnm="cadd_phred", verbose=FALSE), a)
    # a subset keeps the file information
    y <- x[2L]
    expect_s4_class(y, "SeqAnnotGDSList")
    expect_equal(S4Vectors::mcols(y)$filename, normalizePath(gnomad))
    expect_error(seqAnnotOpen(character()), "No annotation")
    expect_error(seqAnnotOpen(tempfile()), "No file")
    expect_error(seqAnnotate(snp, 1L, varnm="cadd_phred", verbose=FALSE),
        "annot_gds")
})

test_that("the chromosomes of an opened file do not depend on its filter", {
    v <- .sample_variants(2L)
    snp <- paste(v$chr, v$pos, v$ref, v$alt, sep="-")
    f <- SeqArray::seqOpen(favor_gds(), allow.duplicate=TRUE)
    on.exit(SeqArray::seqClose(f))
    seqSetFilter(f, variant.sel=integer(), verbose=FALSE)
    # an opened GDS file instead of a list
    a <- seqAnnotate(snp, f, varnm="cadd_phred", verbose=FALSE)
    expect_equal(as.numeric(a$cadd_phred), as.numeric(v$cadd), tolerance=1e-5)
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
