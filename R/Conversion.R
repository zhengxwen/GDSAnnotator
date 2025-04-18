#######################################################################
#
# Package Name: GDSAnnotator
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Data import & Format conversion
#

# Package-wide variables
.packageEnv <- new.env()

# Internal functions
.cat <- function(...) cat(..., "\n", sep="")


# create a new SeqArray GDS file
.gds_new <- function(out_fn, compress, var_id_st="int32")
{
    outfile <- createfn.gds(out_fn)
    on.exit(closefn.gds(outfile))
    put.attr.gdsn(outfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(outfile$root, "FileVersion", "v1.0")
    addfolder.gdsn(outfile, "description")
    # add sample.id
    add.gdsn(outfile, "sample.id", integer())
    # add variant.id
    add.gdsn(outfile, "variant.id", storage=var_id_st, compress=compress)
    # add position
    add.gdsn(outfile, "position", storage="int32", compress=compress)
    # add chromosome
    add.gdsn(outfile, "chromosome", storage="string", compress=compress)
    # add allele
    add.gdsn(outfile, "allele", storage="string", compress=compress)
    # add a folder for genotypes
    nd <- addfolder.gdsn(outfile, "genotype")
    put.attr.gdsn(nd, "VariableName", "GT")
    put.attr.gdsn(nd, "Description", "Genotype")
    # add phase folder
    addfolder.gdsn(outfile, "phase")
    # add annotation folder
    nd_annot <- addfolder.gdsn(outfile, "annotation")
    # add id
    add.gdsn(nd_annot, "id", storage="string", compress=compress)
    # add qual
    add.gdsn(nd_annot, "qual", storage="float32", compress=compress)
    # add filter
    nd <- add.gdsn(nd_annot, "filter", storage="int32", compress=compress)
    put.attr.gdsn(nd, "R.class", "factor")
    put.attr.gdsn(nd, "R.levels", "PASS")
    put.attr.gdsn(nd, "Description", "All filters passed")
    # VCF INFO
    addfolder.gdsn(nd_annot, "info")
    # output
    on.exit()
    outfile
}



#######################################################################
# Convert gnomAD VCF to a SeqArray GDS file
#

seqToGDS_gnomAD <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="gnomAD", verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))

    # compression algorithm
    map_compress <- c(LZMA="LZMA_RA", ZIP="ZIP_RA", none="")
    compress1 <- compress2 <- map_compress[compress]
    if (compress=="LZMA") compress1 <- "ZIP_RA"  # reduce memory usage

    # create the gds file
    if (verbose)
    {
        .cat("gnomAD VCF => GDS (", date(), "):")
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }

    # import from VCF
    seqVCF2GDS(vcf_fn, out_fn, storage.option=compress1, verbose=verbose)

    # recompress?
    if (compress1 != compress2)
    {
        if (verbose)
            .cat("Recompressing (", date(), ") ...")
        seqRecompress(out_fn, compress=compress, verbose=verbose)
    }

    if (verbose) .cat("Done (", date(), ")")
    # output
    invisible(normalizePath(out_fn))
}



#######################################################################
# Convert Ensembl-VEP to a SeqArray GDS file
#

.vep_vcf <- function(vcf_fn, out_fn, compress, root="CSQ", verbose=TRUE)
{
    # vcf => gds
    seqVCF2GDS(vcf_fn, out_fn, storage.option=compress, optimize=FALSE,
        verbose=verbose)
    # split CSQ (Consequence annotations from Ensembl VEP)    
    f <- seqOpen(out_fn, readonly=FALSE)
    on.exit(seqClose(f))
    # need CSQ
    nm_root <- paste0("annotation/info/", root)
    nm_root2 <- paste0("annotation/info/", root, ".list")
    nd <- index.gdsn(f, nm_root)
    desp <- get.attr.gdsn(nd)$Description
    if (!is.character(desp))
        stop(root, " information is not found!")
    desp <- gsub("^.*Format:", "", desp)
    # sub fields in CSQ
    nm_lst <- trimws(unlist(strsplit(desp, "|", fixed=TRUE)))
    # new directory
    seqAddValue(f, nm_root2, NULL, verbose=verbose)
    if (verbose)
        .cat("    ", paste(nm_lst, collapse=","))
    csq <- seqGetData(f, nm_root, .tolist=TRUE)
    for (i in seq_along(nm_lst))
    {
        if (verbose)
            cat("    ", nm_lst[i], " ...\n    ", sep="")
        v <- lapply(csq, function(s)
            vapply(strsplit(s, "|", fixed=TRUE), `[`, "", i=i))
        suppressWarnings(
            seqAddValue(f, paste0(nm_root2, "/", nm_lst[i]), v,
                compress=compress, verbose=verbose, verbose.attr=FALSE))
    }
    invisible()
}

seqToGDS_VEP <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="CSQ", verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))

    # compression algorithm
    map_compress <- c(LZMA="LZMA_RA", ZIP="ZIP_RA", none="")
    compress1 <- compress2 <- map_compress[compress]
    if (compress=="LZMA") compress1 <- "ZIP_RA"  # reduce memory usage

    # create the gds file
    if (verbose)
    {
        .cat("Ensembl-VEP VCF => GDS (", date(), "):")
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }

    # load variable names and types
    if (is.null(.packageEnv$vep))
    {
        .packageEnv$vep <- read.csv(system.file(, package="GDSAnnotator", mustWork=TRUE))
    }

    # import from VCF
    .vep_vcf(vcf_fn, out_fn, compress1, root, verbose)

    # recompress?
    if (compress1 != compress2)
    {
        if (verbose)
            .cat("Recompressing (", date(), ") ...")
        seqRecompress(out_fn, compress=compress, verbose=verbose)
    }

    if (verbose) .cat("Done (", date(), ")")
    # output
    invisible(normalizePath(out_fn))
}



#######################################################################
# Convert SnpEff to a SeqArray GDS file
#

.snpeff_vcf <- function(vcf_fn, out_fn, compress, root_lst=c("ANN", "LOF", "NMD"),
    verbose=TRUE)
{
    # vcf => gds
    seqVCF2GDS(vcf_fn, out_fn, storage.option=compress, optimize=FALSE,
        verbose=verbose)
    # split CSQ (Consequence annotations from Ensembl VEP)    
    f <- seqOpen(out_fn, readonly=FALSE)
    on.exit(seqClose(f))
    # for ANN, LOF & NMD
    for (root in root_lst)
    {
        nm_root <- paste0("annotation/info/", root)
        nm_root2 <- paste0("annotation/info/", root, ".list")
        nd <- index.gdsn(f, nm_root)
        desp <- get.attr.gdsn(nd)$Description
        if (!is.character(desp))
            stop(root, " information is not found!")
        # columns
        s <- regmatches(desp, regexpr("'.*'", desp))
        if (!length(s))
            stop("No column names!")
        s <- gsub("'|\\s", "", s)
        s <- gsub("/", "-", s, fixed=TRUE)
        nm_lst <- trimws(unlist(strsplit(s, "|", fixed=TRUE)))
        # new directory
        seqAddValue(f, nm_root2, NULL, verbose=verbose)
        if (verbose)
            .cat("    ", paste(nm_lst, collapse=","))
        ann <- seqGetData(f, nm_root, .tolist=TRUE)
        for (i in seq_along(nm_lst))
        {
            if (verbose)
                cat("    ", nm_lst[i], " ...\n    ", sep="")
            v <- lapply(ann, function(s)
                vapply(strsplit(s, "|", fixed=TRUE), `[`, "", i=i))
            suppressWarnings(
                seqAddValue(f, paste0(nm_root2, "/", nm_lst[i]), v,
                    compress=compress, verbose=verbose, verbose.attr=FALSE))
        }
    }
    invisible()
}

seqToGDS_SnpEff <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root=c("ANN", "LOF", "NMD"), verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root) > 0)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))

    # compression algorithm
    map_compress <- c(LZMA="LZMA_RA", ZIP="ZIP_RA", none="")
    compress1 <- compress2 <- map_compress[compress]
    if (compress=="LZMA") compress1 <- "ZIP_RA"  # reduce memory usage

    # create the gds file
    if (verbose)
    {
        .cat("SnpEff VCF => GDS (", date(), "):")
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }

    # import from VCF
    .snpeff_vcf(vcf_fn, out_fn, compress1, root, verbose)

    # recompress?
    if (compress1 != compress2)
    {
        if (verbose)
            .cat("Recompressing (", date(), ") ...")
        seqRecompress(out_fn, compress=compress, verbose=verbose)
    }

    if (verbose) .cat("Done (", date(), ")")
    # output
    invisible(normalizePath(out_fn))
}

