#######################################################################
#
# Package Name: GDSAnnotator
# Copyright (C) 2025-2026    Xiuwen Zheng
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Data import & Format conversion
#

# Package-wide variables
.packageEnv <- new.env(parent=emptyenv())

# Internal functions
.cat <- function(...) cat(..., "\n", sep="")

tm <- function() strftime(Sys.time(), "%Y-%m-%d %H:%M:%S")


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
    verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))

    # compression algorithm
    map_compress <- c(LZMA="LZMA_RA", ZIP="ZIP_RA", none="")
    compress1 <- compress2 <- map_compress[compress]
    if (compress=="LZMA") compress1 <- "ZIP_RA"  # reduce memory usage

    # create the gds file
    if (verbose)
    {
        .cat("##< ", tm())
        .cat("gnomAD VCF => GDS")
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }

    # import from VCF
    attr(verbose, "header_no_time") <- TRUE
    seqVCF2GDS(vcf_fn, out_fn, storage.option=compress1, verbose=verbose)

    # recompress?
    if (compress1 != compress2)
    {
        if (verbose)
            .cat("Recompressing (", tm(), ") ...")
        seqRecompress(out_fn, compress=compress, verbose=verbose)
    }

    if (verbose) .cat("##> ", tm())
    # output
    invisible(normalizePath(out_fn))
}



#######################################################################
# Internal: split a pipe-delimited annotation field into sub-fields
# using block-based processing to reduce memory usage
#

.split_annot_blocks <- function(f, nm_root, nm_root2, nm_lst, nm_desp,
    compress, bsize, type_fn=NULL, verbose=TRUE)
{
    # create output folder
    seqAddValue(f, nm_root2, NULL, verbose=verbose)
    if (verbose)
        .cat("    ", paste(nm_lst, collapse=","))
    n_fields <- length(nm_lst)
    if (is.null(nm_desp)) nm_desp <- rep("", n_fields)    
    # pre-create empty GDS nodes for each sub-field
    nd_folder <- index.gdsn(f, nm_root2)
    data_nodes <- vector("list", n_fields)
    idx_nodes <- vector("list", n_fields)
    for (i in seq_len(n_fields))
    {
        nm <- nm_lst[i]
        v <- character()
        if (!is.null(type_fn)) v <- type_fn(nm, v)
        tp <- storage.mode(v)
        if (tp=="double") tp <- "float32"
        data_nodes[[i]] <- add.gdsn(nd_folder, nm,
            storage=tp, valdim=0L, compress=compress)
        if (nzchar(nm_desp[i]))
            put.attr.gdsn(data_nodes[[i]], "Description", nm_desp[i])
        idx_nodes[[i]] <- add.gdsn(nd_folder, paste0("@", nm),
            storage="int32", valdim=0L, compress=compress, visible=FALSE)
    }
    # block-by-block processing
    if (verbose) cat("Processing:\n")
    seqBlockApply(f, nm_root, function(bk)
    {
        if (inherits(bk, "SeqVarDataList"))
        {
            ss <- strsplit(bk$data, "|", fixed=TRUE)
            ns <- bk$length
        } else {
            ss <- strsplit(bk, "|", fixed=TRUE)
            ns <- rep(1L, length(bk))
        }
        for (i in seq_len(n_fields))
        {
            # extract sub-field i from each variant's annotation(s)
            v <- vapply(ss, `[`, "", i=i)
            # apply type conversion if provided
            if (!is.null(type_fn)) v <- type_fn(nm_lst[i], v)
            # append data and index
            suppressWarnings(append.gdsn(data_nodes[[i]], v))
            append.gdsn(idx_nodes[[i]], ns)
        }
        NULL  # return 
    }, as.is="none", bsize=bsize, .progress=verbose)
    # finalize all nodes
    for (i in seq_len(n_fields))
    {
        readmode.gdsn(data_nodes[[i]])
        readmode.gdsn(idx_nodes[[i]])
    }
    # return
    invisible()
}



#######################################################################
# Convert Ensembl-VEP to a SeqArray GDS file
#

# VEP type conversion function
.vep_type_fn <- function(nm, v)
{
    # check
    stopifnot(is.character(nm), length(nm)==1L)
    stopifnot(!is.null(.packageEnv$vep))
    # variable names and types
    j <- match(nm, .packageEnv$vep$Name)
    if (!is.na(j))
    {
        lg_set <- c("1", "T", "TRUE", "YES")
        tp <- .packageEnv$vep$Type[j]
        if (is.atomic(v))
        {
            v <- switch(tp,
                integer = as.integer(v),
                numeric = as.numeric(v),
                logical = is.element(v, set=lg_set),
                v)
        } else {
            v <- switch(tp,
                integer = lapply(v, as.integer),
                numeric = lapply(v, as.numeric),
                logical = lapply(v, is.element, set=lg_set),
                v)
        }
    }
    v
}

.vep_vcf <- function(vcf_fn, out_fn, compress, root, bsize, verbose)
{
    # vcf => gds
    attr(verbose, "header_no_time") <- TRUE
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
    # get descriptions for each field
    nm_desp <- vapply(nm_lst, function(nm) {
        j <- match(nm, .packageEnv$vep$Name)
        if (!is.na(j)) .packageEnv$vep$Description[j] else ""
    }, "")
    # split annotation into sub-fields using block processing
    .split_annot_blocks(f, nm_root, nm_root2, nm_lst, nm_desp, compress,
        bsize, type_fn=.vep_type_fn, verbose=verbose)
    invisible()
}

seqToGDS_VEP <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="CSQ", bsize=100000L, verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root)==1L)
    stopifnot(is.numeric(bsize), length(bsize)==1L, bsize>=1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))

    # compression algorithm
    map_compress <- c(LZMA="LZMA_RA", ZIP="ZIP_RA", none="")
    compress1 <- compress2 <- map_compress[compress]
    if (compress=="LZMA") compress1 <- "ZIP_RA"  # reduce memory usage

    # create the gds file
    if (verbose)
    {
        .cat("##< ", tm())
        .cat("Ensembl-VEP VCF => GDS")
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }

    # load variable names and types
    if (is.null(.packageEnv$vep))
    {
        .packageEnv$vep <- read.csv(system.file("extdata",
            "vep_output_format.csv", package="GDSAnnotator", mustWork=TRUE))
        nm <- c("Name", "Description", "Type")
        if (!all(nm %in% colnames(.packageEnv$vep)))
        {
            stop("The internal 'vep_output_format.csv' should have ",
                "the following columns: ", paste(nm, collapse=","), ".")
        }
    }

    # import from VCF
    .vep_vcf(vcf_fn, out_fn, compress1, root, as.integer(bsize), verbose)

    # recompress?
    if (compress1 != compress2)
    {
        if (verbose)
            .cat("Recompressing (", tm(), ") ...")
        # seqRecompress(out_fn, compress=compress, optimize=TRUE,
        #     verbose=verbose)
    } else {
        # cleanup.gds(out_fn, verbose=verbose)
    }

    if (verbose) .cat("##> ", tm())
    # output
    invisible(normalizePath(out_fn))
}



#######################################################################
# Convert SnpEff to a SeqArray GDS file
#

.snpeff_vcf <- function(vcf_fn, out_fn, compress, root_lst, bsize, verbose)
{
    # vcf => gds
    attr(verbose, "header_no_time") <- TRUE
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
        # split annotation into sub-fields using block processing
        .split_annot_blocks(f, nm_root, nm_root2, nm_lst, nm_desp=NULL,
            compress, bsize, type_fn=NULL, verbose=verbose)
    }
    invisible()
}

seqToGDS_SnpEff <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root=c("ANN", "LOF", "NMD"), bsize=100000L, verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root) > 0)
    stopifnot(is.numeric(bsize), length(bsize)==1L, bsize>=1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))

    # compression algorithm
    map_compress <- c(LZMA="LZMA_RA", ZIP="ZIP_RA", none="")
    compress1 <- compress2 <- map_compress[compress]
    if (compress=="LZMA") compress1 <- "ZIP_RA"  # reduce memory usage

    # create the gds file
    if (verbose)
    {
        .cat("##< ", tm())
        .cat("SnpEff VCF => GDS")
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }

    # import from VCF
    .snpeff_vcf(vcf_fn, out_fn, compress1, root, as.integer(bsize), verbose)

    # recompress?
    if (compress1 != compress2)
    {
        if (verbose)
            .cat("Recompressing (", tm(), ") ...")
        seqRecompress(out_fn, compress=compress, optimize=TRUE,
            verbose=verbose)
    } else {
        cleanup.gds(out_fn, verbose=verbose)
    }
    if (verbose) .cat("##> ", tm())
    # output
    invisible(normalizePath(out_fn))
}
