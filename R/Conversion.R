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
    nm_uniform, compress, bsize, type_fn=NULL, verbose=TRUE)
{
    # create output folder
    seqAddValue(f, nm_root2, NULL, verbose=verbose)
    if (verbose)
        .cat("    ", paste(nm_lst, collapse=","))
    n_fields <- length(nm_lst)
    if (is.null(nm_desp)) nm_desp <- rep("", n_fields)
    if (is.null(nm_uniform)) nm_uniform <- rep(FALSE, n_fields)
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
        if (!nm_uniform[i])
        {
            idx_nodes[[i]] <- add.gdsn(nd_folder, paste0("@", nm),
                storage="int32", valdim=0L, compress=compress, visible=FALSE)
        }
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
            if (nm_uniform[i])
            {
                # uniform: one value per variant (take first of each group)
                idx <- cumsum(c(1L, head(ns, -1L)))
                suppressWarnings(append.gdsn(data_nodes[[i]], v[idx]))
            } else {
                # variable-length: append all data and index
                suppressWarnings(append.gdsn(data_nodes[[i]], v))
                append.gdsn(idx_nodes[[i]], ns)
            }
        }
        NULL  # return 
    }, as.is="none", bsize=bsize, .progress=verbose)
    # finalize all nodes
    for (i in seq_len(n_fields))
    {
        readmode.gdsn(data_nodes[[i]])
        if (!nm_uniform[i])
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

.vep_vcf <- function(vcf_fn, out_fn, compress, root, keep, bsize, verbose)
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
    # get uniform flags for each field
    nm_uniform <- vapply(nm_lst, function(nm) {
        j <- match(nm, .packageEnv$vep$Name)
        if (!is.na(j)) isTRUE(.packageEnv$vep$Uniform[j]) else FALSE
    }, FALSE)
    # split annotation into sub-fields using block processing
    .split_annot_blocks(f, nm_root, nm_root2, nm_lst, nm_desp, nm_uniform,
        compress, bsize, type_fn=.vep_type_fn, verbose=verbose)
    # remove the original root node if keep=FALSE
    if (isFALSE(keep))
    {
        delete.gdsn(index.gdsn(f, nm_root))
        if (verbose)
            .cat("    removed '", nm_root, "'")
    }
    invisible()
}

seqToGDS_VEP <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="CSQ", keep=TRUE, bsize=100000L, verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root)==1L)
    stopifnot(is.logical(keep), length(keep)==1L)
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
    .vep_vcf(vcf_fn, out_fn, compress1, root, keep, as.integer(bsize), verbose)

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



#######################################################################
# Convert SnpEff to a SeqArray GDS file
#

# SnpEff type conversion function
.snpeff_type_fn <- function(nm, v)
{
    # check
    stopifnot(is.character(nm), length(nm)==1L)
    stopifnot(!is.null(.packageEnv$snpeff_sub))
    # variable names and types
    j <- match(nm, .packageEnv$snpeff_sub$Name)
    if (!is.na(j))
    {
        tp <- .packageEnv$snpeff_sub$Type[j]
        if (is.atomic(v))
        {
            v <- switch(tp,
                integer = as.integer(v),
                numeric = as.numeric(v),
                v)
        } else {
            v <- switch(tp,
                integer = lapply(v, as.integer),
                numeric = lapply(v, as.numeric),
                v)
        }
    }
    v
}

.snpeff_vcf <- function(vcf_fn, out_fn, compress, root_lst, keep, bsize,
    verbose)
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
        # look up descriptions, types and uniform flags from CSV
        .packageEnv$snpeff_sub <- .packageEnv$snpeff[
            .packageEnv$snpeff$Field == root, , drop=FALSE]
        nm_desp <- vapply(nm_lst, function(nm) {
            j <- match(nm, .packageEnv$snpeff_sub$Name)
            if (!is.na(j)) .packageEnv$snpeff_sub$Description[j] else ""
        }, "")
        nm_uniform <- vapply(nm_lst, function(nm) {
            j <- match(nm, .packageEnv$snpeff_sub$Name)
            if (!is.na(j)) isTRUE(.packageEnv$snpeff_sub$Uniform[j]) else FALSE
        }, FALSE)
        # determine type_fn: use it if any field has a Type defined
        type_fn <- NULL
        if (any(nzchar(.packageEnv$snpeff_sub$Type)))
            type_fn <- .snpeff_type_fn
        # split annotation into sub-fields using block processing
        .split_annot_blocks(f, nm_root, nm_root2, nm_lst, nm_desp,
            nm_uniform, compress, bsize, type_fn=type_fn, verbose=verbose)
        # remove the original root node if keep=FALSE
        if (isFALSE(keep))
        {
            delete.gdsn(index.gdsn(f, nm_root))
            if (verbose)
                .cat("    removed '", nm_root, "'")
        }
    }
    invisible()
}

seqToGDS_SnpEff <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root=c("ANN", "LOF", "NMD"), keep=TRUE, bsize=100000L, verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root) > 0)
    stopifnot(is.logical(keep), length(keep)==1L)
    stopifnot(is.numeric(bsize), length(bsize)==1L, bsize>=1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))

    # compression algorithm
    map_compress <- c(LZMA="LZMA_RA", ZIP="ZIP_RA", none="")
    compress1 <- compress2 <- map_compress[compress]
    if (compress=="LZMA") compress1 <- "ZIP_RA"  # reduce memory usage

    # load variable names and types
    if (is.null(.packageEnv$snpeff))
    {
        .packageEnv$snpeff <- read.csv(system.file("extdata",
            "snpeff_output_format.csv", package="GDSAnnotator",
            mustWork=TRUE))
        nm <- c("Field", "Name", "Description", "Type")
        if (!all(nm %in% colnames(.packageEnv$snpeff)))
        {
            stop("The internal 'snpeff_output_format.csv' should have ",
                "the following columns: ", paste(nm, collapse=","), ".")
        }
    }

    # create the gds file
    if (verbose)
    {
        .cat("##< ", tm())
        .cat("SnpEff VCF => GDS")
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }

    # import from VCF
    .snpeff_vcf(vcf_fn, out_fn, compress1, root, keep, as.integer(bsize),
        verbose)

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



#######################################################################
# Convert ANNOVAR-annotated VCF to a SeqArray GDS file
#

# ANNOVAR type conversion function
.annovar_type_fn <- function(nm, v)
{
    # check
    stopifnot(is.character(nm), length(nm)==1L)
    stopifnot(!is.null(.packageEnv$annovar))
    # variable names and types
    j <- match(nm, .packageEnv$annovar$Name)
    if (!is.na(j))
    {
        tp <- .packageEnv$annovar$Type[j]
        if (nzchar(tp))
        {
            v <- switch(tp,
                integer = as.integer(v),
                numeric = as.numeric(v),
                v)
        }
    }
    v
}

# Decode ANNOVAR \xHH escape sequences in character vectors
.annovar_decode_hex <- function(x)
{
    idx <- grep("\\\\x[0-9a-fA-F]{2}", x)
    if (length(idx))
    {
        x[idx] <- vapply(x[idx], function(s) {
            while (grepl("\\\\x[0-9a-fA-F]{2}", s))
            {
                m <- regexpr("\\\\x[0-9a-fA-F]{2}", s)
                hex <- substring(s, m + 2, m + 3)
                ch <- rawToChar(as.raw(strtoi(hex, base=16L)))
                s <- sub("\\\\x[0-9a-fA-F]{2}", ch, s, fixed=FALSE)
            }
            s
        }, "", USE.NAMES=FALSE)
    }
    x
}

.annovar_vcf <- function(vcf_fn, out_fn, compress, keep, verbose)
{
    # vcf => gds
    attr(verbose, "header_no_time") <- TRUE
    seqVCF2GDS(vcf_fn, out_fn, storage.option=compress, optimize=FALSE,
        verbose=verbose)
    # post-process: remove ANNOVAR bookkeeping nodes
    f <- seqOpen(out_fn, readonly=FALSE)
    on.exit(seqClose(f))
    nd_info <- index.gdsn(f, "annotation/info")
    nm_all <- ls.gdsn(nd_info)
    # remove ANNOVAR_DATE and ALLELE_END pseudo-nodes
    for (nm in c("ANNOVAR_DATE", "ALLELE_END"))
    {
        if (nm %in% nm_all)
        {
            # also remove the index node if present
            idx_nm <- paste0("@", nm)
            if (idx_nm %in% nm_all)
                delete.gdsn(index.gdsn(nd_info, idx_nm))
            delete.gdsn(index.gdsn(nd_info, nm))
            if (verbose)
                .cat("    removed '", nm, "'")
        }
    }
    # apply type conversion and add descriptions for known ANNOVAR fields
    # also decode \xHH escape sequences in string fields
    nm_all <- ls.gdsn(nd_info)
    annovar_db <- .packageEnv$annovar
    for (nm in nm_all)
    {
        if (startsWith(nm, "@")) next
        nd <- index.gdsn(nd_info, nm)
        dp <- objdesp.gdsn(nd)
        # decode \xHH in character/string nodes
        if (dp$type == "String" || dp$type == "Factor")
        {
            v <- read.gdsn(nd)
            if (any(grepl("\\\\x[0-9a-fA-F]{2}", v)))
            {
                v <- .annovar_decode_hex(v)
                # rewrite the node
                compress_nd <- dp$compress
                delete.gdsn(nd)
                nd <- add.gdsn(nd_info, nm, v, storage="string",
                    compress=compress_nd, replace=TRUE)
                if (verbose)
                    .cat("    decoded '", nm, "'")
            }
        }
        j <- match(nm, annovar_db$Name)
        if (!is.na(j))
        {
            nd <- index.gdsn(nd_info, nm)
            # add description if available
            desp <- annovar_db$Description[j]
            if (nzchar(desp))
            {
                a <- get.attr.gdsn(nd)
                if (is.null(a$Description) || !nzchar(a$Description))
                    put.attr.gdsn(nd, "Description", desp)
            }
        }
    }
    # remove raw nodes if keep=FALSE
    if (isFALSE(keep))
    {
        # remove ANNOVAR-specific annotation nodes (Func.*, Gene.*, etc.)
        nm_all <- ls.gdsn(nd_info)
        annovar_prefixes <- c("Func.", "Gene.", "GeneDetail.",
            "ExonicFunc.", "AAChange.")
        for (nm in nm_all)
        {
            if (startsWith(nm, "@")) next
            is_annovar <- any(vapply(annovar_prefixes,
                function(p) startsWith(nm, p), FALSE))
            if (is_annovar)
            {
                idx_nm <- paste0("@", nm)
                if (idx_nm %in% nm_all)
                    delete.gdsn(index.gdsn(nd_info, idx_nm))
                delete.gdsn(index.gdsn(nd_info, nm))
                if (verbose)
                    .cat("    removed '", nm, "'")
            }
        }
    }
    invisible()
}

seqToGDS_ANNOVAR <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    keep=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.logical(keep), length(keep)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))

    # compression algorithm
    map_compress <- c(LZMA="LZMA_RA", ZIP="ZIP_RA", none="")
    compress1 <- compress2 <- map_compress[compress]
    if (compress=="LZMA") compress1 <- "ZIP_RA"  # reduce memory usage

    # load variable names and types
    if (is.null(.packageEnv$annovar))
    {
        .packageEnv$annovar <- read.csv(system.file("extdata",
            "annovar_output_format.csv", package="GDSAnnotator",
            mustWork=TRUE))
        nm <- c("Name", "Description", "Type")
        if (!all(nm %in% colnames(.packageEnv$annovar)))
        {
            stop("The internal 'annovar_output_format.csv' should have ",
                "the following columns: ", paste(nm, collapse=","), ".")
        }
    }

    # create the gds file
    if (verbose)
    {
        .cat("##< ", tm())
        .cat("ANNOVAR VCF => GDS")
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }

    # import from VCF
    .annovar_vcf(vcf_fn, out_fn, compress1, keep, verbose)

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
