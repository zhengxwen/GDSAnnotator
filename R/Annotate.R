#######################################################################
#
# Package Name: GDSAnnotator
# Copyright (C) 2025-2026    Xiuwen Zheng
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Annotate variants
#


# Define new generic functions
setGeneric("seqAnnotate", function(object, annot_gds, varnm, ..., verbose=TRUE)
        standardGeneric("seqAnnotate"))


# Open the GDS file(s) with variant annotation
.open_annot_gds <- function(gds_fn, verbose=TRUE)
{
    ans <- gds_fn
    if (is.character(gds_fn))
    {
        # check
        stopifnot(length(gds_fn) > 0L)
        # open the GDS file(s)
        ans <- vector("list", length(gds_fn))
        on.exit({ for (f in ans) if (!is.null(f)) seqClose(f) })
        nm <- character(length(gds_fn))
        for (i in seq_along(gds_fn))
        {
            if (isTRUE(verbose))
                cat("Open", sQuote(basename(gds_fn[i])))
            # open the file
            f <- ans[[i]] <- seqOpen(gds_fn[i])
            # check chromosome
            v <- seqGetData(f, "$chromosome")
            if (isTRUE(verbose))
            {
                .cat(" [",
                    prettyNum(length(v), big.mark=",", scientific=FALSE),
                    " variants]")
            }
            if (nrun(v) > 1L)
                stop(gds_fn[i], " should only contain one chromosome.")
            nm[i] <- paste0("chr",
                paste(as.character(unique(v)), collapse="&"))
        }
        names(ans) <- nm
        # output
        on.exit()
    }
    ans
}

# Close the GDS file(s)
.close_annot_gds <- function(annot_gds)
{
    lapply(annot_gds, seqClose)
    invisible()
}

# Check annotation gds list
.check_annot_gds <- function(annot_gds)
{
    if (is.character(annot_gds))
    {
        if (length(annot_gds) == 0L)
            stop("No annotation gds file!")
        fn <- annot_gds[!file.exists(annot_gds)]
        if (length(fn))
            stop("No file(s): ", paste(fn, collapse=", "), ".")
    } else if (is.list(annot_gds))
    {
        if (length(annot_gds) == 0L)
            stop("No annotation gds file!")
        for (i in seq_along(annot_gds))
        {
            if (!inherits(annot_gds[[i]], "SeqVarGDSClass"))
                stop("annot_gds[[", i, "]] is not a GDS object.")
        }
    }
    invisible()
}

# GDS nodes for the INFO field
.gds_varnm_cvt <- c(
    ':chromosome'="chromosome", ':position'="position", ':allele'="allele",
    ':id'="annotation/id", ':qual'="annotation/qual",
    ':filter'="annotation/filter")

.gds_varnm <- function(nm)
{
    if (length(nm))
    {
        if (anyNA(nm))
            stop("'varnm' should not have NA_character_.")
        i <- match(nm, names(.gds_varnm_cvt))
        if (anyNA(i))
            nm[is.na(i)] <- paste0("annotation/info/", nm[is.na(i)])
        j <- which(!is.na(i))
        if (length(j))
            nm[j] <- .gds_varnm_cvt[i[j]]
    }
    nm
}

# get all annotation names in the INFO field
.annot_list <- function(gdsfile)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # process
    nd_info <- index.gdsn(gdsfile, "annotation/info")
    nm <- ls.gdsn(nd_info, recursive=TRUE, include.dirs=FALSE)
}


# Helper function for ann_chr_pos_allele
ann_pos_allele <- function(gds, chr, pos, ref, alt, varnm, verbose=TRUE)
{
    # check
    stopifnot(is.integer(pos))
    stopifnot(is.character(ref))
    stopifnot(is.character(alt))
    stopifnot(length(pos)==length(ref))
    stopifnot(length(pos)==length(alt))
    # GDS node names
    colnm <- names(varnm)
    if (is.null(colnm)) colnm <- varnm
    varnm <- .gds_varnm(varnm)
    # process
    # find matching variants using seqSetFilterPos
    ii <- seqSetFilterPos(gds, chr, pos, ref=ref, alt=alt,
        ret.idx=TRUE, verbose=FALSE)
    # get
    if (length(varnm))
    {
        if (verbose)
        {
            cat("[", basename(gds$filename), "] ", sep="")
            n <- seqSummary(gds, "genotype", verbose=FALSE)$seldim[3L]
            .cat(" # of variants found: ", n)
        }
        l_verbose <- isTRUE(verbose) && (length(varnm)>1L)
        if (l_verbose)
            cat("[", length(varnm), "] ", sep="")
        ans <- lapply(varnm, function(nm)
        {
            if (l_verbose) cat(".")
            seqGetData(gds, nm, .tolist=NA)
        })
        names(ans) <- colnm
        ans <- DataFrame(ans)
        if (l_verbose) cat("\n")
        # ii maps each input to the filtered set row (NA = not found)
        ans[ii, ]
    } else {
        # return the file index and variant index
        DataFrame(variant_idx = seqGetData(gds, "$variant_index")[ii])
    }
}


# Annotate with chromosome, position, reference & alternative alleles
#   dispatch per chromosome, combine results in original input order
ann_chr_pos_allele <- function(chr, pos, ref, alt, annot_gds, varnm,
    verbose=TRUE)
{
    # check
    stopifnot(length(chr)==1L || length(chr)==length(pos))
    stopifnot(length(pos) == length(ref))
    if (!is.null(alt))
        stopifnot(length(pos) == length(alt))
    stopifnot(is.logical(verbose), length(verbose)==1L)
    # check & process
    stopifnot(is.character(varnm))
    # process pos & each chromosome
    if (!is.integer(pos)) pos <- as.integer(pos)
    chr_lst <- unique(chr)
    if (length(chr) == 1L) chr <- rep(chr, length(pos))
    # match chromosome names in the annotation GDS files
    s <- strsplit(gsub("^chr", "", names(annot_gds)), "&", fixed=TRUE)
    gds_idx <- rep(seq_along(annot_gds), each=lengths(s))
    names(gds_idx) <- unlist(s)    
    # process each chromosome subset
    ans_lst <- vector("list", length(chr_lst))
    for (i in seq_along(chr_lst))
    {
        ch <- chr_lst[i]
        # find the GDS file(s) for this chromosome
        gds_ii <- which(names(gds_idx) == ch)
        j <- which(chr == ch)
        for (k in gds_ii)
        {
            d <- ann_pos_allele(annot_gds[[k]], ch,
                pos[j], ref[j], alt[j], varnm, verbose)
            if (length(varnm))
                d$..idx <- j
            else
                d$file_idx <- Rle(k, nrow(d))
            ans_lst[[i]] <- d
        }
    }
    # combine results in original input order
    if (length(chr_lst) == 1L)
    {
        ans_lst[[1L]]
    } else {
        ans <- do.call(rbind, ans_lst)
        remove(ans_lst)
        ans <- ans[order(ans$..idx), ]
        ans$..idx <- NULL  # remove the temporary index column
        ans
    }
}


# Annotate with a data frame including the columns for
#   chromosome, position, reference & alternative alleles
ann_dataframe <- function(object, annot_gds, varnm, col_chr="chr",
    col_pos="pos", col_ref="ref", col_alt="alt", ..., verbose=TRUE)
{
    # check
    stopifnot(is.character(varnm))
    chr <- object[[col_chr]]
    if (is.null(chr)) stop("No 'chr' column.")
    pos <- object[[col_pos]]
    if (is.null(pos)) stop("No 'pos' column.")
    ref <- object[[col_ref]]
    if (is.null(ref)) stop("No 'ref' column.")
    alt <- object[[col_alt]]
    # check & open annotated gds
    .check_annot_gds(annot_gds)
    if_close_gds <- is.character(annot_gds)
    annot_gds <- .open_annot_gds(annot_gds, verbose)
    if (if_close_gds)
        on.exit(.close_annot_gds(annot_gds))
    if (missing(varnm))
        varnm <- .annot_list(annot_gds[[1L]])
    if (is.null(varnm)) varnm <- character()
    # process
    ann_chr_pos_allele(chr, pos, ref, alt, annot_gds, varnm, verbose=verbose)
}


# Annotate a GDS file
ann_gdsfile <- function(object, annot_gds, varnm, add_to_gds=FALSE,
    no_sample=TRUE, root="", ..., verbose=TRUE)
{
    # check
    stopifnot(inherits(object, "SeqVarGDSClass"))
    stopifnot(is.logical(add_to_gds) || is.character(add_to_gds),
        length(add_to_gds)==1L, !is.na(add_to_gds))
    if (isTRUE(add_to_gds))
    {
        if (isTRUE(object$readonly))
            stop("GDS file should not be readonly when 'add_to_gds=TRUE'.")
        z <- seqSummary(object, "genotype", verbose=TRUE)
        if (z$dim[3L] != z$seldim[3L])
            stop("All variants in GDS should be used when 'add_to_gds=TRUE'.")
    }
    stopifnot(is.logical(no_sample), length(no_sample)==1L)
    stopifnot(is.character(root), length(root)==1L, !is.na(root))
    # check & process
    if (missing(varnm))
        varnm <- .annot_list(annot_gds[[1L]])
    if (is.null(varnm)) varnm <- character()
    stopifnot(is.character(varnm))
    # check & open annotated gds
    .check_annot_gds(annot_gds)
    if_close_gds <- is.character(annot_gds)
    annot_gds <- .open_annot_gds(annot_gds, verbose)
    if (if_close_gds)
        on.exit(.close_annot_gds(annot_gds))
    # process
    chr <- seqGetData(object, "$chromosome")
    if (nrun(chr) > 1L)
        stop(basename(object$filename), " should only contain one chromosome.")
    pos <- seqGetData(object, "position")
    allele <- seqGetData(object, "allele")
    ref <- sub(",.*", "", allele)
    alt <- ifelse(grepl(",", allele, fixed=TRUE),
        sub("^[^,]*,", "", allele), NA_character_)
    if (isFALSE(add_to_gds) || length(varnm)==0L)
    {
        # return DataFrame
        ann_chr_pos_allele(chr, pos, ref, alt, annot_gds, varnm,
            verbose=verbose)
    } else {
        # when add_to_gds=TRUE or it is a file name
        # sort position
        idx <- NULL
        if (is.unsorted(pos) || pos[1L] > pos[length(pos)])
        {
            idx <- order(pos)
            pos <- pos[idx]; ref <- ref[idx]; alt <- alt[idx]
            idx <- order(idx)
        }
        map <- ann_chr_pos_allele(chr, pos, ref, alt, annot_gds, varnm="",
            verbose=verbose)
        remove(chr, pos, ref, alt)
        # check
        if (anyNA(map$file_idx) || length(unique(map$file_idx))>1L)
        {
            stop("The GDS file should have only one chromosome.")
        }
        if (verbose)
        {
            n <- sum(is.na(map$variant_idx))
            .cat("Found #: ", nrow(map)-n, " variant(s), Not found #: ", n)
        }
        # output gds file
        outgds <- object
        if (is.character(add_to_gds))
        {
            # add_to_gds is character
            fmt.var <- NULL
            if (isTRUE(no_sample))
            {
                seqSetFilter(object, sample.sel=integer(), action="push+set",
                    verbose=FALSE)
                on.exit(seqFilterPop(object))
                fmt.var <- character()
            }
            seqExport(object, add_to_gds, fmt.var=fmt.var, optimize=FALSE,
                verbose=verbose)
            outgds <- seqOpen(add_to_gds, readonly=FALSE)
            on.exit({
                if (!is.null(outgds)) seqClose(outgds)
            }, add=TRUE)
        }
        if (verbose)
            .cat("Add new data to ", sQuote(basename(outgds$filename)), ":")
        rootnm <- "annotation/info"
        if (root != "")
        {
            rootnm <- paste0(rootnm, "/", root)
            seqAddValue(outgds, rootnm, NULL, replace=TRUE, verbose=verbose)
        }
        # GDS nodes
        colnm <- gsub("^\\:", ".", varnm)
        varnm <- .gds_varnm(varnm)
        f <- annot_gds[[as.integer(map$file_idx[1])]]
        if (verbose)
            cat("    [", basename(f$filename), "] ", sep="")
        ii <- seqSetFilter(f, variant.sel=map$variant_idx, ret.idx=TRUE,
            warn=FALSE, verbose=verbose)$variant_idx
        rerow <- anyNA(ii) || is.unsorted(ii) || ii[1L] > ii[length(ii)]
        if (rerow && !is.null(idx)) ii <- ii[idx]
        if (!rerow && !is.null(idx)) { rerow <- TRUE; ii <- idx }
        # add each node
        for (i in seq_along(varnm))
        {
            if (verbose)
            {
                cat("    adding ", sQuote(colnm[i]), " (", tm(), ") ...\n    ",
                    sep="")
            }
            # set a filter to get data
            v <- seqGetData(f, varnm[i], .tolist=TRUE)
            if (rerow) v <- v[ii]
            seqAddValue(outgds, paste0(rootnm, "/", colnm[i]), v,
                replace=TRUE, verbose=verbose, verbose.attr=FALSE)
            remove(v)
        }
        # optimize ...
        if (is.character(add_to_gds))
        {
            if (verbose)
                cat("Optimize the access efficiency ...\n")
            seqClose(outgds)
            outgds <- NULL
            cleanup.gds(add_to_gds, verbose=verbose)
        }
        # return
        invisible()
    }
}


# Annotate "chr-pos-ref-alt"
ann_variant <- function(object, annot_gds, varnm, split="[-_:]", ...,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(object), length(object)>0L)
    ss <- strsplit(object, split)
    ns <- lengths(ss)
    if (any(ns != 4L))
        stop("Invalid input: ", object[which(ns != 4L)[1L]])
    chr <- vapply(ss, `[`, "", i=1L)
    pos <- as.integer(vapply(ss, `[`, "", i=2L))
    if (anyNA(pos))
        stop("Invalid position: ", object[which(is.na(pos))[1L]])
    ref <- vapply(ss, `[`, "", i=3L)
    alt <- vapply(ss, `[`, "", i=4L)
    # check & open annotated gds
    .check_annot_gds(annot_gds)
    if_close_gds <- is.character(annot_gds)
    annot_gds <- .open_annot_gds(annot_gds, verbose)
    if (if_close_gds)
        on.exit(.close_annot_gds(annot_gds))
    if (missing(varnm))
        varnm <- .annot_list(annot_gds[[1L]])
    if (is.null(varnm)) varnm <- character()
    stopifnot(is.character(varnm))
    # process
    ann_chr_pos_allele(chr, pos, ref, alt, annot_gds, varnm, verbose=verbose)
}


# Annotate a GRanges/GRangesList object
ann_GRanges <- function(object, annot_gds, varnm, ..., verbose=TRUE)
{
    # check & open annotated gds
    .check_annot_gds(annot_gds)
    if_close_gds <- is.character(annot_gds)
    annot_gds <- .open_annot_gds(annot_gds, verbose)
    if (if_close_gds)
        on.exit(.close_annot_gds(annot_gds))
    if (missing(varnm))
        varnm <- .annot_list(annot_gds[[1L]])
    if (is.null(varnm)) varnm <- character()
    stopifnot(is.character(varnm))
    # verbose
    l_verbose <- isTRUE(verbose) && (length(varnm)>1L)
    if (l_verbose)
        cat("[", length(varnm), "] ", sep="")
    # GDS node names
    colnm <- names(varnm)
    if (is.null(colnm)) colnm <- varnm
    varnm <- .gds_varnm(varnm)
    # process
    ans <- lapply(annot_gds, function(gds)
    {
        seqSetFilter(gds, object, verbose=FALSE)
        if (seqSummary(gds, "genotype", verbose=FALSE)$seldim[3L])
        {
            DataFrame(lapply(varnm, function(nm)
            {
                if (l_verbose) cat(".")
                seqGetData(gds, nm, .tolist=NA)
            }))
        } else
            NULL
    })
    if (l_verbose) cat("\n")
    # output
    ans <- do.call(rbind, ans)
    if (!is.null(ans)) names(ans) <- colnm
    ans
}


# Annotate an IRanges object
ann_IRanges <- function(object, annot_gds, varnm, chr, ..., verbose=TRUE)
{
    # check
    stopifnot(is.character(varnm))
    # check & open annotated gds
    .check_annot_gds(annot_gds)
    if_close_gds <- is.character(annot_gds)
    annot_gds <- .open_annot_gds(annot_gds, verbose)
    if (if_close_gds)
        on.exit(.close_annot_gds(annot_gds))
    if (missing(varnm))
        varnm <- .annot_list(annot_gds[[1L]])
    if (is.null(varnm)) varnm <- character()
    # verbose
    l_verbose <- isTRUE(verbose) && (length(varnm)>1L)
    if (l_verbose)
        cat("[", length(varnm), "] ", sep="")
    # GDS node names
    colnm <- names(varnm)
    if (is.null(colnm)) colnm <- varnm
    varnm <- .gds_varnm(varnm)
    # process
    ans <- lapply(annot_gds, function(gds)
    {
        seqSetFilter(gds, object, chr=chr, verbose=FALSE)
        if (seqSummary(gds, "genotype", verbose=FALSE)$seldim[3L])
        {
            DataFrame(lapply(varnm, function(nm)
            {
                if (l_verbose) cat(".")
                seqGetData(gds, nm, .tolist=NA)
            }))
        } else
            NULL
    })
    if (l_verbose) cat("\n")
    # output
    ans <- do.call(rbind, ans)
    if (!is.null(ans)) names(ans) <- colnm
    ans
}


# Set methods
setMethod("seqAnnotate", signature(object="data.frame"), ann_dataframe)
setMethod("seqAnnotate", signature(object="DataFrame"), ann_dataframe)
setMethod("seqAnnotate", signature(object="SeqVarGDSClass"), ann_gdsfile)
setMethod("seqAnnotate", signature(object="character"), ann_variant)
setMethod("seqAnnotate", signature(object="GRanges"), ann_GRanges)
setMethod("seqAnnotate", signature(object="GRangesList"), ann_GRanges)
setMethod("seqAnnotate", signature(object="IRanges"), ann_IRanges)


# Annotate a GDS file with a file name input
seqAnnotateFile <- function(gds_fn, annot_gds, varnm, add_to_gds=FALSE,
    no_sample=TRUE, root="", ..., verbose=TRUE)
{
    # check
    stopifnot(is.character(gds_fn), length(gds_fn)==1L)
    if (isTRUE(verbose)) .cat("Open ", sQuote(gds_fn))
    gds <- seqOpen(gds_fn)
    on.exit(seqClose(gds))
    # process
    ann_gdsfile(gds, annot_gds, varnm, add_to_gds, no_sample, root,
        ..., verbose=verbose)
}


# Annotate a VCF file
seqAnnotateVCF <- function(vcf_fn, annot_gds, varnm, ..., verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)==1L)
    if (!requireNamespace("VariantAnnotation", quietly=TRUE))
        stop("The package 'VariantAnnotation' should be installed.")
    if (!requireNamespace("SummarizedExperiment", quietly=TRUE))
        stop("The package 'SummarizedExperiment' should be installed.")
    # open the file
    if (isTRUE(verbose)) .cat("Open ", sQuote(vcf_fn))
    vcf <- VariantAnnotation::readVcf(vcf_fn)
    gr <- SummarizedExperiment::rowRanges(vcf)
    # output
    seqAnnotate(gr, annot_gds, varnm, ..., verbose=verbose)
}

