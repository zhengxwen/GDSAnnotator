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
    stopifnot(inherits(gds, "SeqVarGDSClass"))
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
        # drop=FALSE: keep a DataFrame even when a single annotation is requested
        ans <- ans[ii, , drop=FALSE]
        ans$..no <- is.na(ii)
        ans
    } else {
        # return the file index and variant index
        DataFrame(
            variant_idx = seqGetData(gds, "$variant_index")[ii],
            ..no = is.na(ii))
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
    gds_idx <- rep(seq_along(annot_gds), times=lengths(s))
    names(gds_idx) <- unlist(s)    
    # process each chromosome subset
    ans <- vector("list", length(chr_lst))
    for (i in seq_along(chr_lst))
    {
        ch <- chr_lst[i]
        # find the GDS file(s) for this chromosome
        gds_ii <- which(names(gds_idx) == ch)
        # indices for this chromosome
        idx <- which(chr == ch)
        # process each GDS file for this chromosome
        v <- lapply(seq_along(gds_ii), function(j)
        {
            k <- gds_ii[j]
            if (!length(idx)) return(NULL)
            # annotate variants in this chromosome subset
            d <- ann_pos_allele(annot_gds[[k]],
                ch, pos[idx], ref[idx], alt[idx], varnm, verbose)
            if (!length(varnm)) d$file_idx <- Rle(k, nrow(d))
            # temporary index for combining results in original input order
            d$..idx <- idx
            if (j < length(gds_ii) && any(d$..no))
            {
                # some variants are not found in this file,
                # so they will be processed in the next file(s)
                idx <<- idx[d$..no]
                d <- d[!d$..no, , drop=FALSE]
            }
            d$..no <- NULL  # remove the temporary column
            # return
            d
        })
        # combine results for this chromosome
        ans[[i]] <- if (length(v) == 1L) v[[1L]] else do.call(rbind, v)
    }
    # combine results in original input order
    ans <- if (length(ans) == 1L) ans[[1L]] else do.call(rbind, ans)
    # return
    ans <- ans[order(ans$..idx), , drop=FALSE]
    ans$..idx <- NULL  # remove the temporary index column
    ans
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
    no_sample=TRUE, root="", cleanup=FALSE, bsize=100000L, ..., verbose=TRUE)
{
    # check
    stopifnot(inherits(object, "SeqVarGDSClass"))
    stopifnot(is.numeric(bsize), length(bsize)==1L, !is.na(bsize), bsize>0L)
    bsize <- as.double(bsize)
    stopifnot(is.logical(add_to_gds) || is.character(add_to_gds),
        length(add_to_gds)==1L, !is.na(add_to_gds))
    if (isTRUE(add_to_gds))
    {
        if (isTRUE(object$readonly))
            stop("GDS file should not be readonly when 'add_to_gds=TRUE'.")
        z <- seqSummary(object, "genotype", verbose=FALSE)
        if (z$dim[3L] != z$seldim[3L])
            stop("All variants in GDS should be used when 'add_to_gds=TRUE'.")
    }
    stopifnot(is.logical(no_sample), length(no_sample)==1L)
    stopifnot(is.character(root), length(root)==1L, !is.na(root))
    # check & open annotated gds
    .check_annot_gds(annot_gds)
    if_close_gds <- is.character(annot_gds)
    annot_gds <- .open_annot_gds(annot_gds, verbose)
    if (if_close_gds)
        on.exit(.close_annot_gds(annot_gds))
    # check annotation variable names
    if (missing(varnm))
        varnm <- .annot_list(annot_gds[[1L]])
    if (is.null(varnm)) varnm <- character()
    stopifnot(is.character(varnm))
    # read data 
    if (isTRUE(verbose)) cat("Reading chromosome")
    chr <- seqGetData(object, "$chromosome")
    if (isTRUE(verbose)) cat(", position")
    pos <- seqGetData(object, "position")
    if (isTRUE(verbose)) cat(", reference allele")
    ref <- seqGetData(object, "$ref")
    if (isTRUE(verbose)) cat(", alternative allele")
    alt <- seqGetData(object, "$alt")
    if (isTRUE(verbose)) cat("\n")
    # process
    if (isTRUE(verbose)) cat("Processing annotation ...\n")
    if (isFALSE(add_to_gds) || length(varnm)==0L)
    {
        # return DataFrame
        ann_chr_pos_allele(chr, pos, ref, alt, annot_gds, varnm,
            verbose=verbose)
    } else {
        # when add_to_gds=TRUE or it is a file name
        map <- ann_chr_pos_allele(chr, pos, ref, alt, annot_gds,
            varnm=character(), verbose=verbose)
        remove(chr, pos, ref, alt)
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
                # 'add=TRUE' to keep the handler closing 'annot_gds'
                on.exit(seqFilterPop(object), add=TRUE, after=FALSE)
                fmt.var <- character()
            }
            seqExport(object, add_to_gds, fmt.var=fmt.var, optimize=FALSE,
                verbose=verbose)
            outgds <- seqOpen(add_to_gds, readonly=FALSE)
            on.exit({
                seqClose(outgds)
                if (isTRUE(cleanup))
                    cleanup.gds(add_to_gds, verbose=verbose)
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
        n_total <- nrow(map)
        # group found variants by annotation file, map$file_idx has no NA
        file_ids <- map$file_idx
        if (!inherits(file_ids, "Rle"))
            file_ids <- Rle(as.integer(file_ids))
        # split the variants into blocks to reduce the memory usage: a list
        #   of n_total elements (as returned by .tolist=TRUE) needs much more
        #   memory than the annotation data itself
        nblock <- max(1L, as.integer(ceiling(n_total / bsize)))
        blk_st <- as.double(seq_len(nblock) - 1L) * bsize + 1
        blk_st <- as.integer(pmin(blk_st, n_total + 1))
        blk_cnt <- diff(c(blk_st, n_total + 1L))
        # get the values of the variable 'varnm[i]' for the k-th block, in
        #   the compact list(length, data) representation
        get_block <- function(i, k)
        {
            fi <- file_ids[seq.int(blk_st[k], length.out=blk_cnt[k])]  # Rle
            vi <- map$variant_idx[seq.int(blk_st[k], length.out=blk_cnt[k])]
            rl <- runLength(fi); rv <- runValue(fi)
            # starting offset of each run (1, rl1+1, rl1+rl2+1, ...)
            ost <- cumsum(c(1L, rl))
            v <- lapply(seq_along(rv), function(j)
            {
                f <- annot_gds[[rv[j]]]
                ii <- vi[seq.int(ost[j], length.out=rl[j])]
                ii <- seqSetFilter(f, variant.sel=ii, ret.idx=TRUE,
                    warn=FALSE, verbose=FALSE)$variant_idx
                seqGetData(f, varnm[i], .tolist=TRUE)[ii]
            })
            v <- if (length(rv) == 1L) v[[1L]] else
                unlist(v, recursive=FALSE, use.names=FALSE)
            v <- lapply(v, function(x) unlist(x, use.names=FALSE))
            structure(list(length=lengths(v),
                data=unlist(v, use.names=FALSE)), class="SeqVarDataList")
        }
        # add each annotation variable
        for (i in seq_along(varnm))
        {
            if (verbose)
            {
                cat("    adding ", sQuote(colnm[i]), " (", tm(),
                    ") ...\n", sep="")
            }
            # write block by block, so that only one block is held in memory
            #   instead of the annotation of all variants
            gen <- function(k, i)
            {
                if (k > nblock) return(NULL)  # no more block
                if (verbose && nblock>1L)
                    .cat("        block ", k, "/", nblock, " (", tm(), ")")
                get_block(i, k)
            }
            if (verbose) cat("    ")
            seqAddValue(outgds, paste0(rootnm, "/", colnm[i]), gen,
                replace=TRUE, param=i, verbose=verbose, verbose.attr=FALSE)
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


# Annotate a GDS file with a file name input
seqAnnotateGDS <- function(gds_fn, annot_gds, varnm, add_to_gds=FALSE,
    no_sample=TRUE, root="", cleanup=FALSE, bsize=100000L, ..., verbose=TRUE)
{
    # check
    stopifnot(is.character(gds_fn), length(gds_fn)==1L)
    if (isTRUE(verbose))
        cat("Open", sQuote(gds_fn))
    gds <- seqOpen(gds_fn, readonly=!isTRUE(add_to_gds))
    if (isTRUE(verbose))
    {
        dm <- seqSummary(gds, "genotype", verbose=FALSE)$dim
        cat(" [", prettyNum(dm[3L], big.mark=",", scientific=FALSE),
            " variants]\n", sep="")
    }
    on.exit({
        seqClose(gds)
        if (isTRUE(add_to_gds) && isTRUE(cleanup))
            cleanup.gds(gds_fn, verbose=verbose)
    })
    # process
    ann_gdsfile(gds, annot_gds, varnm, add_to_gds, no_sample, root,
        cleanup=cleanup, bsize=bsize, ..., verbose=verbose)
}

