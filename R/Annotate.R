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
#   object: the variants to be annotated, see the methods below
#   annot_gds: the annotation GDS file names or a list of opened GDS files
#   varnm: the annotation names, see .gds_varnm()
#   verbose: if TRUE, show the progress information
setGeneric("seqAnnotate", function(object, annot_gds, varnm, ..., verbose=TRUE)
        standardGeneric("seqAnnotate"))


# Normalize chromosome names for matching, ignoring the "chr" prefix
#   (e.g., "chr22" in a VCF file vs "22" in a GDS file)
.norm_chr <- function(x) sub("^chr", "", as.character(x), ignore.case=TRUE)

# Open the GDS file(s) with variant annotation; the returned list is named
#   by the chromosome of each file ("chr" + the chromosome label used in
#   the file), for dispatching the query variants to the files
.open_annot_gds <- function(gds_fn, verbose=TRUE)
{
    ans <- gds_fn
    if (is.character(gds_fn))
    {
        # check
        stopifnot(length(gds_fn) > 0L)
        # open the GDS file(s), which are closed if an error occurs
        ans <- vector("list", length(gds_fn))
        on.exit({ for (f in ans) if (!is.null(f)) seqClose(f) })
        for (i in seq_along(gds_fn))
        {
            if (isTRUE(verbose))
                message("Open ", sQuote(basename(gds_fn[i])), appendLF=FALSE)
            # open the file
            f <- ans[[i]] <- seqOpen(gds_fn[i])
            if (isTRUE(verbose))
            {
                n <- seqSummary(f, "genotype", verbose=FALSE)$dim[3L]
                .cat(" [", prettyNum(n, big.mark=",", scientific=FALSE),
                    " variants]")
            }
        }
    }
    # the chromosome of each file (the files opened above are closed by the
    #   on.exit() handler if an error occurs here)
    nm <- vapply(seq_along(ans), function(i)
    {
        v <- seqGetData(ans[[i]], "$chromosome")
        if (nrun(v) > 1L)
        {
            fn <- if (is.character(gds_fn)) gds_fn[i] else
                paste0("annot_gds[[", i, "]]")
            stop(fn, " should only contain one chromosome.")
        }
        paste(as.character(unique(v)), collapse="&")
    }, "")
    names(ans) <- paste0("chr", nm)
    # output
    on.exit()
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
        ok <- vapply(annot_gds, inherits, FALSE, what="SeqVarGDSClass")
        if (!all(ok))
            stop("annot_gds[[", which(!ok)[1L], "]] is not a GDS object.")
    }
    invisible()
}

# Check and open the annotation GDS files, and set the annotation names;
#   'close' in the returned list is TRUE if the files are opened here (to be
#   closed by the caller)
.annot_open <- function(annot_gds, varnm, all_varnm, verbose)
{
    .check_annot_gds(annot_gds)
    close <- is.character(annot_gds)
    gds <- .open_annot_gds(annot_gds, verbose)
    if (close) on.exit(.close_annot_gds(gds))
    if (all_varnm)
        varnm <- .annot_list(gds[[1L]])
    if (is.null(varnm)) varnm <- character()
    stopifnot(is.character(varnm))
    on.exit()
    list(gds=gds, varnm=varnm, close=close)
}

# GDS nodes for the INFO field
.gds_varnm_cvt <- c(
    ':chromosome'="chromosome", ':position'="position", ':allele'="allele",
    ':id'="annotation/id", ':qual'="annotation/qual",
    ':filter'="annotation/filter")

# the GDS node names of the annotation names 'nm': a name starting with ":"
#   is a variant node (see .gds_varnm_cvt), any other name is under
#   "annotation/info"
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


# Helper function for .ann_chr_pos_allele(): annotate with one file
.ann_pos_allele <- function(gds, chr, pos, ref, alt, varnm, verbose=TRUE)
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
            message("[", basename(gds$filename), "] ", appendLF=FALSE)
            n <- seqSummary(gds, "genotype", verbose=FALSE)$seldim[3L]
            .cat(" # of variants found: ", n)
        }
        l_verbose <- isTRUE(verbose) && (length(varnm)>1L)
        if (l_verbose)
            message("[", length(varnm), "] ", appendLF=FALSE)
        ans <- lapply(varnm, function(nm)
        {
            if (l_verbose) message(".", appendLF=FALSE)
            seqGetData(gds, nm, .tolist=NA)
        })
        names(ans) <- colnm
        ans <- DataFrame(ans)
        if (l_verbose) message()
        # ii maps each input to the filtered set row (NA = not found)
        # drop=FALSE: keep a DataFrame even when a single annotation is
        #   requested
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


# Helper function for .ann_chr_pos_allele(): annotate one chromosome
.ann_one_chr <- function(ch, idx, pos, ref, alt, annot_gds, gds_idx, gds_chr,
    varnm, verbose)
{
    # find the GDS file(s) for this chromosome
    gds_ii <- which(names(gds_idx) == .norm_chr(ch))
    if (!length(gds_ii))
    {
        # no file for this chromosome: all the variants are not found
        d <- .ann_pos_allele(annot_gds[[1L]], ch, pos[idx], ref[idx],
            alt[idx], varnm, verbose=FALSE)
        if (!length(varnm)) d$file_idx <- Rle(1L, nrow(d))
        d$..idx <- idx
        d$..no <- NULL
        return(d)
    }
    # process each GDS file for this chromosome: the variants not found
    #   in a file are looked up in the next file(s), so the files are
    #   processed one after another
    v <- vector("list", length(gds_ii))
    for (j in seq_along(gds_ii))
    {
        if (!length(idx)) break
        k <- gds_idx[gds_ii[j]]
        # annotate variants in this chromosome subset, using the
        #   chromosome label of the file
        d <- .ann_pos_allele(annot_gds[[k]], gds_chr[gds_ii[j]],
            pos[idx], ref[idx], alt[idx], varnm, verbose)
        if (!length(varnm)) d$file_idx <- Rle(k, nrow(d))
        # temporary index for combining results in original input order
        d$..idx <- idx
        if (j < length(gds_ii) && any(d$..no))
        {
            # some variants are not found in this file,
            # so they will be processed in the next file(s)
            idx <- idx[d$..no]
            d <- d[!d$..no, , drop=FALSE]
        }
        d$..no <- NULL  # remove the temporary column
        v[[j]] <- d
    }
    # combine results for this chromosome
    v <- v[!vapply(v, is.null, FALSE)]
    if (length(v) == 1L) v[[1L]] else do.call(rbind, v)
}

# Annotate with chromosome, position, reference & alternative alleles
#   dispatch per chromosome, combine results in original input order
.ann_chr_pos_allele <- function(chr, pos, ref, alt, annot_gds, varnm,
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
    # match chromosome names in the annotation GDS files, ignoring the "chr"
    #   prefix; 'gds_chr' is the chromosome label as stored in each file
    s <- strsplit(sub("^chr", "", names(annot_gds)), "&", fixed=TRUE)
    gds_idx <- rep(seq_along(annot_gds), times=lengths(s))
    gds_chr <- unlist(s)
    names(gds_idx) <- .norm_chr(gds_chr)
    # process each chromosome subset
    ans <- lapply(chr_lst, function(ch)
        .ann_one_chr(ch, which(chr == ch), pos, ref, alt, annot_gds, gds_idx,
            gds_chr, varnm, verbose))
    # combine results in original input order
    ans <- if (length(ans) == 1L) ans[[1L]] else do.call(rbind, ans)
    # return
    ans <- ans[order(ans$..idx), , drop=FALSE]
    ans$..idx <- NULL  # remove the temporary index column
    ans
}


# Annotate with a data frame including the columns for
#   chromosome, position, reference & alternative alleles
.ann_dataframe <- function(object, annot_gds, varnm, col_chr="chr",
    col_pos="pos", col_ref="ref", col_alt="alt", ..., verbose=TRUE)
{
    # check
    chr <- object[[col_chr]]
    if (is.null(chr)) stop("No 'chr' column.")
    pos <- object[[col_pos]]
    if (is.null(pos)) stop("No 'pos' column.")
    ref <- object[[col_ref]]
    if (is.null(ref)) stop("No 'ref' column.")
    alt <- object[[col_alt]]
    # check & open annotated gds
    a <- .annot_open(annot_gds, if (!missing(varnm)) varnm, missing(varnm),
        verbose)
    if (a$close) on.exit(.close_annot_gds(a$gds))
    # process
    .ann_chr_pos_allele(chr, pos, ref, alt, a$gds, a$varnm, verbose=verbose)
}


#######################################################################
# Annotate a GDS file
#

# read the chromosomes, positions, reference & alternative alleles
.read_alleles <- function(object, verbose)
{
    if (isTRUE(verbose)) message("Reading chromosome", appendLF=FALSE)
    chr <- seqGetData(object, "$chromosome")
    if (isTRUE(verbose)) message(", position", appendLF=FALSE)
    pos <- seqGetData(object, "position")
    if (isTRUE(verbose)) message(", reference allele", appendLF=FALSE)
    ref <- seqGetData(object, "$ref")
    if (isTRUE(verbose)) message(", alternative allele", appendLF=FALSE)
    alt <- seqGetData(object, "$alt")
    if (isTRUE(verbose)) message()
    list(chr=chr, pos=pos, ref=ref, alt=alt)
}

# export a GDS file to a new file, opened in the writing mode
.ann_export <- function(object, out_fn, no_sample, verbose)
{
    fmt.var <- NULL
    if (isTRUE(no_sample))
    {
        seqSetFilter(object, sample.sel=integer(), action="push+set",
            verbose=FALSE)
        on.exit(seqFilterPop(object))
        fmt.var <- character()
    }
    seqExport(object, out_fn, fmt.var=fmt.var, optimize=FALSE,
        verbose=verbose)
    seqOpen(out_fn, readonly=FALSE)
}

# create the parent folders for the nested variable names (e.g.,
#   "CSQ.list/Allele" as returned by seqAnnotList()), since seqAddValue()
#   requires the parent node to exist
.ann_add_folders <- function(outgds, rootnm, colnm)
{
    for (s in setdiff(unique(dirname(colnm)), "."))
    {
        ss <- strsplit(s, "/", fixed=TRUE)[[1L]]
        for (j in seq_along(ss))
        {
            nd <- paste0(rootnm, "/", paste(ss[seq_len(j)], collapse="/"))
            if (!exist.gdsn(outgds, nd))
                seqAddValue(outgds, nd, NULL, verbose=FALSE)
        }
    }
    invisible()
}

# split 'n' variants into blocks of 'bsize' variants: the index of the first
#   variant and the number of variants of each block
.blocks <- function(n, bsize)
{
    nblock <- max(1L, as.integer(ceiling(n / bsize)))
    st <- as.double(seq_len(nblock) - 1L) * bsize + 1
    st <- as.integer(pmin(st, n + 1))
    list(st=st, cnt=diff(c(st, n + 1L)))
}

# the values of a GDS variable for a block of variants, in the compact
#   list(length, data) representation
.ann_get_block <- function(annot_gds, nm, fi, vi)
{
    rl <- runLength(fi); rv <- runValue(fi)
    # starting offset of each run (1, rl1+1, rl1+rl2+1, ...)
    ost <- cumsum(c(1L, rl))
    v <- lapply(seq_along(rv), function(j)
    {
        f <- annot_gds[[rv[j]]]
        ii <- vi[seq.int(ost[j], length.out=rl[j])]
        ii <- seqSetFilter(f, variant.sel=ii, ret.idx=TRUE,
            warn=FALSE, verbose=FALSE)$variant_idx
        seqGetData(f, nm, .tolist=TRUE)[ii]
    })
    v <- if (length(rv) == 1L) v[[1L]] else
        unlist(v, recursive=FALSE, use.names=FALSE)
    v <- lapply(v, function(x) unlist(x, use.names=FALSE))
    structure(list(length=lengths(v),
        data=unlist(v, use.names=FALSE)), class="SeqVarDataList")
}

# evaluate 'expr', muffling only the warning of gdsfmt on the missing
#   characters (stored as "" since NA is not allowed in a GDS string variable)
.muffle_na_string <- function(expr)
{
    withCallingHandlers(expr, warning=function(w)
    {
        if (startsWith(conditionMessage(w),
            "Missing characters are converted to"))
        {
            invokeRestart("muffleWarning")
        }
    })
}

# write the annotations of the variants to a GDS file block by block
.ann_write <- function(outgds, root, varnm, map, annot_gds, bsize, verbose)
{
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
    .ann_add_folders(outgds, rootnm, colnm)
    # group found variants by annotation file, map$file_idx has no NA
    file_ids <- map$file_idx
    if (!inherits(file_ids, "Rle"))
        file_ids <- Rle(as.integer(file_ids))
    # split the variants into blocks to reduce the memory usage: a list
    #   of n_total elements (as returned by .tolist=TRUE) needs much more
    #   memory than the annotation data itself
    blk <- .blocks(nrow(map), bsize)
    nblock <- length(blk$st)
    dot_block <- if (nblock > 1L)
        unique(as.integer(round(seq(1L, nblock,
            length.out=min(10L, nblock))))) else integer()
    # the k-th block of the variable 'varnm[i]', NULL if no more block; the
    #   data are written block by block, so that only one block is held in
    #   memory instead of the annotation of all variants
    gen <- function(k, i)
    {
        if (k > nblock) return(NULL)  # no more block
        if (verbose && (k %in% dot_block))
            message(".", appendLF=FALSE)
        s <- seq.int(blk$st[k], length.out=blk$cnt[k])
        .ann_get_block(annot_gds, varnm[i], file_ids[s], map$variant_idx[s])
    }
    # add each annotation variable to the file
    for (i in seq_along(varnm))
    {
        if (verbose)
        {
            message("    adding ", i, "/", length(varnm), " ",
                sQuote(colnm[i]), " (", .tm(), ")\n    ", appendLF=FALSE)
        }
        .muffle_na_string(seqAddValue(outgds, paste0(rootnm, "/", colnm[i]),
            gen, replace=TRUE, param=i, verbose=verbose, verbose.attr=FALSE))
    }
    invisible()
}

# check 'add_to_gds'
.check_add_to_gds <- function(object, add_to_gds)
{
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
    invisible()
}

# Annotate a GDS file
.ann_gdsfile <- function(object, annot_gds, varnm, add_to_gds=FALSE,
    no_sample=TRUE, root="", cleanup=FALSE, bsize=1000000L, ..., verbose=TRUE)
{
    # check
    stopifnot(inherits(object, "SeqVarGDSClass"))
    stopifnot(is.numeric(bsize), length(bsize)==1L, !is.na(bsize), bsize>0L)
    bsize <- as.double(bsize)
    .check_add_to_gds(object, add_to_gds)
    stopifnot(is.logical(no_sample), length(no_sample)==1L)
    stopifnot(is.character(root), length(root)==1L, !is.na(root))
    # check & open annotated gds
    a <- .annot_open(annot_gds, if (!missing(varnm)) varnm, missing(varnm),
        verbose)
    if (a$close) on.exit(.close_annot_gds(a$gds))
    # read data
    v <- .read_alleles(object, verbose)
    # process
    if (isTRUE(verbose)) message("Processing annotation ...")
    if (isFALSE(add_to_gds) || length(a$varnm)==0L)
    {
        # return DataFrame
        return(.ann_chr_pos_allele(v$chr, v$pos, v$ref, v$alt, a$gds,
            a$varnm, verbose=verbose))
    }
    # when add_to_gds=TRUE or it is a file name
    map <- .ann_chr_pos_allele(v$chr, v$pos, v$ref, v$alt, a$gds,
        varnm=character(), verbose=verbose)
    remove(v)
    if (verbose)
    {
        n <- sum(is.na(map$variant_idx))
        .cat("Found #: ", nrow(map)-n, " variant(s), Not found #: ", n)
    }
    # output gds file
    outgds <- object
    if (is.character(add_to_gds))
    {
        outgds <- .ann_export(object, add_to_gds, no_sample, verbose)
        on.exit({
            seqClose(outgds)
            if (isTRUE(cleanup))
                cleanup.gds(add_to_gds, verbose=verbose)
        }, add=TRUE)
    }
    .ann_write(outgds, root, a$varnm, map, a$gds, bsize, verbose)
}


#######################################################################
# Annotate the other classes
#

# Annotate "chr-pos-ref-alt"
.ann_variant <- function(object, annot_gds, varnm, split="[-_:]", ...,
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
    a <- .annot_open(annot_gds, if (!missing(varnm)) varnm, missing(varnm),
        verbose)
    if (a$close) on.exit(.close_annot_gds(a$gds))
    # process
    .ann_chr_pos_allele(chr, pos, ref, alt, a$gds, a$varnm, verbose=verbose)
}


# Annotate all the variants located in genomic ranges
.ann_ranges <- function(object, annot_gds, varnm, all_varnm, chr, verbose)
{
    # check & open annotated gds
    a <- .annot_open(annot_gds, varnm, all_varnm, verbose)
    if (a$close) on.exit(.close_annot_gds(a$gds))
    varnm <- a$varnm
    # verbose
    l_verbose <- isTRUE(verbose) && (length(varnm)>1L)
    if (l_verbose)
        message("[", length(varnm), "] ", appendLF=FALSE)
    # GDS node names
    colnm <- names(varnm)
    if (is.null(colnm)) colnm <- varnm
    varnm <- .gds_varnm(varnm)
    # process
    ans <- lapply(a$gds, function(gds)
    {
        if (is.null(chr))
            seqSetFilter(gds, object, verbose=FALSE)
        else
            seqSetFilter(gds, object, chr=chr, verbose=FALSE)
        if (seqSummary(gds, "genotype", verbose=FALSE)$seldim[3L])
        {
            DataFrame(lapply(varnm, function(nm)
            {
                if (l_verbose) message(".", appendLF=FALSE)
                seqGetData(gds, nm, .tolist=NA)
            }))
        } else
            NULL
    })
    if (l_verbose) message()
    # output
    ans <- do.call(rbind, ans)
    if (!is.null(ans)) names(ans) <- colnm
    ans
}

# Annotate a GRanges/GRangesList object
.ann_GRanges <- function(object, annot_gds, varnm, ..., verbose=TRUE)
{
    .ann_ranges(object, annot_gds, if (!missing(varnm)) varnm,
        missing(varnm), NULL, verbose)
}

# Annotate an IRanges object
.ann_IRanges <- function(object, annot_gds, varnm, chr, ..., verbose=TRUE)
{
    # check
    if (missing(chr))
        stop("'chr' should be specified for an IRanges object.")
    stopifnot(is.character(chr), length(chr)>0L)
    .ann_ranges(object, annot_gds, if (!missing(varnm)) varnm,
        missing(varnm), chr, verbose)
}


# Annotate a VCF object (from the package VariantAnnotation): one row per
#   alternative allele, so a CollapsedVCF and its expand()-ed version give
#   the same output; a row without any alternative allele is skipped.
#   VariantAnnotation is a suggested package, loaded on demand, so the class
#   VCF is not used in a method signature: seqAnnotate() reaches this
#   function via the method for "ANY" below. seqnames() and start() are the
#   generics imported from GenomicRanges and IRanges, whose methods for
#   VCF objects are registered by SummarizedExperiment (a dependency of
#   VariantAnnotation), and ref() and alt() are called with '::' since
#   SeqArray defines its own generics of these names
.ann_VCF <- function(object, annot_gds, varnm, ..., verbose=TRUE)
{
    if (!requireNamespace("VariantAnnotation", quietly=TRUE))
        stop("The package 'VariantAnnotation' should be installed.")
    ref_al <- as.character(VariantAnnotation::ref(object))
    alt_al <- VariantAnnotation::alt(object)
    if (is.character(alt_al) || is(alt_al, "XStringSet"))
    {
        # ExpandedVCF: one alternative allele per row (note that a
        #   DNAStringSet is a List of DNAString, so lengths() would give
        #   the number of bases)
        nalt <- rep.int(1L, length(alt_al))
        alt_al <- as.character(alt_al)
    } else {
        # CollapsedVCF: a DNAStringSetList or CharacterList with zero, one
        #   or more alternative alleles per row
        nalt <- lengths(alt_al)
        alt_al <- as.character(unlist(alt_al, use.names=FALSE))
    }
    variants <- DataFrame(
        chr=rep(as.character(seqnames(object)), nalt),
        pos=rep(start(object), nalt),
        ref=rep(ref_al, nalt), alt=unname(alt_al))
    # skip the rows without alternative allele (ALT is "." in the VCF file)
    i <- !is.na(variants$alt) & nzchar(variants$alt) & variants$alt != "."
    if (!all(i)) variants <- variants[i, , drop=FALSE]
    # process
    seqAnnotate(variants, annot_gds, varnm, ..., verbose=verbose)
}


# Fallback for the classes without a specific method: a VCF object of the
#   package VariantAnnotation (CollapsedVCF or ExpandedVCF) is annotated by
#   .ann_VCF(), any other class is an error
.ann_any <- function(object, annot_gds, varnm, ..., verbose=TRUE)
{
    if (inherits(object, "VCF"))
    {
        .ann_VCF(object, annot_gds, varnm, ..., verbose=verbose)
    } else {
        stop("No seqAnnotate() method for an object of class ",
            sQuote(class(object)[1L]), ".")
    }
}


# Set methods
setMethod("seqAnnotate", signature(object="data.frame"), .ann_dataframe)
setMethod("seqAnnotate", signature(object="DataFrame"), .ann_dataframe)
setMethod("seqAnnotate", signature(object="SeqVarGDSClass"), .ann_gdsfile)
setMethod("seqAnnotate", signature(object="character"), .ann_variant)
setMethod("seqAnnotate", signature(object="GRanges"), .ann_GRanges)
setMethod("seqAnnotate", signature(object="GRangesList"), .ann_GRanges)
setMethod("seqAnnotate", signature(object="IRanges"), .ann_IRanges)
setMethod("seqAnnotate", signature(object="ANY"), .ann_any)


# Annotate a VCF file: read it with VariantAnnotation::readVcf() and
#   annotate the VCF object
#   vcf_fn: the VCF or BCF file name
#   annot_gds, varnm, verbose: see seqAnnotate()
seqAnnotateVCF <- function(vcf_fn, annot_gds, varnm, ..., verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)==1L)
    if (!requireNamespace("VariantAnnotation", quietly=TRUE))
        stop("The package 'VariantAnnotation' should be installed.")
    # open the file
    if (isTRUE(verbose)) .cat("Open ", sQuote(vcf_fn))
    vcf <- VariantAnnotation::readVcf(vcf_fn)
    # output
    .ann_VCF(vcf, annot_gds, varnm, ..., verbose=verbose)
}


# Annotate a GDS file with a file name input
#   gds_fn: the SeqArray GDS file name
#   annot_gds, varnm, verbose: see seqAnnotate()
#   add_to_gds: FALSE to return a DataFrame, TRUE to add the annotations to
#       'gds_fn', or the name of a new GDS file
#   no_sample: if TRUE, the new GDS file has no per-sample data
#   root: the folder under "annotation/info" for the annotations, "" for none
#   cleanup: if TRUE, clean up the fragments of the output GDS file
#   bsize: the number of variants in a block when writing the annotations
seqAnnotateGDS <- function(gds_fn, annot_gds, varnm, add_to_gds=FALSE,
    no_sample=TRUE, root="", cleanup=TRUE, bsize=1000000L, ..., verbose=TRUE)
{
    # check
    stopifnot(is.character(gds_fn), length(gds_fn)==1L)
    if (isTRUE(verbose))
        message("Open ", sQuote(gds_fn), appendLF=FALSE)
    gds <- seqOpen(gds_fn, readonly=!isTRUE(add_to_gds))
    on.exit(seqClose(gds))
    if (isTRUE(verbose))
    {
        dm <- seqSummary(gds, "genotype", verbose=FALSE)$dim
        message(" [", prettyNum(dm[3L], big.mark=",", scientific=FALSE),
            " variants]")
    }
    if (isTRUE(add_to_gds) && isTRUE(cleanup))
    {
        on.exit(cleanup.gds(gds_fn, verbose=verbose), add=TRUE)
    }

    # process
    .ann_gdsfile(gds, annot_gds, varnm, add_to_gds, no_sample, root,
        cleanup=cleanup, bsize=bsize, ..., verbose=verbose)
}
