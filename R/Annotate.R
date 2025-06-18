#######################################################################
#
# Package Name: GDSAnnotator
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Data import & Format conversion
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
        # open the GDS file(s)
        ans <- vector("list", length(gds_fn))
        on.exit({ for (f in ans) if (!is.null(f)) seqClose(f) })
        nm <- character(length(gds_fn))
        for (i in seq_along(gds_fn))
        {
            if (isTRUE(verbose))
                .cat("Open ", sQuote(gds_fn[i]))
            # open the file
            f <- ans[[i]] <- seqOpen(gds_fn[i])
            # check chromosome
            v <- seqGetData(f, "$chromosome")
            if (nrun(v) > 1L)
                stop(gds_fn[i], " should only contain one chromosome.")
            nm[i] <- as.character(v[1L])
            # check position (should in in increasing order)
            v <- seqGetData(f, "position")
            if (anyNA(v))
            {
                stop("Positions in ", gds_fn[i],
                    " should not have any misssing value.")
            }
            if (is.unsorted(v) || v[1L] > v[length(v)])
            {
                stop("Positions in ", gds_fn[i],
                    " should be in increasing order.")
            }
        }
        # re-order
        names(ans) <- nm
        ans <- ans[order(suppressWarnings(as.numeric(nm)), nm)]
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
        if (anyNA(nm)) stop("'varnm' should not have NA_character_.")
        i <- match(nm, names(.gds_varnm_cvt))
        if (anyNA(i))
            nm[is.na(i)] <- paste0("annotation/info/", nm[is.na(i)])
        j <- which(!is.na(i))
        if (length(j))
            nm[j] <- .gds_varnm_cvt[i[j]]
    }
    nm
}


# Annotate with the same chromosome
ann_pos_allele <- function(annot_gds, annot_gds_idx, pos, allele, varnm,
    verbose=TRUE)
{
    # check
    stopifnot(length(annot_gds)==length(annot_gds_idx))
    stopifnot(is.integer(pos))
    stopifnot(is.character(allele))
    stopifnot(length(pos)==length(allele))
    # position
    idx <- NULL
    if (is.unsorted(pos) || pos[1L] > pos[length(pos)])
    {
        idx <- order(pos)
        pos <- pos[idx]; allele <- allele[idx]
        idx <- order(idx)
    }
    # GDS node names
    colnm <- varnm
    varnm <- .gds_varnm(varnm)
    # process
    if (length(annot_gds))
    {
        # use the first gds
        f <- annot_gds[[1L]]
        ptr <- SeqArray:::.buffer_position(f)  # pointer to positions
        allele_nd <- index.gdsn(f, "allele")
        # call C to find variant indices
        ii <- .Call(SEQ_Find_Position, ptr, allele_nd, pos, allele)
        # clear position buffer
        SeqArray:::.buffer_position(f, clear=TRUE)
        # get
        if (length(varnm))
        {
            # set a filter to get data
            if (verbose)
                cat("[", basename(f$filename), "] ", sep="")
            ii <- seqSetFilter(f, variant.sel=ii, ret.idx=TRUE,
                verbose=verbose)$variant_idx
            ans <- seqGetData(f, varnm, .tolist=NA)
            ans <- DataFrame(ans)
            # check
            rerow <- anyNA(ii) || is.unsorted(ii) || ii[1L] > ii[length(ii)]
            if (rerow)
            {
                if (!is.null(idx))
                    ans <- ans[ii[idx], ]
                else
                    ans <- ans[ii, ]
            } else {
                if (!is.null(idx)) ans <- ans[idx, ]
            }
            names(ans) <- colnm
        } else {
            ans <- DataFrame(file_idx=Rle(annot_gds_idx[1L], length(ii)),
                variant_idx=ii)
            if (!is.null(idx)) ans <- ans[idx, ]
        }
    } else {
        # no annotation
        v <- rep(NA, length(pos))
        ans <- lapply(seq_along(varnm), function(i) v)
        names(ans) <- colnm
    }
    # output
    ans
}

# Annotate with chromosome, position, reference & alternative alleles
#   where 'alt' can be NULL to use 'ref' only
ann_chr_pos_allele <- function(chr, pos, ref, alt, annot_gds, varnm,
    verbose=TRUE)
{
    # check
    stopifnot(length(chr)==1L || length(chr)==length(pos))
    stopifnot(length(pos) == length(ref))
    if (!missing(alt) && !is.null(alt))
    {
        stopifnot(length(pos) == length(alt))
        s <- ref
        ref <- paste0(ref, ",", alt)  # TODO: optimize
        if (anyNA(alt))
        {
            x <- is.na(alt)
            ref[x] <- s[x]
            remove(x)
        }
        remove(s)
    }
    stopifnot(is.logical(verbose), length(verbose)==1L)
    # check & process
    stopifnot(is.character(varnm))
    if (identical(varnm, "")) varnm <- character()
    # process pos & each chromosome
    if (!is.integer(pos)) pos <- as.integer(pos)
    chr_lst <- unique(chr)
    ans <- lapply(chr_lst, function(ch)
    {
        ii <- which(names(annot_gds) == ch)
        a <- annot_gds[ii]
        ann_pos_allele(a, ii, pos, ref, varnm, verbose)
    })
    # output
    ans[[1L]]
}


# Annotate with a data frame including the columns for
#   chromosome, position, reference & alternative alleles
ann_dataframe <- function(object, annot_gds, varnm, col_chr="chr",
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
    .check_annot_gds(annot_gds)
    if_close_gds <- is.character(annot_gds)
    annot_gds <- .open_annot_gds(annot_gds)
    if (if_close_gds)
        on.exit(.close_annot_gds(annot_gds))
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
    stopifnot(is.character(varnm))
    if (identical(varnm, "")) varnm <- character()
    # check & open annotated gds
    .check_annot_gds(annot_gds)
    if_close_gds <- is.character(annot_gds)
    annot_gds <- .open_annot_gds(annot_gds)
    if (if_close_gds)
        on.exit(.close_annot_gds(annot_gds))
    # process
    chr <- seqGetData(object, "$chromosome")
    if (nrun(chr) > 1L)
        stop(basename(object$filename), " should only contain one chromosome.")
    pos <- seqGetData(object, "position")
    allele <- seqGetData(object, "allele")
    if (isFALSE(add_to_gds) || length(varnm)==0L)
    {
        # return DataFrame
        ann_chr_pos_allele(chr, pos, allele, NULL, annot_gds, varnm,
            verbose=verbose)
    } else {
        # when add_to_gds=TRUE or it is a file name
        # sort position
        idx <- NULL
        if (is.unsorted(pos) || pos[1L] > pos[length(pos)])
        {
            idx <- order(pos)
            pos <- pos[idx]; allele <- allele[idx]
            idx <- order(idx)
        }
        map <- ann_chr_pos_allele(chr, pos, allele, NULL, annot_gds, varnm="",
            verbose=verbose)
        remove(chr, pos, allele)
        # check
        if (anyNA(map$file_idx) || length(unique(map$file_idx))>1L)
        {
            stop("The GDS file should have only one chromosome.")
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
        ii <- seqSetFilter(f, variant.sel=map$variant_idx, ret.idx=TRUE,
            verbose=verbose)$variant_idx
        rerow <- anyNA(ii) || is.unsorted(ii) || ii[1L] > ii[length(ii)]
        if (rerow && !is.null(idx)) ii <- ii[idx]
        if (!rerow && !is.null(idx)) { rerow <- TRUE; ii <- idx }
        # add each node
        for (i in seq_along(varnm))
        {
            if (verbose)
                cat("    adding", sQuote(colnm[i]), "...\n    ")
            # set a filter to get data
            v <- seqGetData(f, varnm[i], .tolist=TRUE)
            if (rerow) v <- v[ii]
            seqAddValue(outgds, paste0(rootnm, "/", colnm[i]), v,
                replace=TRUE, verbose=verbose, verbose.attr=FALSE)
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


# Annotate a GDS file with a file name input
ann_filename <- function(object, annot_gds, varnm, ..., verbose=TRUE)
{
    # check
    stopifnot(is.character(object), length(object)==1L)
    if (isTRUE(verbose)) .cat("Open ", sQuote(object))
    object <- seqOpen(object)
    on.exit(seqClose(object))
    # process
    ann_gdsfile(object, annot_gds, varnm, ..., verbose=TRUE)
}


# Set methods
setMethod("seqAnnotate", signature(object="data.frame"), ann_dataframe)
setMethod("seqAnnotate", signature(object="DataFrame"), ann_dataframe)
setMethod("seqAnnotate", signature(object="SeqVarGDSClass"), ann_gdsfile)
setMethod("seqAnnotate", signature(object="character"), ann_filename)

