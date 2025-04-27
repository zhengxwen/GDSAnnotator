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


# Open the GDS files with variant annotation
.open_annot_gds <- function(gds_fn, verbose=TRUE)
{
    # check
    stopifnot(is.character(gds_fn), length(gds_fn)>0L)
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
            stop("Positions in ", gds_fn[i], " should not have any misssing value.")
        if (is.unsorted(v) || v[1L] > v[length(v)])
            stop("Positions in ", gds_fn[i], " should be in increasing order.")
    }
    # re-order
    names(ans) <- nm
    ans <- ans[order(suppressWarnings(as.numeric(nm)), nm)]
    # output
    on.exit()
    ans
}

.close_annot_gds <- function(annot_gds)
{
    lapply(annot_gds, seqClose)
}


# Annotate with the same chromosome
ann_pos_allele <- function(annot_gds, pos, allele, varnm)
{
    # check
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
    # process
    if (length(annot_gds))
    {
        # use the first gds
        f <- annot_gds[[1L]]
        ptr <- SeqArray:::.buffer_position(f)  # pointer to positions
        allele_nd <- index.gdsn(f, "allele")
        # call C
        ii <- .Call(SEQ_Find_Position, ptr, allele_nd, pos, allele)
        # clear position buffer
        SeqArray:::.buffer_position(f, clear=TRUE)
        # set a filter to get data
        ii <- seqSetFilter(f, variant.sel=ii, ret.idx=TRUE)$variant_idx
        ans <- seqGetData(f, paste0("annotation/info/", varnm),
            .tolist=TRUE)
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
    } else {
        # no annotation
        v <- rep(NA, length(pos))
        ans <- lapply(seq_along(varnm), function(i) v)
    }
    # output
    names(ans) <- varnm
    ans
}

# Annotate with chromosome, position, reference & alternative alleles
ann_chr_pos_allele <- function(chr, pos, ref, alt, annot_gds, varnm,
    verbose=TRUE)
{
    # check
    stopifnot(length(chr)==1L || length(chr)==length(pos))
    stopifnot(length(pos) == length(ref))
    if (!is.null(alt))
    {
        stopifnot(length(pos) == length(alt))
        ref <- paste0(ref, ",", alt)
    }
    stopifnot(is.character(varnm))
    stopifnot(is.logical(verbose), length(verbose)==1L)
    # check annotation gds file(s)
    stopifnot(is.character(annot_gds), length(annot_gds)>0L)
    # open gds file(s)
    annot_gds <- .open_annot_gds(annot_gds)
    on.exit(.close_annot_gds(annot_gds))
    if (length(varnm) == 0L) return(NULL)
    # prepare
    if (!is.integer(pos)) pos <- as.integer(pos)
    # process each chromosome
    chr_lst <- unique(chr)
    ans <- lapply(chr_lst, function(ch)
    {
        a <- annot_gds[names(annot_gds) == ch]
        ann_pos_allele(a, pos, ref, varnm)
    })
    # output
    ans[[1L]]
}


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
    if (is.null(alt)) stop("No 'alt' column.")
    # process
    ann_chr_pos_allele(chr, pos, ref, alt, annot_gds, varnm, verbose=verbose)
}


# Set methods
setMethod("seqAnnotate", signature(object="data.frame"), ann_dataframe)

