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


# Check whether in increasing order or not, return false if any missing value
.test_increasing <- function(x)
{
    .Call(gds_test_increasing, x)
}



# Open the GDS files with variant annotation
.load_annot_gds <- function(gds_fn, verbose=TRUE)
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
        # check position
        v <- seqGetData(f, "position")
        if (anyNA(v))
            stop("Positions in ", gds_fn[i], " should not have any misssing value.")
        if (is.unsorted(v) || v[1L] > v[length(v)])
            stop("Positions in ", gds_fn[i], " should be in an increasing order.")
    }
    # re-order
    names(ans) <- nm
    ans <- ans[order(suppressWarnings(as.numeric(nm)), nm)]
    # output
    on.exit()
    ans
}


# Annotate with the same chromosome
ann_pos_allele <- function(annot_gds, pos, allele, varnm)
{
    if (length(annot_gds))
    {
        # i <- .Call
        ans <- seqGetData(annot_gds[[1L]], paste0("annotation/info/", varnm),
            .tolist=TRUE)
        if (anyNA(i) || duplicated(i) || is.unsorted(i) || i[1L]>i[length(i)])
            ans <- lapply(ans, function(v) v[i])
    } else {
        # no annotation
        v <- rep(NA, length(pos))
        ans <- vector("list", length(nm))
        for (i in seq_along(ans)) ans[[i]] <- v
    }
    # output
    names(ans) <- varnm
    ans
}

# Annotate with chromosome, position, reference & alternative alleles
ann_chr_pos_allele <- function(chr, pos, ref, alt, annot_gds, varnm, verbose=TRUE)
{
    # check
    stopifnot(length(chr) == length(pos))
    stopifnot(length(chr) == length(ref))
    if (!is.null(alt))
    {
        stopifnot(length(chr) == length(alt))
        ref <- paste0(ref, ",", alt)
    }
    stopifnot(is.character(varnm))
    stopifnot(is.logical(verbose), length(verbose)==1L)
    # check annotation gds file(s)
    stopifnot(is.character(annot_gds), length(annot_gds)>0L)
    annot_gds <- .load_annot_gds(annot_gds)
    if (length(varnm) == 0L) return(NULL)

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


ann_dataframe <- function(object, annot_gds, varnm, ..., verbose=TRUE)
{

}

setMethod("seqAnnotate", signature(object="data.frame"), ann_dataframe)

