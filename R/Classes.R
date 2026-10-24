#######################################################################
#
# Package Name: GDSAnnotator
# Copyright (C) 2025-2026    Xiuwen Zheng
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Classes: the opened annotation GDS files
#


# A list of opened annotation GDS files, with the file name, the chromosomes
#   (the labels used in the file) and the number of variants of each file in
#   mcols(), which are kept in sync with the files by subsetting
setClass("SeqAnnotGDSList", contains="SimpleList",
    prototype=prototype(elementType="SeqVarGDSClass"))

setValidity("SeqAnnotGDSList", function(object)
{
    nm <- c("filename", "chromosome", "n_variant")
    if (!all(nm %in% names(mcols(object))))
    {
        return(paste0("mcols() should have the columns ",
            paste(sQuote(nm, FALSE), collapse=", "), "."))
    }
    TRUE
})

# the number of variants per chromosome of an opened GDS file, using all the
#   variants instead of the selected ones (the filter is restored afterwards)
.gds_chr_count <- function(f)
{
    seqFilterPush(f)
    on.exit(seqFilterPop(f))
    seqResetFilter(f, verbose=FALSE)
    seqSummary(f, "chromosome", verbose=FALSE)
}

# a SeqAnnotGDSList object from a list of opened GDS files
.new_annot_gds <- function(lst)
{
    cnt <- lapply(lst, .gds_chr_count)
    m <- DataFrame(
        filename = vapply(lst, function(f)
            normalizePath(f$filename, mustWork=FALSE), "", USE.NAMES=FALSE),
        chromosome = CharacterList(lapply(cnt, names)),
        n_variant = vapply(cnt, sum, 0, USE.NAMES=FALSE))
    new("SeqAnnotGDSList", listData=lst, elementMetadata=m)
}

# Open the annotation GDS file(s)
seqAnnotOpen <- function(gds_fn, allow.duplicate=FALSE)
{
    # check
    stopifnot(is.character(gds_fn))
    if (length(gds_fn) == 0L)
        stop("No annotation gds file!")
    if (anyNA(gds_fn))
        stop("'gds_fn' should not contain NA.")
    fn <- gds_fn[!file.exists(gds_fn)]
    if (length(fn))
        stop("No file(s): ", paste(fn, collapse=", "), ".")
    stopifnot(is.logical(allow.duplicate), length(allow.duplicate)==1L)
    # open the GDS file(s), which are closed if an error occurs
    lst <- vector("list", length(gds_fn))
    on.exit({ for (f in lst) if (!is.null(f)) seqClose(f) })
    for (i in seq_along(gds_fn))
        lst[[i]] <- seqOpen(gds_fn[i], allow.duplicate=allow.duplicate)
    ans <- .new_annot_gds(lst)
    # output
    on.exit()
    ans
}

# Close all the files
setMethod("seqClose", "SeqAnnotGDSList", function(object)
{
    for (f in as.list(object)) seqClose(f)
    invisible()
})

# Show the files
setMethod("show", "SeqAnnotGDSList", function(object)
{
    m <- mcols(object)
    cat("SeqAnnotGDSList object with ", length(object),
        " annotation GDS file(s)\n", sep="")
    idx <- format(paste0("[", seq_along(object), "]"), justify="right")
    for (i in seq_along(object))
    {
        chr <- m$chromosome[[i]]
        if (length(chr) > 6L)
            chr <- c(chr[seq_len(5L)], "...", chr[length(chr)])
        cat("  ", idx[i], " ", basename(m$filename[i]), " (",
            prettyNum(m$n_variant[i], big.mark=",", scientific=FALSE),
            " variants, chromosome ", paste(chr, collapse=","), ")\n", sep="")
    }
    invisible()
})

# The annotation GDS files as a SeqAnnotGDSList object: the files given by
#   names are opened (to be closed by the caller), and an opened GDS file or
#   a list of opened GDS files is wrapped
.as_annot_gds <- function(annot_gds, verbose=TRUE)
{
    if (is.character(annot_gds))
    {
        ans <- seqAnnotOpen(annot_gds)
        if (isTRUE(verbose))
        {
            m <- mcols(ans)
            for (i in seq_along(ans))
            {
                .cat("Open ", sQuote(basename(m$filename[i])), " [",
                    prettyNum(m$n_variant[i], big.mark=",", scientific=FALSE),
                    " variants]")
            }
        }
        return(ans)
    }
    # check
    if (inherits(annot_gds, "SeqVarGDSClass"))
        annot_gds <- list(annot_gds)
    if (!is(annot_gds, "SeqAnnotGDSList"))
    {
        if (!is.list(annot_gds))
        {
            stop("'annot_gds' should be GDS file names, a SeqAnnotGDSList ",
                "object or a list of opened GDS files.")
        }
        for (i in seq_along(annot_gds))
        {
            if (!inherits(annot_gds[[i]], "SeqVarGDSClass"))
                stop("annot_gds[[", i, "]] is not a GDS object.")
        }
        annot_gds <- .new_annot_gds(annot_gds)
    }
    if (length(annot_gds) == 0L)
        stop("No annotation gds file!")
    annot_gds
}
