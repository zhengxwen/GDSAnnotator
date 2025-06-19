#######################################################################
#
# Package Name: GDSAnnotator
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Miscellaneous functions
#


# Call base::table
.table_var <- function(x, ...)
{
    if (!is.list(x) || inherits(x, "SeqVarDataList"))
        x <- list(x)
    for (i in seq_along(x))
    {
        z <- x[[i]]
        if (inherits(z, "SeqVarDataList")) x[[i]] <- z$data
    }
    x <- append(x, list(exclude=NULL))
    do.call(base::table, x)
}

# Return the counts of unique values in a GDS node
seqValueCounts <- function(gdsfile, varnm, FUN=NULL, ..., bsize=100000L,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(varnm), length(varnm)>0L)
    stopifnot(is.null(FUN) || is.function(FUN))
    stopifnot(is.numeric(bsize), length(bsize)==1L, bsize>0L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    # check gdsfile
    if (is.character(gdsfile))
    {
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    }

    # user-defined function?
    varnm_name <- names(varnm)
    if (is.null(FUN))
    {
        FUN <- .table_var
        if (!is.null(names(varnm))) names(varnm) <- NULL
    }

    # process
    lst <- seqBlockApply(gdsfile, varnm, FUN,
        as.is="list", bsize=bsize, .tolist=FALSE, .progress=verbose, ...)

    # merge
    nm_lst <- lapply(seq_along(dimnames(lst[[1L]])), function(i)
    {
        sort(unique(unlist(lapply(lst, function(x) dimnames(x)[[i]]))),
            na.last=TRUE)
    })
    i_na <- vapply(nm_lst, function(nm) which(is.na(nm))[1L], 0L)
    cnt <- array(0L, dim=lengths(nm_lst))
    for (i in seq_along(lst))
    {
        ss <- dimnames(lst[[i]])
        ii <- lapply(seq_along(ss), function(j) {
            k <- match(ss[[j]], nm_lst[[j]])
            if (anyNA(k)) k[is.na(k)] <- i_na[j]  # check NA
            k
        })
        v <- do.call(`[`, c(list(cnt), ii)) + unname(lst[[i]])
        cnt <- do.call(`[<-`, c(list(cnt), ii, list(v)))
    }

    # output
    dimnm <- nm_lst
    names(dimnm) <- varnm_name
    ans <- array(cnt, dim=dim(cnt), dimnames=dimnm)
    class(ans) <- "table"
    ans
}



# Group the variants and return a SeqUnitListClass object
seqGroup <- function(gdsfile, verbose=TRUE)
{
    # check
    stopifnot(is.logical(verbose), length(verbose)==1L)
    if (is.character(gdsfile))
    {
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    }
}

