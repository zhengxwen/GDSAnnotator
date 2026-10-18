#######################################################################
#
# Package Name: GDSAnnotator
# Copyright (C) 2025-2026    Xiuwen Zheng
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Miscellaneous functions
#


# Return all annotation in the INFO field
#   gdsfile: a GDS file name or an opened GDS file
# return a DataFrame with the name, type, trait and description of each node
seqAnnotList <- function(gdsfile)
{
    # check gdsfile
    if (is.character(gdsfile))
    {
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    }
    # process
    nd_info <- index.gdsn(gdsfile, "annotation/info")
    nm <- ls.gdsn(nd_info, recursive=TRUE, include.dirs=FALSE)
    v <- lapply(nm, function(s) {
        nd <- index.gdsn(nd_info, s)
        dp <- objdesp.gdsn(nd)
        s <- get.attr.gdsn(nd)$Description
        if (is.null(s)) s <- NA_character_
        data.frame(type=dp$type, trait=dp$trait, description=s[1L])
    })
    df <- do.call(rbind, v)
    # output
    DataFrame(name=nm, df)
}


# Call base::table
.table_var <- function(x, ...)
{
    if (!is.list(x) || inherits(x, "SeqVarDataList"))
        x <- list(x)
    x <- lapply(x, function(z)
        if (inherits(z, "SeqVarDataList")) z$data else z)
    x <- append(x, list(exclude=NULL))
    do.call(base::table, x)
}

# Call base::table for per-variant counting
.table_var2 <- function(x, ...)
{
    if (!is.list(x) || inherits(x, "SeqVarDataList"))
        x <- list(x)
    # if SeqVarDataList
    a <- vapply(x, inherits, FALSE, what="SeqVarDataList")
    if (any(a))
    {
        ns <- x[[which(a)[1L]]]$length
        b <- data.frame(.index=rep.int(seq_along(ns), times=ns))
        b <- cbind(b, as.data.frame(lapply(x[a], function(z) z$data)))
        a <- !duplicated(b)
        x <- lapply(x, function(z)
            if (inherits(z, "SeqVarDataList")) z$data[a] else z)
    }
    x <- append(x, list(exclude=NULL))
    do.call(base::table, x)
}

# Merge tables
.merge_table <- function(lst, var_name)
{
    # merge objects of class 'table'
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
    names(dimnm) <- var_name
    ans <- array(cnt, dim=dim(cnt), dimnames=dimnm)
    class(ans) <- "table"
    ans
}


# check the GDS file names
.check_gds_fn <- function(fn, dup=TRUE)
{
    if (length(fn) == 0L)
        stop("'gdsfile' should be a file name.")
    if (anyNA(fn))
        stop("'gdsfile' should not contain NA.")
    if (!dup && anyDuplicated(fn))
        stop("'gdsfile' should not contain duplicated file names.")
    invisible()
}

# check that the tables can be merged by .merge_table()
.check_tables <- function(lst)
{
    z <- vapply(lst, is.table, FALSE)
    if (!all(z)) stop("All internal lst[[...]] should be a 'table' object.")
    z <- vapply(lst, function(x) length(dimnames(x)), 0L)
    if (anyNA(z) || any(z!=z[1L]))
        stop("All internal lst[[...]] should have the same dimension.")
    invisible()
}

# the tables of the values in a GDS file, one per block
.value_counts_blocks <- function(gdsfile, varnm, FUN, per_variant, parallel,
    bsize, verbose, ...)
{
    if (is.null(FUN))
    {
        if (length(list(...)) > 0L)
            stop("User-defined parameters should be none when 'FUN=NULL'.")
        FUN <- if (isTRUE(per_variant)) .table_var2 else .table_var
    } else {
        user_fun <- FUN
        FUN <- function(x, ...)
        {
            v <- user_fun(x, ...)
            if (!inherits(v, "table")) v <- .table_var(v, ...)
            v
        }
    }
    # blocking process
    seqBlockApply(gdsfile, varnm, FUN, as.is="list", parallel=parallel,
        bsize=bsize, .tolist=FALSE, .progress=verbose, ...)
}

# Return the counts of unique values in a GDS node
#   gdsfile: the GDS file name(s), or an opened GDS file
#   varnm: the GDS node names, with the names of the dimensions if named
#   use_info: if TRUE, the nodes are under "annotation/info"
#   FUN: NULL or a user-defined function returning a 'table' object (or the
#       values to be tabulated) for each block
#   per_variant: if TRUE, the same values of a variant are counted once
#   parallel, bsize: passed to seqBlockApply()
#   verbose: if TRUE, show the progress information
#   ...: passed to 'FUN'
seqValueCounts <- function(gdsfile, varnm, use_info=TRUE, FUN=NULL,
    per_variant=FALSE, parallel=FALSE, bsize=100000L, verbose=TRUE, ...)
{
    # check
    stopifnot(is.character(varnm), length(varnm)>0L)
    stopifnot(is.logical(use_info), length(use_info)==1L)
    stopifnot(is.null(FUN) || is.function(FUN))
    stopifnot(is.logical(per_variant), length(per_variant)==1L)
    stopifnot(is.numeric(bsize), length(bsize)==1L, bsize>0L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    var_name <- names(varnm)
    if (is.null(var_name)) var_name <- basename(varnm)
    # check gdsfile
    if (is.character(gdsfile))
    {
        .check_gds_fn(gdsfile, dup=FALSE)
        if (length(gdsfile) > 1L)
        {
            # count each file and merge the counts
            lst <- lapply(seq_along(gdsfile), function(i)
            {
                if (isTRUE(verbose))
                    message("[", i, "/", length(gdsfile), "] ", appendLF=FALSE)
                seqValueCounts(gdsfile[i], varnm, use_info=use_info, FUN=FUN,
                    per_variant=per_variant, parallel=parallel, bsize=bsize,
                    verbose=verbose, ...)
            })
            return(.merge_table(lst, var_name))
        }
        # when length(gdsfile)==1
        if (isTRUE(verbose))
            .cat("Open ", sQuote(basename(gdsfile)))
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    }
    names(varnm) <- NULL
    if (isTRUE(use_info))
        varnm <- paste0("annotation/info/", varnm)
    # the table of each block
    lst <- .value_counts_blocks(gdsfile, varnm, FUN, per_variant, parallel,
        bsize, verbose, ...)
    lst <- lst[!vapply(lst, is.null, FALSE)]
    if (length(lst)==0L) return(NULL)
    .check_tables(lst)
    # merge objects of class 'table'
    .merge_table(lst, var_name)
}
