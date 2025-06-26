#######################################################################
#
# Package Name: GDSAnnotator
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Miscellaneous functions
#


# Return all annotation in the INFO field
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
    for (i in seq_along(x))
    {
        z <- x[[i]]
        if (inherits(z, "SeqVarDataList")) x[[i]] <- z$data
    }
    x <- append(x, list(exclude=NULL))
    do.call(base::table, x)
}

# Return the counts of unique values in a GDS node
seqValueCounts <- function(gdsfile, varnm, use_info=TRUE, FUN=NULL,
    bsize=100000L, verbose=TRUE, ...)
{
    # check
    stopifnot(is.character(varnm), length(varnm)>0L)
    stopifnot(is.logical(use_info), length(use_info)==1L)
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
    var_name <- names(varnm)
    if (isTRUE(use_info))
        varnm <- paste0("annotation/info/", varnm)
    if (!is.null(names(varnm))) names(varnm) <- NULL
    if (is.null(FUN))
    {
        FUN <- .table_var
        lst <- seqBlockApply(gdsfile, varnm, FUN,
            as.is="list", bsize=bsize, .tolist=FALSE, .progress=verbose, ...)
    } else if (is.function(FUN))
    {
        lst <- seqBlockApply(gdsfile, varnm,
            FUN = function(x, ...)
            {
                v <- FUN(x, ...)
                if (!inherits(v, "table")) v <- .table_var(v, ...)
                v
            }, as.is="list", bsize=bsize, .tolist=FALSE, .progress=verbose, ...)
    }

    # process
    i <- vapply(lst, is.null, FALSE)
    if (any(i)) lst <- lst[!i]
    if (length(lst)==0L) return(NULL)

    # check
    z <- vapply(lst, is.table, FALSE)
    if (!all(z)) stop("All internal lst[[...]] should be a 'table' object.")
    z <- vapply(lst, function(x) length(dimnames(x)), 0L)
    if (anyNA(z) || any(z!=z[1L]))
        stop("All internal lst[[...]] should have the same dimension.")

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


# Group the variants by annotation and return a SeqUnitListClass object
seqUnitGroupAnnot <- function(gdsfile, varnm, by=1L, cond=NULL,
    bsize=100000L, verbose=TRUE, ...)
{
    # check
    stopifnot(is.character(varnm), length(varnm)>0L)
    stopifnot(is.character(by) || is.numeric(by), length(by)>0L)
    if (is.character(by)) by <- match(by, varnm)
    if (anyNA(by) || !all(1<=by | by<=length(varnm)))
    {
        stop("'by' should be one or more of ", paste(varnm, collapse=", "),
            ", or the index/indices.")
    }
    cond_err <- "'cond' should be NULL or a list according to 'varnm'."
    if (is.list(cond))
    {
        if (length(cond) != length(varnm))
            stop(cond_err)
        a <- vapply(cond, function(x)
            is.null(x) || is.vector(x) || is.function(x), FALSE)
        if (!all(a))
        {
            a <- which(!a)[1L]
            stop("'cond[[", a, "]]' should be NULL, a vector or a function.")
        }
    } else if (!is.null(cond))
        stop(cond_err)

    stopifnot(is.logical(verbose), length(verbose)==1L)
    # check gdsfile
    if (is.character(gdsfile))
    {
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    }

    # initialize
    by2 <- setdiff(seq_along(varnm), by)

    # grouping function
    group <- function(x, ...)
    {
        # condition
        for (i in seq_along(cond))
        {
            v <- cond[[i]]
            if (is.vector(v))
            {
                if (inherits(x[[i]], "SeqVarDataList"))
                    x[[i]]$data <- x[[i]]$data %in% v
                else
                    x[[i]] <- x[[i]] %in% v
            } else if (is.function(v))
            {
                x[[i]] <- v(x[[i]], ...)
            }
        }
        # filter based on condition
        if (length(by2))
        {
            y <- lapply(x[by2], function(z) z$data)  # assuming SeqVarDataList
            y <- Reduce(`&`, y)
            for (i in by) x[[i]]$data[!y] <- NA  # exclude based on 'cond'
        }
        # group variable
        grp <- lapply(by, function(i) {
            s <- x[[i]]$data; s[s==""] <- NA; s
        })
        # indices
        ii <- rep(x[[length(x)]], times=x[[by[1L]]]$length)
        # grouping
        ans <- base::split(ii, grp, sep="\xFF")
        lapply(ans, function(z) sort(unique(z)))
    }

    # process
    lst <- seqBlockApply(gdsfile, c(varnm, "$variant_index"), FUN=group,
        as.is="list", bsize=bsize, .tolist=FALSE, .progress=verbose, ...)
    i <- vapply(lst, is.null, FALSE)
    if (any(i)) lst <- lst[!i]
    if (length(lst)==0L) return(NULL)

    # merge
    ans <- lst[[1L]]
    for (i in seq_along(lst)[-1L])
    {
        v <- lst[[i]]
        ss <- intersect(names(ans), names(v))
        for (s in ss)
            ans[[s]] <- sort(unique(c(ans[[s]], v[[s]])))
        if (length(ss))
            v <- v[-match(ss, names(v))]
        ans <- c(ans, v)
    }
    # sort by starting index/position
    ans <- ans[order(vapply(ans, `[`, 0L, i=1L))]

    # data.frame
    if (length(by) > 1L)
    {
        v <- strsplit(names(ans), '\xFF', fixed=TRUE)
        nc <- max(lengths(v))
        df <- lapply(seq_len(nc), function(i)
            vapply(v, `[`, "", i=i))
    } else {
        df <- list(names(ans))
    }
    # column names
    nm <- names(by)
    if (is.null(nm))
    {
        nm <- names(varnm)[by]
        if (is.null(nm)) nm <- rep("", length(by))
    }
    if (anyNA(nm)) nm[is.na(nm)] <- ""
    if (any(nm==""))
        nm[nm==""] <- paste0("group", seq_len(sum(nm=="")))
    names(df) <- nm
    # chr, position
    i1 <- vapply(ans, `[`, 0L, i=1L)
    i2 <- vapply(ans, function(z) z[max(length(z), 1L)], 0L)
    ii <- seqSetFilter(gdsfile, variant.sel=c(i1, i2), action="push+set",
        ret.idx=TRUE, warn=FALSE, verbose=FALSE)$variant_idx
    on.exit(seqFilterPop(gdsfile), add=TRUE)
    i1 <- seq_along(ans)
    df$chr <- seqGetData(gdsfile, "chromosome")[ii[i1]]
    pos <- seqGetData(gdsfile, "position")[ii]
    df$start <- pos[i1]
    df$end <- pos[seq.int(length(ans)+1L, length.out=length(ans))]

    # output
    ans <- list(desp=as.data.frame(df), index=ans)
    class(ans) <- "SeqUnitListClass"
    ans
}

