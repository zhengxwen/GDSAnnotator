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
        # open the GDS file(s)
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
    # match chromosome names in the annotation GDS files, ignoring the "chr"
    #   prefix; 'gds_chr' is the chromosome label as stored in each file
    s <- strsplit(sub("^chr", "", names(annot_gds)), "&", fixed=TRUE)
    gds_idx <- rep(seq_along(annot_gds), times=lengths(s))
    gds_chr <- unlist(s)
    names(gds_idx) <- .norm_chr(gds_chr)
    # process each chromosome subset
    ans <- vector("list", length(chr_lst))
    for (i in seq_along(chr_lst))
    {
        ch <- chr_lst[i]
        # find the GDS file(s) for this chromosome
        gds_ii <- which(names(gds_idx) == .norm_chr(ch))
        # indices for this chromosome
        idx <- which(chr == ch)
        if (!length(gds_ii))
        {
            # no file for this chromosome: all the variants are not found
            d <- ann_pos_allele(annot_gds[[1L]], ch, pos[idx], ref[idx],
                alt[idx], varnm, verbose=FALSE)
            if (!length(varnm)) d$file_idx <- Rle(1L, nrow(d))
            d$..idx <- idx
            d$..no <- NULL
            ans[[i]] <- d
            next
        }
        # process each GDS file for this chromosome: the variants not found
        #   in a file are looked up in the next file(s)
        v <- vector("list", length(gds_ii))
        for (j in seq_along(gds_ii))
        {
            if (!length(idx)) break
            k <- gds_idx[gds_ii[j]]
            # annotate variants in this chromosome subset, using the
            #   chromosome label of the file
            d <- ann_pos_allele(annot_gds[[k]], gds_chr[gds_ii[j]],
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
    stopifnot(is.character(varnm))
    # process
    ann_chr_pos_allele(chr, pos, ref, alt, annot_gds, varnm, verbose=verbose)
}


#######################################################################
# Write the annotations of the variants to a GDS file block by block

# split 'n' variants into blocks of 'bsize' variants: the index of the first
#   variant and the number of variants of each block
.blocks <- function(n, bsize)
{
    nblock <- max(1L, as.integer(ceiling(n / bsize)))
    st <- as.double(seq_len(nblock) - 1L) * bsize + 1
    st <- as.integer(pmin(st, n + 1))
    list(st=st, cnt=diff(c(st, n + 1L)))
}

# whether 'x' returned by seqGetData() is a vector with one value per variant
#   or a SeqVarDataList object with a vector of data
.is_vec <- function(x) is.atomic(x) && is.null(dim(x))
.is_dl <- function(x) inherits(x, "SeqVarDataList") && is.null(dim(x$data))

# the values for the variants 'k' (the indices in 'x', NA for a variant not
#   found) in the compact list(length, data) representation, where 'x' is a
#   vector or a SeqVarDataList object returned by seqGetData(, .padNA=FALSE):
#   the values of a SeqVarDataList object are kept as they are stored, and a
#   variant not found has no value; a vector has one value per variant, and a
#   variant not found has a missing value
.annot_subset <- function(x, k)
{
    in_order <- !anyNA(k) && !is.unsorted(k, strictly=TRUE)
    if (.is_dl(x))
    {
        len <- x$length
        if (in_order && length(k)==length(len))
            return(list(length=len, data=x$data))
        n <- len[k]
        n[is.na(n)] <- 0L
        j <- which(n > 0L)
        # the offset of the values of each variant in x$data
        off <- c(0, cumsum(as.double(len)))[k[j]]
        list(length=n, data=x$data[rep.int(off, n[j]) + sequence(n[j])])
    } else {
        if (!in_order || length(k)!=length(x)) x <- x[k]
        list(length=rep.int(1L, length(k)), data=x)
    }
}

# the values of the GDS variable 'nm' in the annotation file 'f' for the
#   variants 'ii' (the indices in 'f', NA for a variant not found, possibly
#   unsorted and repeated)
.annot_get <- function(f, nm, ii)
{
    k <- seqSetFilter(f, variant.sel=ii, ret.idx=TRUE, warn=FALSE,
        verbose=FALSE)$variant_idx
    x <- NULL
    if (length(objdesp.gdsn(index.gdsn(f, nm))$dim) == 1L)
        x <- seqGetData(f, nm, .padNA=FALSE)
    if (.is_vec(x) || .is_dl(x))
        return(.annot_subset(x, k))
    # other forms (e.g., a matrix): one list element per variant
    v <- seqGetData(f, nm, .tolist=TRUE)[k]
    v <- lapply(v, function(x) unlist(x, use.names=FALSE))
    list(length=lengths(v), data=unlist(v, use.names=FALSE))
}

# the value of each variable for no variant as returned by seqGetData(), which
#   gives the type of the values of the variants not found; NULL if
#   it is neither a vector nor a SeqVarDataList object (e.g., a matrix), or if
#   no file has the variable
.annot_template <- function(annot_gds, varnm)
{
    ans <- vector("list", length(varnm))
    empty <- logical(length(annot_gds))  # TRUE if no variant is selected
    for (i in seq_along(varnm))
    {
        for (j in seq_along(annot_gds))
        {
            f <- annot_gds[[j]]
            if (is.null(index.gdsn(f, varnm[i], silent=TRUE))) next
            if (!empty[j])
            {
                seqSetFilter(f, variant.sel=integer(), warn=FALSE,
                    verbose=FALSE)
                empty[j] <- TRUE
            }
            x <- seqGetData(f, varnm[i], .padNA=FALSE)
            if (.is_vec(x) || .is_dl(x)) ans[i] <- list(x)
            break
        }
    }
    ans
}

# the levels of each factor variable (e.g., "annotation/filter") in all the
#   annotation files, in the order of the files; NULL if it is not a factor
.annot_levels <- function(annot_gds, varnm)
{
    lapply(varnm, function(nm)
    {
        lv <- NULL
        for (f in annot_gds)
        {
            n <- index.gdsn(f, nm, silent=TRUE)
            if (is.null(n)) next
            a <- get.attr.gdsn(n)
            if ("factor" %in% a$R.class)
                lv <- union(lv, as.character(a$R.levels))
        }
        lv
    })
}

# the values of the GDS variable 'nm' for a block of variants in the compact
#   list(length, data) representation, where 'fi' and 'vi' are the indices of
#   the annotation file and of the variant in that file (NA if not found),
#   get_file(i) returns the i-th annotation file, and 'tmpl' is the value of
#   'nm' for no variant (see .annot_template())
.annot_block <- function(get_file, nm, fi, vi, tmpl)
{
    found <- !is.na(vi)
    if (!any(found) && (!is.null(tmpl) || !length(vi)))
    {
        # no variant found, no file to read
        v <- if (length(vi)) .annot_subset(tmpl, vi) else
            list(length=integer(), data=NULL)
    } else {
        # a variant not found joins the run of the previous variant found (or
        #   of the first one), so that a file is read once per run
        if (any(found) && !all(found))
            fi <- fi[which(found)][pmax(cumsum(found), 1L)]
        # the runs of consecutive variants in the same file
        n <- length(fi)
        ed <- c(which(fi[-1L] != fi[-n]), n)
        st <- c(1L, ed[-length(ed)] + 1L)
        v <- lapply(seq_along(st), function(j)
            .annot_get(get_file(fi[st[j]]), nm, vi[seq.int(st[j], ed[j])]))
        if (length(v) == 1L)
        {
            v <- v[[1L]]
        } else {
            v <- list(
                length=unlist(lapply(v, `[[`, "length"), use.names=FALSE),
                data=unlist(lapply(v, `[[`, "data"), use.names=FALSE))
        }
    }
    structure(v, class="SeqVarDataList")
}


#######################################################################
# Write the annotation variables in parallel: each process writes one
#   variable at a time to a temporary GDS file, and the main process copies
#   the finished variables to the output file in the order of 'varnm'

# the number of processes for 'parallel'
.num_parallel <- function(parallel)
{
    if (is.null(parallel) || isFALSE(parallel)) return(1L)
    if (inherits(parallel, "cluster")) return(length(parallel))
    if (isTRUE(parallel))
    {
        n <- detectCores() - 1L
        return(if (is.na(n) || n < 2L) 2L else n)
    }
    if (is.numeric(parallel) && length(parallel)==1L && !is.na(parallel) &&
        parallel >= 1)
    {
        return(as.integer(parallel))
    }
    stop("'parallel' should be FALSE, TRUE, a positive integer or ",
        "a cluster object.")
}

# the total compressed size of each variable in the annotation files, so that
#   the processes start with the largest variables
.annot_size <- function(annot_gds, varnm)
{
    vapply(varnm, function(nm)
    {
        s <- 0
        for (f in annot_gds)
        {
            for (p in c(nm, paste0(dirname(nm), "/@", basename(nm))))
            {
                n <- index.gdsn(f, p, silent=TRUE)
                if (!is.null(n))
                {
                    z <- objdesp.gdsn(n)$size
                    if (isTRUE(is.finite(z))) s <- s + z
                }
            }
        }
        s
    }, 0, USE.NAMES=FALSE)
}

# write the file and variant indices of the variants to a temporary GDS file,
#   which each process reads block by block, and return the last block in
#   which each annotation file has a variant found (0 if none)
.ann_map_write <- function(fn, file_ids, variant_idx, blk, nfile)
{
    f <- createfn.gds(fn)
    on.exit(closefn.gds(f))
    n1 <- add.gdsn(f, "file_idx", storage="int32", valdim=0L,
        compress="LZ4_RA")
    n2 <- add.gdsn(f, "variant_idx", storage="int32", valdim=0L,
        compress="LZ4_RA")
    last <- integer(nfile)
    for (k in seq_along(blk$st))
    {
        s <- seq.int(blk$st[k], length.out=blk$cnt[k])
        fi <- as.integer(file_ids[s])
        vi <- variant_idx[s]
        append.gdsn(n1, fi)
        append.gdsn(n2, vi)
        last[unique(fi[!is.na(vi)])] <- k
    }
    readmode.gdsn(n1)
    readmode.gdsn(n2)
    last
}

# the VCF type of a GDS node
.vcf_type <- function(node)
{
    switch(as.character(objdesp.gdsn(node)$type[1L]),
        Raw="Integer", Integer="Integer", Real="Float", Factor="String",
        String="String", Logical="Flag", ".")
}

# make the integer node 'n' a factor with the levels 'levels' (if not NULL),
#   in the same way as gdsfmt stores a factor
.put_levels <- function(n, levels)
{
    if (!is.null(levels))
    {
        put.attr.gdsn(n, "R.class", "factor")
        put.attr.gdsn(n, "R.levels", levels)
    }
    invisible()
}

# write the blocks of values returned by gen(k) (NULL if no more block) to
#   the node 'nm' and its index node '@nm' in the folder 'node', in the same
#   way as SeqArray::seqAddValue() for "annotation/info/nm" with its default
#   arguments, except that a factor is stored as the integer codes of
#   'levels' (the levels in all the annotation files) instead of strings, and
#   that the type of 'tmpl' (see .annot_template()) is kept when there is no
#   value at all; it is used with and without 'parallel', so the output does
#   not depend on it
.ann_add_stream <- function(node, nm, gen, nvar, tmpl=NULL, levels=NULL,
    compress="LZMA_RA")
{
    # the storage mode of the values 'x'
    stmode <- function(x)
    {
        s <- if (is.factor(x) || is.character(x)) "string" else storage.mode(x)
        if (identical(s, "double")) "float32" else s
    }
    n <- nidx <- st <- NULL
    cnt <- 0
    pending <- list()  # the lengths of the leading blocks without value
    k <- 0L
    repeat {
        k <- k + 1L
        v <- gen(k)
        if (is.null(v)) break
        len <- v$length
        len[is.na(len) | len<0L] <- 0L
        dat <- v$data
        if (is.null(dat)) dat <- integer()
        if (!is.null(levels)) dat <- match(as.character(dat), levels)
        remove(v)
        cnt <- cnt + length(len)
        if (is.null(n))
        {
            if (!length(dat))
            {
                pending[[length(pending) + 1L]] <- len
                next
            }
            # the first block with values defines the storage mode
            st <- stmode(dat)
            n <- add.gdsn(node, nm, dat, storage=st, compress=compress,
                replace=TRUE)
            .put_levels(n, levels)
            nidx <- add.gdsn(node, paste0("@", nm),
                unlist(pending, use.names=FALSE), storage="int",
                compress=compress, replace=TRUE, visible=FALSE)
            pending <- list()
            append.gdsn(nidx, len)
        } else {
            if (is.factor(dat)) dat <- as.character(dat)
            if (length(dat) && is.character(dat) != (st == "string"))
            {
                stop("The data block ", k,
                    " does not have the same data type as the previous ones.")
            }
            if (length(dat)) append.gdsn(n, dat)
            append.gdsn(nidx, len)
        }
        remove(len, dat)
    }
    if (is.null(n))
    {
        # no value at all in any block: the type of 'tmpl'
        dat <- integer()
        if (is.null(levels) && .is_dl(tmpl))
            dat <- tmpl$data[0L]
        n <- add.gdsn(node, nm, dat, storage=stmode(dat), compress=compress,
            replace=TRUE)
        .put_levels(n, levels)
        nidx <- add.gdsn(node, paste0("@", nm),
            unlist(pending, use.names=FALSE), storage="int",
            compress=compress, replace=TRUE, visible=FALSE)
    }
    readmode.gdsn(n)
    readmode.gdsn(nidx)
    if (cnt != nvar)
        stop("The total number of variants should be ", nvar, ", but ", cnt)
    put.attr.gdsn(n, "Number", ".")
    put.attr.gdsn(n, "Type", .vcf_type(n))
    put.attr.gdsn(n, "Description", "")
    digest.gdsn(n, algo="md5", action="add")
    digest.gdsn(nidx, algo="md5", action="add")
    invisible(n)
}

# run in a process for the j-th job: write the i-th variable (i = ord[j]) to
#   the temporary file tmp_fn[i]; the variants are read from the file 'map_fn'
#   block by block, and an annotation file is opened when needed and closed
#   after the last block using it, so that the memory usage depends neither
#   on the total number of variants nor on the number of annotation files
.ann_par_var <- function(j, ord, map_fn, annot_fn, last_blk, varnm, tmpl,
    lv, bsize, tmp_fn)
{
    i <- ord[j]
    mf <- openfn.gds(map_fn, allow.duplicate=TRUE)
    on.exit(closefn.gds(mf))
    n_fi <- index.gdsn(mf, "file_idx")
    n_vi <- index.gdsn(mf, "variant_idx")
    nvar <- objdesp.gdsn(n_vi)$dim
    blk <- .blocks(nvar, bsize)
    # the opened annotation files, named by their index
    e <- new.env(parent=emptyenv())
    on.exit(for (s in ls(e)) seqClose(e[[s]]), add=TRUE)
    get_file <- function(k)
    {
        s <- as.character(k)
        if (is.null(e[[s]]))
            assign(s, seqOpen(annot_fn[k], allow.duplicate=TRUE), envir=e)
        e[[s]]
    }
    gen <- function(k)
    {
        if (k > length(blk$st)) return(NULL)  # no more block
        st <- blk$st[k]; cnt <- blk$cnt[k]
        v <- .annot_block(get_file, varnm[i],
            read.gdsn(n_fi, start=st, count=cnt),
            read.gdsn(n_vi, start=st, count=cnt), tmpl[[i]])
        # close the files not needed by the remaining blocks
        for (s in ls(e))
        {
            if (last_blk[as.integer(s)] <= k)
            {
                seqClose(e[[s]])
                rm(list=s, envir=e)
            }
        }
        v
    }
    f <- createfn.gds(tmp_fn[i])
    on.exit(closefn.gds(f), add=TRUE)
    # suppress the warnings "Missing characters are converted to" as in the
    #   serial processing, since NA is not allowed in a GDS string variable
    suppressWarnings(.ann_add_stream(f$root, "data", gen, nvar, tmpl[[i]],
        lv[[i]]))
    i
}

# copy the variable written by .ann_par_var() in the file 'fn' to 'path' in
#   the output file, replacing the existing node at the same position as
#   add.gdsn(, replace=TRUE) does in .ann_add_stream()
.ann_copy_node <- function(gdsfile, path, fn)
{
    f <- openfn.gds(fn, allow.duplicate=TRUE)
    on.exit(closefn.gds(f))
    folder <- index.gdsn(gdsfile, dirname(path))
    nm <- basename(path)
    for (s in c("", "@"))
    {
        src <- index.gdsn(f, paste0(s, "data"))
        old <- index.gdsn(folder, paste0(s, nm), silent=TRUE)
        if (is.null(old))
        {
            copyto.gdsn(folder, src, paste0(s, nm))
        } else {
            tmp <- paste0(s, nm, ".tmp")
            while (!is.null(index.gdsn(folder, tmp, silent=TRUE)))
                tmp <- paste0(tmp, "_")
            copyto.gdsn(folder, src, tmp)
            moveto.gdsn(index.gdsn(folder, tmp), old, relpos="replace+rename")
        }
    }
    invisible()
}

# write the variables 'varnm' to 'rootnm/colnm' in the output file with
#   'npar' processes, starting with the largest variables (see 'size')
.ann_par_write <- function(outgds, rootnm, colnm, varnm, tmpl, lv, size,
    map_fn, annot_fn, last_blk, bsize, parallel, npar, verbose)
{
    # the temporary files of the variables, next to the output file
    fn <- outgds$filename
    tmp_fn <- tempfile(paste0(basename(fn), ".tmp", seq_along(varnm), "_"),
        dirname(fn), ".gds")
    on.exit(unlink(c(map_fn, tmp_fn), force=TRUE))
    # the processes are started with makeCluster() instead of forking, so
    #   they do not share the memory pages of the current process
    cl <- parallel
    if (!inherits(cl, "cluster"))
    {
        cl <- makeCluster(npar)
        on.exit(stopCluster(cl), add=TRUE)
        clusterCall(cl, .libPaths, .libPaths())
    }
    if (verbose)
        .cat("    ", length(varnm), " variables, ", npar, " processes")
    # copy each finished variable as soon as the previous ones are copied
    st <- new.env(parent=emptyenv())
    st$done <- logical(length(varnm))
    st$i <- 1L
    combine <- function(i)
    {
        st$done[i] <- TRUE
        while (st$i <= length(varnm) && st$done[st$i])
        {
            j <- st$i
            .ann_copy_node(outgds, paste0(rootnm, "/", colnm[j]), tmp_fn[j])
            unlink(tmp_fn[j], force=TRUE)
            if (verbose)
            {
                .cat("    added ", j, "/", length(varnm), " ",
                    sQuote(colnm[j]), " (", tm(), ")")
            }
            st$i <- j + 1L
        }
        invisible()
    }
    seqParallel(cl, length(varnm), FUN=.ann_par_var, split="none",
        .combine=combine, ord=order(size, decreasing=TRUE), map_fn=map_fn,
        annot_fn=annot_fn, last_blk=last_blk, varnm=varnm, tmpl=tmpl, lv=lv,
        bsize=bsize, tmp_fn=tmp_fn)
    invisible()
}


# Annotate a GDS file
ann_gdsfile <- function(object, annot_gds, varnm, add_to_gds=FALSE,
    no_sample=TRUE, root="", cleanup=FALSE, bsize=1000000L, parallel=FALSE,
    ..., verbose=TRUE)
{
    # check
    stopifnot(inherits(object, "SeqVarGDSClass"))
    stopifnot(is.numeric(bsize), length(bsize)==1L, !is.na(bsize), bsize>0L)
    bsize <- as.double(bsize)
    npar <- .num_parallel(parallel)
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
    if (isTRUE(verbose)) message("Reading chromosome", appendLF=FALSE)
    chr <- seqGetData(object, "$chromosome")
    if (isTRUE(verbose)) message(", position", appendLF=FALSE)
    pos <- seqGetData(object, "position")
    if (isTRUE(verbose)) message(", reference allele", appendLF=FALSE)
    ref <- seqGetData(object, "$ref")
    if (isTRUE(verbose)) message(", alternative allele", appendLF=FALSE)
    alt <- seqGetData(object, "$alt")
    if (isTRUE(verbose)) message()
    # process
    if (isTRUE(verbose)) message("Processing annotation ...")
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
        # create the parent folders for the nested variable names (e.g.,
        #   "CSQ.list/Allele" as returned by seqAnnotList()), since
        #   seqAddValue() requires the parent node to exist
        for (s in setdiff(unique(dirname(colnm)), "."))
        {
            ss <- strsplit(s, "/", fixed=TRUE)[[1L]]
            for (j in seq_along(ss))
            {
                nd <- paste0(rootnm, "/",
                    paste(ss[seq_len(j)], collapse="/"))
                if (!exist.gdsn(outgds, nd))
                    seqAddValue(outgds, nd, NULL, verbose=FALSE)
            }
        }
        n_total <- nrow(map)
        # group found variants by annotation file, map$file_idx has no NA
        file_ids <- map$file_idx
        if (!inherits(file_ids, "Rle"))
            file_ids <- Rle(as.integer(file_ids))
        # split the variants into blocks to reduce the memory usage, and the
        #   values of a block are kept as vectors instead of a list with one
        #   element per variant (which needs much more memory)
        blk <- .blocks(n_total, bsize)
        nblock <- length(blk$st)
        # the value of each variable for no variant, and the factor levels
        tmpl <- .annot_template(annot_gds, varnm)
        lv <- .annot_levels(annot_gds, varnm)
        npar <- min(npar, length(varnm))
        if (npar > 1L)
        {
            # each process opens the annotation files again
            annot_fn <- normalizePath(vapply(annot_gds,
                function(f) f$filename, ""))
            size <- .annot_size(annot_gds, varnm)
            if (if_close_gds)
            {
                .close_annot_gds(annot_gds)
                annot_gds <- list()  # not to be closed again on exit
            }
            # the processes read the file and variant indices from a file
            fn <- outgds$filename
            map_fn <- tempfile(paste0(basename(fn), ".map"), dirname(fn),
                ".gds")
            on.exit(unlink(map_fn, force=TRUE), add=TRUE)
            last_blk <- .ann_map_write(map_fn, file_ids, map$variant_idx,
                blk, length(annot_fn))
            remove(map, file_ids)
            .ann_par_write(outgds, rootnm, colnm, varnm, tmpl, lv, size,
                map_fn, annot_fn, last_blk, bsize, parallel, npar, verbose)
            return(invisible())
        }
        dot_block <- if (nblock > 1L)
            unique(as.integer(round(seq(1L, nblock,
                length.out=min(10L, nblock))))) else integer()
        # get the values of the variable 'varnm[i]' for the k-th block, in
        #   the compact list(length, data) representation
        get_file <- function(j) annot_gds[[j]]
        get_block <- function(i, k)
        {
            s <- seq.int(blk$st[k], length.out=blk$cnt[k])
            .annot_block(get_file, varnm[i], as.integer(file_ids[s]),
                map$variant_idx[s], tmpl[[i]])
        }
        # add each annotation variable
        for (i in seq_along(varnm))
        {
            if (verbose)
            {
                message("    adding ", i, "/", length(varnm), " ",
                    sQuote(colnm[i]), " (", tm(), ")\n    ", appendLF=FALSE)
            }
            # write block by block, so that only one block is held in memory
            #   instead of the annotation of all variants
            gen <- function(k)
            {
                if (k > nblock) return(NULL)  # no more block
                if (verbose && (k %in% dot_block))
                    message(".", appendLF=FALSE)
                get_block(i, k)
            }
            # suppress the warnings from writing the annotation data, e.g.,
            #   "Missing characters are converted to" since NA is not allowed
            #   in a GDS string variable
            p <- paste0(rootnm, "/", colnm[i])
            n <- suppressWarnings(.ann_add_stream(
                index.gdsn(outgds, dirname(p)), basename(p), gen, n_total,
                tmpl[[i]], lv[[i]]))
            if (verbose)
                .cat(capture.output(print(n, attribute=FALSE)))
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
        message("[", length(varnm), "] ", appendLF=FALSE)
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


# Annotate an IRanges object
ann_IRanges <- function(object, annot_gds, varnm, chr, ..., verbose=TRUE)
{
    # check
    if (missing(chr))
        stop("'chr' should be specified for an IRanges object.")
    stopifnot(is.character(chr), length(chr)>0L)
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
        message("[", length(varnm), "] ", appendLF=FALSE)
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
ann_VCF <- function(object, annot_gds, varnm, ..., verbose=TRUE)
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
#   ann_VCF(), any other class is an error
ann_any <- function(object, annot_gds, varnm, ..., verbose=TRUE)
{
    if (inherits(object, "VCF"))
    {
        ann_VCF(object, annot_gds, varnm, ..., verbose=verbose)
    } else {
        stop("No seqAnnotate() method for an object of class ",
            sQuote(class(object)[1L]), ".")
    }
}


# Set methods
setMethod("seqAnnotate", signature(object="data.frame"), ann_dataframe)
setMethod("seqAnnotate", signature(object="DataFrame"), ann_dataframe)
setMethod("seqAnnotate", signature(object="SeqVarGDSClass"), ann_gdsfile)
setMethod("seqAnnotate", signature(object="character"), ann_variant)
setMethod("seqAnnotate", signature(object="GRanges"), ann_GRanges)
setMethod("seqAnnotate", signature(object="GRangesList"), ann_GRanges)
setMethod("seqAnnotate", signature(object="IRanges"), ann_IRanges)
setMethod("seqAnnotate", signature(object="ANY"), ann_any)


# Annotate a VCF file: read it with VariantAnnotation::readVcf() and
#   annotate the VCF object
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
    ann_VCF(vcf, annot_gds, varnm, ..., verbose=verbose)
}


# Annotate a GDS file with a file name input
seqAnnotateGDS <- function(gds_fn, annot_gds, varnm, add_to_gds=FALSE,
    no_sample=TRUE, root="", cleanup=TRUE, bsize=1000000L, parallel=FALSE,
    ..., verbose=TRUE)
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
    ann_gdsfile(gds, annot_gds, varnm, add_to_gds, no_sample, root,
        cleanup=cleanup, bsize=bsize, parallel=parallel, ..., verbose=verbose)
}
