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

# Close the GDS file(s)
.close_annot_gds <- function(annot_gds)
{
    lapply(annot_gds, seqClose)
    invisible()
}

# GDS nodes for the INFO field
.gds_varnm <- function(nm)
{
    if (identical(nm, "")) nm <- character()
    if (length(nm))
    {
        if (anyNA(nm)) stop("'varnm' should not have NA_character_.")
        s <- nm
        nm <- paste0("annotation/info/", nm)
        nm[s==":id"] <- "annotation/id"
        nm[s==":qual"] <- "annotation/qual"
        nm[s==":filter"] <- "annotation/filter"
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
            names(ans) <- varnm
        } else {
            ans <- DataFrame(file_idx=Rle(annot_gds_idx[1L], length(ii)),
                variant_idx=ii)
            if (!is.null(idx)) ans <- ans[idx, ]
        }
    } else {
        # no annotation
        v <- rep(NA, length(pos))
        ans <- lapply(seq_along(varnm), function(i) v)
        names(ans) <- varnm
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
    varnm <- .gds_varnm(varnm)
    # check annotation gds file(s)
    stopifnot(is.character(annot_gds), length(annot_gds)>0L)
    # open gds file(s)
    annot_gds <- .open_annot_gds(annot_gds)
    on.exit(.close_annot_gds(annot_gds))
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
    # process
    ann_chr_pos_allele(chr, pos, ref, alt, annot_gds, varnm, verbose=verbose)
}


# Annotate a GDS file
ann_gdsfile <- function(object, annot_gds, varnm, ..., verbose=TRUE)
{
    # check
    stopifnot(inherits(object, "SeqVarGDSClass"))
    # process
    map <- ann_chr_pos_allele(
        seqGetData(object, "$chromosome"),
        seqGetData(object, "position"),
        seqGetData(object, "allele"), NULL,
        annot_gds, "", verbose=verbose)
    map
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

