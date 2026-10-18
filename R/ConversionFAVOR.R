#######################################################################
#
# Package Name: GDSAnnotator
# Copyright (C) 2025-2026    Xiuwen Zheng
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Data import & Format conversion from FAVOR (https://favor.genohub.org)
#


#######################################################################
# Convert FAVOR CSV files to a SeqArray GDS file
#

# the columns of a FAVOR CSV file identifying the variants
.favor_col_csv <- c("position", "chromosome", "ref_vcf", "alt_vcf",
    "variant_vcf")

# the GDS nodes of the variants, as used by .gds_new()
.gds_var_nodes <- c("variant.id", "position", "chromosome", "allele",
    "annotation/id", "annotation/qual", "annotation/filter")

# check that all the CSV files have the same header line
.check_csv_header <- function(fn)
{
    if (length(fn) < 2L) return(invisible())
    hd <- lapply(fn, readLines, n=1L)
    if (length(hd[[1L]]) != 1L)
        stop("No header line in the first file of 'csv_fn'.")
    ok <- vapply(hd, identical, FALSE, y=hd[[1L]])
    if (!all(ok))
    {
        stop(sprintf("The header of 'csv_fn[%d]' differs from ",
            which(!ok)[1L]), "the header of the first file.")
    }
    invisible()
}

# append the variants and the annotations of a FAVOR CSV file
.favor_append <- function(outfile, nd_root, df, nm_others, first,
    use_float32, compress, verbose)
{
    if (!all(.favor_col_csv %in% colnames(df)))
    {
        stop("The csv header should contain all of ",
            paste(.favor_col_csv, collapse=", "), ".")
    }
    # basic
    n <- nrow(df)
    val <- list(
        if ("vid" %in% colnames(df)) df$vid else seq_len(n),
        df$position, df$chromosome, paste0(df$ref_vcf, ",", df$alt_vcf),
        df$variant_vcf, rep(NaN, n), rep(1L, n))
    for (i in seq_along(val))
        append.gdsn(index.gdsn(outfile, .gds_var_nodes[i]), val[[i]])
    # annotation
    for (nm in nm_others)
    {
        if (verbose) .cat("    [adding ", nm, "]\t")
        v <- .na_to_empty(df[[nm]])
        if (first)
        {
            st <- typeof(v)
            if (st=="double" && isTRUE(use_float32)) st <- "float32"
            nd <- add.gdsn(nd_root, nm, v, storage=st, compress=compress)
        } else {
            nd <- append.gdsn(index.gdsn(nd_root, nm), v)
        }
        if (verbose)
            message("    ", paste(capture.output(nd), collapse="\n    "))
    }
    invisible()
}

# switch the GDS nodes to the read mode, and add the MD5 digests
.finish_nodes <- function(nodes)
{
    for (nd in nodes)
    {
        readmode.gdsn(nd)
        digest.gdsn(nd, "md5", action="add")
    }
    invisible()
}

# Convert the FAVOR CSV files to a SeqArray GDS file
#   csv_fn: the file name(s) of the FAVOR CSV files (possibly gzip-compressed)
#       with the same header
#   out_fn: the output GDS file name
#   compress: the compression method of the output file, "LZMA", "ZIP" or
#       "none"
#   root: the folder under "annotation/info" for the annotations, "" for none
#   use_float32: if TRUE, store the floating-point numbers as 32-bit float
#   verbose: if TRUE, show the progress information
# return the normalized output file name invisibly
seqFAVOR2GDS <- function(csv_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="FAVOR", use_float32=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(csv_fn), length(csv_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root)==1L)
    stopifnot(is.logical(use_float32), length(use_float32)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))
    .check_csv_header(csv_fn)
    opt <- .compress_opt(compress)
    # create the gds file
    .conv_begin("FAVOR CSV => GDS", out_fn, compress, verbose)
    outfile <- .gds_new(out_fn, opt$first, var_id_st="double")
    on.exit(closefn.gds(outfile))
    nd_root <- index.gdsn(outfile, "annotation/info")
    if (!is.na(root) && root!="")
        nd_root <- addfolder.gdsn(nd_root, root)
    # load the csv files one by one, to limit the memory usage
    nm_others <- character()
    for (i in seq_along(csv_fn))
    {
        if (verbose) .cat("Reading ", csv_fn[i], " ...")
        df <- read.csv(csv_fn[i])
        if (verbose) .cat("    ", nrow(df), " x ", ncol(df))
        if (i == 1L)
            nm_others <- setdiff(colnames(df), .favor_col_csv)
        .favor_append(outfile, nd_root, df, nm_others, i==1L, use_float32,
            opt$first, verbose)
        remove(df)
    }
    # close the nodes
    .finish_nodes(c(lapply(.gds_var_nodes, index.gdsn, node=outfile),
        lapply(nm_others, index.gdsn, node=nd_root)))
    # recompress?
    on.exit()
    .conv_end(out_fn, compress, opt$recompress, verbose, outfile=outfile)
}



#######################################################################
# Convert FAVOR tar file (containing CSV) to a SeqArray GDS file
#

# the columns of a FAVOR CSV file in a tar file identifying the variants
.favor_col_tar <- c("variant_vcf", "chromosome", "position", "ref_vcf",
    "alt_vcf")

# open a file in the tar.gz file
.tar_open <- function(tar_fn, sub_fn, tar_cmd="tar")
{
    cmd <- paste(tar_cmd, "-xOzf", shQuote(tar_fn), shQuote(sub_fn))
    pipe(cmd, "rt")
}

# the CSV files to be read in the tar file, ordered by the number after "_"
#   in the file names (e.g., "chr22_1.csv", "chr22_2.csv", ...)
.favor_tar_files <- function(tar_fn, fn_in_tar, tar_cmd)
{
    if (!is.null(fn_in_tar))
        return(unique(fn_in_tar))
    fn <- untar(tar_fn, list=TRUE, tar=tar_cmd)
    fn <- fn[grepl("\\.csv$", fn, ignore.case=TRUE)]
    s <- basename(fn)
    m <- regexpr("_\\d+", s)
    n <- rep(NA_integer_, length(s))
    n[m > 0L] <- as.integer(substring(regmatches(s, m), 2L))
    fn[order(n)]
}

# create the annotation nodes with the types and the descriptions in
#   'favor_csv_header.csv'; each node has an index node, since a variant
#   without a value is not stored
.favor_tar_nodes <- function(nd_info, nm, favor_head, use_float32)
{
    i <- match(nm, favor_head$column)
    tp <- favor_head$type[i]
    tp_vcf <- c(character="String", integer="Integer", numeric="Float")[tp]
    tp[tp == "numeric"] <- if (use_float32) "float32" else "float64"
    desp <- favor_head$description[i]
    desp[is.na(desp)] <- ""
    # the data node and its index node, created one after the other
    nd <- lapply(seq_along(nm), function(k)
    {
        n <- add.gdsn(nd_info, nm[k], storage=tp[k], compress="ZIP_ra")
        put.attr.gdsn(n, "Number", ".")
        put.attr.gdsn(n, "Type", tp_vcf[[k]])
        put.attr.gdsn(n, "Description", desp[k])
        list(n, add.gdsn(nd_info, paste0("@", nm[k]), storage="bit1",
            compress="ZIP_ra", visible=FALSE))
    })
    names(nd) <- nm
    list(data=lapply(nd, `[[`, 1L), index=lapply(nd, `[[`, 2L))
}

# append the variants and the annotations of a block of a FAVOR CSV file
.favor_tar_append <- function(outfile, nodes, df, nline)
{
    # basic
    n <- nrow(df)
    val <- list(nline + seq_len(n), df$position, df$chromosome,
        paste0(df$ref_vcf, ",", df$alt_vcf),
        if (is.null(df$rsid)) rep("", n) else df$rsid,
        rep(NaN, n), rep(NA_integer_, n))
    for (i in seq_along(val))
        append.gdsn(index.gdsn(outfile, .gds_var_nodes[i]), val[[i]])
    # annotation: only the values are stored, with an index node
    for (nm in names(nodes$data))
    {
        v <- df[[nm]]
        if (!is.null(v))
        {
            x <- if (is.character(v)) !is.na(v) & v!="" else is.finite(v)
            if (any(x))
                append.gdsn(nodes$data[[nm]], v[x])
        } else {
            x <- rep(0L, n)
        }
        append.gdsn(nodes$index[[nm]], x)
    }
    invisible()
}

# read and check the header line of a FAVOR CSV file
.favor_tar_header <- function(con, favor_head, verbose)
{
    hr <- unlist(strsplit(readLines(con, n=1L), ",", fixed=TRUE))
    if (verbose)
        .cat("  (", length(hr), " columns)")
    if (anyDuplicated(hr))
        stop("Duplicated column names!")
    if (!all(.favor_col_tar %in% hr))
    {
        stop("The csv file should have the columns: ",
            paste(.favor_col_tar, collapse=", "))
    }
    if (!all(hr %in% favor_head$column))
    {
        stop("No type defined for the column(s): ",
            paste(setdiff(hr, favor_head$column), collapse=","))
    }
    hr
}

# read a FAVOR CSV file in the tar file block by block
.favor_tar_read <- function(tar_fn, fn_in_tar, i, tar_cmd, outfile, nd_info,
    favor_head, nodes, nline, use_float32, block_size, verbose)
{
    if (verbose)
    {
        message("Loading file (", i, "):\n    ", fn_in_tar[i],
            appendLF=FALSE)
    }
    con <- .tar_open(tar_fn, fn_in_tar[i], tar_cmd)
    on.exit(close(con))
    hr <- .favor_tar_header(con, favor_head, verbose)
    # first file? create new nodes
    if (is.null(nodes))
    {
        nodes <- .favor_tar_nodes(nd_info, setdiff(hr, .favor_col_tar),
            favor_head, use_float32)
    }
    # column types
    hr_type <- favor_head$type[match(hr, favor_head$column)]
    ii <- 0L
    while (NROW(df <- read.csv(con, header=FALSE, nrows=block_size,
        colClasses=hr_type)))
    {
        names(df) <- hr
        .favor_tar_append(outfile, nodes, df, nline)
        # show progress
        nline <- nline + nrow(df)
        remove(df)
        ii <- ii + 1L
        if (ii >= 100L)
        {
            if (verbose)
            {
                .cat("    ", prettyNum(nline, big.mark=",", scientific=FALSE),
                    "\t", .tm())
            }
            ii <- 0L
            gc(verbose=FALSE, reset=TRUE, full=TRUE)
        }
        # check EOF
        s <- readLines(con, 1L)
        if (length(s)) pushBack(s, con) else break
    }
    if (verbose)
    {
        .cat("    ", prettyNum(nline, big.mark=",", scientific=FALSE),
            "\t", .tm())
    }
    list(nodes=nodes, nline=nline)
}

# finish writing the variants: RLE-coded chromosomes (for faster loading),
#   and remove the index nodes not needed (a value for every variant)
.favor_tar_finish <- function(outfile, verbose)
{
    gc(verbose=FALSE, reset=TRUE, full=TRUE)
    sync.gds(outfile)
    seqOptimize(outfile, "chromosome", verbose=FALSE)
    nvar <- objdesp.gdsn(index.gdsn(outfile, "variant.id"))$dim
    nm <- ls.gdsn(index.gdsn(outfile, "annotation/info"))
    for (s in nm[!grepl("^@", nm)])
    {
        nd <- index.gdsn(outfile, paste0("annotation/info/", s))
        if (identical(objdesp.gdsn(nd)$dim, nvar))
            delete.gdsn(index.gdsn(outfile, paste0("annotation/info/@", s)))
    }
    if (verbose) message("Loading Done.")
    invisible()
}

# create the output file of seqFAVORTar2GDS()
.favor_tar_new <- function(out_fn, root)
{
    outfile <- .gds_new(out_fn, "ZIP_ra")
    on.exit(closefn.gds(outfile))
    add.gdsn(index.gdsn(outfile, "description"), "reference", "GRCh38",
        visible=FALSE)
    if (root != "")
        addfolder.gdsn(index.gdsn(outfile, "annotation/info"), root)
    on.exit()
    outfile
}

# Convert a FAVOR tar file of CSV files to a SeqArray GDS file
#   tar_fn: the file name of the tar file (gzip-compressed) of the FAVOR CSV
#       files
#   out_fn: the output GDS file name
#   fn_in_tar: the CSV file names in the tar file, or NULL for all
#   compress: the compression method of the output file, "LZMA", "ZIP" or
#       "none"
#   root: the folder under "annotation/info" for the annotations, "" for none
#   use_float32: if TRUE, store the floating-point numbers as 32-bit float
#   block_size: the number of lines read from a CSV file at a time
#   tar_cmd: the tar program, "tar" if it is ""
#   verbose: if TRUE, show the progress information
# return the normalized output file name invisibly
seqFAVORTar2GDS <- function(tar_fn, out_fn, fn_in_tar=NULL,
    compress=c("LZMA", "ZIP", "none"), root="", use_float32=TRUE,
    block_size=100000L, tar_cmd=Sys.getenv("TAR"), verbose=TRUE)
{
    # check
    stopifnot(is.character(tar_fn), length(tar_fn)==1L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    if (!is.null(fn_in_tar))
        stopifnot(is.character(fn_in_tar), length(fn_in_tar)>0L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root)==1L)
    stopifnot(is.logical(use_float32), length(use_float32)==1L)
    stopifnot(is.numeric(block_size), length(block_size)==1L, block_size>=1000)
    stopifnot(is.character(tar_cmd), length(tar_cmd)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))
    if (verbose)
        .cat("##< ", .tm(), "\nFAVOR csv tar => GDS\n    input: ", tar_fn,
            "\n    using the program ", shQuote(tar_cmd))
    # list the file in tar
    if (tar_cmd == "") tar_cmd <- "tar"
    fn_in_tar <- .favor_tar_files(tar_fn, fn_in_tar, tar_cmd)
    if (length(fn_in_tar) == 0L)
        stop("No csv files in the input tar file.")
    if (verbose)
    {
        .cat("    files: ", paste(basename(fn_in_tar), collapse=", "))
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }
    # load FAVOR CSV header
    favor_head <- read.csv(system.file("extdata", "favor_csv_header.csv",
            package="GDSAnnotator", mustWork=TRUE), header=TRUE)
    # create the output
    outfile <- .favor_tar_new(out_fn, root)
    on.exit(closefn.gds(outfile))
    nd_info <- index.gdsn(outfile, "annotation/info")
    if (root != "") nd_info <- index.gdsn(nd_info, root)
    # read the files one by one (the nodes are created with the first file)
    r <- list(nodes=NULL, nline=0L)
    for (i in seq_along(fn_in_tar))
    {
        r <- .favor_tar_read(tar_fn, fn_in_tar, i, tar_cmd, outfile, nd_info,
            favor_head, r$nodes, r$nline, use_float32, block_size, verbose)
    }
    .favor_tar_finish(outfile, verbose)
    # recompress?
    on.exit()
    .conv_end(out_fn, compress, compress=="LZMA", verbose, outfile=outfile)
}
