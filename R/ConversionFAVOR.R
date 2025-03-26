#######################################################################
#
# Package Name: GDSAnnotator
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Data import & Format conversion from FAVOR (https://favor.genohub.org)
#


#######################################################################
# Convert FAVOR CSV files to a SeqArray GDS file
#

seqToGDS_FAVOR <- function(csv_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="FAVOR", use_float32=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(csv_fn), length(csv_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root)==1L)
    stopifnot(is.logical(use_float32), length(use_float32)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))
    if (length(csv_fn) > 1L)
    {
        # check the header of each file, should be the same!
        s0 <- readLines(csv_fn[1L], n=1L)
        if (length(s0) != 1L)
            stop("No header line in the first file of 'csv_fn'.")
        for (i in 2L:length(csv_fn))
        {
            s <- readLines(csv_fn[i], n=1L)
            if (length(s)!=1L || s!=s0)
                stop(sprintf("Error in the header of 'csv_fn[%d]'", i))
        }
    }

    # compression algorithm
    map_compress <- c(LZMA="LZMA_ra", ZIP="ZIP_RA", none="")
    compress1 <- compress2 <- map_compress[compress]
    if (compress=="LZMA") compress1 <- "ZIP_RA"  # reduce memory usage

    # create the gds file
    if (verbose)
    {
        .cat("FAVOR => GDS (", date(), "):")
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }
    outfile <- .gds_new(out_fn, compress1, var_id_st="double")
    on.exit(closefn.gds(outfile))
    nm_lst <- c("vid", "position", "chromosome", "ref_vcf", "alt_vcf", "variant_vcf")
    nd_root <- index.gdsn(outfile, "annotation/info")
    if (!is.na(root) && root!="")
        nd_root <- addfolder.gdsn(nd_root, root)

    # load csv
    for (i in seq_along(csv_fn))
    {
        .cat("Reading ", csv_fn[i], " ...")
        df <- read_csv(csv_fn[i], progress=verbose, show_col_types=FALSE)
        .cat("    ", nrow(df), " x ", ncol(df))
        if (!all(nm_lst %in% colnames(df)))
        {
            stop("The csv header should contain all of ",
                paste(nm_lst, collapse=", "), ".")
        }
        # basic
        append.gdsn(index.gdsn(outfile, "variant.id"), df$vid)
        append.gdsn(index.gdsn(outfile, "position"), df$position)
        append.gdsn(index.gdsn(outfile, "chromosome"), df$chromosome)
        append.gdsn(index.gdsn(outfile, "allele"), paste0(df$ref_vcf, ",", df$alt_vcf))
        append.gdsn(index.gdsn(outfile, "annotation/id"), df$variant_vcf)
        append.gdsn(index.gdsn(outfile, "annotation/qual"), rep(NaN, nrow(df)))
        append.gdsn(index.gdsn(outfile, "annotation/filter"), rep(1L, nrow(df)))
        # annotation
        for (nm in setdiff(colnames(df), nm_lst))
        {
            if (verbose) .cat("    [adding ", nm, "]\t")
            v <- df[[nm]]
            if (i == 1L)
            {
                st <- typeof(v)
                if (st=="double" && isTRUE(use_float32)) st <- "float32"
                nd <- suppressWarnings(
                    add.gdsn(nd_root, nm, v, storage=st, compress=compress1))
            } else {
                nd <- suppressWarnings(append.gdsn(index.gdsn(nd_root, nm), v))
            }
            if (verbose) { cat("    "); print(nd) }
        }
        remove(df)
    }

    # close the nodes
    for (nm in c("variant.id", "position", "chromosome", "allele",
        "annotation/id", "annotation/qual", "annotation/filter"))
    {
        nd <- index.gdsn(outfile, nm)
        readmode.gdsn(nd)
        SeqArray:::.DigestCode(nd, verbose=FALSE)
    }
    for (nm in setdiff(colnames(df), nm_lst))
    {
        nd <- index.gdsn(nd_root, nm)
        readmode.gdsn(nd)
        SeqArray:::.DigestCode(nd, verbose=FALSE)
    }

    # recompress?
    on.exit()
    if (compress1 != compress2)
    {
        if (verbose)
            .cat("Recompressing (", date(), ") ...")
        closefn.gds(outfile)
        cleanup.gds(out_fn, verbose=FALSE)
        seqRecompress(out_fn, compress=compress, verbose=verbose)
    } else {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        closefn.gds(outfile)
        cleanup.gds(out_fn, verbose=verbose)
    }

    if (verbose) .cat("Done (", date(), ")")
    # output
    invisible(normalizePath(out_fn))
}

