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
        .cat("##< ", tm())
        .cat("FAVOR CSV => GDS")
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }
    outfile <- .gds_new(out_fn, compress1, var_id_st="double")
    on.exit(closefn.gds(outfile))
    nm_lst <- c("position", "chromosome", "ref_vcf", "alt_vcf", "variant_vcf")
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
        if ("vid" %in% colnames(df))
        {
            append.gdsn(index.gdsn(outfile, "variant.id"), df$vid)
        } else {
            append.gdsn(index.gdsn(outfile, "variant.id"), seq_len(nrow(df)))
        }
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
            .cat("Recompressing (", tm(), ") ...")
        closefn.gds(outfile)
        cleanup.gds(out_fn, verbose=FALSE)
        seqRecompress(out_fn, compress=compress, verbose=verbose)
    } else {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        closefn.gds(outfile)
        cleanup.gds(out_fn, verbose=verbose)
    }

    if (verbose) .cat("##> ", tm())
    # output
    invisible(normalizePath(out_fn))
}



#######################################################################
# Convert FAVOR tar file (containing CSV) to a SeqArray GDS file
#

# open a file in the tar.gz file
tar_open <- function(tar_fn, sub_fn, tar_cmd="tar")
{
	cmd <- paste(tar_cmd, "-xOzf", tar_fn, sub_fn)
	pipe(cmd, "rt")
}


seqToGDS_FAVOR_tar <- function(tar_fn, out_fn, fn_in_tar=NULL,
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

    # create the gds file
    if (verbose)
    {
        .cat("##< ", tm())
        .cat("FAVOR csv tar => GDS")
        .cat("    input: ", tar_fn)
        .cat("    using the program ", shQuote(tar_cmd))
    }

    # list the file in tar
    if (tar_cmd == "") tar_cmd <- "tar"
    if (is.null(fn_in_tar))
    {
        fn_in_tar <- untar(tar_fn, list=TRUE, tar=tar_cmd)
        fn_in_tar <- fn_in_tar[grepl("\\.csv$", fn_in_tar, ignore.case=TRUE)]
        if (length(fn_in_tar))
        {
            n <- vapply(basename(fn_in_tar), function(s)
            {
                i <- regexpr("_\\d+", s)
                if (i > 0L)
                    as.integer(gsub("^_", "", regmatches(s, i)))
                else
                    NA_integer_
            }, 0L)
            fn_in_tar <- fn_in_tar[order(n)]
        }
	} else {
	    fn_in_tar <- unique(fn_in_tar)
	}
	if (length(fn_in_tar) == 0L)
	    stop("No csv files in the input tar file.")
    if (verbose)
    {
        .cat("        ", paste(basename(fn_in_tar), collapse=", "))
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }

    # load FAVOR CSV header
    favor_head <- read.csv(system.file("extdata", "favor_csv_header.csv",
        package="GDSAnnotator", mustWork=TRUE))

    # create the output
    outfile <- createfn.gds(out_fn)
    put.attr.gdsn(outfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(outfile$root, "FileVersion", "v1.0")
    addfolder.gdsn(outfile, "description")
    # add sample.id
    add.gdsn(outfile, "sample.id", integer())
    # add variant.id
    nd_varid <- add.gdsn(outfile, "variant.id", storage="int32", compress="ZIP_ra")
    # add position
    nd_pos <- add.gdsn(outfile, "position", storage="int32", compress="ZIP_ra")
    # add chromosome
    nd_chr <- add.gdsn(outfile, "chromosome", storage="string", compress="ZIP_ra")
    # add allele
    nd_allele <- add.gdsn(outfile, "allele", storage="string", compress="ZIP_ra")
    # add a folder for genotypes
    nd_geno <- addfolder.gdsn(outfile, "genotype")
    put.attr.gdsn(nd_geno, "VariableName", "GT")
    put.attr.gdsn(nd_geno, "Description", "Genotype")
    # add phase folder
    nd_phase <- addfolder.gdsn(outfile, "phase")
    # add annotation folder
    nd_annot <- addfolder.gdsn(outfile, "annotation")
    # add rsID
    nd_id <- add.gdsn(nd_annot, "id", storage="string", compress="ZIP_ra")
    # add qual
    nd_qual <- add.gdsn(nd_annot, "qual", storage="float", compress="ZIP_ra")
    # add filter
    nd_flt <- add.gdsn(nd_annot, "filter", storage="int32", compress="ZIP_ra")
    put.attr.gdsn(nd_flt, "R.class", "factor")
    put.attr.gdsn(nd_flt, "R.levels", "PASS")
    put.attr.gdsn(nd_flt, "Description", "All filters passed")
    # VCF INFO
    nd_info <- addfolder.gdsn(nd_annot, "info")
    if (root != "")
        nd_info <- addfolder.gdsn(nd_info, root)
    nd_lst <- nd_idx_lst <- list()

    nm_head <- c("variant_vcf", "chromosome", "position", "ref_vcf", "alt_vcf")
    nline <- 0L

    # for-loop for each file
    for (i_fn in seq_along(fn_in_tar))
    {
        fn <- fn_in_tar[i_fn]
        if (verbose)
            cat("Loading file (", i_fn, "):\n    ", fn, sep="")
        InF <- tar_open(tar_fn, fn, tar_cmd)
        # check the header
        hr <- readLines(InF, n=1L)
        hr <- unlist(strsplit(hr, ",", fixed=TRUE))
        if (verbose)
            .cat("  (", length(hr), " columns)")
        if (anyDuplicated(hr))
            stop("Duplicated column names!")
        if (!all(nm_head %in% hr))
        {
            stop("The csv file should have the columns: ",
                paste(nm_head, collapse=", "))
        }
        if (!all(hr %in% favor_head$column))
        {
            stop("No type defined for the column(s): ",
                paste(setdiff(hr, favor_head$column), collapse=","))
        }
        # first time? create new nodes
        if (i_fn == 1L)
        {
            nm_others <- setdiff(hr, nm_head)
            for (nm in nm_others)
            {
                i <- which(nm == favor_head$column)
                tp <- favor_head$type[i]
                if (tp == "numeric")
                    tp <- ifelse(use_float32, "float32", "float64")
                nd <- add.gdsn(nd_info, nm, storage=tp, compress="ZIP_ra")
                put.attr.gdsn(nd, "Number", ".")
                s <- favor_head$description[i]
                if (is.na(s)) s <- ""
                put.attr.gdsn(nd, "Description", s)
                nd_lst[[nm]] <- nd
                nd_idx_lst[[nm]] <- add.gdsn(nd_info, paste0("@", nm),
                    storage="bit1", compress="ZIP_ra", visible=FALSE)
            }
        }
        # column types
        hr_type <- favor_head$type[match(hr, favor_head$column)]
        ii <- 0L
        # read
        while(length(txt <- readLines(InF, n=block_size)))
        {
            df <- read.csv(text=txt, header=FALSE, nrows=block_size,
                colClasses=hr_type)
            names(df) <- hr
            # basic
            append.gdsn(nd_varid, nline + seq_len(nrow(df)))
            append.gdsn(nd_pos, df$position)
            append.gdsn(nd_chr, df$chromosome)
            append.gdsn(nd_allele, paste0(df$ref_vcf, ",", df$alt_vcf))
            if (!is.null(df$rsid))
                append.gdsn(nd_id, df$rsid)
            else
                append.gdsn(nd_id, rep("", nrow(df)))
            append.gdsn(nd_qual, rep(NaN, nrow(df)))
            append.gdsn(nd_flt, rep(NA_integer_, nrow(df)))
            # annotation
            for (nm in nm_others)
            {
                v <- df[[nm]]
                if (!is.null(v))
                {
                    if (is.character(v))
                        x <- !is.na(v) & v!=""
                    else
                        x <- is.finite(v)
                    # if (any(x))
                        append.gdsn(nd_lst[[nm]], v[x])
                } else {
                    x <- rep(0L, nrow(df))
                }
                append.gdsn(nd_idx_lst[[nm]], x)
            }
            # show progress
            nline <- nline + nrow(df)
            ii <- ii + 1L
            if (ii>=100L && verbose)
            {
                .cat("    ", prettyNum(nline, big.mark=",", scientific=FALSE),
                    "\t", tm())
                ii <- 0L
            }
        }
        # close the pipe
        close(InF)
        if (verbose)
        {
            .cat("    ", prettyNum(nline, big.mark=",", scientific=FALSE),
                "\t", tm())
        }
    }

    # get the number of variants
    sync.gds(outfile)
    nvar <- objdesp.gdsn(index.gdsn(outfile, "variant.id"))$dim
    # RLE-coded chromosome (for faster loading)
    SeqArray:::.optim_chrom(outfile)

    # whether has missing values or not
    nm_lst <- ls.gdsn(index.gdsn(outfile, "annotation/info"))
    nm_lst <- nm_lst[!grepl("^@", nm_lst)]
    for (nm in nm_lst)
    {
        nd <- index.gdsn(outfile, paste0("annotation/info/", nm))
        if (identical(objdesp.gdsn(nd)$dim, nvar))
        {
            # delete the indexing node
            delete.gdsn(index.gdsn(outfile, paste0("annotation/info/@", nm)))
        }
    }
    if (verbose) cat("Loading Done.\n")

    # recompress?
    on.exit()
    if (compress == "LZMA")
    {
        if (verbose)
            .cat("Recompressing (", tm(), ") ...")
        closefn.gds(outfile)
        cleanup.gds(out_fn, verbose=FALSE)
        seqRecompress(out_fn, compress=compress, verbose=verbose)
    } else {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        closefn.gds(outfile)
        cleanup.gds(out_fn, verbose=verbose)
    }

    if (verbose) .cat("##> ", tm())
    # output
    invisible(normalizePath(out_fn))
}
