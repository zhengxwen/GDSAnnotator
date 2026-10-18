#######################################################################
#
# Package Name: GDSAnnotator
# Copyright (C) 2025-2026    Xiuwen Zheng
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Data import & Format conversion
#

# Package-wide variables
.packageEnv <- new.env(parent=emptyenv())

# Internal functions: show a line of progress information, as a message
#   (i.e., on stderr, and it can be suppressed by suppressMessages())
.cat <- function(...) message(...)

# the current date and time, shown in the progress information
.tm <- function() strftime(Sys.time(), "%Y-%m-%d %H:%M:%S")

# replace NA by "" in a character vector, since NA is not allowed in a GDS
#   string variable
.na_to_empty <- function(x)
{
    if (is.character(x)) x[is.na(x)] <- ""
    x
}


#######################################################################
# Internal: the steps shared by the conversion functions
#

# the compression methods of a conversion: 'first' when writing the data
#   (ZIP_RA for LZMA, which needs less memory), and whether to recompress the
#   output file afterwards
.compress_opt <- function(compress)
{
    final <- c(LZMA="LZMA_RA", ZIP="ZIP_RA", none="")[[compress]]
    first <- if (compress == "LZMA") "ZIP_RA" else final
    list(first=first, recompress=(first != final))
}

# show the beginning of a conversion
.conv_begin <- function(title, out_fn, compress, verbose)
{
    if (verbose)
    {
        .cat("##< ", .tm())
        .cat(title)
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
    }
    invisible()
}

# finish a conversion: close 'outfile' if it is open, recompress the output
#   file or clean up its fragments, and return the file name
.conv_end <- function(out_fn, compress, recompress, verbose, outfile=NULL,
    cleanup=TRUE)
{
    if (recompress)
    {
        if (verbose)
            .cat("Recompressing (", .tm(), ") ...")
        if (!is.null(outfile))
        {
            closefn.gds(outfile)
            cleanup.gds(out_fn, verbose=FALSE)
        }
        seqRecompress(out_fn, compress=compress, verbose=verbose)
    } else {
        if (!is.null(outfile))
        {
            if (verbose)
                message("Optimize the access efficiency ...")
            closefn.gds(outfile)
        }
        if (cleanup)
            cleanup.gds(out_fn, verbose=verbose)
    }
    if (verbose) .cat("##> ", .tm())
    # output
    invisible(normalizePath(out_fn))
}

# load an internal table of the output format of an annotator (cached)
.format_table <- function(key, fn, cols)
{
    if (is.null(.packageEnv[[key]]))
    {
        tab <- read.csv(system.file("extdata", fn, package="GDSAnnotator",
            mustWork=TRUE))
        if (!all(cols %in% colnames(tab)))
        {
            stop("The internal '", fn, "' should have the following columns: ",
                paste(cols, collapse=","), ".")
        }
        .packageEnv[[key]] <- tab
    }
    .packageEnv[[key]]
}

# the descriptions and the 'uniform' flags of the sub-fields
.field_info <- function(tab, nm)
{
    j <- match(nm, tab$Name)
    desp <- tab$Description[j]
    desp[is.na(j)] <- ""
    uniform <- if (is.null(tab$Uniform)) logical(length(nm)) else
        tab$Uniform[j] %in% TRUE
    list(desp=desp, uniform=uniform)
}

# create a new SeqArray GDS file
.gds_new <- function(out_fn, compress, var_id_st="int32")
{
    outfile <- createfn.gds(out_fn)
    on.exit(closefn.gds(outfile))
    put.attr.gdsn(outfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(outfile$root, "FileVersion", "v1.0")
    addfolder.gdsn(outfile, "description")
    # add sample.id
    add.gdsn(outfile, "sample.id", integer())
    # add variant.id
    add.gdsn(outfile, "variant.id", storage=var_id_st, compress=compress)
    # add position
    add.gdsn(outfile, "position", storage="int32", compress=compress)
    # add chromosome
    add.gdsn(outfile, "chromosome", storage="string", compress=compress)
    # add allele
    add.gdsn(outfile, "allele", storage="string", compress=compress)
    # add a folder for genotypes
    nd <- addfolder.gdsn(outfile, "genotype")
    put.attr.gdsn(nd, "VariableName", "GT")
    put.attr.gdsn(nd, "Description", "Genotype")
    # add phase folder
    addfolder.gdsn(outfile, "phase")
    # add annotation folder
    nd_annot <- addfolder.gdsn(outfile, "annotation")
    # add id
    add.gdsn(nd_annot, "id", storage="string", compress=compress)
    # add qual
    add.gdsn(nd_annot, "qual", storage="float32", compress=compress)
    # add filter
    nd <- add.gdsn(nd_annot, "filter", storage="int32", compress=compress)
    put.attr.gdsn(nd, "R.class", "factor")
    put.attr.gdsn(nd, "R.levels", "PASS")
    put.attr.gdsn(nd, "Description", "All filters passed")
    # VCF INFO
    addfolder.gdsn(nd_annot, "info")
    # output
    on.exit()
    outfile
}

# delete the INFO nodes and their index nodes if they exist
.delete_info <- function(nd_info, nm, verbose)
{
    all_nm <- ls.gdsn(nd_info)
    for (s in intersect(nm, all_nm))
    {
        if (paste0("@", s) %in% all_nm)
            delete.gdsn(index.gdsn(nd_info, paste0("@", s)))
        delete.gdsn(index.gdsn(nd_info, s))
        if (verbose)
            .cat("    removed '", s, "'")
    }
    invisible()
}


#######################################################################
# Convert gnomAD VCF to a SeqArray GDS file
#   vcf_fn: the file name(s) of the gnomAD VCF files
#   out_fn: the output GDS file name
#   compress: the compression method of the output file, "LZMA", "ZIP" or
#       "none"
#   verbose: if TRUE, show the progress information
# return the normalized output file name invisibly
#

seqGnomAD2GDS <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))
    opt <- .compress_opt(compress)
    .conv_begin("gnomAD VCF => GDS", out_fn, compress, verbose)
    # import from VCF
    attr(verbose, "header_no_time") <- TRUE
    seqVCF2GDS(vcf_fn, out_fn, storage.option=opt$first, verbose=verbose)
    # recompress?
    .conv_end(out_fn, compress, opt$recompress, verbose, cleanup=FALSE)
}



#######################################################################
# Internal: split a pipe-delimited annotation field into sub-fields
# using block-based processing to reduce memory usage
#

# Target number of raw annotation bytes processed at any one time.
# gregexpr() + substring() avoids the full strsplit() expansion but still
# allocates an integer boundary matrix; chunking by byte size caps the
# working set independently of 'bsize'.
.annot_chunk_bytes <- 8L * 1048576L  # 8 MB

# split the annotation field 'nm_root' into the sub-fields stored in the
#   folder 'nm_root2'
.split_annot_blocks <- function(f, nm_root, nm_root2, nm_lst, nm_desp,
    nm_uniform, compress, bsize, type_fn=NULL, paren=FALSE, verbose=TRUE)
{
    # create output folder
    seqAddValue(f, nm_root2, NULL, verbose=verbose)
    if (verbose)
        .cat("    ", paste(nm_lst, collapse=","))
    n_fields <- length(nm_lst)
    if (is.null(nm_desp)) nm_desp <- rep("", n_fields)
    if (is.null(nm_uniform)) nm_uniform <- rep(FALSE, n_fields)
    # pre-create empty GDS nodes for each sub-field
    nd_folder <- index.gdsn(f, nm_root2)
    nodes <- .split_new_nodes(nd_folder, nm_lst, nm_desp, type_fn, compress)
    ctx <- list(nodes=nodes, nm_lst=nm_lst, nm_uniform=nm_uniform,
        type_fn=type_fn, paren=paren)
    # block-by-block processing
    if (verbose) message("Processing:")
    # .padNA=FALSE: without it, a field whose variants never carry more than
    #   one entry (e.g. SnpEff NMD) comes back padded to one value per variant
    #   with NA for the absent ones, so the sub-fields would get more entries
    #   than nm_root's index node accounts for.
    seqBlockApply(f, nm_root, function(bk) .split_block(bk, ctx),
        as.is="none", bsize=bsize, .padNA=FALSE, .progress=verbose)
    # finalize all data nodes
    for (nd in nodes) readmode.gdsn(nd)
    # The per-variant index for non-uniform fields is identical to nm_root's
    # index node; copy it directly instead of rebuilding during the loop.
    nd_idx <- index.gdsn(f, paste0(dirname(nm_root), "/@", basename(nm_root)),
        silent=TRUE)
    if (!is.null(nd_idx))
    {
        for (nm in nm_lst[!nm_uniform])
            copyto.gdsn(nd_folder, nd_idx, paste0("@", nm))
    }
    # return
    invisible()
}

# create the empty GDS nodes of the sub-fields
.split_new_nodes <- function(nd_folder, nm_lst, nm_desp, type_fn, compress)
{
    lapply(seq_along(nm_lst), function(i)
    {
        v <- character()
        if (!is.null(type_fn)) v <- type_fn(nm_lst[i], v)
        tp <- storage.mode(v)
        if (tp=="double") tp <- "float32"
        nd <- add.gdsn(nd_folder, nm_lst[i], storage=tp, valdim=0L,
            compress=compress)
        if (nzchar(nm_desp[i]))
            put.attr.gdsn(nd, "Description", nm_desp[i])
        nd
    })
}

# the sub-chunks of a block, cut at variant boundaries
.split_chunks <- function(dat, ns)
{
    # Cut the block into sub-chunks of ~.annot_chunk_bytes of raw text, so
    # that only a small part of it is expanded at a time. Cuts are made at
    # variant boundaries to keep each variant's entries (and hence the index
    # and 'uniform' selection) intact. Sizing by bytes rather than by entry
    # count adapts to gene-dense regions where a single variant carries many
    # long transcript annotations.
    cs <- cumsum(ns)  # index of each variant's last entry
    # as.double() is required here: this accumulates *bytes*, and an
    # integer cumsum() silently yields NA past 2^31 (2.1 GB). The NAs would
    # drop every chunk boundary after that point and collapse the rest of
    # the block into a single huge sub-chunk -- i.e. it would fail exactly
    # on the large blocks this chunking exists for.
    # The leading 0 and [cs + 1L] keep the result aligned with 'ns' when a
    # variant has no entries at all (cs == 0), as is common for SnpEff
    # LOF/NMD; a plain [cs] drops those silently and shifts every
    # subsequent chunk boundary.
    cw <- c(0, cumsum(as.double(nchar(dat, type="bytes",
        keepNA=FALSE))))[cs + 1L]
    g <- floor(cw / .annot_chunk_bytes)
    st <- c(1L, which(g[-1L] != g[-length(g)]) + 1L)  # first variant
    # A variant with no entries must stay in the same sub-chunk as the
    # variant that follows it: for a 'uniform' field the original code
    # resolves it to the *next* variant's first entry, so a cut may not
    # fall on a run of empty variants. Boundaries land on non-empty
    # variants, so walk each one back over the preceding empty run.
    if (length(st) > 1L)
    {
        d <- which(ns > 0L)
        j <- match(st[-1L], d) - 1L
        st <- unique(c(1L, ifelse(j >= 1L, d[pmax(j, 1L)] + 1L, 1L)))
    }
    en <- c(st[-1L] - 1L, length(ns))  # last variant
    list(st=st, en=en, cs=cs)
}

# the boundaries of the sub-fields in each entry, using the positions of the
#   pipes instead of strsplit(), which avoids a full copy of all the field
#   strings at once
.split_bounds <- function(sub, n_fields, paren)
{
    # no entries in this sub-chunk at all
    if (!length(sub)) return(matrix(0L, 0L, n_fields + 1L))
    n_expect <- n_fields - 1L  # pipes per well-formed entry
    pipe_pos <- gregexpr("|", sub, fixed=TRUE)
    nc <- nchar(sub)
    np <- lengths(pipe_pos)
    # gregexpr() reports a single -1 when there is no match at all;
    # which() keeps this safe if an entry is NA
    np[which(np == 1L & vapply(pipe_pos, `[`, 0L, 1L) < 0L)] <- 0L
    if (paren)
    {
        # SnpEff wraps each LOF/NMD entry in parentheses, which are not part
        # of the first and last field. Skip them by moving the outer
        # boundaries instead of rewriting the strings, and test each entry
        # so data without them is left alone.
        lo <- as.integer(startsWith(sub, "("))
        hi <- nc - as.integer(endsWith(sub, ")"))
    } else {
        lo <- 0L; hi <- nc
    }
    # n_expect == 0 must take the slow path: with no pipes at all gregexpr()
    # yields -1, which would become a bogus column
    if (n_expect > 0L && all(np == n_expect))
        return(cbind(lo, do.call(rbind, pipe_pos), hi + 1L))
    # pad short entries; for an over-long entry keep the pipe that ends the
    # last field, so it is cut where strsplit() would cut it rather than
    # running to the end of the string
    pipe_pos <- lapply(seq_along(pipe_pos), function(j) {
        p <- pipe_pos[[j]]
        if (np[j] >= n_fields) p[seq_len(n_fields)]
        else c(p[seq_len(np[j])], rep.int(hi[j] + 1L, n_fields - np[j]))
    })
    cbind(lo, do.call(rbind, pipe_pos))
}

# split the entries of a block returned by seqBlockApply() and append the
#   sub-fields to their nodes
.split_block <- function(bk, ctx)
{
    if (inherits(bk, "SeqVarDataList"))
    {
        dat <- bk$data
        ns <- bk$length
    } else {
        dat <- bk
        ns <- rep(1L, length(bk))
    }
    ch <- .split_chunks(dat, ns)
    # the sub-chunks are processed one after another to limit the memory usage
    for (k in seq_along(ch$st))
    {
        if (length(ch$st) > 1L)
        {
            # entries of the variants st[k]:en[k]; a run of variants with no
            # entries at all gives an empty sub-chunk
            i1 <- ch$cs[ch$st[k]] - ns[ch$st[k]] + 1L
            i2 <- ch$cs[ch$en[k]]
            sub <- if (i2 >= i1) dat[seq.int(i1, i2)] else dat[0L]
            nsub <- ns[seq.int(ch$st[k], ch$en[k])]
        } else {
            sub <- dat; nsub <- ns  # single chunk: avoid a copy
        }
        .split_write(ctx, sub, nsub,
            .split_bounds(sub, length(ctx$nm_lst), ctx$paren))
    }
    invisible()
}

# write each sub-field of the entries to the corresponding data node
.split_write <- function(ctx, sub, nsub, bnd)
{
    # the first entry of each variant, for the 'uniform' sub-fields
    first <- cumsum(c(1L, head(nsub, -1L)))
    for (i in seq_along(ctx$nodes))
    {
        # A field missing from an entry (a trailing '|', or fewer fields
        # than the header declares) yields "". NOTE: this is a deliberate
        # departure from strsplit(), which would give NA. An entry that is
        # itself NA is stored as "" as well.
        v <- substring(sub, bnd[, i] + 1L, bnd[, i + 1L] - 1L)
        if (!is.null(ctx$type_fn)) v <- ctx$type_fn(ctx$nm_lst[i], v)
        if (ctx$nm_uniform[i]) v <- v[first]
        append.gdsn(ctx$nodes[[i]], .na_to_empty(v))
    }
    invisible()
}



#######################################################################
# Convert Ensembl-VEP to a SeqArray GDS file
#

# VEP type conversion function
.vep_type_fn <- function(nm, v)
{
    # check
    stopifnot(is.character(nm), length(nm)==1L)
    stopifnot(!is.null(.packageEnv$vep))
    # variable names and types
    j <- match(nm, .packageEnv$vep$Name)
    if (!is.na(j))
    {
        lg_set <- c("1", "T", "TRUE", "YES")
        tp <- .packageEnv$vep$Type[j]
        if (is.atomic(v))
        {
            v <- switch(tp,
                integer = as.integer(v),
                numeric = as.numeric(v),
                logical = is.element(v, set=lg_set),
                v)
        } else {
            v <- switch(tp,
                integer = lapply(v, as.integer),
                numeric = lapply(v, as.numeric),
                logical = lapply(v, is.element, set=lg_set),
                v)
        }
    }
    v
}

# import a VEP-annotated VCF file and split the CSQ field
.vep_vcf <- function(vcf_fn, out_fn, opt, root, keep, bsize, verbose)
{
    # vcf => gds
    attr(verbose, "header_no_time") <- TRUE
    seqVCF2GDS(vcf_fn, out_fn, storage.option=opt$first, optimize=FALSE,
        verbose=verbose)
    # split CSQ (Consequence annotations from Ensembl VEP)
    f <- seqOpen(out_fn, readonly=FALSE)
    on.exit(seqClose(f))
    # need CSQ
    nm_root <- paste0("annotation/info/", root)
    nm_root2 <- paste0("annotation/info/", root, ".list")
    nd <- index.gdsn(f, nm_root)
    desp <- get.attr.gdsn(nd)$Description
    if (!is.character(desp))
        stop(root, " information is not found!")
    desp <- gsub("^.*Format:", "", desp)
    # sub fields in CSQ, with their descriptions and uniform flags
    nm_lst <- trimws(unlist(strsplit(desp, "|", fixed=TRUE)))
    info <- .field_info(.packageEnv$vep, nm_lst)
    # split annotation into sub-fields using block processing
    .split_annot_blocks(f, nm_root, nm_root2, nm_lst, info$desp, info$uniform,
        opt$first, bsize, type_fn=.vep_type_fn, verbose=verbose)
    # remove the original root node if keep is not TRUE
    if (!isTRUE(keep))
    {
        delete.gdsn(index.gdsn(f, nm_root))
        if (verbose)
            .cat("    removed '", nm_root, "'")
    }
    invisible()
}

# Convert a VEP-annotated VCF file to a SeqArray GDS file
#   vcf_fn: the file name(s) of the VCF files annotated by Ensembl VEP
#   out_fn: the output GDS file name
#   compress: the compression method of the output file, "LZMA", "ZIP" or
#       "none"
#   root: the INFO field of the VEP annotation, split into the sub-fields
#       stored in the folder "annotation/info/<root>.list"
#   keep: if TRUE, keep the original field 'root'
#   bsize: the number of variants in a block when splitting the field
#   verbose: if TRUE, show the progress information
# return the normalized output file name invisibly
seqVEP2GDS <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="CSQ", keep=FALSE, bsize=50000L, verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root)==1L)
    stopifnot(is.logical(keep), length(keep)==1L)
    stopifnot(is.numeric(bsize), length(bsize)==1L, bsize>=1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))
    opt <- .compress_opt(compress)
    .conv_begin("Ensembl-VEP VCF => GDS", out_fn, compress, verbose)
    # load variable names and types
    .format_table("vep", "vep_output_format.csv",
        c("Name", "Description", "Type"))
    # import from VCF
    .vep_vcf(vcf_fn, out_fn, opt, root, keep, as.integer(bsize), verbose)
    # recompress?
    .conv_end(out_fn, compress, opt$recompress, verbose)
}



#######################################################################
# Convert SnpEff to a SeqArray GDS file
#

# SnpEff type conversion function
.snpeff_type_fn <- function(nm, v)
{
    # check
    stopifnot(is.character(nm), length(nm)==1L)
    stopifnot(!is.null(.packageEnv$snpeff_sub))
    # variable names and types
    j <- match(nm, .packageEnv$snpeff_sub$Name)
    if (!is.na(j))
    {
        tp <- .packageEnv$snpeff_sub$Type[j]
        if (is.atomic(v))
        {
            v <- switch(tp,
                integer = as.integer(v),
                numeric = as.numeric(v),
                v)
        } else {
            v <- switch(tp,
                integer = lapply(v, as.integer),
                numeric = lapply(v, as.numeric),
                v)
        }
    }
    v
}

# split a SnpEff field (ANN, LOF or NMD)
.snpeff_split <- function(f, root, opt, keep, bsize, verbose)
{
    nm_root <- paste0("annotation/info/", root)
    nm_root2 <- paste0("annotation/info/", root, ".list")
    nd <- index.gdsn(f, nm_root)
    desp <- get.attr.gdsn(nd)$Description
    if (!is.character(desp))
        stop(root, " information is not found!")
    # columns
    s <- regmatches(desp, regexpr("'.*'", desp))
    if (!length(s))
        stop("No column names!")
    s <- gsub("'|\\s", "", s)
    s <- gsub("/", "-", s, fixed=TRUE)
    nm_lst <- trimws(unlist(strsplit(s, "|", fixed=TRUE)))
    # look up descriptions, types and uniform flags from CSV
    sub <- .packageEnv$snpeff[.packageEnv$snpeff$Field == root, , drop=FALSE]
    .packageEnv$snpeff_sub <- sub
    info <- .field_info(sub, nm_lst)
    # determine type_fn: use it if any field has a Type defined
    type_fn <- if (any(nzchar(sub$Type))) .snpeff_type_fn else NULL
    # split annotation into sub-fields using block processing
    .split_annot_blocks(f, nm_root, nm_root2, nm_lst, info$desp, info$uniform,
        opt$first, bsize, type_fn=type_fn,
        paren=(root %in% c("LOF", "NMD")), verbose=verbose)
    # remove the original root node if keep is not TRUE
    if (!isTRUE(keep))
    {
        delete.gdsn(index.gdsn(f, nm_root))
        if (verbose)
            .cat("    removed '", nm_root, "'")
    }
    invisible()
}

# import a SnpEff-annotated VCF file and split the fields
.snpeff_vcf <- function(vcf_fn, out_fn, opt, root_lst, keep, bsize, verbose)
{
    # vcf => gds
    attr(verbose, "header_no_time") <- TRUE
    seqVCF2GDS(vcf_fn, out_fn, storage.option=opt$first, optimize=FALSE,
        verbose=verbose)
    # split the fields (each modifies the GDS file)
    f <- seqOpen(out_fn, readonly=FALSE)
    on.exit(seqClose(f))
    for (root in root_lst)
        .snpeff_split(f, root, opt, keep, bsize, verbose)
    invisible()
}

# Convert a SnpEff-annotated VCF file to a SeqArray GDS file
#   vcf_fn: the file name(s) of the VCF files annotated by SnpEff
#   out_fn: the output GDS file name
#   compress: the compression method of the output file, "LZMA", "ZIP" or
#       "none"
#   root: the INFO fields of the SnpEff annotation, each split into the
#       sub-fields stored in the folder "annotation/info/<root>.list"
#   keep: if TRUE, keep the original fields 'root'
#   bsize: the number of variants in a block when splitting the fields
#   verbose: if TRUE, show the progress information
# return the normalized output file name invisibly
seqSnpEff2GDS <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root=c("ANN", "LOF", "NMD"), keep=FALSE, bsize=50000L, verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root) > 0)
    stopifnot(is.logical(keep), length(keep)==1L)
    stopifnot(is.numeric(bsize), length(bsize)==1L, bsize>=1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))
    opt <- .compress_opt(compress)
    # load variable names and types
    .format_table("snpeff", "snpeff_output_format.csv",
        c("Field", "Name", "Description", "Type"))
    .conv_begin("SnpEff VCF => GDS", out_fn, compress, verbose)
    # import from VCF
    .snpeff_vcf(vcf_fn, out_fn, opt, root, keep, as.integer(bsize), verbose)
    # recompress?
    .conv_end(out_fn, compress, opt$recompress, verbose)
}



#######################################################################
# Convert ANNOVAR-annotated VCF to a SeqArray GDS file
#

# Decode ANNOVAR \xHH escape sequences in character vectors
.annovar_decode_hex <- function(x)
{
    idx <- grep("\\\\x[0-9a-fA-F]{2}", x)
    if (length(idx))
    {
        x[idx] <- vapply(x[idx], function(s) {
            while (grepl("\\\\x[0-9a-fA-F]{2}", s))
            {
                m <- regexpr("\\\\x[0-9a-fA-F]{2}", s)
                hex <- substring(s, m + 2, m + 3)
                ch <- rawToChar(as.raw(strtoi(hex, base=16L)))
                s <- sub("\\\\x[0-9a-fA-F]{2}", ch, s, fixed=FALSE)
            }
            s
        }, "", USE.NAMES=FALSE)
    }
    x
}

# decode the \xHH escape sequences in an INFO node, and add the description
#   from 'annovar_output_format.csv' if the node has none
.annovar_node <- function(nd_info, nm, verbose)
{
    nd <- index.gdsn(nd_info, nm)
    dp <- objdesp.gdsn(nd)
    # decode \xHH in character/string nodes
    if (dp$type == "String" || dp$type == "Factor")
    {
        v <- read.gdsn(nd)
        if (any(grepl("\\\\x[0-9a-fA-F]{2}", v)))
        {
            # rewrite the node
            delete.gdsn(nd)
            nd <- add.gdsn(nd_info, nm, .annovar_decode_hex(v),
                storage="string", compress=dp$compress, replace=TRUE)
            if (verbose)
                .cat("    decoded '", nm, "'")
        }
    }
    # add description if available
    j <- match(nm, .packageEnv$annovar$Name)
    desp <- if (is.na(j)) "" else .packageEnv$annovar$Description[j]
    if (nzchar(desp))
    {
        a <- get.attr.gdsn(nd)
        if (is.null(a$Description) || !nzchar(a$Description))
            put.attr.gdsn(nd, "Description", desp)
    }
    invisible()
}

# import an ANNOVAR-annotated VCF file and post-process the INFO nodes
.annovar_vcf <- function(vcf_fn, out_fn, opt, keep, verbose)
{
    # vcf => gds
    attr(verbose, "header_no_time") <- TRUE
    seqVCF2GDS(vcf_fn, out_fn, storage.option=opt$first, optimize=FALSE,
        verbose=verbose)
    # post-process: remove ANNOVAR bookkeeping nodes
    f <- seqOpen(out_fn, readonly=FALSE)
    on.exit(seqClose(f))
    nd_info <- index.gdsn(f, "annotation/info")
    .delete_info(nd_info, c("ANNOVAR_DATE", "ALLELE_END"), verbose)
    # decode \xHH escape sequences in string fields, and add descriptions for
    #   known ANNOVAR fields (each call may rewrite a node)
    nm <- ls.gdsn(nd_info)
    for (s in nm[!startsWith(nm, "@")])
        .annovar_node(nd_info, s, verbose)
    # remove ANNOVAR-specific annotation nodes (Func.*, Gene.*, etc.)
    if (isFALSE(keep))
    {
        nm <- ls.gdsn(nd_info)
        .delete_info(nd_info,
            nm[grepl("^(Func|Gene|GeneDetail|ExonicFunc|AAChange)\\.", nm)],
            verbose)
    }
    invisible()
}

# Convert an ANNOVAR-annotated VCF file to a SeqArray GDS file
#   vcf_fn: the file name(s) of the VCF files written by table_annovar.pl
#       with -vcfinput
#   out_fn: the output GDS file name
#   compress: the compression method of the output file, "LZMA", "ZIP" or
#       "none"
#   keep: if FALSE, remove the gene-based annotation fields (Func.*, Gene.*,
#       GeneDetail.*, ExonicFunc.* and AAChange.*)
#   verbose: if TRUE, show the progress information
# return the normalized output file name invisibly
seqANNOVAR2GDS <- function(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    keep=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf_fn), length(vcf_fn)>0L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    compress <- match.arg(compress)
    stopifnot(is.logical(keep), length(keep)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))
    opt <- .compress_opt(compress)
    # load variable names and descriptions
    .format_table("annovar", "annovar_output_format.csv",
        c("Name", "Description", "Type"))
    .conv_begin("ANNOVAR VCF => GDS", out_fn, compress, verbose)
    # import from VCF
    .annovar_vcf(vcf_fn, out_fn, opt, keep, verbose)
    # recompress?
    .conv_end(out_fn, compress, opt$recompress, verbose)
}
