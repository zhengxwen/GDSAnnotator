#######################################################################
#
# Package Name: GDSAnnotator
# Copyright (C) 2025-2026    Xiuwen Zheng
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Data import & Format conversion from Illumina Nirvana JSON output
#         (https://github.com/Illumina/Nirvana)
#


#######################################################################
# Convert Illumina Nirvana JSON output to a SeqArray GDS file
#
# Nirvana writes a bgzip-compressed JSON file with a streaming layout:
#     {"header":{...},"positions":[
#     {position_1},
#     {position_2},
#     ...
#     ],"genes":[
#     {gene_1},
#     ...
#     ]}
# i.e. the header is on the first line, then one position object per line
# (comma-terminated), then the gene section. This importer reads the file
# line-by-line in blocks so that whole-genome outputs can be converted with
# a small, constant memory footprint.
#

# Load the field catalogue 'nirvana_output_format.csv' (columns Name, Path,
# Type and optionally Default & Description) describing which Nirvana fields
# can be imported and how they are typed. Cached in the package environment.
.nirvana_load_format <- function()
{
    if (is.null(.packageEnv$nirvana))
    {
        fn <- system.file("extdata", "nirvana_output_format.csv",
            package="GDSAnnotator", mustWork=TRUE)
        d <- read.csv(fn, stringsAsFactors=FALSE)
        req <- c("Name", "Path", "Type")
        if (!all(req %in% colnames(d)))
        {
            stop("'nirvana_output_format.csv' should have the columns: ",
                paste(req, collapse=", "), ".")
        }
        if (is.null(d$Default)) d$Default <- TRUE
        if (is.null(d$Description)) d$Description <- ""
        .packageEnv$nirvana <- d
    }
    .packageEnv$nirvana
}

# Resolve the requested 'fields' against the catalogue, returning a data.frame
# with columns Name, Path, Type (Type=NA means "infer from the data").
#   fields=NULL              => all catalogue rows flagged Default=TRUE
#   a catalogue 'Name'       => that row (rename via a named element)
#   any other string         => treated as a JSON dot-path, type inferred
.nirvana_select_fields <- function(fields, fmt)
{
    if (is.null(fields))
    {
        keep <- !is.na(as.logical(fmt$Default)) & as.logical(fmt$Default)
        d <- fmt[keep, c("Name", "Path", "Type", "Description"), drop=FALSE]
    } else {
        nm_in <- names(fields)
        if (is.null(nm_in)) nm_in <- rep("", length(fields))
        rows <- lapply(seq_along(fields), function(i)
        {
            f <- unname(fields[i]); nmi <- nm_in[i]
            j <- match(f, fmt$Name)
            if (!is.na(j))
            {
                r <- fmt[j, c("Name", "Path", "Type", "Description"),
                    drop=FALSE]
                if (nzchar(nmi)) r$Name <- nmi
            } else {
                r <- data.frame(
                    Name = if (nzchar(nmi)) nmi else gsub("[.]", "_", f),
                    Path = f, Type = NA_character_, Description = "",
                    stringsAsFactors = FALSE)
            }
            r
        })
        d <- do.call(rbind, rows)
    }
    if (anyDuplicated(d$Name)) d$Name <- make.unique(d$Name)
    rownames(d) <- NULL
    d
}

# Map a catalogue 'Type' to a GDS storage mode; NA => "" (infer from data).
.nirvana_storage <- function(type, float32)
{
    if (length(type) != 1L || is.na(type)) return("")
    t <- tolower(trimws(as.character(type)))
    if (t %in% c("integer", "int", "int32")) "int32"
    else if (t %in% c("float", "float32", "float64", "numeric", "double",
        "real")) (if (isTRUE(float32)) "float32" else "float64")
    else if (t %in% c("logical", "bool", "boolean", "flag")) "int32"
    else "string"
}

# Walk a parsed Nirvana variant object by a vector of keys and return a flat
# vector of values (or NULL). An array-of-objects extracts 'key' from each
# element, so e.g. c("clinvar","significance") collects the significance of
# every ClinVar record for the variant.
.nirvana_extract <- function(x, parts)
{
    for (key in parts)
    {
        if (is.null(x)) return(NULL)
        nm <- names(x)
        if (is.list(x) && !is.null(nm))
        {
            if (!(key %in% nm)) return(NULL)
            x <- x[[key]]                       # object field
        } else if (is.list(x))
        {
            # array of objects: pull 'key' from each element
            x <- lapply(x, function(e)
                if (is.list(e) && !is.null(e[[key]])) e[[key]] else NULL)
            x <- x[!vapply(x, is.null, FALSE)]
            if (!length(x)) return(NULL)
            x <- unlist(x, recursive=TRUE, use.names=FALSE)
        } else {
            return(NULL)                        # scalar, cannot descend
        }
    }
    unlist(x, recursive=TRUE, use.names=FALSE)
}

# Collapse an extracted value to a single value matching the declared type
# ("string" => unique values joined by 'sep'; numeric => first finite value).
.nirvana_value <- function(x, type, sep)
{
    if (is.null(x) || !length(x))
        return(if (type == "string") NA_character_ else NA_real_)
    if (type == "string")
    {
        x <- as.character(x)
        x <- x[!is.na(x) & nzchar(x)]
        if (!length(x)) NA_character_ else paste(unique(x), collapse=sep)
    } else {
        x <- suppressWarnings(as.numeric(x))
        x <- x[is.finite(x)]
        if (!length(x)) NA_real_ else x[1L]
    }
}

# Infer a GDS storage type from sample values.
.nirvana_infer <- function(vals, float32)
{
    flat <- unlist(vals, use.names=FALSE)
    if (is.logical(flat)) flat <- as.integer(flat)
    flat <- flat[!is.na(flat)]
    if (is.character(flat)) flat <- flat[nzchar(flat)]
    if (!length(flat)) return("string")
    num_t <- if (isTRUE(float32)) "float32" else "float64"
    if (is.numeric(flat))
        return(if (all(flat == floor(flat))) "int32" else num_t)
    sup <- suppressWarnings(as.numeric(flat))
    if (!anyNA(sup))
        return(if (all(sup == floor(sup))) "int32" else num_t)
    "string"
}

# Read & parse the first (header) line; return the parsed header list.
.nirvana_header <- function(con)
{
    line <- readLines(con, n=1L, warn=FALSE)
    if (!length(line) || !grepl('^\\{"header"', line))
    {
        stop("Not a Nirvana JSON file (the file should start with ",
            "'{\"header\":...').")
    }
    # first line: {"header":{...},"positions":[
    s <- sub(',"positions":\\[[[:space:]]*$', '}', line)
    h <- tryCatch(jsonlite::fromJSON(s, simplifyVector=FALSE),
        error=function(e) NULL)
    if (is.null(h)) list() else h$header
}


seqToGDS_Nirvana <- function(json_fn, out_fn, fields=NULL,
    compress=c("LZMA", "ZIP", "none"), root="", rm_chr_prefix=TRUE,
    float32=TRUE, sep=";", bsize=10000L, verbose=TRUE)
{
    # check
    stopifnot(is.character(json_fn), length(json_fn)==1L)
    stopifnot(is.character(out_fn), length(out_fn)==1L)
    if (!is.null(fields))
        stopifnot(is.character(fields), length(fields)>0L)
    compress <- match.arg(compress)
    stopifnot(is.character(root), length(root)==1L, !is.na(root))
    stopifnot(is.logical(rm_chr_prefix), length(rm_chr_prefix)==1L)
    stopifnot(is.logical(float32), length(float32)==1L)
    stopifnot(is.character(sep), length(sep)==1L)
    stopifnot(is.numeric(bsize), length(bsize)==1L, bsize>=1L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))
    if (!requireNamespace("jsonlite", quietly=TRUE))
        stop("The package 'jsonlite' is required: install.packages('jsonlite').")

    # compression algorithm
    map_compress <- c(LZMA="LZMA_RA", ZIP="ZIP_RA", none="")
    compress1 <- compress2 <- map_compress[compress]
    if (compress == "LZMA") compress1 <- "ZIP_RA"  # reduce memory usage

    # resolve fields against the 'nirvana_output_format.csv' catalogue
    fmt <- .nirvana_load_format()
    sel <- .nirvana_select_fields(fields, fmt)
    nf <- nrow(sel)
    if (nf == 0L)
        stop("No fields selected for import.")
    fld_nm <- sel$Name
    fld_desp <- sel$Description
    parts <- strsplit(sel$Path, ".", fixed=TRUE)
    # storage type from the catalogue ("" => infer from the first data block)
    storage_str <- vapply(sel$Type, .nirvana_storage, "", float32=float32)

    if (verbose)
    {
        .cat("##< ", tm())
        .cat("Illumina Nirvana JSON => GDS")
        .cat("    input: ", json_fn)
        .cat("    output: ", out_fn)
        .cat("    compression: ", compress)
        .cat("    fields: ", paste(fld_nm, collapse=", "))
    }

    # open the JSON file (bgzip is gzip-compatible)
    con <- gzfile(json_fn, "rt")
    on.exit(close(con))
    header <- .nirvana_header(con)

    # ---- block reader: read up to 'n' position objects -------------------
    finished <- FALSE
    read_block <- function(n)
    {
        lines <- readLines(con, n=n, warn=FALSE)
        if (!length(lines)) { finished <<- TRUE; return(list()) }
        lines <- trimws(lines)
        # the positions section ends at a line starting with ']'
        endi <- which(startsWith(lines, "]"))
        if (length(endi))
        {
            finished <<- TRUE
            lines <- if (endi[1L] > 1L) lines[seq_len(endi[1L]-1L)] else
                character()
        }
        lines <- sub(",$", "", lines)
        lines <- lines[nzchar(lines)]
        lapply(lines, function(t)
            tryCatch(jsonlite::fromJSON(t, simplifyVector=FALSE),
                error=function(e) NULL))
    }

    # ---- turn a block of positions into a list of per-variant records ----
    make_rows <- function(pos_list, types)
    {
        rows <- lapply(pos_list, function(p)
        {
            if (is.null(p)) return(NULL)
            vars <- p$variants
            if (is.null(vars) || !length(vars)) return(NULL)
            chrom <- p$chromosome
            if (is.null(chrom)) chrom <- NA_character_
            else if (isTRUE(rm_chr_prefix)) chrom <- sub("^chr", "", chrom)
            posn <- if (is.null(p$position)) NA_integer_ else p$position
            refA <- if (is.null(p$refAllele)) NA_character_ else p$refAllele
            alts <- p$altAlleles
            lapply(seq_along(vars), function(i)
            {
                v <- vars[[i]]
                alt <- if (i <= length(alts)) alts[[i]] else v$altAllele
                if (is.null(alt)) alt <- NA_character_
                fv <- lapply(seq_len(nf), function(j)
                    .nirvana_value(.nirvana_extract(v, parts[[j]]),
                        types[j], sep))
                list(chr=as.character(chrom), pos=as.integer(posn),
                    ref=as.character(refA), alt=as.character(alt),
                    id=.nirvana_value(.nirvana_extract(v, "dbsnp"),
                        "string", sep),
                    fld=fv)
            })
        })
        unlist(rows, recursive=FALSE)
    }

    # ---- read the first block; infer types only for custom (unknown) paths
    block1 <- read_block(bsize)
    block1 <- block1[!vapply(block1, is.null, FALSE)]
    if (length(block1) == 0L)
        stop("No variant positions found in the Nirvana JSON file.")
    types <- storage_str
    for (jj in which(types == ""))
    {
        s <- unlist(lapply(block1, function(p)
        {
            vars <- p$variants
            if (is.null(vars)) return(NULL)
            unlist(lapply(vars, function(v) .nirvana_extract(v, parts[[jj]])),
                use.names=FALSE)
        }), use.names=FALSE)
        types[jj] <- .nirvana_infer(s, float32)
    }
    if (verbose)
        .cat("    types: ",
            paste0(fld_nm, "(", types, ")", collapse=", "))

    # ---- create the output GDS skeleton ----------------------------------
    outfile <- .gds_new(out_fn, compress1, var_id_st="int32")
    on.exit(closefn.gds(outfile), add=TRUE)
    # provenance from the Nirvana header (stored on 'description')
    nd_desp <- index.gdsn(outfile, "description")
    put1 <- function(nm, val)
        if (!is.null(val) && length(val))
            put.attr.gdsn(nd_desp, nm, as.character(val)[1L])
    put1("annotator", header$annotator)
    put1("genomeAssembly", header$genomeAssembly)
    put1("schemaVersion", header$schemaVersion)
    put1("dataVersion", header$dataVersion)
    if (!is.null(header$dataSources))
    {
        ds <- vapply(header$dataSources, function(d)
            paste0(d$name, ":", if (is.null(d$version)) "" else d$version), "")
        put.attr.gdsn(nd_desp, "dataSources", paste(ds, collapse="; "))
    }
    # the INFO root (optionally a sub-folder)
    nd_info <- index.gdsn(outfile, "annotation/info")
    if (nzchar(root)) nd_info <- addfolder.gdsn(nd_info, root)
    # one node per requested field
    vcf_type <- c(int32="Integer", float32="Float", float64="Float",
        string="String")
    fld_nodes <- vector("list", nf)
    for (j in seq_len(nf))
    {
        nd <- add.gdsn(nd_info, fld_nm[j], storage=types[j], valdim=0L,
            compress=compress1)
        put.attr.gdsn(nd, "Number", "1")
        put.attr.gdsn(nd, "Type", vcf_type[types[j]])
        desp <- fld_desp[j]
        if (length(desp) != 1L || is.na(desp) || !nzchar(desp))
            desp <- paste0("Nirvana: ", sel$Path[j])
        put.attr.gdsn(nd, "Description", desp)
        fld_nodes[[j]] <- nd
    }
    nd_varid  <- index.gdsn(outfile, "variant.id")
    nd_pos    <- index.gdsn(outfile, "position")
    nd_chr    <- index.gdsn(outfile, "chromosome")
    nd_allele <- index.gdsn(outfile, "allele")
    nd_id     <- index.gdsn(outfile, "annotation/id")
    nd_qual   <- index.gdsn(outfile, "annotation/qual")
    nd_filter <- index.gdsn(outfile, "annotation/filter")

    # ---- append one block of variant records -----------------------------
    nline <- 0L
    field_seen <- logical(nf)
    append_rows <- function(rows)
    {
        rows <- rows[!vapply(rows, is.null, FALSE)]
        n <- length(rows)
        if (!n) return(invisible())
        chr <- vapply(rows, `[[`, "", "chr")
        pos <- vapply(rows, `[[`, 0L, "pos")
        ref <- vapply(rows, `[[`, "", "ref")
        alt <- vapply(rows, `[[`, "", "alt")
        id  <- vapply(rows, function(r)
            if (is.na(r$id)) "" else r$id, "")
        append.gdsn(nd_varid, nline + seq_len(n))
        append.gdsn(nd_pos, as.integer(pos))
        append.gdsn(nd_chr, chr)
        append.gdsn(nd_allele, paste0(ref, ",", alt))
        append.gdsn(nd_id, id)
        append.gdsn(nd_qual, rep(NaN, n))
        append.gdsn(nd_filter, rep(NA_integer_, n))
        for (j in seq_len(nf))
        {
            if (types[j] == "string")
            {
                v <- vapply(rows, function(r)
                {
                    z <- r$fld[[j]]
                    if (is.null(z) || is.na(z)) "" else as.character(z)
                }, "")
                if (any(nzchar(v))) field_seen[j] <<- TRUE
            } else {
                v <- vapply(rows, function(r)
                {
                    z <- r$fld[[j]]
                    if (is.null(z)) NA_real_ else as.numeric(z)
                }, 0)
                if (any(is.finite(v))) field_seen[j] <<- TRUE
            }
            suppressWarnings(append.gdsn(fld_nodes[[j]], v))
        }
        nline <<- nline + n
        invisible()
    }

    # ---- process all blocks ----------------------------------------------
    if (verbose) cat("Processing:\n")
    append_rows(make_rows(block1, types))
    repeat
    {
        if (finished) break
        blk <- read_block(bsize)
        if (length(blk)) append_rows(make_rows(blk, types))
        if (verbose && nline > 0L)
            .cat("    ", prettyNum(nline, big.mark=",", scientific=FALSE),
                " variants\t", tm())
    }

    # ---- finalize --------------------------------------------------------
    for (nm in c("variant.id", "position", "chromosome", "allele",
        "annotation/id", "annotation/qual", "annotation/filter"))
    {
        nd <- index.gdsn(outfile, nm)
        readmode.gdsn(nd)
        digest.gdsn(nd, "md5", action="add")
    }
    # drop fields that never carried a value; finalize the rest
    for (j in seq_len(nf))
    {
        if (!field_seen[j])
        {
            delete.gdsn(fld_nodes[[j]])
        } else {
            readmode.gdsn(fld_nodes[[j]])
            digest.gdsn(fld_nodes[[j]], "md5", action="add")
        }
    }
    if (verbose)
    {
        .cat("Imported ", prettyNum(nline, big.mark=",", scientific=FALSE),
            " variant(s); kept ", sum(field_seen), "/", nf, " field(s).")
    }

    # close input & output, then (re)compress on disk
    on.exit()
    close(con)
    closefn.gds(outfile)
    if (compress1 != compress2)
    {
        if (verbose) .cat("Recompressing (", tm(), ") ...")
        cleanup.gds(out_fn, verbose=FALSE)
        seqRecompress(out_fn, compress=compress, verbose=verbose)
    } else {
        if (verbose) cat("Optimizing the access efficiency ...\n")
        cleanup.gds(out_fn, verbose=verbose)
    }

    if (verbose) .cat("##> ", tm())
    # output
    invisible(normalizePath(out_fn))
}
