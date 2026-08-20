#######################################################################
#
# Package Name: GDSAnnotator
# Copyright (C) 2026    Xiuwen Zheng
#
# Description:
#     Variant annotation data manipulation using GDS files
#     Summary statistics & HTML report (VEP / SnpEff-like)
#


#######################################################################
# Internal: count accumulators
#
# Every statistic is stored as a named numeric vector (category => count),
# so that the counters can be added up across blocks, across GDS files and
# (potentially) across parallel workers with the same piece of code.
#

.acc_new <- function() new.env(parent=emptyenv())

# add a named count vector to the accumulator slot 'key'
.acc_add <- function(e, key, tab)
{
    if (!length(tab)) return(invisible())
    tab <- tab[!is.na(names(tab))]
    if (!length(tab)) return(invisible())
    old <- e[[key]]
    if (is.null(old))
    {
        e[[key]] <- tab
        return(invisible())
    }
    i <- match(names(tab), names(old))
    j <- !is.na(i)
    if (any(j)) old[i[j]] <- old[i[j]] + tab[j]
    if (any(!j)) old <- c(old, tab[!j])
    e[[key]] <- old
    invisible()
}

# tabulate 'x' and add the counts to the accumulator slot 'key'
.acc_tab <- function(e, key, x)
{
    x <- x[!is.na(x)]
    if (length(x)) .acc_add(e, key, table(x, useNA="no"))
    invisible()
}

# element-wise min/max accumulator (e.g., the position range per chromosome)
.acc_minmax <- function(e, key, v, FUN)
{
    if (!length(v)) return(invisible())
    old <- e[[key]]
    if (is.null(old))
    {
        e[[key]] <- v
        return(invisible())
    }
    i <- match(names(v), names(old))
    j <- !is.na(i)
    if (any(j)) old[i[j]] <- FUN(old[i[j]], v[j])
    if (any(!j)) old <- c(old, v[!j])
    e[[key]] <- old
    invisible()
}

# normalize an INFO field to the list(data=, length=) form, so that scalar
# (Number=1) and variable-length annotation fields are handled alike
.as_datalist <- function(v, n)
{
    if (inherits(v, "SeqVarDataList")) v else
        structure(list(length=rep_len(1L, n), data=v), class="SeqVarDataList")
}


#######################################################################
# Internal: the Sequence Ontology consequence table
#
# A single table drives the consequence ranking (most severe consequence
# per variant), the impact classification, the SnpEff-style functional
# class (MISSENSE / NONSENSE / SILENT) and the genomic region. Ensembl-VEP
# and SnpEff both report SO terms, so the two sources share the code.
#

.so_table <- function()
{
    if (is.null(.packageEnv$so))
    {
        so <- read.csv(system.file("extdata", "so_consequence.csv",
            package="GDSAnnotator", mustWork=TRUE), colClasses="character")
        nm <- c("Term", "Rank", "Impact", "Class", "Region")
        if (!all(nm %in% colnames(so)))
        {
            stop("The internal 'so_consequence.csv' should have the ",
                "following columns: ", paste(nm, collapse=","), ".")
        }
        so$Rank <- as.integer(so$Rank)
        .packageEnv$so <- so
    }
    .packageEnv$so
}

# minimum SO rank of each (possibly '&'-joined) consequence string
.so_rank <- function(x)
{
    so <- .so_table()
    r <- so$Rank; names(r) <- so$Term
    ans <- unname(r[x])  # fast path: a single term, no splitting needed
    i <- which(is.na(ans) & grepl("&", x, fixed=TRUE))
    if (length(i))
    {
        ans[i] <- vapply(strsplit(x[i], "&", fixed=TRUE), function(z)
            {
                v <- r[z]
                if (all(is.na(v))) NA_integer_ else min(v, na.rm=TRUE)
            }, 0L)
    }
    ans
}

# split the '&'-joined terms, only for the entries that need it
.so_split <- function(x)
{
    i <- grepl("&", x, fixed=TRUE)
    if (any(i))
        c(x[!i], unlist(strsplit(x[i], "&", fixed=TRUE), use.names=FALSE))
    else
        x
}


#######################################################################
# Internal: annotation field profiles
#
# Map the source-specific sub-field names to a common set of slots, so
# that the same counting code works for Ensembl-VEP and SnpEff (and for
# any other source added later).
#

.annot_slots <- c("consequence", "impact", "gene", "gene_id", "feature",
    "feature_type", "biotype", "sift", "polyphen", "errors")

.annot_profile <- function(f)
{
    hit <- function(s) !is.null(index.gdsn(f, s, silent=TRUE))
    pf <- if (hit("annotation/info/CSQ.list"))
    {
        list(source="Ensembl VEP", root="annotation/info/CSQ.list",
            consequence="Consequence", impact="IMPACT", gene="SYMBOL",
            gene_id="Gene", feature="Feature", feature_type="Feature_type",
            biotype="BIOTYPE", sift="SIFT", polyphen="PolyPhen",
            errors=NA_character_)
    } else if (hit("annotation/info/ANN.list"))
    {
        list(source="SnpEff", root="annotation/info/ANN.list",
            consequence="Annotation", impact="Annotation_Impact",
            gene="Gene_Name", gene_id="Gene_ID", feature="Feature_ID",
            feature_type="Feature_Type", biotype="Transcript_BioType",
            sift=NA_character_, polyphen=NA_character_,
            errors="ERRORS-WARNINGS-INFO")
    } else
        return(NULL)
    # keep only the sub-fields that really exist in the GDS file
    for (s in .annot_slots)
    {
        v <- pf[[s]]
        if (is.na(v) || !hit(paste0(pf$root, "/", v)))
            pf[[s]] <- NA_character_
    }
    pf
}


#######################################################################
# Internal: statistics derived from the 'allele' node
#

.stat_allele <- function(e, al)
{
    # split "REF,ALT1,ALT2,..." without expanding the bi-allelic majority
    ref <- sub(",.*$", "", al)
    rest <- substring(al, nchar(ref) + 2L)
    i <- grepl(",", rest, fixed=TRUE)
    n_alt <- rep.int(1L, length(rest))
    if (any(i))
        n_alt[i] <- lengths(gregexpr(",", rest[i], fixed=TRUE)) + 1L
    .acc_tab(e, "n_allele", n_alt + 1L)  # 2 = bi-allelic
    if (any(i))
    {
        alt <- c(rest[!i], unlist(strsplit(rest[i], ",", fixed=TRUE),
            use.names=FALSE))
        ref <- c(ref[!i], rep.int(ref[i], n_alt[i]))
    } else
        alt <- rest
    # variant class, per alternate allele
    nr <- nchar(ref); na <- nchar(alt)
    cls <- rep.int("SNV", length(alt))
    cls[nr==na & nr>1L] <- "substitution (MNV)"
    cls[nr < na] <- "insertion"
    cls[nr > na] <- "deletion"
    cls[startsWith(alt, "<") | alt=="*" | grepl("[^ACGTNacgtn]", alt)] <- "other"
    .acc_tab(e, "var_class", cls)
    # base changes (SNVs only)
    i <- which(cls == "SNV")
    if (length(i))
        .acc_tab(e, "base_change", paste0(toupper(ref[i]), ">", toupper(alt[i])))
    # InDel length, truncated at +/- 20bp
    i <- which(cls=="insertion" | cls=="deletion")
    if (length(i))
    {
        d <- na[i] - nr[i]
        d[d >  20L] <-  21L
        d[d < -20L] <- -21L
        .acc_tab(e, "indel_len", d)
    }
    invisible()
}


#######################################################################
# Calculate the summary statistics from GDS file(s)
#

# QUAL histogram breaks
.qual_break <- c(-Inf, 0, 10, 20, 30, 40, 50, 100, 200, 500, 1000, Inf)

seqAnnotStat <- function(gdsfile, bsize=100000L, verbose=TRUE)
{
    # check
    stopifnot(is.numeric(bsize), length(bsize)==1L, bsize>0L)
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))
    # multiple files: calculate one by one and merge the counters
    if (is.character(gdsfile))
    {
        if (length(gdsfile) == 0L)
            stop("'gdsfile' should be a file name.")
        if (anyNA(gdsfile))
            stop("'gdsfile' should not contain NA.")
        if (length(gdsfile) > 1L)
        {
            lst <- lapply(seq_along(gdsfile), function(i)
            {
                if (isTRUE(verbose))
                    cat("[", i, "/", length(gdsfile), "] ", sep="")
                seqAnnotStat(gdsfile[i], bsize=bsize, verbose=verbose)
            })
            return(Reduce(.stat_merge, lst))
        }
        # when length(gdsfile)==1
        if (isTRUE(verbose))
            .cat("Open ", sQuote(basename(gdsfile)))
        fn <- normalizePath(gdsfile, mustWork=FALSE)
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
        fn <- gdsfile$filename
    }
    f <- gdsfile

    # the annotation source (VEP, SnpEff, or none)
    pf <- .annot_profile(f)
    # the GDS nodes to be read, in one single pass
    varnm <- c(chr="chromosome", pos="position", allele="allele")
    v <- c(id="annotation/id", qual="annotation/qual",
        filter="annotation/filter")
    for (i in seq_along(v))
    {
        if (!is.null(index.gdsn(f, v[i], silent=TRUE)))
            varnm[names(v)[i]] <- v[i]
    }
    if (!is.null(pf))
    {
        for (s in .annot_slots)
        {
            if (!is.na(pf[[s]]))
                varnm[s] <- paste0(pf$root, "/", pf[[s]])
        }
    }

    # initialize
    e <- .acc_new()
    n_var <- 0L
    so <- .so_table()

    # block-by-block processing, all the counters are updated at once
    seqBlockApply(f, unname(varnm), function(x)
    {
        names(x) <- names(varnm)
        n <- length(x$chr)
        n_var <<- n_var + n

        # ---- variant-level statistics
        .acc_tab(e, "chrom", x$chr)
        .acc_minmax(e, "pos_min", vapply(split(x$pos, x$chr), min, 0L), pmin)
        .acc_minmax(e, "pos_max", vapply(split(x$pos, x$chr), max, 0L), pmax)
        .stat_allele(e, x$allele)
        if (!is.null(x$id))
        {
            v <- !is.na(x$id) & nzchar(x$id) & x$id!="."
            .acc_add(e, "known", c(known=sum(v), novel=sum(!v)))
        }
        if (!is.null(x$qual))
        {
            v <- is.na(x$qual)
            if (any(v)) .acc_add(e, "qual_na", c(`NA`=sum(v)))
            if (!all(v))
                .acc_add(e, "qual", table(cut(x$qual[!v], .qual_break)))
        }
        if (!is.null(x$filter))
            .acc_tab(e, "filter", as.character(x$filter))

        # ---- annotation records
        for (s in .annot_slots)
        {
            if (!is.null(x[[s]])) x[[s]] <- .as_datalist(x[[s]], n)
        }
        if (!is.null(x$consequence))
        {
            ns <- x$consequence$length
            d <- x$consequence$data
            .acc_add(e, "n_annot", c(n=length(d)))
            # (a) all consequence terms, with '&' expanded
            trm <- .so_split(d)
            .acc_tab(e, "cons_all", trm)
            # (b) genomic region & functional class of every term
            i <- match(trm, so$Term)
            .acc_tab(e, "region", so$Region[i])
            v <- so$Class[i]
            .acc_tab(e, "func_class", v[!is.na(v) & nzchar(v)])
            # (c) the most severe consequence of each variant
            rk <- .so_rank(d)
            iv <- rep.int(seq_along(ns), ns)
            k <- order(iv, rk, na.last=TRUE)
            k <- k[!duplicated(iv[k])]
            i <- match(rk[k], so$Rank)
            .acc_tab(e, "cons_severe", so$Term[i])
            .acc_tab(e, "impact_severe", so$Impact[i])
        }
        for (s in c("impact", "biotype", "feature_type", "feature", "sift",
            "polyphen"))
        {
            v <- x[[s]]
            if (!is.null(v))
            {
                v <- v$data
                v <- v[!is.na(v) & nzchar(v)]
                # e.g., SIFT is stored as "deleterious(0.02)"
                if (s=="sift" || s=="polyphen") v <- sub("\\(.*$", "", v)
                .acc_tab(e, s, v)
            }
        }
        if (!is.null(x$errors))
        {
            v <- x$errors$data
            .acc_tab(e, "errors", .so_split(v[!is.na(v) & nzchar(v)]))
        }
        # gene x impact counts (the SnpEff 'genes.txt' equivalent)
        if (!is.null(x$gene))
        {
            g <- x$gene$data
            im <- if (!is.null(x$impact)) x$impact$data else rep.int("", length(g))
            i <- !is.na(g) & nzchar(g)
            if (any(i))
                .acc_tab(e, "gene_impact", paste0(g[i], "\r", im[i]))
        }
        NULL  # return
    }, as.is="none", bsize=bsize, .progress=verbose)

    # the VCF header, for the annotator version & command line
    hd <- NULL
    nd <- index.gdsn(f, "description/vcf.header", silent=TRUE)
    if (!is.null(nd)) hd <- read.gdsn(nd)

    # output
    ans <- list(file=fn, filesize=file.size(fn), n_variant=n_var,
        source=if (is.null(pf)) "unknown" else pf$source, profile=pf,
        header=hd, counts=as.list(e), created=Sys.time())
    class(ans) <- "SeqAnnotStat"
    ans
}


# Merge two SeqAnnotStat objects
.stat_merge <- function(x, y)
{
    e <- .acc_new()
    for (s in union(names(x$counts), names(y$counts)))
    {
        if (s == "pos_min")
        {
            .acc_minmax(e, s, x$counts[[s]], pmin)
            .acc_minmax(e, s, y$counts[[s]], pmin)
        } else if (s == "pos_max")
        {
            .acc_minmax(e, s, x$counts[[s]], pmax)
            .acc_minmax(e, s, y$counts[[s]], pmax)
        } else {
            .acc_add(e, s, x$counts[[s]])
            .acc_add(e, s, y$counts[[s]])
        }
    }
    x$counts <- as.list(e)
    x$n_variant <- x$n_variant + y$n_variant
    x$file <- c(x$file, y$file)
    x$filesize <- c(x$filesize, y$filesize)
    if (is.null(x$header)) x$header <- y$header
    x
}


print.SeqAnnotStat <- function(x, ...)
{
    .cat("SeqAnnotStat [", x$source, "]")
    .cat("    file: ", paste(basename(x$file), collapse=", "))
    .cat("    # of variants: ", format(x$n_variant, big.mark=","))
    if (!is.null(x$counts$n_annot))
    {
        .cat("    # of annotation records: ",
            format(x$counts$n_annot[[1L]], big.mark=","))
    }
    .cat("    statistics: ", paste(names(x$counts), collapse=", "))
    invisible(x)
}


# Return the per-gene counts as a data.frame (the SnpEff 'genes.txt'
# equivalent): one row per gene, one column per impact category
seqAnnotGeneTable <- function(stat)
{
    stopifnot(inherits(stat, "SeqAnnotStat"))
    v <- stat$counts$gene_impact
    if (is.null(v)) return(NULL)
    s <- strsplit(names(v), "\r", fixed=TRUE)
    gene <- vapply(s, `[`, "", i=1L)
    im <- vapply(s, function(z) if (length(z)>1L) z[2L] else "", "")
    im[!nzchar(im)] <- "unknown"
    ug <- sort(unique(gene))
    ui <- unique(im)
    ui <- c(intersect(c("HIGH","MODERATE","LOW","MODIFIER"), ui),
        setdiff(ui, c("HIGH","MODERATE","LOW","MODIFIER")))
    cnt <- matrix(0L, length(ug), length(ui), dimnames=list(NULL, ui))
    cnt[cbind(match(gene, ug), match(im, ui))] <- as.integer(v)
    df <- data.frame(gene=ug, total=rowSums(cnt), stringsAsFactors=FALSE)
    df <- cbind(df, as.data.frame(cnt))
    df <- df[order(df$total, decreasing=TRUE), ]
    rownames(df) <- NULL
    df
}


#######################################################################
# Internal: build the report content
#
# The report is first described as a list of sections, independently of
# the output format; each section is a title plus a list of items. The
# format-specific renderers below turn the items into HTML or Markdown,
# so that the two output formats always show the same content.
#
# Item types:
#     kv         key-value table (the general information)
#     bar        horizontal bar chart of a named count vector
#     col        column chart, e.g., the counts per chromosome
#     table      count / percentage table
#     basechange the 4x4 REF-by-ALT matrix with the Ts/Tv ratio
#     note       a short remark
#

.item <- function(type, ...) c(list(type=type), list(...))

# format a count with thousands separators
.fmt <- function(x)
    format(x, big.mark=",", trim=TRUE, scientific=FALSE)

# get a value from the VCF header data.frame
.hdr_val <- function(hd, id)
{
    if (is.null(hd)) return(NA_character_)
    i <- match(id, hd$id)
    if (is.na(i)) NA_character_ else hd$value[i]
}

# order the chromosomes numerically, then alphabetically
.chrom_order <- function(s)
{
    v <- suppressWarnings(as.integer(sub("^chr", "", s)))
    i <- is.na(v)
    if (any(i)) v[i] <- 1000L + order(s[i])
    order(v)
}

.rep_sections <- function(stat)
{
    cnt <- stat$counts
    hd <- stat$header
    sec <- list()
    add <- function(title, ...)
        sec[[length(sec)+1L]] <<- list(title=title, items=list(...))

    # ---- general information
    gi <- c("Report generated"=format(stat$created, "%Y-%m-%d %H:%M:%S"),
        "GDS file"=paste(basename(stat$file), collapse=", "),
        "File size"=paste(format(round(sum(stat$filesize)/1024^2, 1)), "MB"),
        "Annotation source"=stat$source)
    for (s in c("VEP", "VEP-command-line", "SnpEffVersion", "SnpEffCmd",
        "source", "reference"))
    {
        v <- .hdr_val(hd, s)
        if (!is.na(v)) gi[s] <- v
    }
    gi["# of variants"] <- .fmt(stat$n_variant)
    if (!is.null(cnt$n_annot))
        gi["# of annotation records"] <- .fmt(cnt$n_annot[[1L]])
    if (!is.null(cnt$gene_impact))
    {
        gi["# of genes with annotation"] <- .fmt(length(unique(
            sub("\r.*$", "", names(cnt$gene_impact)))))
    }
    if (!is.null(cnt$feature))
        gi["# of features (transcripts)"] <- .fmt(length(cnt$feature))
    if (!is.null(cnt$pos_min) && !is.null(cnt$pos_max))
    {
        v <- sum(as.double(cnt$pos_max) - as.double(cnt$pos_min) + 1)
        gi["Genomic span covered"] <- paste(.fmt(round(v/1e6, 1)), "Mb")
        gi["Variant rate"] <- paste("one variant every",
            .fmt(round(v/stat$n_variant)), "bases")
    }
    add("General information", .item("kv", key=names(gi), value=unname(gi)))

    # ---- variant-level sections
    v <- cnt$chrom[.chrom_order(names(cnt$chrom))]
    add("Variants by chromosome",
        .item("col", tab=v, sorted=FALSE,
            colnm=c("Chromosome", "Count", "Percentage")),
        .item("table", tab=v, sorted=FALSE,
            colnm=c("Chromosome", "Count", "Percentage")))
    add("Variant classes", .item("bar", tab=cnt$var_class),
        .item("table", tab=cnt$var_class))
    add("Base changes (SNVs)", .item("basechange", tab=cnt$base_change))
    if (!is.null(cnt$indel_len))
    {
        v <- cnt$indel_len
        v <- v[order(as.integer(names(v)))]
        # the length is truncated at +/- 20bp in seqAnnotStat()
        names(v)[names(v) == "-21"] <- "<=-21"
        names(v)[names(v) ==  "21"] <- ">=21"
        add("InDel length distribution",
            .item("col", tab=v, sorted=FALSE,
                color="#e8913a",
                colnm=c("InDel length (bp)", "Count", "Percentage")))
    }
    add("Number of alleles per variant",
        .item("table", tab=cnt$n_allele, sorted=FALSE,
            colnm=c("# of alleles", "Count", "Percentage")))
    if (!is.null(cnt$known))
        add("Known vs. novel (the ID field)", .item("table", tab=cnt$known))
    if (!is.null(cnt$filter))
        add("FILTER", .item("table", tab=cnt$filter))
    if (!is.null(cnt$qual))
        add("QUAL distribution", .item("bar", tab=cnt$qual, sorted=FALSE))

    # ---- annotation sections
    if (!is.null(cnt$cons_all))
    {
        add("Consequences (all annotation records)",
            .item("bar", tab=cnt$cons_all), .item("table", tab=cnt$cons_all))
        add("Consequences (most severe per variant)",
            .item("bar", tab=cnt$cons_severe),
            .item("table", tab=cnt$cons_severe))
        add("Variants by region",
            .item("bar", tab=cnt$region, color="#7a5ea8"))
    }
    if (!is.null(cnt$impact))
    {
        add("Effects by impact (all records)",
            .item("bar", tab=cnt$impact, color=.impact_color, sorted=FALSE),
            .item("table", tab=cnt$impact))
    }
    if (!is.null(cnt$impact_severe))
    {
        add("Impact (most severe per variant)",
            .item("bar", tab=cnt$impact_severe, color=.impact_color,
                sorted=FALSE))
    }
    if (!is.null(cnt$func_class))
    {
        v <- cnt$func_class
        it <- list(.item("bar", tab=v, color="#4a9d5b"), .item("table", tab=v))
        if (all(c("MISSENSE", "SILENT") %in% names(v)))
        {
            it <- c(it, list(.item("note", text=sprintf(
                "Missense / Silent ratio = %.4f", v[["MISSENSE"]]/v[["SILENT"]]),
                bold=TRUE)))
        }
        sec[[length(sec)+1L]] <- list(title="Effects by functional class",
            items=it)
    }
    if (!is.null(cnt$biotype))
        add("Transcript biotypes", .item("bar", tab=cnt$biotype))
    if (!is.null(cnt$feature_type))
        add("Feature types", .item("table", tab=cnt$feature_type))
    if (!is.null(cnt$sift))
        add("SIFT predictions", .item("table", tab=cnt$sift))
    if (!is.null(cnt$polyphen))
        add("PolyPhen predictions", .item("table", tab=cnt$polyphen))
    if (!is.null(cnt$errors))
    {
        add("Warnings & errors reported by the annotator",
            .item("table", tab=cnt$errors))
    }
    if (!is.null(cnt$gene_impact))
    {
        v <- cnt$gene_impact
        v <- vapply(split(unname(v), sub("\r.*$", "", names(v))), sum, 0)
        add("Top annotated genes", .item("bar", tab=v, top=20L),
            .item("note", text=paste("The full gene-by-impact table is",
                "available from seqAnnotGeneTable().")))
    }

    sec
}

# sort and truncate a count vector, as requested by a 'bar' item
.bar_tab <- function(it)
{
    tab <- it$tab
    if (is.null(it$sorted) || isTRUE(it$sorted))
        tab <- sort(tab, decreasing=TRUE)
    top <- if (is.null(it$top)) 25L else it$top
    n_more <- max(0L, length(tab) - top)
    if (n_more > 0L) tab <- tab[seq_len(top)]
    list(tab=tab, n_more=n_more)
}

# the 4x4 REF-by-ALT matrix and the Ts/Tv ratio
.base_matrix <- function(tab)
{
    b <- c("A", "C", "G", "T")
    m <- matrix(0, 4L, 4L, dimnames=list(b, b))
    s <- strsplit(names(tab), ">", fixed=TRUE)
    for (i in seq_along(tab))
    {
        r <- s[[i]][1L]; a <- s[[i]][2L]
        if (r %in% b && a %in% b) m[r, a] <- m[r, a] + tab[[i]]
    }
    n_ts <- m["A","G"] + m["G","A"] + m["C","T"] + m["T","C"]
    list(m=m, n_ts=n_ts, n_tv=sum(m) - n_ts - sum(diag(m)))
}


#######################################################################
# Internal: the HTML renderer
#
# Plain HTML tables and inline SVG charts only: no JavaScript, no pandoc,
# and no additional package is needed, so the output is a single file
# that can be viewed offline.
#

.esc <- function(s)
{
    s <- gsub("&", "&amp;", s, fixed=TRUE)
    s <- gsub("<", "&lt;", s, fixed=TRUE)
    gsub(">", "&gt;", s, fixed=TRUE)
}

.impact_color <- c(
    HIGH="#d64550", MODERATE="#e8913a", LOW="#4a9d5b", MODIFIER="#6c8ebf")

# horizontal bar chart with counts & percentages, as inline SVG
.svg_bar <- function(it)
{
    if (is.null(it$tab) || !length(it$tab))
        return("<p><i>not available</i></p>")
    z <- .bar_tab(it)
    tab <- z$tab
    color <- if (is.null(it$color)) "#4a76b8" else it$color
    n <- length(tab); tot <- sum(tab)
    rh <- 20L; wd <- 760L; lab_w <- 250L; bar_w <- 380L
    hg <- n*rh + 6L
    y <- (seq_len(n)-1L)*rh + 14L
    len <- as.integer(pmax(1, round(bar_w * tab / max(tab))))
    col <- if (length(color) == 1L) rep_len(color, n) else
        ifelse(is.na(color[names(tab)]), "#4a76b8", color[names(tab)])
    s <- sprintf(paste0(
        '<text x="%d" y="%d" class="lab">%s</text>',
        '<rect x="%d" y="%d" width="%d" height="12" fill="%s"/>',
        '<text x="%d" y="%d" class="val">%s (%.2f%%)</text>'),
        lab_w-4L, y, .esc(names(tab)), lab_w+2L, y-10L, len, col,
        as.integer(lab_w+8L+len), y, .fmt(tab), 100*tab/tot)
    paste0('<svg class="chart" viewBox="0 0 ', wd, ' ', hg,
        '" width="100%" height="', hg, '">', paste(s, collapse=""), '</svg>',
        if (z$n_more > 0L)
            paste0('<p class="note">... ', z$n_more,
                ' more categories not shown</p>')
        else "")
}

# vertical column chart (e.g., the number of variants per chromosome)
.svg_col <- function(it)
{
    tab <- it$tab
    if (is.null(tab) || !length(tab)) return("<p><i>not available</i></p>")
    color <- if (is.null(it$color)) "#4a76b8" else it$color
    n <- length(tab); wd <- max(760L, n*26L); hg <- 200L; bar_h <- 150L
    bw <- as.integer((wd - 40L) %/% n)
    x <- 30L + (seq_len(n)-1L)*bw
    h <- as.integer(pmax(1, round(bar_h * tab / max(tab))))
    s <- sprintf(paste0(
        '<rect x="%d" y="%d" width="%d" height="%d" fill="%s">',
        '<title>%s: %s</title></rect>',
        '<text x="%d" y="%d" class="ax" text-anchor="middle">%s</text>'),
        x, 20L+bar_h-h, max(bw-3L, 2L), h, color, .esc(names(tab)), .fmt(tab),
        as.integer(x+bw/2), 20L+bar_h+14L, .esc(names(tab)))
    paste0('<svg class="chart" viewBox="0 0 ', wd, ' ', hg, '" width="100%">',
        sprintf('<text x="0" y="24" class="ax">%s</text>', .fmt(max(tab))),
        paste(s, collapse=""), '</svg>')
}

.html_table <- function(it)
{
    tab <- it$tab
    if (is.null(tab) || !length(tab)) return("<p><i>not available</i></p>")
    if (is.null(it$sorted) || isTRUE(it$sorted))
        tab <- sort(tab, decreasing=TRUE)
    colnm <- if (is.null(it$colnm)) c("Category", "Count", "Percentage") else
        it$colnm
    tot <- sum(tab)
    paste0('<table><tr>', paste0('<th>', colnm, '</th>', collapse=""),
        '</tr>', paste0(sprintf(
            '<tr><td>%s</td><td class="n">%s</td><td class="n">%.2f%%</td></tr>',
            .esc(names(tab)), .fmt(tab), 100*tab/tot), collapse=""),
        '</table>')
}

.html_keyval <- function(it)
{
    paste0('<table class="kv">', paste0(sprintf(
        '<tr><th>%s</th><td>%s</td></tr>', .esc(it$key), .esc(it$value)),
        collapse=""), '</table>')
}

.html_base_change <- function(it)
{
    if (is.null(it$tab) || !length(it$tab))
        return("<p><i>not available</i></p>")
    z <- .base_matrix(it$tab)
    m <- z$m; b <- rownames(m); mx <- max(m)
    cell <- function(i, j)
    {
        if (i == j) return('<td class="diag">-</td>')
        sprintf('<td class="n" style="background:rgba(74,118,184,%.2f)">%s</td>',
            if (mx>0) 0.75*m[i,j]/mx else 0, .fmt(m[i,j]))
    }
    rows <- vapply(b, function(i) paste0('<tr><th>', i, '</th>',
        paste0(vapply(b, function(j) cell(i, j), ""), collapse=""), '</tr>'), "")
    paste0('<table><tr><th>REF \\ ALT</th>',
        paste0('<th>', b, '</th>', collapse=""), '</tr>',
        paste(rows, collapse=""), '</table>',
        sprintf(paste0('<p class="note">Transitions: %s &nbsp;|&nbsp; ',
            'Transversions: %s &nbsp;|&nbsp; <b>Ts/Tv = %.3f</b></p>'),
            .fmt(z$n_ts), .fmt(z$n_tv), z$n_ts/z$n_tv))
}

.html_note <- function(it)
{
    s <- .esc(it$text)
    if (isTRUE(it$bold)) s <- paste0("<b>", s, "</b>")
    paste0('<p class="note">', s, '</p>')
}

.html_css <- 'body{font-family:sans-serif;
margin:0 auto;max-width:1000px;padding:24px;color:#222;line-height:1.45}
h1{border-bottom:3px solid #4a76b8;padding-bottom:8px}
h2{margin-top:34px;font-size:1.15em;background:#eef2f8;padding:6px 10px;
border-left:4px solid #4a76b8}
table{border-collapse:collapse;margin:10px 0;font-size:.9em}
th,td{border:1px solid #ccd;padding:3px 10px;text-align:left}
th{background:#f2f5fa}
td.n{text-align:right;font-variant-numeric:tabular-nums}
td.diag{background:#eee;text-align:center}
table.kv th{width:230px}
table.kv td{font-family:monospace;font-size:.85em;
word-break:break-all}
svg.chart{margin:8px 0}
svg text.lab{font-size:11px;text-anchor:end;fill:#333}
svg text.val{font-size:11px;fill:#555}
svg text.ax{font-size:10px;fill:#666}
p.note{color:#666;font-size:.85em}
code{background:#f4f4f4;padding:1px 4px}'

.render_html <- function(sec, title)
{
    body <- vapply(sec, function(s)
    {
        v <- vapply(s$items, function(it) switch(it$type,
            kv = .html_keyval(it),
            bar = .svg_bar(it),
            col = .svg_col(it),
            table = .html_table(it),
            basechange = .html_base_change(it),
            note = .html_note(it), ""), "")
        paste0('<h2>', .esc(s$title), '</h2>', paste(v, collapse=""))
    }, "")
    paste0('<!DOCTYPE html>\n<html><head><meta charset="utf-8">',
        '<title>', .esc(title), '</title>\n<style>', .html_css,
        '</style></head>\n<body><h1>', .esc(title), '</h1>\n',
        paste(body, collapse="\n"),
        '\n<hr><p class="note">Generated by GDSAnnotator::seqAnnotReport()',
        '</p></body></html>')
}


#######################################################################
# Internal: the Markdown renderer
#
# GitHub-flavored Markdown with pipe tables; the bar charts are drawn
# with block characters, so the output is readable as plain text as well.
#

.md_esc <- function(s)
    gsub("|", "\\|", s, fixed=TRUE)

.md_bar <- function(x, width=24L)
{
    n <- as.integer(round(width * x / max(x)))
    # "\u2588" is the full block character
    vapply(n, function(i) strrep("\u2588", max(i, 1L)), "")
}

.md_table <- function(colnm, ...)
{
    col <- list(...)
    s <- vapply(col, function(z) paste0("| ", z, " "), character(length(col[[1L]])))
    if (!is.matrix(s)) s <- matrix(s, nrow=1L)
    c(paste0("| ", paste(colnm, collapse=" | "), " |"),
        paste0("|", paste(rep(" --- ", length(colnm)), collapse="|"), "|"),
        paste0(apply(s, 1L, paste, collapse=""), "|"), "")
}

.md_count_table <- function(it, with_bar=FALSE)
{
    tab <- it$tab
    if (is.null(tab) || !length(tab)) return(c("*not available*", ""))
    n_more <- 0L
    if (isTRUE(with_bar))
    {
        z <- .bar_tab(it); tab <- z$tab; n_more <- z$n_more
    } else if (is.null(it$sorted) || isTRUE(it$sorted))
        tab <- sort(tab, decreasing=TRUE)
    colnm <- if (is.null(it$colnm)) c("Category", "Count", "Percentage") else
        it$colnm
    tot <- sum(tab)
    s <- if (isTRUE(with_bar))
    {
        .md_table(c(colnm, ""), .md_esc(names(tab)), .fmt(tab),
            sprintf("%.2f%%", 100*tab/tot), .md_bar(tab))
    } else {
        .md_table(colnm, .md_esc(names(tab)), .fmt(tab),
            sprintf("%.2f%%", 100*tab/tot))
    }
    if (n_more > 0L)
        s <- c(head(s, -1L), paste0("*... ", n_more,
            " more categories not shown*"), "")
    s
}

.md_base_change <- function(it)
{
    if (is.null(it$tab) || !length(it$tab)) return(c("*not available*", ""))
    z <- .base_matrix(it$tab)
    m <- z$m; b <- rownames(m)
    s <- c(paste0("| REF \\ ALT | ", paste(b, collapse=" | "), " |"),
        paste0("|", paste(rep(" --- ", length(b)+1L), collapse="|"), "|"),
        vapply(b, function(i) paste0("| **", i, "** | ", paste(vapply(b,
            function(j) if (i==j) "-" else .fmt(m[i,j]), ""),
            collapse=" | "), " |"), ""), "")
    c(unname(s), sprintf(
        "Transitions: %s | Transversions: %s | **Ts/Tv = %.3f**",
        .fmt(z$n_ts), .fmt(z$n_tv), z$n_ts/z$n_tv), "")
}

# a signature of the data shown by an item, to detect duplicated content
.md_key <- function(tab)
{
    tab <- sort(tab, decreasing=TRUE)
    paste(names(tab), tab, sep="=", collapse=",")
}

.render_md <- function(sec, title)
{
    ans <- c(paste("#", title), "")
    for (s in sec)
    {
        ans <- c(ans, paste("##", s$title), "")
        shown <- character()
        for (it in s$items)
        {
            if (it$type=="bar" || it$type=="col")
            {
                z <- .bar_tab(it)
                # a truncated chart is still followed by the full table
                if (z$n_more == 0L) shown <- c(shown, .md_key(it$tab))
            } else if (it$type == "table")
            {
                # in Markdown a chart is rendered as a table, so a plain
                # table of the same counts would just repeat it
                if (.md_key(it$tab) %in% shown) next
            }
            ans <- c(ans, switch(it$type,
                kv = .md_table(c("Item", "Value"), .md_esc(it$key),
                    paste0("`", .md_esc(it$value), "`")),
                bar = ,
                col = .md_count_table(it, with_bar=TRUE),
                table = .md_count_table(it),
                basechange = .md_base_change(it),
                note = c(if (isTRUE(it$bold)) paste0("**", it$text, "**")
                    else paste0("*", it$text, "*"), ""),
                NULL))
        }
    }
    c(ans, "---", "", "Generated by `GDSAnnotator::seqAnnotReport()`")
}


#######################################################################
# Internal: the R Markdown renderer
#
# The counts are embedded in the document as R code, and the charts are
# drawn by the code chunks when the document is knitted. The output is
# still a single file with no companion image directory, but the charts
# are real plots that the user can restyle.
#

# embed a vector as R code, wrapped over several lines if needed
.rmd_data <- function(x, quote=FALSE)
{
    x <- if (isTRUE(quote)) as.character(x) else as.numeric(x)
    paste(deparse(unname(x), width.cutoff=60L), collapse="\n    ")
}

# embed a named count vector as R code
.rmd_count <- function(x)
    paste0("structure(", .rmd_data(x), ",\n    names = ",
        .rmd_data(names(x), quote=TRUE), ")")

.rmd_chunk <- function(label, opt, code)
{
    c(paste0("```{r ", label, if (nzchar(opt)) paste0(", ", opt) else "", "}"),
        code, "```", "")
}

# the left margin needed by the category labels of a horizontal bar chart
.rmd_margin <- function(nm)
    max(6, min(24, round(0.55 * max(nchar(nm)) + 2)))

.rmd_bar <- function(it, label)
{
    if (is.null(it$tab) || !length(it$tab)) return(c("*not available*", ""))
    z <- .bar_tab(it)
    tab <- rev(z$tab)  # barplot() draws the first element at the bottom
    color <- if (is.null(it$color)) "#4a76b8" else it$color
    col <- if (length(color) == 1L) color else
        ifelse(is.na(color[names(tab)]), "#4a76b8", color[names(tab)])
    c(.rmd_chunk(label,
        sprintf("fig.height=%.1f, fig.width=7",
            max(2, 0.24*length(tab) + 1)),
        c(paste0("x <- ", .rmd_count(tab)),
            sprintf("op <- par(mar=c(4, %d, 1, 2), cex=0.8)",
                .rmd_margin(names(tab))),
            sprintf("barplot(x, horiz=TRUE, las=1, xlab=\"Count\", col=%s)",
                if (length(col) == 1L) paste0('"', col, '"')
                else .rmd_count_col(col)),
            "par(op)")),
        if (z$n_more > 0L)
            c(paste0("*... ", z$n_more, " more categories not shown*"), "")
        else NULL)
}

# a character vector of colors, as R code
.rmd_count_col <- function(x)
    .rmd_data(unname(x), quote=TRUE)

.rmd_col <- function(it, label)
{
    tab <- it$tab
    if (is.null(tab) || !length(tab)) return(c("*not available*", ""))
    color <- if (is.null(it$color)) "#4a76b8" else it$color
    .rmd_chunk(label, "fig.height=3.2, fig.width=7",
        c(paste0("x <- ", .rmd_count(tab)),
            "op <- par(mar=c(4, 4, 1, 2), cex=0.8)",
            sprintf("barplot(x, las=2, ylab=\"Count\", col=\"%s\")",
                color[1L]),
            "par(op)"))
}

.rmd_table <- function(it, label)
{
    tab <- it$tab
    if (is.null(tab) || !length(tab)) return(c("*not available*", ""))
    if (is.null(it$sorted) || isTRUE(it$sorted))
        tab <- sort(tab, decreasing=TRUE)
    colnm <- if (is.null(it$colnm)) c("Category", "Count", "Percentage") else
        it$colnm
    .rmd_chunk(label, "",
        c(paste0("x <- ", .rmd_count(tab)),
            "knitr::kable(data.frame(names(x), format(x, big.mark=\",\"),",
            "        sprintf(\"%.2f%%\", 100*x/sum(x)), row.names=NULL),",
            paste0("    col.names = ", .rmd_data(colnm, quote=TRUE), ")")))
}

.rmd_keyval <- function(it, label)
{
    .rmd_chunk(label, "",
        c(paste0("x <- data.frame(Item = ", .rmd_data(it$key, quote=TRUE), ","),
            paste0("    Value = ", .rmd_data(it$value, quote=TRUE), ")"),
            "knitr::kable(x)"))
}

.rmd_base_change <- function(it, label)
{
    if (is.null(it$tab) || !length(it$tab)) return(c("*not available*", ""))
    z <- .base_matrix(it$tab)
    b <- rownames(z$m)
    c(.rmd_chunk(label, "",
        c(paste0("m <- matrix(", .rmd_data(as.vector(z$m)), ", 4L, 4L,"),
            paste0("    dimnames = list(", .rmd_data(b, quote=TRUE), ", ",
                .rmd_data(b, quote=TRUE), "))"),
            "knitr::kable(m, caption=\"REF (row) by ALT (column)\")")),
        sprintf("Transitions: %s | Transversions: %s | **Ts/Tv = %.3f**",
            .fmt(z$n_ts), .fmt(z$n_tv), z$n_ts/z$n_tv), "")
}

.render_rmd <- function(sec, title)
{
    ans <- c("---",
        paste0('title: "', gsub('"', "'", title, fixed=TRUE), '"'),
        "output:", "  html_document:", "    toc: true", "    toc_depth: 2",
        "---", "",
        .rmd_chunk("setup", "include=FALSE",
            paste0("# set 'echo=TRUE' to show the code of each chart\n",
                "knitr::opts_chunk$set(echo=FALSE, comment=NA)")))
    for (i in seq_along(sec))
    {
        s <- sec[[i]]
        ans <- c(ans, paste("##", s$title), "")
        for (j in seq_along(s$items))
        {
            it <- s$items[[j]]
            lab <- paste0("s", i, "_", j)
            ans <- c(ans, switch(it$type,
                kv = .rmd_keyval(it, lab),
                bar = .rmd_bar(it, lab),
                col = .rmd_col(it, lab),
                table = .rmd_table(it, lab),
                basechange = .rmd_base_change(it, lab),
                note = c(if (isTRUE(it$bold)) paste0("**", it$text, "**")
                    else paste0("*", it$text, "*"), ""),
                NULL))
        }
    }
    c(ans, "---", "", "Generated by `GDSAnnotator::seqAnnotReport()`")
}


#######################################################################
# Generate a summary report from a GDS file
#

seqAnnotReport <- function(stat, out_fn="annot_summary.html",
    format=c("auto", "html", "markdown", "rmd"),
    title="GDSAnnotator Summary Report", bsize=100000L, verbose=TRUE)
{
    # check
    stopifnot(is.character(out_fn), length(out_fn)==1L, !is.na(out_fn))
    format <- match.arg(format)
    stopifnot(is.character(title), length(title)==1L, !is.na(title))
    stopifnot(is.logical(verbose), length(verbose)==1L, !is.na(verbose))
    # the output format, guessed from the file extension if "auto"
    if (format == "auto")
    {
        s <- tolower(sub("^.*\\.", "", basename(out_fn)))
        format <- switch(s, rmd="rmd", md=, markdown="markdown", "html")
    }
    # the summary statistics
    if (!inherits(stat, "SeqAnnotStat"))
        stat <- seqAnnotStat(stat, bsize=bsize, verbose=verbose)

    # build the format-independent content, then render it
    sec <- .rep_sections(stat)
    s <- switch(format,
        html = .render_html(sec, title),
        markdown = .render_md(sec, title),
        rmd = .render_rmd(sec, title))
    writeLines(s, out_fn)

    if (verbose)
        .cat("Report [", format, "]: ", normalizePath(out_fn))
    # output
    invisible(normalizePath(out_fn))
}
