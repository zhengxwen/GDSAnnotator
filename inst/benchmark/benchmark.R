#######################################################################
#
# GDSAnnotator benchmark
#
# Measures (1) on-disk storage of an annotation GDS file versus an
# equivalent flat text table, and (2) random-access query latency for a
# range of query sizes. Designed to be run on a real, whole-genome
# annotation GDS (e.g. the full FAVOR database), where the advantages of
# the columnar, block-compressed GDS layout are most pronounced.
#
# Usage:
#   Rscript benchmark.R  /path/to/annotation.gds  [n_annotations]
#
# If no GDS file is supplied, the small packaged FAVOR example is used.
#

suppressMessages({
    library(gdsfmt)
    library(SeqArray)
    library(GDSAnnotator)
})

args <- commandArgs(trailingOnly = TRUE)
gds <- if (length(args) >= 1L) args[1L] else
    system.file("extdata", "favor_chr22_sub.gds", package = "GDSAnnotator")
max_ann <- if (length(args) >= 2L) as.integer(args[2L]) else Inf
stopifnot(file.exists(gds))

cat("Annotation GDS :", gds, "\n")
ann <- seqAnnotList(gds)
varnm <- head(ann$name, max_ann)
cat("Annotations    :", length(varnm), "of", nrow(ann), "\n")

f <- seqOpen(gds, allow.duplicate = TRUE)
n_var <- objdesp.gdsn(index.gdsn(f, "variant.id"))$dim
chr <- seqGetData(f, "$chromosome"); pos <- seqGetData(f, "position")
ref <- seqGetData(f, "$ref");        alt <- seqGetData(f, "$alt")
seqClose(f)
cat("Variants       :", prettyNum(n_var, big.mark = ","), "\n\n")

## ---- 1. storage comparison -----------------------------------------------
# Stream a flat text dump so that very large files do not blow up memory.
cat("== Storage ==\n")
sz_gds <- file.size(gds)
cat(sprintf("  GDS                : %10.2f MB\n", sz_gds / 1024^2))

# Estimate flat-text size from a sample of variants (to keep this cheap on
# whole-genome inputs), then scale up.
samp <- sort(sample(n_var, min(n_var, 5000L)))
f <- seqOpen(gds, allow.duplicate = TRUE)
seqSetFilter(f, variant.id = seqGetData(f, "variant.id")[samp], verbose = FALSE)
tab <- data.frame(chr = seqGetData(f, "$chromosome"),
    pos = seqGetData(f, "position"), ref = seqGetData(f, "$ref"),
    alt = seqGetData(f, "$alt"))
for (nm in varnm)
    tab[[nm]] <- seqGetData(f, paste0("annotation/info/", nm), .tolist = NA)
seqClose(f)
tmp <- tempfile(fileext = ".tsv")
write.table(tab, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
bytes_per_var <- file.size(tmp) / nrow(tab)
raw <- readBin(tmp, "raw", file.size(tmp))
gz_per_var <- length(memCompress(raw, "gzip")) / nrow(tab)
unlink(tmp, force = TRUE)
est_txt <- bytes_per_var * n_var
est_gz  <- gz_per_var * n_var
cat(sprintf("  flat text (est.)   : %10.2f MB  (%.1fx GDS)\n",
    est_txt / 1024^2, est_txt / sz_gds))
cat(sprintf("  flat text+gzip est : %10.2f MB  (%.1fx GDS)\n\n",
    est_gz / 1024^2, est_gz / sz_gds))

## ---- 2. random-access query latency --------------------------------------
cat("== Random-access query latency (full annotation set) ==\n")
for (nq in c(10L, 100L, 1000L, 10000L)) {
    if (nq > n_var) next
    k <- sample(n_var, nq)
    q <- paste(chr[k], pos[k], ref[k], alt[k], sep = "-")
    el <- system.time(seqAnnotate(q, gds, varnm = varnm,
        verbose = FALSE))[["elapsed"]]
    cat(sprintf("  %6d variants : %7.3f s  (%.3f ms/variant)\n",
        nq, el, 1000 * el / nq))
}
cat("\nDone.\n")
