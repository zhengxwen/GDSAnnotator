\name{seqAnnotReport}
\alias{seqAnnotReport}
\title{Variant annotation summary report}
\description{
    Generate a summary report of the variants and their annotations stored
in a SeqArray GDS file, in HTML or Markdown format, similar to the reports
of Ensembl-VEP and SnpEff.
}
\usage{
seqAnnotReport(stat, out_fn="annot_summary.html",
    format=c("auto", "html", "markdown", "rmd"),
    title="GDSAnnotator Summary Report", parallel=FALSE, bsize=100000L,
    verbose=TRUE)
}
\arguments{
    \item{stat}{a \code{SeqAnnotStat} object returned from
        \code{\link{seqAnnotStat}}, or a \code{SeqVarGDSClass} object, a GDS
        file name or a character vector for multiple file names, for which
        \code{\link{seqAnnotStat}} is called first}
    \item{out_fn}{the file name of the output report}
    \item{format}{the output format: \code{"html"} and \code{"markdown"}
        are finished reports, ready to be read as they are, and their charts
        are drawn by the renderer (inline SVG in HTML, block characters in
        Markdown); \code{"rmd"} is instead the R Markdown \emph{source} of a
        report, in which the counts are embedded as R code and the charts are
        drawn by \code{barplot()} only when the document is knitted, so it
        can be edited before rendering; \code{"auto"} (by default) takes the
        format from the file extension of \code{out_fn}: ".md" and
        ".markdown" give \code{"markdown"}, ".Rmd" gives \code{"rmd"}, and
        anything else gives \code{"html"} (the extension is matched
        case-insensitively)}
    \item{title}{the title of the report}
    \item{parallel}{the number of processes, passed to
        \code{\link{seqAnnotStat}} if \code{stat} is not a
        \code{SeqAnnotStat} object; see \code{\link{seqAnnotStat}}}
    \item{bsize}{block size, passed to \code{\link{seqAnnotStat}} if
        \code{stat} is not a \code{SeqAnnotStat} object}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    The three output formats show the same content: the report is first built
as a format-independent list of sections, then rendered. None of them writes
a companion image directory, i.e., the report is always a single file.

    The HTML output is self-contained, using plain HTML tables and inline SVG
charts, so it requires no JavaScript, no pandoc and no additional package,
and it can be viewed offline.

    The Markdown output uses GitHub-flavored pipe tables, with the bar charts
drawn using block characters, so that it is also readable as plain text.

    The R Markdown output embeds the counts in the document as R code, and
the charts are drawn by \code{barplot()} when the document is knitted, e.g.,
by \code{rmarkdown::render()}. Use this format to restyle the charts, to add
text of your own, or to include the report in a larger document. The code
chunks are hidden by \code{echo=FALSE} in the setup chunk.

    The report contains a general information section (the annotator version
and command line taken from the VCF header, the number of variants,
annotation records, genes and transcripts, and the variant rate), the
variant-level sections (variants per chromosome, variant classes, base
changes with the Ts/Tv ratio, InDel length distribution, alleles per variant,
known vs. novel, FILTER and QUAL) and, if the file was annotated by
Ensembl-VEP or SnpEff, the annotation sections (consequences for all records
and for the most severe consequence per variant, region, impact, functional
class with the Missense/Silent ratio, transcript biotypes, feature types,
SIFT and PolyPhen predictions, the warnings reported by the annotator, and
the most frequently annotated genes).
}
\value{
    Return the file name of the report with an absolute path.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAnnotStat}}, \code{\link{seqAnnotGeneTable}}
}

\examples{
# Ensembl-VEP annotated variants of chromosome 22 (subset)
gds_fn <- system.file("extdata", "example_wgs_sites_chr22_vep.gds",
    package="GDSAnnotator")

# HTML report
html_fn <- tempfile(fileext=".html")
seqAnnotReport(gds_fn, html_fn, title="Ensembl-VEP, chromosome 22",
    verbose=FALSE)

# Markdown report, reusing the statistics
st <- seqAnnotStat(gds_fn, verbose=FALSE)
md_fn <- tempfile(fileext=".md")
seqAnnotReport(st, md_fn, title="Ensembl-VEP, chromosome 22", verbose=FALSE)

cat(head(readLines(md_fn), 20), sep="\n")

# R Markdown source, to be knitted by rmarkdown::render()
rmd_fn <- tempfile(fileext=".Rmd")
seqAnnotReport(st, rmd_fn, title="Ensembl-VEP, chromosome 22", verbose=FALSE)

# remove the temporary files
unlink(c(html_fn, md_fn, rmd_fn), force=TRUE)
}

\keyword{gds}
\keyword{annotation}
