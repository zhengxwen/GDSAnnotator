\name{seqAnnotGeneTable}
\alias{seqAnnotGeneTable}
\title{Per-gene annotation counts}
\description{
    Return the number of annotation records of each gene, broken down by
impact category.
}
\usage{
seqAnnotGeneTable(stat)
}
\arguments{
    \item{stat}{a \code{SeqAnnotStat} object returned from
        \code{\link{seqAnnotStat}}}
}
\details{
    The output is the equivalent of the "genes.txt" file of SnpEff. It is
\code{NULL} if the GDS file has no gene annotation, e.g., if the file was
not annotated by Ensembl-VEP or SnpEff.
}
\value{
    Returns a \code{data.frame} with one row per gene, sorted by the total
count in decreasing order: the column \code{gene} for the gene symbol,
\code{total} for the total number of annotation records, and one column for
each impact category (HIGH, MODERATE, LOW and MODIFIER).
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAnnotStat}}, \code{\link{seqAnnotReport}}
}

\examples{
# SnpEff-annotated variants of chromosome 22 (subset)
gds_fn <- system.file("extdata", "example_wgs_sites_chr22_snpeff.gds",
    package="GDSAnnotator")

st <- seqAnnotStat(gds_fn, verbose=FALSE)
head(seqAnnotGeneTable(st))
}

\keyword{gds}
\keyword{annotation}
