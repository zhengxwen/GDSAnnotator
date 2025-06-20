\name{seqValueCounts}
\alias{seqValueCounts}
\title{Counts of unique values}
\description{
    Calculate the counts of unique values.
}
\usage{
seqValueCounts(gdsfile, varnm, FUN=NULL, bsize=100000L, verbose=TRUE, ...)
}
\arguments{
    \item{gdsfile}{a \code{SeqVarGDSClass} object or a GDS file name}
    \item{varnm}{a character vector containing GDS node name(s)}
    \item{FUN}{the function to be applied}
    \item{bsize}{block size}
    \item{verbose}{if \code{TRUE}, show information}
    \item{...}{additional arguments passed to specific methods}
}
\value{
    Returns a contingency table, an object of class \code{table}, an array of
integer values.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqOpen}}
}

\examples{
# Ensembl-VEP annotation
vcf_fn <- system.file("extdata", "example_wgs_sites_chr22_vep.vcf.gz",
    package="GDSAnnotator")

seqToGDS_VEP(vcf_fn, "tmp.gds")

# delete the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{annotation}
