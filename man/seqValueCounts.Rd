\name{seqValueCounts}
\alias{seqValueCounts}
\title{Counts of unique values}
\description{
    Calculate the counts of unique values.
}
\usage{
seqValueCounts(gdsfile, varnm, use_info=TRUE, FUN=NULL, bsize=100000L,
    verbose=TRUE, ...)
}
\arguments{
    \item{gdsfile}{a \code{SeqVarGDSClass} object or a GDS file name}
    \item{varnm}{a character vector containing GDS node name(s)}
    \item{use_info}{\code{TRUE} or \code{FALSE}; if \code{TRUE}, to use the
        node(s) in the folder "annotation/info"; otherwise, full node name(s)
        should be specified in \code{varnm}}
    \item{FUN}{\code{NULL} to use the default \code{table}, or a user-defined
        function to be applied}
    \item{bsize}{block size}
    \item{verbose}{if \code{TRUE}, show information}
    \item{...}{additional arguments passed to \code{FUN} if it is a function}
}
\value{
    Returns a contingency table, an object of class \code{table}, an array of
integer values.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAnnotList}}
}

\examples{
# gnomAD v4.1.0 dataset of chromosome 22
gds_fn <- system.file("extdata", "gnomad.genomes.v4.chr22_sub.gds",
    package="GDSAnnotator")

# list the annotations in the gds file
seqAnnotList(gds_fn)

seqValueCounts(gds_fn, "variant_type")

seqValueCounts(gds_fn, c("allele_type", "variant_type"))
}

\keyword{gds}
\keyword{annotation}
