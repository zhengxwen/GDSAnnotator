\name{seqUnitGroupAnnot}
\alias{seqUnitGroupAnnot}
\title{Grouping the variants}
\description{
    Group the variants based on the variant annotations.
}
\usage{
seqUnitGroupAnnot(gdsfile, varnm, by=1L, cond=NULL, bsize=100000L,
    verbose=TRUE, ...)
}
\arguments{
    \item{gdsfile}{a \code{SeqVarGDSClass} object or a GDS file name}
    \item{varnm}{a character vector containing GDS node name(s)}
    \item{by}{}
    \item{cond}{}
    \item{FUN}{the function to be applied}
    \item{bsize}{block size}
    \item{verbose}{if \code{TRUE}, show information}
    \item{...}{additional arguments passed to the user-defined function(s) in
        \code{cond}}
}
\value{
    Returns a contingency table, an object of class \code{table}, an array of
integer values.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAnnotList}}, \code{\link{seqAnnotate}}
}

\examples{
library(GDSAnnotator)

# FAVOR essential dataset of chromosome 22 (subset)
annot_gds <- system.file("extdata", "favor_chr22_sub.gds",
    package="GDSAnnotator")

# list the annotations in the gds file
(ann <- seqAnnotList(annot_gds))

# seqUnitGroupAnnot(annot_gds, "")
}

\keyword{gds}
\keyword{annotation}
\keyword{WGS}
