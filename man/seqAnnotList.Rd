\name{seqAnnotList}
\alias{seqAnnotList}
\title{List annotations}
\description{
    List the annotations stored in the GDS file.
}
\usage{
seqAnnotList(gdsfile)
}
\arguments{
    \item{gdsfile}{a \code{SeqVarGDSClass} object or a GDS file name}
}
\value{
    Returns a \code{DataFrame} with the columns 'name', 'type', 'trait' and
'description'.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAnnotate}}
}

\examples{
# gnomAD v4.1.0 dataset of chromosome 22
gds_fn <- system.file("extdata", "gnomad.genomes.v4.chr22_sub.gds",
    package="GDSAnnotator")

# list the annotations in the gds file
seqAnnotList(gds_fn)
}

\keyword{gds}
\keyword{annotation}
