\name{seqToGDS_FAVOR}
\alias{seqToGDS_FAVOR}
\title{Reformat FAVOR files}
\description{
    Convert FAVOR CSV files to a SeqArray GDS file.
}
\usage{
seqToGDS_FAVOR(csv_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="FAVOR", use_float32=TRUE, verbose=TRUE)
}
\arguments{
    \item{csv_fn}{the path(es) of FAVOR csv file}
    \item{out_fn}{the file name of output SeqArray GDS file}
    \item{compress}{}
    \item{root}{}
    \item{use_float32}{}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\references{
    Zhou H, Arapoglou T, ..., Lin X. FAVOR: Functional Annotation of Variants Online
Resource and Annotator for Variation across the Human Genome . NAR 2022 Nov 9; gkac966.
PMID: 36350676. DOI: 10.1093/nar/gkac966.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
# delete the temporary file
unlink(c("tmp.gds", "tmp_p2.gds"), force=TRUE)
}

\keyword{gds}
\keyword{annotation}
\keyword{WGS}
