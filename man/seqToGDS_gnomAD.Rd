\name{seqToGDS_gnomAD}
\alias{seqToGDS_gnomAD}
\title{Reformat gnomAD VCF}
\description{
    Convert gnomAD VCF files to a SeqArray GDS file.
}
\usage{
seqToGDS_gnomAD(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="gnomAD", verbose=TRUE)
}
\arguments{
    \item{vcf_fn}{the path(es) of gnomAD VCF file}
    \item{out_fn}{the file name of output SeqArray GDS file}
    \item{compress}{}
    \item{root}{}
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
\keyword{gnomAD}
