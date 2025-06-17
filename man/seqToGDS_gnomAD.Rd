\name{seqToGDS_gnomAD}
\alias{seqToGDS_gnomAD}
\title{Reformat gnomAD VCF}
\description{
    Convert gnomAD VCF files to a SeqArray GDS file.
}
\usage{
seqToGDS_gnomAD(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    verbose=TRUE)
}
\arguments{
    \item{vcf_fn}{the path(es) of gnomAD VCF file}
    \item{out_fn}{the file name of output SeqArray GDS file}
    \item{compress}{compression method (LZMA by default), \code{"LZMA"} has
        a higher compression ratio than \code{"ZIP"}}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\references{
    Chen S, Francioli LC, Goodrich, JK, Collins, RL, etc.,
Nature. 625, 92-100 (2024). doi.org/10.1038/s41586-023-06045-0.
PMID: 38057664
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
# gnomAD v4.1.0 dataset of chromosome 22
vcf_fn <- system.file("extdata", "gnomad.genomes.v4.chr22_sub.vcf.gz",
    package="GDSAnnotator")

seqToGDS_gnomAD(vcf_fn, "tmp.gds")

# delete the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{annotation}
\keyword{WGS}
\keyword{gnomAD}
