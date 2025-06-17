\name{seqToGDS_SnpEff}
\alias{seqToGDS_SnpEff}
\title{Reformat SnpEff VCF}
\description{
    Convert SnpEff VCF output to a SeqArray GDS file.
}
\usage{
seqToGDS_SnpEff(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root=c("ANN", "LOF", "NMD"), verbose=TRUE)
}
\arguments{
    \item{vcf_fn}{the path(es) of VCF file from SnpEff output}
    \item{out_fn}{the file name of output SeqArray GDS file}
    \item{compress}{compression method (LZMA by default), \code{"LZMA"} has
        a higher compression ratio than \code{"ZIP"}}
    \item{root}{the name(s) of INFO field in the input VCF file to be parsed}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\references{
    Cingolani, P., Platts, A., Wang le L., Coon, M., Nguyen, T., Wang, L.,
Land, S.J., Lu, X., & Ruden, D.M. (2012). A program for annotating and
predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the
genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly, 6(2), 80-92.
PMID: 22728672.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
# SnpEff annotation
vcf_fn <- system.file("extdata", "example_wgs_sites_chr22_snpeff.vcf.gz",
    package="GDSAnnotator")

seqToGDS_SnpEff(vcf_fn, "tmp.gds")

# delete the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{annotation}
\keyword{WGS}
\keyword{SnpEff}
