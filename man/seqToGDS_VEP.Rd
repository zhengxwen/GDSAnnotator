\name{seqToGDS_VEP}
\alias{seqToGDS_VEP}
\title{Reformat Ensembl-VEP VCF}
\description{
    Convert Ensembl-VEP VCF output to a SeqArray GDS file.
}
\usage{
seqToGDS_VEP(vcf_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="CSQ", verbose=TRUE)
}
\arguments{
    \item{vcf_fn}{the path(es) of VCF file from Ensembl-VEP output}
    \item{out_fn}{the file name of output SeqArray GDS file}
    \item{compress}{compression method (LZMA by default), \code{"LZMA"} has
        a higher compression ratio than \code{"ZIP"}}
    \item{root}{the name of INFO field in the input VCF file, the default is
        \code{"CSQ"}}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\references{
    McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F.
The Ensembl Variant Effect Predictor. Genome Biology Jun 6;17(1):122. (2016)
DOI:10.1186/s13059-016-0974-4.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
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
\keyword{WGS}
\keyword{VEP}
