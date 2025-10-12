\name{seqToGDS_FAVOR}
\alias{seqToGDS_FAVOR}
\alias{seqToGDS_FAVOR_tar}
\title{Reformat FAVOR files}
\description{
    Convert FAVOR CSV files to a SeqArray GDS file.
}
\usage{
seqToGDS_FAVOR(csv_fn, out_fn, compress=c("LZMA", "ZIP", "none"),
    root="FAVOR", use_float32=TRUE, verbose=TRUE)
seqToGDS_FAVOR_tar(tar_fn, out_fn, fn_in_tar=NULL,
    compress=c("LZMA", "ZIP", "none"), root="", use_float32=TRUE,
    block_size=10000L, tar_cmd=Sys.getenv("TAR"), verbose=TRUE)
}
\arguments{
    \item{csv_fn}{the path(es) of FAVOR csv file}
    \item{tar_fn}{the path of a FAVOR tar file containing csv files}
    \item{out_fn}{the file name of output SeqArray GDS file}
    \item{compress}{compression method (LZMA by default), \code{"LZMA"} has
        a higher compression ratio than \code{"ZIP"}}
    \item{root}{specify the GDS node, the annotation will be stored in
        "annotation/info/root"; if \code{root=""}, "annotation/info" is used}
    \item{use_float32}{if \code{TRUE}, use 32-bit float to store floating-point
        numbers instead of 64-bit float}
    \item{fn_in_tar}{a character vector for the csv file names in the tar file}
    \item{block_size}{the block size (# of lines) for reading a csv file}
    \item{tar_cmd}{the path of TAR to the command to be used}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\references{
    Zhou H, Arapoglou T, ..., Lin X. FAVOR: Functional Annotation of Variants
Online Resource and Annotator for Variation across the Human Genome.
NAR 2022 Nov 9; gkac966. PMID: 36350676. DOI: 10.1093/nar/gkac966.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
# FAVOR essential dataset of chromosome 22 (subset) in CSV format
csv_fn <- system.file("extdata", "favor_chr22_sub.csv.gz",
    package="GDSAnnotator")

seqToGDS_FAVOR(csv_fn, "tmp.gds")

# delete the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{annotation}
\keyword{WGS}
\keyword{FAVOR}
