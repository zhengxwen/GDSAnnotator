\name{seqAnnotate-methods}
\alias{seqAnnotate}
\alias{seqAnnotate,data.frame-method}
\alias{seqAnnotate,DataFrame-method}
\alias{seqAnnotate,SeqVarGDSClass-method}
\alias{seqAnnotate,character-method}
\title{Annotate a list of variants}
\description{
    Annotate a list of variants using the annotation stored in GDS files.
}

\usage{
\S4method{seqAnnotate}{data.frame}(object, annot_gds, varnm, col_chr="chr",
    col_pos="pos", col_ref="ref", col_alt="alt", ..., verbose=TRUE)
\S4method{seqAnnotate}{DataFrame}(object, annot_gds, varnm, col_chr="chr",
    col_pos="pos", col_ref="ref", col_alt="alt", ..., verbose=TRUE)
}
\arguments{
    \item{object}{a data.frame}
    \item{annot_gds}{a character vector for GDS file names or a list of opened
        GDS files}
    \item{varnm}{a character vector for GDS node names without the prefix
        "annotation/info/"; ":id" for accessing "annotation/id", ":qual" for
        "annotation/qual", ":filter" for "annotation/filter"; if
        \code{varnm=character()} or \code{varnm=""}, see details}
    \item{col_chr}{the column name for chromosome}
    \item{col_pos}{the column name for position}
    \item{col_ref}{the column name for the reference alleles}
    \item{col_alt}{the column name for the alternative alleles}
    \item{verbose}{if \code{TRUE}, show information}
    \item{...}{additional arguments passed to specific methods}
}
\value{
    Returns a DataFrame or NULL.
}
\details{
    If \code{varnm=character()} or \code{varnm=""}, \code{seqAnnotate} returns
a DataFrame object with two columns 'file_idx' and 'variant_idx' to indicate
where the variants are found, and NA in 'variant_idx' if the variant is not
found.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqToGDS_VEP}}, \code{\link{seqToGDS_SnpEff}}
}

\examples{
# FAVOR essential dataset of chromosome 22 in CSV format
csv_fn <- system.file("extdata", "favor_chr22_sub.csv.gz",
    package="GDSAnnotator")

seqToGDS_FAVOR(csv_fn, "tmp.gds")

# delete the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{annotation}
\keyword{WGS}
