\name{seqAnnotate-methods}
\alias{seqAnnotate}
\alias{seqAnnotate,data.frame-method}
\alias{seqAnnotate,DataFrame-method}
\alias{seqAnnotate,SeqVarGDSClass-method}
\alias{seqAnnotate,character-method}
\alias{seqAnnotate,GRanges-method}
\alias{seqAnnotate,GRangesList-method}
\alias{seqAnnotate,IRanges-method}
\alias{seqAnnotateFile}
\title{Annotate a list of variants}
\description{
    Annotate a list of variants using the annotation stored in GDS files.
}

\usage{
\S4method{seqAnnotate}{data.frame}(object, annot_gds, varnm, col_chr="chr",
    col_pos="pos", col_ref="ref", col_alt="alt", ..., verbose=TRUE)
\S4method{seqAnnotate}{DataFrame}(object, annot_gds, varnm, col_chr="chr",
    col_pos="pos", col_ref="ref", col_alt="alt", ..., verbose=TRUE)
\S4method{seqAnnotate}{SeqVarGDSClass}(object, annot_gds, varnm,
    add_to_gds=FALSE, no_sample=TRUE, root="", ..., verbose=TRUE)
\S4method{seqAnnotate}{character}(object, annot_gds, varnm, split="[-_:]", ...,
    verbose=TRUE)
\S4method{seqAnnotate}{GRanges}(object, annot_gds, varnm, ..., verbose=TRUE)
\S4method{seqAnnotate}{GRangesList}(object, annot_gds, varnm, ..., verbose=TRUE)
seqAnnotateFile(gds_fn, annot_gds, varnm, add_to_gds=FALSE, no_sample=TRUE,
    root="", ..., verbose=TRUE)
}
\arguments{
    \item{object}{a \code{data.frame}, \code{DataFrame}, a SeqArray GDS object
        \code{SeqVarGDSClass}, character vector; if \code{character}, it should
        be a combination of chromosome, position, reference and alternative
        alleles, e.g., "22-23413241-A-G"}
    \item{gds_fn}{character for a GDS file name}
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
    \item{add_to_gds}{logical or character for an output file name; if
        \code{FALSE}, return DataFrame for annotation; if \code{TRUE}, add
        annotation to the input GDS file; if \code{add_to_gds} is a file name,
        output the annotation to a new GDS file}
    \item{no_sample}{logical, when \code{add_to_gds} is a file name and
        \code{no_sample=TRUE}, the output GDS file excludes per-sample data}
    \item{root}{character, to save the annotated GDS nodes to a \code{root}
        node (using "annotation/info/root"); if \code{root=""}, use
        "annotation/info" instead}
    \item{split}{regular expression to split the character vector in
        \code{object} by using \code{strsplit}, using "-", "_" or ":" by
        default}
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
    \code{\link{seqAnnotList}}, \code{\link{seqValueCounts}}
}

\examples{
# FAVOR essential dataset of chromosome 22 (subset)
annot_gds <- system.file("extdata", "favor_chr22_sub.gds",
    package="GDSAnnotator")

# list the annotations in the gds file
(ann <- seqAnnotList(annot_gds))

varnm <- ann$name  # use all annotation
snp <- c("22-10510007-T-C", "22-10510038-T-C", "22-10510282-G-C",
    "22-10510303-G-C"  # a fake SNP, return NA
    )
seqAnnotate(snp, annot_gds, varnm)

gr <- GenomicRanges::GRanges(seqnames = c("chr22", "chr22"),
	ranges = IRanges::IRanges(start = c(10510100, 10510225),
	    end = c(10510125, 10510250)))
seqAnnotate(gr, annot_gds, varnm)
}

\keyword{gds}
\keyword{annotation}
\keyword{WES}
\keyword{WGS}
