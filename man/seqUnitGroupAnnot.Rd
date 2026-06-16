\name{seqUnitGroupAnnot}
\alias{seqUnitGroupAnnot}
\title{Grouping the variants}
\description{
    Group the variants based on the variant annotations.
}
\usage{
seqUnitGroupAnnot(gdsfile, varnm, by=1L, cond=NULL, use_info=TRUE,
    bsize=100000L, verbose=TRUE, ...)
}
\arguments{
    \item{gdsfile}{a \code{SeqVarGDSClass} object or a GDS file name}
    \item{varnm}{a character vector of annotation field name(s); when
        \code{use_info=TRUE} (the default) the prefix "annotation/info/" is
        added automatically, so e.g. \code{"cadd_phred"} refers to
        "annotation/info/cadd_phred"; set \code{use_info=FALSE} to pass full
        GDS node paths instead}
    \item{use_info}{logical; if \code{TRUE}, prepend "annotation/info/" to each
        name in \code{varnm} (consistent with \code{\link{seqValueCounts}})}
    \item{by}{a numeric or character vector specifying which annotation
        variable(s) to group by; if numeric, it indicates the index/indices
        in \code{varnm}; if character, it should be variable name(s) from
        \code{varnm}; default is \code{1L} (group by the first variable)}
    \item{cond}{an optional list with the same length as \code{varnm} for
        filtering annotation values before grouping; each element can be
        \code{NULL} (no filtering), a vector (replace the values by whether
        they match the vector, i.e. \code{value \%in\% vector}), or a function
        applied to the value vector to transform it (e.g. binning a numeric
        score into categories)}
    \item{bsize}{block size}
    \item{verbose}{if \code{TRUE}, show information}
    \item{...}{additional arguments passed to the user-defined function(s) in
        \code{cond}}
}
\value{
    Returns an object of class \code{SeqUnitListClass}, a list with two
components: \code{desp}, a \code{data.frame} describing each group (one column
per grouping variable, plus \code{chr}, \code{start} and \code{end} giving the
chromosome and the position range spanned by the group); and \code{index}, a
list of integer vectors, each holding the variant indices belonging to that
group. The result can be passed to SeqArray aggregation functions such as
\code{seqUnitApply}.
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

# group the variants by binning the CADD-PHRED score via a 'cond' function
# ('cadd_phred' resolves to 'annotation/info/cadd_phred' as use_info=TRUE)
unit <- seqUnitGroupAnnot(annot_gds, varnm="cadd_phred",
    cond=list(function(x) ifelse(x >= 10, "damaging", "tolerated")))
unit$desp                # one row per group
lengths(unit$index)      # number of variants per group
}

\keyword{gds}
\keyword{annotation}
\keyword{WGS}
