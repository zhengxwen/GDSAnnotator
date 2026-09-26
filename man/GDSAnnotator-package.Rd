\name{GDSAnnotator-package}
\alias{GDSAnnotator-package}
\alias{GDSAnnotator}
\docType{package}
\title{Whole-genome variant annotation data manipulation using GDS files}
\description{
    Fast and memory-efficient manipulation of whole-genome variant annotation
data stored in Genomic Data Structure (GDS) files.
}
\details{
    \tabular{ll}{
        Package: \tab GDSAnnotator \cr
        Type: \tab Package \cr
        License: \tab GPL-3 \cr
        Depends: \tab gdsfmt, SeqArray \cr
    }

    The package stores per-variant annotations in SeqArray GDS files, one
chromosome per file, under the node "annotation/info". The columnar, typed
and block-compressed layout is compact and supports random access, so a small
set of variants can be annotated without reading the whole database.

    \itemize{
        \item Import: \code{\link{seqToGDS_FAVOR}},
            \code{\link{seqToGDS_gnomAD}}, \code{\link{seqToGDS_VEP}},
            \code{\link{seqToGDS_SnpEff}} and \code{\link{seqToGDS_ANNOVAR}}
            convert the annotation resources and the outputs of the
            annotators to GDS files.
        \item Annotate: \code{\link{seqAnnotate}} annotates variants given as
            "chr-pos-ref-alt" strings, a \code{data.frame} or
            \code{DataFrame}, a \code{GRanges} object, a \code{VCF} object
            or file (\code{\link{seqAnnotateVCF}}) or a SeqArray GDS file
            (\code{\link{seqAnnotateGDS}}), and can write the annotations
            back into a GDS file block by block.
        \item Explore & summarize: \code{\link{seqAnnotList}} lists the
            annotations in a file, \code{\link{seqValueCounts}} tabulates
            their values, and \code{\link{seqAnnotStat}},
            \code{\link{seqAnnotReport}} and \code{\link{seqAnnotGeneTable}}
            produce VEP/SnpEff-like summary statistics and reports in a
            single pass over the file.
    }
}
\author{
    Xiuwen Zheng \email{xiuwen.zheng@abbvie.com}
}
\seealso{
    \code{\link[SeqArray]{seqOpen}}, \code{\link[SeqArray]{seqGetData}},
    \code{\link[SeqArray]{seqVCF2GDS}}
}

\examples{
# FAVOR essential dataset of chromosome 22 (subset)
annot_gds <- system.file("extdata", "favor_chr22_sub.gds",
    package="GDSAnnotator")

# the annotations in the file
seqAnnotList(annot_gds)

# annotate two variants
seqAnnotate(c("22-10510007-T-C", "22-10510038-T-C"), annot_gds,
    varnm=c("rsid", "cadd_phred", "linsight"))
}

\keyword{package}
\keyword{gds}
\keyword{annotation}
