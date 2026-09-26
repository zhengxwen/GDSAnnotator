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

    Storing variant annotations in a GDS file was introduced by FAVOR with
the annotated GDS (aGDS) format, in which the FAVOR functional annotations of
the variants of a study are added to the SeqArray GDS file holding the
genotypes. This package generalizes the format to annotation GDS files of any
source, with or without genotypes: a whole-genome annotation database is kept
as one GDS file per chromosome and queried by position and alleles, and the
annotations of the variants of a study are written into its genotype GDS file
(with the FAVOR annotations, this produces an aGDS file).

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
\references{
    Zhou H, Arapoglou T, Li X, et al. FAVOR: functional annotation of
variants online resource and annotator for variation across the human
genome. Nucleic Acids Res. 2023;51(D1):D1300-D1311.
DOI: 10.1093/nar/gkac966.

    Zheng X, Gogarten SM, Lawrence M, et al. SeqArray: a storage-efficient
high-performance data format for WGS variant calls. Bioinformatics.
2017;33(15):2251-2257. DOI: 10.1093/bioinformatics/btx145.
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
