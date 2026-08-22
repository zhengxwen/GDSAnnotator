\name{seqAnnotStat}
\alias{seqAnnotStat}
\alias{print.SeqAnnotStat}
\title{Summary statistics of variant annotation}
\description{
    Calculate the summary statistics of the variants and their annotations
stored in a SeqArray GDS file, in a single pass over the file.
}
\usage{
seqAnnotStat(gdsfile, parallel=FALSE, bsize=100000L, verbose=TRUE)
\method{print}{SeqAnnotStat}(x, ...)
}
\arguments{
    \item{gdsfile}{a \code{SeqVarGDSClass} object, a GDS file name or a
        character vector for multiple file names}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details}
    \item{bsize}{block size, i.e., the number of variants to be processed
        at a time}
    \item{verbose}{if \code{TRUE}, show information}
    \item{x}{a \code{SeqAnnotStat} object}
    \item{...}{unused arguments}
}
\details{
    All statistics are collected in a single call to
\code{\link{seqBlockApply}}, so that the file is scanned only once no matter
how many statistics are reported. Each statistic is stored as a named vector
of counts, and the counts are additive: when \code{gdsfile} is a character
vector of multiple file names (e.g., one GDS file per chromosome), the files
are summarized one by one and the counters are merged into a single object.
The memory usage is proportional to the number of distinct categories (e.g.,
the number of genes), not to the number of variants.

    Since the counters are additive, \code{parallel} splits the variants
between the processes, each process calculates its own counters, and the
results are merged. The output does not depend on \code{parallel} or on
\code{bsize}.

    The variant-level statistics (the number of variants per chromosome, the
variant classes, the base changes and the Ts/Tv ratio, the InDel length
distribution, the number of alleles, known vs. novel, FILTER and QUAL) are
derived from the nodes "chromosome", "position", "allele", "annotation/id",
"annotation/qual" and "annotation/filter".

    The annotation statistics are calculated when the GDS file contains the
sub-fields of Ensembl-VEP ("annotation/info/CSQ.list", see
\code{\link{seqToGDS_VEP}}) or SnpEff ("annotation/info/ANN.list", see
\code{\link{seqToGDS_SnpEff}}). Both annotators report Sequence Ontology
terms, and an internal table ("so_consequence.csv" in the "extdata" folder)
is used to rank the consequences by severity, to classify them into impact
categories (HIGH, MODERATE, LOW and MODIFIER), into functional classes
(MISSENSE, NONSENSE and SILENT) and into genomic regions. The consequences
are counted twice: once for every annotation record (e.g., one record per
transcript), and once per variant using only the most severe consequence.
}
\value{
    Returns an object of class \code{SeqAnnotStat}, a list with the
components: \code{file}, \code{filesize}, \code{n_variant}, \code{source}
(the annotator), \code{profile} (the annotation sub-fields in use),
\code{header} (the VCF header), \code{counts} (a list of named count
vectors) and \code{created}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAnnotReport}}, \code{\link{seqAnnotGeneTable}},
    \code{\link{seqValueCounts}}, \code{\link{seqAnnotList}}
}

\examples{
# SnpEff-annotated variants of chromosome 22 (subset)
gds_fn <- system.file("extdata", "example_wgs_sites_chr22_snpeff.gds",
    package="GDSAnnotator")

(st <- seqAnnotStat(gds_fn, verbose=FALSE))

# the most severe consequence of each variant
st$counts$cons_severe

# the impact of every annotation record
st$counts$impact

# transitions & transversions
st$counts$base_change

# using two processes
st2 <- seqAnnotStat(gds_fn, parallel=2L, verbose=FALSE)
identical(st$counts, st2$counts)
}

\keyword{gds}
\keyword{annotation}
