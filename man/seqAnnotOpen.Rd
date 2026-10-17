\name{seqAnnotOpen}
\alias{seqAnnotOpen}
\alias{SeqAnnotGDSList}
\alias{SeqAnnotGDSList-class}
\alias{seqClose,SeqAnnotGDSList-method}
\alias{show,SeqAnnotGDSList-method}
\title{Open annotation GDS files}
\description{
    Open one or more annotation GDS files once, so that they can be used in
more than one call of \code{\link{seqAnnotate}} without opening the files
again.
}

\usage{
seqAnnotOpen(gds_fn, allow.duplicate=FALSE)
\S4method{seqClose}{SeqAnnotGDSList}(object)
\S4method{show}{SeqAnnotGDSList}(object)
}
\arguments{
    \item{gds_fn}{a character vector for annotation GDS file names, e.g., one
        file per chromosome; a file can also contain more than one
        chromosome}
    \item{allow.duplicate}{if \code{TRUE}, a file can be opened more than
        once, see \code{\link[SeqArray]{seqOpen}}}
    \item{object}{a \code{SeqAnnotGDSList} object}
}
\value{
    \code{seqAnnotOpen()} returns a \code{SeqAnnotGDSList} object.
\code{seqClose()} closes all the files and returns \code{NULL} invisibly.
}
\details{
    A \code{SeqAnnotGDSList} object is a list of opened SeqArray GDS files
(\code{SeqVarGDSClass} objects), extending the class
\code{\link[S4Vectors:SimpleList-class]{SimpleList}}: \code{x[[i]]} is the
i-th file, and \code{x[i]} is a \code{SeqAnnotGDSList} object with a subset of
the files. \code{\link[S4Vectors:Vector-class]{mcols}(x)} has one row per
file with the columns 'filename' (the normalized file name), 'chromosome' (a
\code{\link[IRanges:AtomicList-class]{CharacterList}} of the chromosomes, as
they are stored in the file) and 'n_variant' (the number of variants). The
chromosomes and the number of variants are obtained from all the variants of
a file, whatever its variant filter is.

    The object can be passed to \code{\link{seqAnnotate}},
\code{\link{seqAnnotateGDS}} and \code{\link{seqAnnotateVCF}} as
\code{annot_gds}: the query variants are dispatched to the files of their
chromosomes (ignoring the "chr" prefix), and the files stay open after the
call. These functions change the variant filter of the files, so the filter
set by the user is not kept. The files are closed by \code{seqClose()}.

    The object holds the connections to the opened files, so it cannot be
saved (e.g., by \code{saveRDS()}) or sent to other R processes; with
\code{parallel} in \code{\link{seqAnnotateGDS}}, each process opens the files
again by their names.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAnnotate}}, \code{\link[SeqArray]{seqOpen}},
    \code{\link[SeqArray]{seqClose}}
}

\examples{
# FAVOR essential dataset of chromosome 22 (subset)
fn <- system.file("extdata", "favor_chr22_sub.gds", package="GDSAnnotator")

# open the file once
annot <- seqAnnotOpen(fn)
annot
S4Vectors::mcols(annot)

# annotate the variants in more than one call
seqAnnotate(c("22-10510007-T-C", "22-10510038-T-C"), annot,
    varnm=c("rsid", "cadd_phred"), verbose=FALSE)
df <- data.frame(chr="22", pos=10510282L, ref="G", alt="C")
seqAnnotate(df, annot, varnm=c("rsid", "cadd_phred"), verbose=FALSE)

# close the file
seqClose(annot)
}

\keyword{gds}
\keyword{annotation}
