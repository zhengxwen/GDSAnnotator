\name{seqToGDS_Nirvana}
\alias{seqToGDS_Nirvana}
\title{Reformat Illumina Nirvana JSON}
\description{
    Convert Illumina Nirvana JSON output to a SeqArray GDS file.
}
\usage{
seqToGDS_Nirvana(json_fn, out_fn, fields=NULL,
    compress=c("LZMA", "ZIP", "none"), root="", rm_chr_prefix=TRUE,
    float32=TRUE, sep=";", bsize=10000L, verbose=TRUE)
}
\arguments{
    \item{json_fn}{the file name of a Nirvana JSON output (a bgzip- or
        gzip-compressed \code{.json.gz} file)}
    \item{out_fn}{the file name of output SeqArray GDS file}
    \item{fields}{the fields to import. If \code{NULL} (the default), all
        fields flagged as default in the field catalogue
        "nirvana_output_format.csv" are used. Otherwise a character vector
        where each element is either a catalogue field \code{Name} (its
        \code{Path} and \code{Type} are taken from the catalogue) or a JSON
        dot-path not in the catalogue (its type is inferred from the data);
        the names of the vector, if any, override the GDS node names. See
        Details.}
    \item{compress}{compression method (LZMA by default), \code{"LZMA"} has a
        higher compression ratio than \code{"ZIP"}}
    \item{root}{character; if not \code{""}, the annotation nodes are written
        under "annotation/info/root" instead of "annotation/info"}
    \item{rm_chr_prefix}{if \code{TRUE}, remove the leading "chr" from the
        chromosome names (e.g. "chr22" becomes "22")}
    \item{float32}{if \code{TRUE}, store non-integer numeric fields as 32-bit
        floats to save space}
    \item{sep}{the separator used to join multiple values of a field for a
        single variant (e.g. multiple consequences)}
    \item{bsize}{the block size (number of variant positions per chunk) for
        streaming; a smaller value reduces peak memory usage}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    Nirvana writes a streaming JSON file: the first line holds the header,
each subsequent line holds one position object, followed by a gene section.
The importer reads the file block-by-block so that whole-genome outputs can
be converted with a constant memory footprint. Each alternate allele becomes
one variant record (\code{allele} stored as "ref,alt"). Provenance from the
Nirvana header (annotator, genome assembly, data version and data sources) is
stored as attributes of the "description" node. Fields that are entirely
missing in the input are dropped from the output.

    The set of importable fields and their data types is defined by the
catalogue file installed at \code{system.file("extdata",
"nirvana_output_format.csv", package="GDSAnnotator")}, with columns
\code{Name} (the GDS node name), \code{Path} (the dot-separated location in
the Nirvana JSON, where an array-of-objects extracts the key from each
element), \code{Type} (\code{integer}, \code{float}, \code{string} or
\code{logical}), \code{Default} (whether the field is imported when
\code{fields=NULL}) and \code{Description}. Editing this CSV is the
recommended way to change which fields are imported by default and how they
are typed. For a field whose type is not given by the catalogue (a custom
dot-path), the type is inferred from the first block of variants. For a field
with several values per variant, string fields are joined by \code{sep} while
numeric fields keep the first finite value.
}
\value{
    Return the file name of GDS format with an absolute path.
}
\references{
    Stromberg M, et al. Nirvana: Clinical Grade Variant Annotator.
ACM-BCB (2017). \url{https://github.com/Illumina/Nirvana}
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqToGDS_VEP}}, \code{\link{seqToGDS_gnomAD}},
    \code{\link{seqAnnotate}}
}

\examples{
# a small example of Nirvana JSON output
json_fn <- system.file("extdata", "example_nirvana.json.gz",
    package="GDSAnnotator")

out_fn <- tempfile(fileext=".gds")
seqToGDS_Nirvana(json_fn, out_fn)

# list the imported annotations
seqAnnotList(out_fn)

# delete the temporary file
unlink(out_fn, force=TRUE)
}

\keyword{gds}
\keyword{annotation}
\keyword{WGS}
\keyword{Nirvana}
