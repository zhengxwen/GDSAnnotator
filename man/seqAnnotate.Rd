\name{seqAnnotate-methods}
\alias{seqAnnotate}
\alias{seqAnnotate,data.frame-method}
\alias{seqAnnotate,DataFrame-method}
\alias{seqAnnotate,SeqVarGDSClass-method}
\alias{seqAnnotate,character-method}
\alias{seqAnnotate,GRanges-method}
\alias{seqAnnotate,GRangesList-method}
\alias{seqAnnotate,IRanges-method}
\alias{seqAnnotate,ANY-method}
\alias{seqAnnotateGDS}
\alias{seqAnnotateVCF}
\title{Annotate a list of variants}
\description{
    Annotate a list of variants using the annotation stored in GDS files.
}

\usage{
\S4method{seqAnnotate}{character}(object, annot_gds, varnm, split="[-_:]",
    ..., verbose=TRUE)
\S4method{seqAnnotate}{data.frame}(object, annot_gds, varnm, col_chr="chr",
    col_pos="pos", col_ref="ref", col_alt="alt", ..., verbose=TRUE)
\S4method{seqAnnotate}{DataFrame}(object, annot_gds, varnm, col_chr="chr",
    col_pos="pos", col_ref="ref", col_alt="alt", ..., verbose=TRUE)
\S4method{seqAnnotate}{GRanges}(object, annot_gds, varnm, ..., verbose=TRUE)
\S4method{seqAnnotate}{GRangesList}(object, annot_gds, varnm, ...,
    verbose=TRUE)
\S4method{seqAnnotate}{IRanges}(object, annot_gds, varnm, chr, ...,
    verbose=TRUE)
\S4method{seqAnnotate}{ANY}(object, annot_gds, varnm, ..., verbose=TRUE)
\S4method{seqAnnotate}{SeqVarGDSClass}(object, annot_gds, varnm,
    add_to_gds=FALSE, no_sample=TRUE, root="", cleanup=FALSE,
    bsize=1000000L, parallel=FALSE, ..., verbose=TRUE)
seqAnnotateGDS(gds_fn, annot_gds, varnm, add_to_gds=FALSE, no_sample=TRUE,
    root="", cleanup=TRUE, bsize=1000000L, parallel=FALSE, ..., verbose=TRUE)
seqAnnotateVCF(vcf_fn, annot_gds, varnm, ..., verbose=TRUE)
}
\arguments{
    \item{object}{the variants to be annotated: a character vector of
        variant identifiers combining chromosome, position, reference and
        alternative alleles (e.g., "22-23413241-A-G", see \code{split}); a
        \code{data.frame} or \code{DataFrame} with the columns
        \code{col_chr}, \code{col_pos}, \code{col_ref} and \code{col_alt}; a
        \code{GRanges}, \code{GRangesList} or \code{IRanges} object, in which
        case all the variants located in the ranges are annotated; a
        \code{VCF} object of the package \pkg{VariantAnnotation}
        (\code{CollapsedVCF} or \code{ExpandedVCF}), annotated with one row
        per alternative allele by the method for \code{ANY}; or a SeqArray
        GDS object (\code{SeqVarGDSClass})}
    \item{gds_fn}{character for a SeqArray GDS file name}
    \item{vcf_fn}{character for a VCF or BCF file name, read by
        \code{\link[VariantAnnotation]{readVcf}}}
    \item{annot_gds}{a character vector for GDS file names, a
        \code{\link{SeqAnnotGDSList}} object returned by
        \code{\link{seqAnnotOpen}}, an opened GDS file or a list of opened GDS
        files (\code{SeqVarGDSClass} objects), e.g., one chromosome per file,
        while a file can also contain more than one chromosome; the files
        given by names are opened and closed in the call, and the opened
        files stay open; if a variant is found in more than one file, the
        annotation is taken from the first file}
    \item{varnm}{a character vector for GDS node names without the prefix
        "annotation/info/"; ":id" for accessing "annotation/id", ":qual" for
        "annotation/qual", ":filter" for "annotation/filter"; if it is missing,
        to use all annotations in the first file of \code{annot_gds}; if
        \code{varnm=character()} or \code{NULL}, return which variants are
        selected; see details}
    \item{col_chr}{the column name for chromosome}
    \item{col_pos}{the column name for position}
    \item{col_ref}{the column name for the reference alleles}
    \item{col_alt}{the column name for the alternative alleles}
    \item{chr}{a character vector of chromosome(s), used with an
        \code{IRanges} object; passed to \code{\link[SeqArray]{seqSetFilter}}}
    \item{add_to_gds}{logical or character for an output file name; if
        \code{FALSE}, return DataFrame for annotation; if \code{TRUE}, add
        annotation to the input GDS file; if \code{add_to_gds} is a file name,
        output the annotation to a new GDS file}
    \item{no_sample}{logical, when \code{add_to_gds} is a file name and
        \code{no_sample=TRUE}, the output GDS file excludes per-sample data}
    \item{root}{character, to save the annotated GDS nodes to a \code{root}
        node (using "annotation/info/root"); if \code{root=""}, use
        "annotation/info" instead}
    \item{cleanup}{if \code{TRUE}, call \code{cleanup.gds} on the output
        GDS file when \code{add_to_gds} is a file name}
    \item{bsize}{the block size, i.e., the number of variants whose
        annotation data are loaded and written to the output GDS file at a
        time, when \code{add_to_gds} is \code{TRUE} or a file name; a smaller
        block size reduces the memory usage of each process}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE}
        (\code{parallel::detectCores() - 1} processes), a positive integer for
        the number of processes, or a cluster object created by
        \code{\link[parallel]{makeCluster}}; only used when
        \code{add_to_gds} is \code{TRUE} or a file name, see details}
    \item{split}{regular expression to split the character vector in
        \code{object} by using \code{strsplit}, using "-", "_" or ":" by
        default, e.g., "22-23413241-A-G" will be split into "22", "23413241",
        "A" and "G"}
    \item{verbose}{if \code{TRUE}, show information}
    \item{...}{additional arguments passed to specific methods, e.g., from
        \code{seqAnnotateGDS()} to the \code{SeqVarGDSClass} method, and
        from the \code{VCF} method to the \code{DataFrame} method}
}
\value{
    Returns a \code{DataFrame} with one row per input variant (one row per
variant located in the ranges for a \code{GRanges}, \code{GRangesList} or
\code{IRanges} input, and one row per alternative allele for a \code{VCF}
input) and one column per annotation in \code{varnm}; the
annotations of a variant not found in \code{annot_gds} are \code{NA}. If
\code{varnm=character()} or \code{NULL}, the \code{DataFrame} has the two
columns 'file_idx' and 'variant_idx' instead, see details. When
\code{add_to_gds} is \code{TRUE} or a file name, the annotations are written
to the GDS file and \code{NULL} is returned invisibly.
}
\details{
    If \code{varnm=character()} or \code{NULL}, \code{seqAnnotate} returns
a DataFrame object with two columns 'file_idx' and 'variant_idx' to indicate
where the variants are found, and NA in 'variant_idx' if the variant is not
found.

    \code{seqAnnotateGDS()} opens the GDS file \code{gds_fn} and calls the
\code{SeqVarGDSClass} method. \code{seqAnnotateVCF()} reads the VCF file with
\code{\link[VariantAnnotation]{readVcf}} and annotates the \code{VCF} object.
The package \pkg{VariantAnnotation} should be installed for the VCF files and
objects: it is a suggested package loaded on demand, so the class \code{VCF}
is not used in a method signature, and a \code{VCF} object is handled by the
method for \code{ANY}, which checks the class of the object.

    A \code{VCF} object is annotated for each alternative allele of each row: a
\code{CollapsedVCF} with multi-allelic sites gives one row per alternative
allele, a row without any alternative allele (\code{ALT} is ".") is skipped,
so the output can have more or fewer rows than the \code{VCF} object, and
\code{\link[VariantAnnotation:VCF-class]{expand}} can be used first to get one
alternative allele per row. The output is the same for a \code{CollapsedVCF}
and for its expanded version.

    When the annotations are written to a GDS file (\code{add_to_gds} is
\code{TRUE} or a file name), the variants are first matched to the annotation
files, and then the annotation data are loaded and written block by block (see
\code{bsize}), so that the memory usage of writing does not depend on the
total number of variants. The genotype GDS file then holds the annotations of
its variants under "annotation/info" (or under "annotation/info/root"); with
the annotations of FAVOR, this produces an annotated GDS (aGDS) file, see
\code{\link{seqFAVOR2GDS}}. The values are written as they are stored in
the annotation files: for a variable with a variable number of values per
variant (e.g., the fields of VEP CSQ or SnpEff ANN), a variant without value
or not found in the annotation files has no value, instead of a missing value;
for a variable with one value per variant, a variant not found has a missing
value. A factor (e.g., \code{":filter"}) is written as a factor with the
levels in all the annotation files.

    With \code{parallel}, the annotation variables in \code{varnm} are
distributed to the processes, and the matching of the variants is not
parallelized. Each process writes one variable at a time to a temporary GDS
file in the directory of the output file, and the variables are copied to the
output file in the order of \code{varnm}, so the output file does not depend
on \code{parallel}. A number or \code{TRUE} starts new R processes with
\code{\link[parallel]{makeCluster}} instead of forking the current process,
since a job scheduler (e.g., SLURM) may count the memory shared by the forked
processes once per process. Each process loads the packages, opens the
annotation files it needs and holds one block of one variable, so the total
memory usage grows with the number of processes. The temporary files need as
much disk space as the new annotations, and there is no gain with a single
variable in \code{varnm}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAnnotOpen}}, \code{\link{seqAnnotList}},
    \code{\link{seqValueCounts}}, \code{\link{seqFAVOR2GDS}},
    \code{\link{seqVEP2GDS}}
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

# a subset of the annotations, from a data.frame
df <- data.frame(chr="22", pos=c(10510007L, 10510038L, 10510282L),
    ref=c("T", "T", "G"), alt=c("C", "C", "C"))
seqAnnotate(df, annot_gds, varnm=c("rsid", "cadd_phred", "linsight"))

# all the variants located in the genomic ranges
gr <- GenomicRanges::GRanges(seqnames = c("chr22", "chr22"),
    ranges = IRanges::IRanges(start = c(10510100, 10510225),
        end = c(10510125, 10510250)))
seqAnnotate(gr, annot_gds, varnm)

# the variants of a VCF file, one row per alternative allele (the package
#   VariantAnnotation should be installed)
vcf_fn <- system.file("extdata", "gnomad.genomes.v4.chr22_sub.vcf.gz",
    package="GDSAnnotator")
if (requireNamespace("VariantAnnotation", quietly=TRUE))
{
    res <- seqAnnotateVCF(vcf_fn, annot_gds, varnm=c("rsid", "cadd_phred"))
    print(head(res[!is.na(res$cadd_phred), ]))

    # or from a VCF object of the package VariantAnnotation
    vcf <- VariantAnnotation::readVcf(vcf_fn)
    res2 <- seqAnnotate(vcf, annot_gds, varnm=c("rsid", "cadd_phred"))
    print(identical(res, res2))
}

# the variants of a SeqArray GDS file (the gnomAD sites of the same region)
gds_fn <- system.file("extdata", "gnomad.genomes.v4.chr22_sub.gds",
    package="GDSAnnotator")
res <- seqAnnotateGDS(gds_fn, annot_gds, varnm=c("rsid", "cadd_phred"))
head(res[!is.na(res$cadd_phred), ])

# write the annotations to a new GDS file
out_fn <- tempfile(fileext=".gds")
seqAnnotateGDS(gds_fn, annot_gds, varnm=c("rsid", "cadd_phred"),
    add_to_gds=out_fn)
f <- seqOpen(out_fn)
head(seqGetData(f, "annotation/info/cadd_phred"))
seqClose(f)

# delete the temporary file
unlink(out_fn, force=TRUE)
}

\keyword{gds}
\keyword{annotation}
\keyword{WES}
\keyword{WGS}
