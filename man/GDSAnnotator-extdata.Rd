\name{GDSAnnotator-extdata}
\alias{GDSAnnotator-extdata}
\docType{data}
\title{Example data files}
\description{
    Small example files in the \code{extdata} folder of the package, used in
the examples, the unit tests and the vignette. All of them are subsets of
chromosome 22 (GRCh38).
}
\format{
    \describe{
        \item{\code{favor_chr22_sub.csv.gz}}{1,000 consecutive variants at
            the start of chromosome 22 (positions 10,510,001-10,510,333) from
            the FAVOR essential database, in the CSV format distributed by
            FAVOR: the columns "variant_vcf", "chromosome", "position",
            "ref_vcf" and "alt_vcf" identify the variants, and the other 30
            columns are FAVOR annotation scores (aPC scores, CADD, LINSIGHT,
            GeneHancer, rDHS, ...).}
        \item{\code{favor_chr22_sub.gds}}{the CSV file above converted with
            \code{\link{seqFAVOR2GDS}(root="")}.}
        \item{\code{favor_csv_header.csv}}{the column names, types and
            descriptions of the FAVOR essential database CSV files, used by
            \code{\link{seqFAVORTar2GDS}}; compiled from the FAVOR
            documentation.}
        \item{\code{gnomad.genomes.v4.chr22_sub.vcf.gz} (and \code{.csi})}{
            150 variants (positions 10,510,033-10,510,756) from the gnomAD
            v4.1 genomes sites VCF file of chromosome 22, with the header and
            the INFO fields unchanged.}
        \item{\code{gnomad.genomes.v4.chr22_sub.gds}}{the VCF file above
            converted with \code{\link{seqGnomAD2GDS}}.}
        \item{\code{example_wgs_sites_chr22_vep.vcf.gz} (and \code{.csi})}{
            250 variant sites annotated with Ensembl VEP v113 (cache
            GRCh38.p14, \code{--everything --flag_pick --nearest gene}, see
            the "##VEP" and "##VEP-command-line" header lines).}
        \item{\code{example_wgs_sites_chr22_vep.gds}}{the VCF file above
            converted with \code{\link{seqVEP2GDS}}, the CSQ field being split
            into the sub-fields under "annotation/info/CSQ.list".}
        \item{\code{example_wgs_sites_chr22_snpeff.vcf.gz} (and
            \code{.csi})}{200 variant sites annotated with SnpEff 5.2f
            (database GRCh38.p14, \code{-lof}, see the "##SnpEffVersion" and
            "##SnpEffCmd" header lines).}
        \item{\code{example_wgs_sites_chr22_snpeff.gds}}{the VCF file above
            converted with \code{\link{seqSnpEff2GDS}}, the ANN, LOF and NMD
            fields being split into the sub-fields under
            "annotation/info/ANN.list", "LOF.list" and "NMD.list".}
        \item{\code{example_wgs_sites_chr22_annovar.vcf.gz}}{10 variant sites
            annotated with ANNOVAR (\code{table_annovar.pl -vcfinput}, with
            the refGeneWithVer, gnomAD exome and dbNSFP-based protocols, see
            the INFO header lines).}
        \item{\code{vep_output_format.csv}, \code{snpeff_output_format.csv},
            \code{annovar_output_format.csv}}{the sub-fields or INFO fields
            written by Ensembl VEP, SnpEff and ANNOVAR, with their
            descriptions and R types, compiled from the documentation of the
            three tools; used by the conversion functions.}
        \item{\code{so_consequence.csv}}{the Sequence Ontology consequence
            terms ranked by severity, with the impact category, the
            functional class and the genomic region, compiled from the
            Ensembl variation documentation and the SnpEff documentation;
            used by \code{\link{seqAnnotStat}}.}
    }
}
\details{
    The VEP, SnpEff and ANNOVAR files are derived from the same sites-only
VCF file of chromosome 22 (a set of imputed variant sites written by Minimac4
v1.0.2, with no sample and no genotype): the variant sites were annotated by
each tool, and a subset of the annotated variants was kept.

    The GDS files are generated from the text files by the script
\code{make_extdata.R} in the \code{scripts} folder of the installed package,
where \code{README.md} describes the provenance of each file in more detail.
}
\source{
    FAVOR: \url{https://favor.genohub.org}.

    gnomAD v4.1: \url{https://gnomad.broadinstitute.org/downloads}, released
under the CC0 1.0 license.

    Ensembl VEP: \url{https://www.ensembl.org/info/docs/tools/vep/};
    SnpEff: \url{https://pcingola.github.io/SnpEff/};
    ANNOVAR: \url{https://annovar.openbioinformatics.org/}.
}
\references{
    Zhou H, Arapoglou T, Li X, et al. FAVOR: functional annotation of
variants online resource and annotator for variation across the human
genome. Nucleic Acids Res. 2023;51(D1):D1300-D1311.
DOI: 10.1093/nar/gkac966.
}
\seealso{
    \code{\link{seqFAVOR2GDS}}, \code{\link{seqGnomAD2GDS}},
    \code{\link{seqVEP2GDS}}, \code{\link{seqSnpEff2GDS}},
    \code{\link{seqANNOVAR2GDS}}
}

\examples{
# the example files
list.files(system.file("extdata", package="GDSAnnotator"))

# the annotations in the FAVOR example file
seqAnnotList(system.file("extdata", "favor_chr22_sub.gds",
    package="GDSAnnotator"))

# the script generating the GDS files
cat(readLines(system.file("scripts", "make_extdata.R",
    package="GDSAnnotator")), sep="\n")
}

\keyword{datasets}
