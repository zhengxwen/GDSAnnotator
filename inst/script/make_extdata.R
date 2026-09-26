#######################################################################
#
# Regenerate the GDS example files in 'inst/extdata' from the text files
# shipped with the package. See 'inst/script/README.md' for the source
# and the provenance of each file.
#
# Usage (from the package root directory):
#     Rscript inst/script/make_extdata.R
#

suppressPackageStartupMessages(library(GDSAnnotator))

extdata <- "inst/extdata"
stopifnot(dir.exists(extdata))
fn <- function(s) file.path(extdata, s)

# FAVOR essential database, chromosome 22 (subset of 1,000 variants)
#   the annotations are stored directly under 'annotation/info' (root="")
seqToGDS_FAVOR(fn("favor_chr22_sub.csv.gz"), fn("favor_chr22_sub.gds"),
    root="", verbose=FALSE)

# gnomAD v4.1 genomes, chromosome 22 (subset of 150 variants)
seqToGDS_gnomAD(fn("gnomad.genomes.v4.chr22_sub.vcf.gz"),
    fn("gnomad.genomes.v4.chr22_sub.gds"), verbose=FALSE)

# Ensembl-VEP output (subset of 250 variants), the CSQ field is split into
#   the sub-fields under 'annotation/info/CSQ.list'
seqToGDS_VEP(fn("example_wgs_sites_chr22_vep.vcf.gz"),
    fn("example_wgs_sites_chr22_vep.gds"), verbose=FALSE)

# SnpEff output (subset of 200 variants), the ANN, LOF and NMD fields are
#   split into the sub-fields under 'annotation/info/{ANN,LOF,NMD}.list'
seqToGDS_SnpEff(fn("example_wgs_sites_chr22_snpeff.vcf.gz"),
    fn("example_wgs_sites_chr22_snpeff.gds"), verbose=FALSE)

for (s in c("favor_chr22_sub.gds", "gnomad.genomes.v4.chr22_sub.gds",
    "example_wgs_sites_chr22_vep.gds", "example_wgs_sites_chr22_snpeff.gds"))
{
    cat(s, ":", nrow(seqAnnotList(fn(s))), "annotations\n")
}
