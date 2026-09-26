# Example data in `inst/extdata`

All the example files are small subsets of chromosome 22 (GRCh38), intended
only for the examples, the unit tests and the vignette. The `.gds` files are
generated from the text files next to them by the package's own conversion
functions, see `make_extdata.R` in this directory
(`Rscript inst/script/make_extdata.R` from the package root).

## Annotation resources

| File | Description |
|:-----|:------------|
| `favor_chr22_sub.csv.gz` | 1,000 consecutive variants at the start of chromosome 22 (positions 10,510,001-10,510,333) from the FAVOR essential database (https://favor.genohub.org, GRCh38), in the CSV format distributed by FAVOR; the columns "variant_vcf", "chromosome", "position", "ref_vcf" and "alt_vcf" identify the variants, and the other 30 columns are the FAVOR annotation scores (aPC scores, CADD, LINSIGHT, GeneHancer, rDHS, ...). |
| `favor_chr22_sub.gds` | the CSV file above converted with `seqToGDS_FAVOR(root="")`. |
| `favor_csv_header.csv` | the column names, types and descriptions of the FAVOR essential database CSV files, used by `seqToGDS_FAVOR_tar()` to read the CSV files with the correct column types; compiled from the FAVOR documentation. |
| `gnomad.genomes.v4.chr22_sub.vcf.gz` (+ `.csi`) | 150 variants (positions 10,510,033-10,510,756) from the gnomAD v4.1 genomes sites VCF of chromosome 22 (https://gnomad.broadinstitute.org/downloads, released under the CC0 license); the header and the INFO fields are unchanged. |
| `gnomad.genomes.v4.chr22_sub.gds` | the VCF file above converted with `seqToGDS_gnomAD()`. |

## Annotator outputs

The three annotated VCF files are derived from the same sites-only VCF of
chromosome 22 (`c22_b0_v1.sites.vcf.gz`, a set of imputed variant sites
written by Minimac4 v1.0.2; the file contains no sample and no genotype).
The variant sites were annotated with the following tools, and a subset of
the annotated variants was kept:

| File | Description |
|:-----|:------------|
| `example_wgs_sites_chr22_vep.vcf.gz` (+ `.csi`) | 250 variants annotated with Ensembl VEP v113 (cache GRCh38.p14, Ensembl 113, `--everything --flag_pick --nearest gene`, see the `##VEP` and `##VEP-command-line` header lines). |
| `example_wgs_sites_chr22_vep.gds` | converted with `seqToGDS_VEP()`: the CSQ field is split into the sub-fields under `annotation/info/CSQ.list`. |
| `example_wgs_sites_chr22_snpeff.vcf.gz` (+ `.csi`) | 200 variants annotated with SnpEff 5.2f (database GRCh38.p14, `-lof`, see the `##SnpEffVersion` and `##SnpEffCmd` header lines). |
| `example_wgs_sites_chr22_snpeff.gds` | converted with `seqToGDS_SnpEff()`: the ANN, LOF and NMD fields are split into the sub-fields under `annotation/info/ANN.list`, `LOF.list` and `NMD.list`. |
| `example_wgs_sites_chr22_annovar.vcf.gz` | 10 variants annotated with ANNOVAR (`table_annovar.pl -vcfinput`, with the refGeneWithVer, gnomAD exome and dbNSFP-based protocols, see the INFO header lines). |

## Format tables

The following CSV files are hand-curated tables used by the conversion and
summary functions:

| File | Description |
|:-----|:------------|
| `vep_output_format.csv` | the CSQ sub-fields of Ensembl VEP with their descriptions and R types, compiled from the VEP documentation (https://www.ensembl.org/info/docs/tools/vep/vep_formats.html); "Uniform" marks the fields that take the same value for all the transcripts of a variant, which are stored once per variant. |
| `snpeff_output_format.csv` | the ANN, LOF and NMD sub-fields of SnpEff, compiled from the SnpEff documentation (https://pcingola.github.io/SnpEff/snpeff/inputoutput/). |
| `annovar_output_format.csv` | the INFO fields written by ANNOVAR's `table_annovar.pl`, with descriptions and R types, compiled from the ANNOVAR documentation (https://annovar.openbioinformatics.org/). |
| `so_consequence.csv` | the Sequence Ontology consequence terms ranked by severity, with the impact category (HIGH, MODERATE, LOW, MODIFIER), the functional class (MISSENSE, NONSENSE, SILENT) and the genomic region, compiled from the Ensembl variation documentation (https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html) and the SnpEff impact classification. |
