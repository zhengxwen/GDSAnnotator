GDSAnnotator: Whole-genome variant annotation data manipulation using GDS files
===

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Introduction

Variant annotation is a critical step in genomic data analysis, facilitating the interpretation of genetic variations and their potential biological significance. As large-scale whole-genome sequencing studies continue to expand, efficient storage, retrieval, and manipulation of variant annotation data have become increasingly important. The Genomic Data Structure ([GDS](https://www.bioconductor.org/packages/SeqArray/)) file format offers a scalable and storage-efficient solution for large genomic datasets, enabling rapid access and computational efficiency.

GDSAnnotator is an R package designed for fast and memory-efficient annotation of variants stored in GDS files. It integrates seamlessly into Bioconductor tools, and utilizes external data resources, e.g, Ensembl Variant Effect Predictor ([VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)) and FAVOR ([Functional Annotation of Variants Online Resource](https://favor.genohub.org)), to provide comprehensive functional annotations, including predicted variant effects and functional characteristics of non-coding variants. Compared to a plain text file, GDS format is 30-40 times more space-efficient when storing whole-genome annotations of FAVOR.


## Bioconductor

Coming soon.


## Installation

* Development version from Github (for developers/testers only):
```R
library("devtools")
install_github("zhengxwen/GDSAnnotator")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](https://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.



## Examples (1)

```R
library(GDSAnnotator)

# FAVOR essential dataset of chromosome 22 (subset)
annot_gds <- system.file("extdata", "favor_chr22_sub.gds", package="GDSAnnotator")

# list the annotations in the gds file (return a DataFrame)
(ann <- seqAnnotList(annot_gds))
## DataFrame with 30 rows and 4 columns
##                       name     type       trait description
##                <character> <factor> <character> <character>
## 1         apc_conservation     Real     Float32          NA
## 2      apc_conservation_v2     Real     Float32          NA
## 3          apc_epigenetics     Real     Float32          NA
## 4   apc_epigenetics_active     Real     Float32          NA
## 5   apc_epigenetics_repr..     Real     Float32          NA
## ...                    ...      ...         ...         ...
## 26  genecode_comprehensi..  Logical       Int32          NA
## 27              genehancer  Logical       Int32          NA
## 28                linsight  Real        Float32          NA
## 29              cadd_phred  Real        Float32          NA
## 30                    rdhs  Logical       Int32          NA

varnm <- ann$name  # use all annotation
snp <- c("22-10510007-T-C", "22-10510038-T-C", "22-10510282-G-C",
    "22-10510303-G-C"  # a fake SNP, return NA
    )  # a list of variants

seqAnnotate(snp, annot_gds, varnm)
## Open ‘favor_chr22_sub.gds’ [1,000 variants]
## [favor_chr22_sub.gds]  # of variants found: 3
## [30] ..............................
## DataFrame with 4 rows and 30 columns
##   apc_conservation apc_conservation_v2 apc_epigenetics apc_epigenetics_active  ...
##          <numeric>           <numeric>       <numeric>              <numeric>
## 1         2.583288            2.232373       0.0874126               0.226559
## 2         4.766045            4.498693       0.0750365               0.226559
## 3         0.393781            0.323713       0.1819110               0.226559
## 4               NA                  NA              NA                     NA
## ......
```

## Examples (2)

https://gds-stat.s3.amazonaws.com/download/favor/Example1_FAVOR_GDS.html


## Key functions in GDSAnnotator

| Function        | Description |
|:----------------|:-------------------------------------------|
| seqToGDS_FAVOR  | Convert FAVOR CSV files to GDS |
| seqToGDS_gnomAD | Convert gnomAD VCF files to GDS |
| seqToGDS_VEP    | Convert Ensembl VEP VCF output to GDS |
| seqToGDS_SnpEff | Convert SnpEff VCF output to GDS |
| seqAnnotate     | Annotate variants using the annotation stored in GDS |
| seqAnnotList    | List the annotations stored in the GDS file |
| seqValueCounts  | Calculate the counts of unique values |
| seqUnitGroupAnnot | Group the variants based on the variant annotations |


## See also

* [SeqArray](https://www.bioconductor.org/packages/SeqArray): Data management of large-scale whole-genome sequence variant calls using GDS files
