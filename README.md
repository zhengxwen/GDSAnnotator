GDSAnnotator: Whole-genome variant annotation data manipulation using GDS files
===

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Introduction

Variant annotation is a critical step in genomic data analysis, facilitating the interpretation of genetic variations and their potential biological significance. As large-scale whole-genome sequencing studies continue to expand, efficient storage, retrieval, and manipulation of variant annotation data have become increasingly important. The Genomic Data Structure ([GDS](https://www.bioconductor.org/packages/SeqArray/)) file format offers a scalable and storage-efficient solution for large genomic datasets, enabling rapid access and computational efficiency.

GDSAnnotator is an R package designed for fast and memory-efficient annotation of variants stored in GDS files. It integrates seamlessly into Bioconductor tools, and utilizes external data resources, e.g, Ensembl Variant Effect Predictor ([VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)) and FAVOR ([Functional Annotation of Variants Online Resource](https://favor.genohub.org)), to provide comprehensive functional annotations, including predicted variant effects and functional characteristics of non-coding variants. Compared to a plain text file, GDS format is 30-40 times more space-efficient when storing whole-genome annotations of FAVOR.


## Bioconductor:

Release version: -

[http://www.bioconductor.org/packages/GDSAnnotator](http://www.bioconductor.org/packages/GDSAnnotator)


## Installation

* Bioconductor repository:
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("GDSAnnotator")
```



## Examples

```R
library(GDSAnnotator)

```


## Key functions in GDSAnnotator

| Function        | Description |
|:----------------|:-------------------------------------------|
| seqToGDS_FAVOR  | Convert FAVOR CSV files to GDS |
| seqToGDS_gnomAD | Convert gnomAD VCF files to GDS |
| seqToGDS_VEP    | Convert Ensembl VEP VCF output to GDS |
| seqToGDS_SnpEff | Convert SnpEff VCF output to GDS |
| seqAnnotate     | Annotate variants using the annotation stored in GDS |


## See also

* [SeqArray](https://www.bioconductor.org/packages/SeqArray): Data management of large-scale whole-genome sequence variant calls using GDS files
