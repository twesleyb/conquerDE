# conquerDE
_conquering DE analysis_

#### Work in progress

* inspired by [csoneson/conquer_comparison](https://github.com/csoneson/conquer_comparison)
 and the work by [Soneson _et al._, 2018](https://pubmed.ncbi.nlm.nih.gov/29481549/)
* exploring and comparing methods for analysis of gene differential expression
	in RNAseq experiments

## Usage
```
devtools::install_github("twesleyb/conquerDE")

data(L) # list of counts, design, contrast, and meta

conquerDE(L, method="edgeRQLF")
conquerDE(L, method="DESeq2")
conquerDE(L, method="glmGamPoi")

```
