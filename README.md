# conquerDE
_conquering DE analysis_

#### Work in progress

* inspired by [csoneson/conquer_comparison](https://github.com/csoneson/conquer_comparison)
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
