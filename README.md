# `GDmat`: Additive and Dominance Relationship Matrices



# How to Install

To install this package, use devtools:

devtools::install_github("camult/GDmat")


## Description


 Calculates both realized additive and dominance relationship matrices.


## Usage

```r
GRM(X, smallValue)
DRM(X)
```


## Arguments

Argument      |Description
------------- |----------------
```X```     |     (numeric, $n$ ) matrix of genotypes with a number of rows identical to the number of genotyped individuals (n) and a number of columns identical to the number of SNPs (m). Values in the matrix are 0 and 2 for homozygous and 1 for heterozygous.
```smallValue```     |     add a small value to the diagonal to ensure the matrix will always be positive definite.

## Value


 Additive and Dominance Relationship Matrices.


## References


 VanRaden, Paul M. Efficient methods to compute genomic predictions. Journal of dairy science 91.11 (2008): 4414-4423.
 
 Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293


## Examples

```r 
 ## Not run:
 
 ## GRM(Markers, smallValue=0.01)
 
 ## RRM(Markers)
 
 ## End(Not run)
 ``` 

