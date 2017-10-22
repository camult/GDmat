# `GRM`: Calculates the realized additive relationship matrix

## Description


 Calculates realized additive (genomic) relationship matrix.


## Usage

```r
GRM(X, smallValue)
```


## Arguments

Argument      |Description
------------- |----------------
```X```     |     is a (numeric, $n$ ) matrix of genotypes with a number of rows identical to the number of SNPs (m) and a number of columns identical to the number of genotyped individuals (n). Values in the matrix are 0 and 2 for homozygous and 1 for heterozygous.
```smallValue```     |     is a numeric valeu that can be added to the diagonal of the genomic relationship matrix (GRM).

## References


 VanRaden, Paul M. Efficient methods to compute genomic predictions. Journal of dairy science 91.11 (2008): 4414-4423.
 
 Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293


## Examples

```r 
 ## Not to run ##
 
 ## GRM(X, smallValue)
 
 ## End(Not run)
 
 ``` 

# `DRM`: Calculates the realized dominance relationship matrix

## Description


 Calculates realized dominance relationship matrix.


## Usage

```r
DRM(X)
```


## Arguments

Argument      |Description
------------- |----------------
```X```     |     is a (numeric, $n$ ) matrix of genotypes with a number of rows identical to the number of SNPs (m) and a number of columns identical to the number of genotyped individuals (n). Values in the matrix are 0 and 2 for homozygous and 1 for heterozygous.

## References


 VanRaden, Paul M. Efficient methods to compute genomic predictions. Journal of dairy science 91.11 (2008): 4414-4423.
 
 Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293


## Examples

```r 
 ## Not to run ##
 
 ## DRM(X)
 
 ## End(Not run)
 
 ``` 

# `NRM`: Calculates the Pedigree-based Additive Relationship Matrix

## Description


 Calculates the the Pedigree-based Additive Relationship Matrix. This is twice the pedigree based kinship matrix.


## Usage

```r
NRM(Pedig)
```


## Arguments

Argument      |Description
------------- |----------------
```Pedig```     |     Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam. Missing parents are coded as NA.

## Value


 Additive relationship matrix.


## Examples

```r 
 ## Not to run ##
 
 ## NRM(Pedig)
 
 ## End(Not run)
 
 ``` 

