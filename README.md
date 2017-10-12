# `GDmat-package`: Additive and Dominance Relationship Matrices

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
```X```     |     (numeric, $n$ ) matrix of genotypes with a number of rows identical to the number of SNPs (m) and a number of columns identical to the number of genotyped individuals (n). Values in the matrix are 0 and 2 for homozygous and 1 for heterozygous.
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

# `NRM`: Calculates the Pedigree-based Additive Relationship Matrix

## Description


 Calculates the the Pedigree-based Additive Relationship Matrix. This is twice the pedigree based kinship matrix.


## Usage

```r
NRM(Pedig, keep.only = NULL, keep = keep.only, AFounder = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
```Pedig```     |     Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam. Missing parents are coded as NA. Both parents must either be missing or present.
```keep.only```     |     If `keep.only` is provided then kinships are computed only for these animals.
```keep```     |     If `keep` is provided then kinships are computed only for these animals and their ancestors.
```AFounder```     |     Additive relationship matrix of the founders. The row names are the ids of the founders. By default, founders are assumed to be unrelated. Founders not included in this matrix are also assumed to be unrelated.

## Value


 Additive relationship matrix.


## References


 Hyun Min Kang, Noah A. Zaitlen, Claire M. Wade, Andrew Kirby, David Heckerman, Mark J. Daly and Eleazar Eskin, 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723. doi:10.1534/genetics.107.080101.
 
 Gianola D, Schon C-C, 2016. Cross-Validation Without Doing Cross-Validation in Genome-Enabled Prediction. G3: Genes|Genomes|Genetics. 6(10):3107-3128. doi:10.1534/g3.116.033381.


## Examples

```r 
 ## Not to run ##
 
 ## NRM(Pedigree, keep.only=NULL, keep=keep.only, AFounder=NULL)
 
 ## End(Not run)
 
 ``` 

