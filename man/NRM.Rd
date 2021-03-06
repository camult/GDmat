% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NRM.R
\name{NRM}
\alias{NRM}
\title{Calculates the Pedigree-based Additive Relationship Matrix}
\usage{
NRM(Pedig)
}
\arguments{
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam. Missing parents are coded as NA.}
}
\value{
Additive relationship matrix.
}
\description{
Calculates the the Pedigree-based Additive Relationship Matrix. This is twice the pedigree based kinship matrix.
}
\examples{
## Not to run ##

## NRM(Pedig)

## End(Not run)

}
