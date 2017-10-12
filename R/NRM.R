#' @title Calculates the Pedigree-based Additive Relationship Matrix 
#' 
#' 
#' @description Calculates the the Pedigree-based Additive Relationship Matrix. This is twice the pedigree based kinship matrix.
#' 
#' 
#' @param Pedig Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam. Missing parents are coded as NA. Both parents must either be missing or present.
#' @param keep If \code{keep} is provided then kinships are computed only for these animals and their ancestors.
#' @param keep.only If \code{keep.only} is provided then kinships are computed only for these animals.
#' @param AFounder Additive relationship matrix of the founders. The row names are the ids of the founders. By default, founders are assumed to be unrelated. Founders not included in this matrix are also assumed to be unrelated. 
#' 
#' 
#' @return Additive relationship matrix.
#' 
#' 
#' @examples
#' ## Not to run ##
#' 
#' ## NRM(Pedigree, keep.only=NULL, keep=keep.only, AFounder=NULL)
#'
#' ## End(Not run)
#' 
#' @export NRM
#' 
#' @import stats
#' @useDynLib GDmat
#' @importFrom Rcpp evalCpp
NRM <- function(Pedig, keep.only=NULL, keep=keep.only, AFounder=NULL){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  if(PedigAsDataTable){
    Pedig <- as.data.frame(Pedig)
    setDF(Pedig)
    }
  
  ids   <- as.character(Pedig[[1]])
  Pedig <- prePed(Pedig[,1:3], keep=keep, addNum=TRUE)

  if(is.null(keep.only)){
    keep.only <- ids
  }else{
    keep.only <- as.character(keep.only)
    keep.only <- setdiff(keep.only, c(NA, "", " ", "0"))
    keep.only <- ids[ids %in% keep.only]
  }

  if(is.null(AFounder)){
    AFounder   <- matrix(1, 0, 0)
    numFounder <- integer(0)
  }else{
    idF        <- Pedig$Indiv[Pedig$Indiv %in% rownames(AFounder)]
    AFounder   <- AFounder[idF, idF, drop=FALSE]
    numFounder <- Pedig[idF, "numIndiv"]
  }
  if(length(keep.only)<0.5*nrow(Pedig)){
    indKeep <- Pedig$Indiv[Pedig$Indiv %in% keep.only]
    numKeep <- Pedig[indKeep,"numIndiv"]
    Pedig$nOff <- 0
    x <- table(Pedig$Sire)
    Pedig[names(x),"nOff"] <- x
    x <- table(Pedig$Dam)
    Pedig[names(x),"nOff"] <- x
    pedKin  <- rcpp_makeA_lowMem(as.integer(Pedig$numSire), as.integer(Pedig$numDam), AFounder, as.integer(numFounder-1), as.character(indKeep), as.integer(numKeep-1), as.integer(Pedig$Indiv %in% indKeep), as.integer(Pedig$nOff))
  }else{
    pedKin <- rcpp_makeA(as.integer(Pedig$numSire), as.integer(Pedig$numDam), AFounder, as.integer(numFounder-1), as.character(Pedig$Indiv))
  }
  if(identical(Pedig$Indiv, keep.only)){return(pedKin)}
  pedKin[keep.only, keep.only]
}
