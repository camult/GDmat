#' @title Calculates the Pedigree-based Additive Relationship Matrix 
#' 
#' 
#' @description Calculates the the Pedigree-based Additive Relationship Matrix. This is twice the pedigree based kinship matrix.
#' 
#' 
#' @param Pedig Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam. Missing parents are coded as NA.
#' 
#' 
#' @return Additive relationship matrix.
#' 
#' 
#' @examples
#' ## Not to run ##
#' 
#' ## NRM(Pedig)
#'
#' ## End(Not run)
#' 
#' @export NRM
NRM <- function(Pedig){
  x <- ped
  n <- nrow(x)
  ret <- matrix(0, n, n)
  diag(ret) <- 1
  rownames(ret) <- colnames(ret) <- x[,1]
  asc1 <- as.character(x[,2])
  asc2 <- as.character(x[,3])
  testAsc1 <- !is.na(x[,2])
  testAsc2 <- !is.na(x[,3])
  testAsc <- testAsc1 & testAsc2
  set <- which(!(!testAsc1 & !testAsc2))
  for (i in set) {
    if (testAsc[i]) { ret[i, i] <- 1 + 0.5 * ret[asc1[i], asc2[i]] }
    j <- 1:(i - 1)
    if (testAsc1[i]) { tmp1 <- 0.5 * ret[asc1[i], j] } else { tmp1 <- 0 }
    if (testAsc2[i]) { tmp2 <- 0.5 * ret[asc2[i], j] } else { tmp2 <- 0 }
    ret[i, j] <- ret[j, i] <- tmp1 + tmp2
  }
  return(ret)
}
