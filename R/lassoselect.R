#' Lasso procedure based on bamlss
#'
#' This function is used to simplify the model based on lasso method with bamlss() function.
#'
#' @usage lasso.select(x, coefficient, threshold)
#'
#' @param x  The independent variables
#' @param coefficent The coefficients estimated by bamlss()
#' @param threshold The threshold controls which coefficients should be set to 0. Default value is 10e-3. Coefficients smaller than threshold would be set to 0.
#'
#' @return A list of results that contains lambda parameters of common and factor variables,
#' lasso index which represents the index of coefficients that are set to be zero and the simplified coefficients.
#' @export

lasso.select <- function(x, coefficent, threshold = 1e-3){
  n <- length(coefficent) #n is the number of the coefficients
  ncommon <- length(x[,-ind]) #the number of the common variables
  nfactor <- length(levels(x[,ind]))  #find the number of levels of factor variables
  #nfactor <- length(unique(x[,ind]))
  #Extract the penalized parameters of common and factor variables separately.
  tau.common <- coefficent[ncommon+1]
  tau.factor <- coefficent[ncommon+nfactor+1]
  #The relation between tau and lambda is an inverse relationship
  lambda.common <- 1/tau.common
  names(lambda.common) <- "mu.s.la(f1).lambda"
  lambda.factor <- 1/tau.factor
  names(lambda.factor) <- "mu.s.la(f2).lambda"
  index.tau <- c(ncommon+1, ncommon+nfactor+1)
  coefi.new <- coefficent[-index.tau] #the new coeffients doesn't contain the penalized parameters
  lasso.index <- NULL
  #Users can set this threshold as whatever they like. The default value is 1e-3
  for(i in 1:length(coefi.new)){
    if(abs(coefi.new[i]) < threshold){
      coefi.new[i] = 0
      lasso.index <- c(lasso.index, i)
    }
  }
  return(list(lambda.common = lambda.common, lambda.factor = lambda.factor,
              lasso.index = lasso.index, modified.coefficients = coefi.new))
}
