#The procedure of lasso selection
set.seed(123)
library(bamlss)
d <- GAMart()
x <- d[,c(9:11,13)]
head(x)
y <- d$bnum
xname <- colnames(x)
ind <- sapply(x, is.factor)
#select the factor variables, which are penalized differently compared to the common variable
ind <- which(ind) #find the index of the factor variables
f1 <- as.formula(paste(" ~ ", paste(xname[-ind], collapse = "+")))
f2 <- as.formula(paste(" ~ ", paste(xname[ind], collapse = "+")))
formula <- y ~ la(f1) + la(f2) #form the formula with two kinds of variables
b <- bamlss(formula, data = d, family = "beta", sampler = FALSE,
                    optimizer = lasso, nlambda = 10, upper = 1e+08, lower = 1e+03, multiple = TRUE)
coefi <- lasso.coef(b)
#The coefficients include the common variable, the lambda parameter(penalized with lasso) of common variables,
#the factor variables(different levels of factor variables have different coefficients),
#the lambda parameter(penalized with lasso) of factor variables and the intercept

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

a <- lasso.select(x, coefi)
a$lambda.common
a$lambda.factor
a$lasso.index
a$modified.coefficients
