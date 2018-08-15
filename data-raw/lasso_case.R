setwd("E:/Master/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
x <- as.data.frame(x, row.names = TRUE)
d <- data.frame(y,x)
d$sole <- as.factor(d$sole)
subdat <- d[y<1,]
suby <- subdat[,1]
subx <- subdat[,-1]
xname <- colnames(subx)
ind <- sapply(subx, is.factor)
ind <- which(ind) #find the index of the factor variables
f1 <- as.formula(paste(" ~ ", paste(xname[-ind], collapse = "+")))
f2 <- as.formula(paste(" ~ ", paste(xname[ind], collapse = "+")))
formula <- y ~ la(f1) + la(f2) #form the formula with two kinds of variables
b <- bamlss(formula, data = subdat, family = "beta", sampler = FALSE,
            optimizer = lasso, nlambda = 100, upper = 1e+03, lower = 1e-01,
            multiple = FALSE)
coefi <- coef(b)

lasso.select <- function(x, coefficent, threshold = 1e-3){
  n <- length(coefficent) #n is the number of the coefficients
  ncommon <- length(x[,-ind]) #the number of the common variables
  nfactor <- length(levels(x[,ind]))  #find the number of levels of factor variables
  #nfactor <- length(unique(x[,ind]))
  #Extract the penalized parameters of common and factor variables separately.
  tau.common <- coefi[ncommon+1]
  tau.factor <- coefi[ncommon+nfactor+1]
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
  return(list(lambda.common = lambda.common, lambda.factor = lambda.factor, lasso.index = lasso.index, modified.coefficients = coefi.new))
}
a <- lasso.select(subx, coefi)
a


gy1 <- betareg::betareg(suby ~ mrate + lemp + lemp2 + age + age2 + sole, data = subx)
gy1$coefficients
b <- bamlss(suby ~ mrate + lemp + lemp2 + age + age2 + sole, data = subdat, family = "beta", sampler = FALSE,
            multiple = FALSE)
coefi <- coef(b)
