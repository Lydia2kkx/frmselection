lasso_select <- function(x, y, family = "beta", nlambda = 100, upper = 1e+08, lower = 1e+03, data = data){
  xname <- colnames(x)
  ind <- sapply(x, is.factor)
  #select the factor variables, which are penalized differently compared to the common variable
  ind <- which(ind)
  f1 <- as.formula(paste(" ~ ", paste(xname[-ind], collapse = "+")))
  f2 <- as.formula(paste(" ~ ", paste(xname[ind], collapse = "+")))

  #browser()
  formula <- y ~ bamlss::la(f1) + bamlss::la(f2)
  b <- bamlss::bamlss(formula, data = data, family = family, sampler = FALSE,
              optimizer = lasso, nlambda = nlambda, upper = upper, lower = lower,
              multiple = FALSE)
  coefi <- coef(b)
  return(list(formula, f1, f2))
}

lasso_select(d1, bnum, data = d)

library("bamlss")
set.seed(123)
d <- GAMart()
d1 <- d[,c(9:11,13)]
bnum <- d$bnum
head(d1)
class(d1)


lasso_select <- function(x, y, family = "beta", nlambda = 100, upper = 1e+08, lower = 1e+03, data = data){
  xname <- colnames(x)
  ind <- sapply(x, is.factor)
  #select the factor variables, which are penalized differently compared to the common variable
  ind <- which(ind)
  f1 <- as.formula(paste(" ~ ", paste(xname, collapse = "+")))
  #f2 <- as.formula(paste(" ~ ", paste(xname[ind], collapse = "+")))
  formula <- y ~ la(f1)
  b <- bamlss(formula, data = data, family = family, sampler = FALSE,
              optimizer = lasso, nlambda = nlambda, upper = upper, lower = lower,
              multiple = FALSE)
  coefi <- coef(b)
  return(list(formula, f1))
}

setwd("E:/Master/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
#colnames(x)
#head(x)
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
x <- as.data.frame(x)
x[,6] <- as.factor(x[,6])
subdat <- dat[x[,1]<=1,]
suby <- subdat[,1]
subx <- x[x[,1]<=1,]


ind <- sapply(subx, is.factor)
ind <- which(ind)
xname <- colnames(subx)
f1 <- as.formula(paste(" ~ ", paste(xname[-ind], collapse = "+")))
f2 <- as.formula(paste(" ~ ", paste(xname[ind], collapse = "+")))
formula <- suby ~ la(f1) + la(f2)
b <- bamlss(formula, data = subdat, family = "beta", sampler = FALSE,
            optimizer = lasso, nlambda = 100, upper = 1e+08, lower = 1e+03,
            multiple = FALSE)



lasso_select <- function(x, y, family = "beta", nlambda = 100, upper = 1e+08, lower = 1e+03, data = data){
  xname <- colnames(x)
  ind <- sapply(x, is.factor)
  #select the factor variables, which are penalized differently compared to the common variable
  ind <- which(ind)
  f1 <- reformulate(paste(xname[-ind], collapse = "+"))
  f2 <- reformulate(paste(xname[ind], collapse = "+"))
  formula <- y ~ la(f1) + la(f2)
  b <- bamlss(formula, data = data, family = family, sampler = FALSE,
              optimizer = lasso, nlambda = nlambda, upper = upper, lower = lower,
              multiple = FALSE)
  coefi <- coef(b)
  return(list(formula, f1, f2))
}
lasso_select(d1, d$bnum, data = d)
