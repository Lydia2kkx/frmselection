lasso_select <- function(x, y, family = "beta", nlambda = 100, upper = 1e+08, lower = 1e+03, data = data){
  xname <- colnames(x)
  ind <- sapply(x, is.factor)
  #select the factor variables, which are penalized differently compared to the common variable
  ind <- which(ind)
  f1 <- as.formula(paste(" ~ ", paste(xname[-ind], collapse = "+")))
  f2 <- as.formula(paste(" ~ ", paste(xname[ind], collapse = "+")))
  formula <- y ~ la(f1) + la(f2)
  b <- bamlss(formula, data = data, family = family, sampler = FALSE,
              optimizer = lasso, nlambda = nlambda, upper = upper, lower = lower,
              multiple = FALSE)
  coefi <- coef(b)
  return(list(formula, f1, f2))
}

set.seed(123)
d <- GAMart()
d1 <- d[,c(9:11,13)]
head(d1)
class(d1)
lasso_select(d1, d$bnum, data = d)



xnam <- paste0("x", 1:3)
f <- as.formula(paste(" ~ ", paste(xnam, collapse = "+")))
ind <- sapply(d1, is.factor)
ind <- which(ind)
fa <- d1[,ind]
head(fa)
ot <- d1[,-ind]
head(ot)
f1 <- bnum ~ la(f)+la(id)
b <- bamlss(f1, data = d, family = "beta", sampler = FALSE,
            optimizer = lasso, nlambda = 100, upper = 1e+08, lower = 1e+03,
            multiple = FALSE)

ind <- sapply(d1, is.factor)
ind <- which(ind)
xname <- colnames(d1)
f1 <- as.formula(paste(" ~ ", paste(xname[-ind], collapse = "+")))
f2 <- as.formula(paste(" ~ ", paste(xname[ind], collapse = "+")))
formula <- bnum ~ la(f1) + la(f2)
b <- bamlss(formula, data = d, family = "beta", sampler = FALSE,
            optimizer = lasso, nlambda = 100, upper = 1e+08, lower = 1e+03,
            multiple = FALSE)
summary(b)

f1 <- as.formula(paste(" ~ ", paste(xname[-ind], collapse = "+")))
f1
f2 <- as.formula(paste(" ~ ", paste(xname[ind], collapse = "+")))
f2
formula <- bnum ~ la(f1) + la(f2)
formula
b <- bamlss(formula, data = d, family = "beta", sampler = FALSE,
            optimizer = lasso, nlambda = 100, upper = 1e+08, lower = 1e+03,
            multiple = FALSE)
summary(b)

e <- quote(`foo bar`)
deparse(e)
deparse(e, backtick = TRUE)
e <- quote(`foo bar`+1)
deparse(e)
deparse(e, control = "all")
f1 <- reformulate(paste(xname[-ind], collapse = "+"))
f2 <- reformulate(paste(xname[ind], collapse = "+"))
formula <- bnum ~ la(f1) + la(f2)
formula
b <- bamlss(formula, data = d, family = "beta", sampler = FALSE,
            optimizer = lasso, nlambda = 100, upper = 1e+08, lower = 1e+03,
            multiple = FALSE)


a<-reformulate(paste(xname[-ind], collapse = "+"), y)
newa <- update(a, NULL ~ .)
