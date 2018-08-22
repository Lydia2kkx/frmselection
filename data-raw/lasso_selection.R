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

a <- lasso.select(x, coefi)
a$lambda.common
a$lambda.factor
a$lasso.index
a$modified.coefficients
