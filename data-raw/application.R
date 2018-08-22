#case
setwd("E:/Master/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
summary(dat)
#The standard deviation of variables
sqrt(diag(var(dat)))
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
#QMLE
library(frm)
f <- frm(y, x, linkfrac = "logit", table = FALSE)
f$p
#The standard deviation of estimated coefficients
sqrt(diag(f$p.var))

x <- as.data.frame(x, row.names = TRUE)
#OLS
ols <- lm(y ~ x$mrate + x$lemp + x$lemp2 + x$age + x$age2 +x$sole)
coef(ols)

#frmselect() to do the model selection
#The four examples
x <- as.matrix(x, row.names = TRUE)
frmselect(x,y, criterion = "AIC",linkfrac = "logit")#The default is forward
#frmselect(x,y, criterion = "AIC",linkfrac = "logit", method = "allsubsets")
#frmselect(x,y, criterion = "AIC",linkfrac = "logit", method = "backward")
#frmselect(x,y, criterion = "AIC",linkfrac = "logit", method = "both")

#betaselect() to do the model selection
#The four examples
newdat <- data.frame(y, x)
subdat <- newdat[y<1,]
suby <- subdat[,1]
subx <- subdat[,-1]
subx <- as.matrix(subdat[,-1])
nrow(subx)
betaselect(subx, suby)
#betaselect(subx, suby, method = "backward")
#betaselect(subx, suby, method = "both")
#betaselect(subx, suby, method = "allsubsets")

#Another way by replacing 1 with 0.9999
y[y==1] <- 0.9999
betaselect(x, y)
#betaselect(x, y, method = "backward")
#betaselect(x, y, method = "both")
#betaselect(x, y, method = "allsubsets")

#Testing whether the modified quasi-binomial family frm_bamlss()function works
x <- as.data.frame(x, row.names = TRUE)
y <- dat[,1]
d <- data.frame(y,x)
formula <- as.formula(paste("y ~ ", paste(colnames(x), collapse = "+")))
library(bamlss)
b <- bamlss(formula, data = d, family = frm_bamlss(link = "logit") ,sampler = FALSE,
            multiple = FALSE)
coef(b)
f <- frm::frm(y, x, linkfrac = "logit", table = FALSE)
f$p

#Lasso Procedure
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
x <- as.data.frame(x, row.names = TRUE)
x$sole <- as.factor(x$sole)
d <- data.frame(y,x)
xname <- colnames(x)
ind <- sapply(x, is.factor)
ind <- which(ind) #find the index of the factor variables
f1 <- as.formula(paste(" ~ ", paste(xname[-ind], collapse = "+")))
f2 <- as.formula(paste(" ~ ", paste(xname[ind], collapse = "+")))
formula <- y ~ la(f1) + la(f2) #form the formula with two kinds of variables
b <- bamlss(formula, data = d, family = frm_bamlss(link = "logit"), sampler = FALSE,
            optimizer = lasso, nlambda = 100, upper = 1e+08, lower = 1e+03, multiple = TRUE)
coefi <- lasso.coef(b)
lasso.select(x, coefi)


subdat <- d[y<1,]
suby <- subdat[,1]
subx <- subdat[,-1]
xname <- colnames(subx)
ind1 <- sapply(subx, is.factor)
ind1 <- which(ind1) #find the index of the factor variables
f11 <- as.formula(paste(" ~ ", paste(xname[-ind], collapse = "+")))
f21 <- as.formula(paste(" ~ ", paste(xname[ind], collapse = "+")))
formula1 <- y ~ la(f11) + la(f21) #form the formula with two kinds of variables
b1 <- bamlss(formula1, data = subdat, family = "beta", sampler = FALSE,
            optimizer = lasso, nlambda = 10, upper = 10e+8, lower = 10e+3,
            multiple = TRUE)
coefi1 <- lasso.coef(b1)
coefi1
lasso.select(subx, coefi1)




gy1 <- betareg::betareg(suby ~ mrate + lemp + lemp2 + age + age2 + sole, data = subx)
gy1$coefficients
b <- bamlss(suby ~ mrate + lemp + lemp2 + age + age2 + sole, data = subdat, family = "beta", sampler = FALSE,
            multiple = FALSE)
coefi <- coef(b)


b <- bamlss(formula, data = d, family = frm_bamlss(link = "logit"), sampler = FALSE,
       optimizer = lasso, multiple = FALSE)

a <- lars(x, y, type = "lasso")
plot(a)

aa <- glmnet(as.matrix(d[,-1]), d[,1], standardize = TRUE, alpha = 1)
plot(aa)
aa$beta

b <- bamlss(y~la(~mrate + lemp + lemp2 + age + age2 + sole), data = d, family = frm_bamlss(link = "logit"), sampler = FALSE,
            optimizer = lasso, nlambda = 100, upper = 400, lower = 0.001, multiple = FALSE)
coefi <- coef(b)
coefi
