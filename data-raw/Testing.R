#Testing whether the modified function works
library(frm)
setwd("E:/Master/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
x <- as.data.frame(x, row.names = TRUE)
d <- data.frame(y,x)
formula <- as.formula(paste("y ~ ", paste(colnames(x), collapse = "+")))

###If the response data could be equal to 1, then the "logit" link works, but the "probit" link doesn't work.
###logit link
library(bamlss)
b <- bamlss(formula, data = d, family = frm_bamlss(link = "logit") ,sampler = FALSE,
            multiple = FALSE)
coef(b)
f <- frm(y, x, linkfrac = "logit", table = FALSE)
f$p
###probit link
b1 <- bamlss(formula, data = d, family = frm_bamlss(link = "probit") ,sampler = FALSE,
            multiple = FALSE)
coef(b1)
b2 <- bamlss(formula, data = d, start = coef(b1), sampler = FALSE,
             multiple = FALSE) #doesn't work!
f1 <- frm(y, x, linkfrac = "probit", table = FALSE)
f1$p
#############If the response data are smaller than 1, then it works
subdat <- d[y<1,]
suby <- subdat[,1]
subx <- subdat[,-1]
#############
###probit link
bsub <- bamlss(formula, data = subdat, family = frm_bamlss(link = "probit") ,sampler = FALSE,
             multiple = FALSE)
coef(bsub)
f2 <- frm(suby, subx, linkfrac = "probit", table = FALSE)
f2$p

###########
b <- bamlss(formula, data = subdat, family = frm_bamlss(link = "logit") ,sampler = FALSE,
            multiple = FALSE)
coef(b)
f11 <- frm(suby, subx, linkfrac = "logit", table = FALSE)
f11$p
