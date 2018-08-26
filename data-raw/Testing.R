#Testing whether the modified quasi-binomial family function works
set.seed(123)
library(bamlss)
d <- GAMart()
head(d)
library(frmselection)
###logit link
b <- bamlss(bnum ~ x1 + x2 + x3, data = d, family = frm_bamlss(link = "logit"),
            sampler = FALSE, multiple = FALSE)
coef(b)
###Using frm() to check, get the same result
x <- d[, 9:11]
f <- frm::frm(d$bnum, x, linkfrac = "logit", table = FALSE)
f$p

###probit link
b1 <- bamlss(bnum ~ x1 + x2 + x3, data = d, family = frm_bamlss(link = "probit"),
            sampler = FALSE, multiple = FALSE)
coef(b1)
###Using frm() to check, get the very similar result
f1 <- frm::frm(d$bnum, x, linkfrac = "probit", table = FALSE)
f1$p


