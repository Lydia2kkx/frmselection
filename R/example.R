setwd("E:/Master/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
head(x)
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
head(x)
b<-frm(y, x, linkfrac="logit")

lb <- lm(y~x)
lb
summary(lb)

p <- ncol(x)
n <- nrow(x)
x <- cbind(rep(1,n), x)




