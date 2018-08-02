setwd("E:/Master/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
#head(x)
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
#head(x)
###验证backward stepwise 结果是对的！
b<-frm(y, x, linkfrac="logit")
-2*ll_logit(b$p, cbind(rep(1,nrow(x)), x), y) + 14
b<-frm(y, x[,1:5], linkfrac="logit")
-2*ll_logit(b$p, cbind(rep(1,nrow(x)), x[,1:5]), y) + 12
head(x[,1:5])
###验证forward stepwise 结果是对的，和backward stepwise一样
b<-frm(y, x[,3:5,drop = FALSE], linkfrac="logit")
-2*ll_logit(b$p, cbind(rep(1,nrow(x)), x[,3:5]), y) + 6
###验证allsubsets 结果和上面2种方法一样

b<-frm(y, x[,1, drop=FALSE], linkfrac="logit")
-2*ll_logit(b$p, cbind(rep(1,nrow(x)), x[,1]), y) + 4



p <- ncol(x)
n <- nrow(x)
x <- cbind(rep(1,n), x[,1:5])



head(x)

