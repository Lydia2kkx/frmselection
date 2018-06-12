setwd("E:/Master/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
x <- as.data.frame(x, col.names = TRUE)
subdat1 <- dat[y<1,]
head(dat)
head(subdat1)
nrow(subdat1)
suby1 <- subdat1[,1]
subx1 <- x[y<1,]
class(subx1)
head(subx1)
library(bamlss)
aa <- bamlss::bamlss(suby1 ~ subx1$mrate + subx1$lemp + subx1$lemp2 + subx1$age +subx1$age2 +subx1$sole, family = "beta")
summary(aa)
aa$model.stats$optimizer$logLik
#############the dependent variables of beta regression must be in (0, 1)
subdat <- dat[x[,1]<=1,]
nrow(subdat)
suby <- subdat[,1]
subx <- x[x[,1]<=1,]
g <- betareg::betareg(suby ~ mrate + lemp + lemp2 + age + age2 + sole, data = subx)
summary(g)
#############################################
f <- suby1 ~ subx1$mrate + subx1$lemp + subx1$lemp2 + subx1$age +subx1$age2 +subx1$sole
f <- list(suby1 ~  subx1$mrate + subx1$lemp + subx1$lemp2 + subx1$age +subx1$age2 +subx1$sole,
          sigma2 ~ subx1$mrate + subx1$lemp + subx1$age +subx1$sole)
f <- suby1 ~ subx1$mrate + subx1$lemp + subx1$lemp2 + subx1$age +subx1$age2

b <- bamlss(f, family = "beta", sampler = FALSE,
            optimizer = lasso, nlambda = 100, upper = 1e+03, lower = 1,
            multiple = FALSE)
summary(b)
class(subx1)

b$parameters
b$model.stats
b$model.stats$optimizer$lasso.stats[1,]
b$model.stats$optimizer$lasso.stats[2,]
b$model.stats$optimizer$lasso.stats[3,]
b$model.stats$optimizer$lasso.stats[20,]
b$model.stats$optimizer$lasso.stats[40,]
b$model.stats$optimizer$lasso.stats[69,]
b$model.stats$optimizer$lasso.stats[80,]

gy <- betareg(suby1 ~ subx1[,1] + subx1[,2], data = subx1)
summary(gy)
gy1 <- betareg::betareg(suby1 ~ mrate + lemp + lemp2 + age + age2 + sole, data = subx1)
summary(gy1)
gy1$coefficients
gy2 <- betareg::betareg(suby1 ~ mrate + lemp + lemp2 + age + age2,  data = subx1)
summary(gy2)
gy1$loglik
gy2$loglik
gy$df.null
gy$df.residual
gy1$df.null
gy1$df.residual
gy2$df.null
gy2$df.residual
gy3 <- betareg::betareg(suby1 ~ mrate, data = subx1)
summary(gy3)
gy3$loglik
gy3$df.null
gy3$df.residual

colnames(subx1)
formula <- cbind(subx1[,1], subx1[,3])
gy <- betareg::betareg(suby1 ~ formula + subx1[,2], data = subx1)
summary(gy)

