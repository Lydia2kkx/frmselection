# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
summary(lm.D9)
##This is a package of model selection in fractional response model
##In the first part I will rely on the frm package to do the model selection
##In the second part I will rely on the bamlss package and let users to apply beta distribution.
##

library(frm)
N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2","X3")

ym <- exp(X[,1]+X[,2]+X[,3]+u)/(1+exp(X[,1]+X[,2]+X[,3]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1
#frm estimation of a logit fractional regression model
a<-frm(y, X, linkfrac="logit")

library(MASS)
stepAIC(a, direction = "forward")
#报错
step(a, direction = "forward")
#报错
library(StepReg)
stepwise(cbind(y,rep(1,250),X),y,selection = "forward",select = "AIC")
#结果很可疑


frmselet<-function(y,x, method = c("both", "backward", "forward", "allset"),
                   criteria = c("AIC", "BIC", ...)){



}




#Define log-likelihood function
#Using the Bernoulli quasi-log-likelihood function
#Logit model
ll_logit<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(exp(x[i,]%*%params)/(1+exp(x[i,]%*%params))) + (1-y[i])*log(1-exp(x[i,]%*%params)/(1+exp(x[i,]%*%params)))
  }
  return(sum(f))
}
x <- cbind(rep(1,nrow(X)), X)
ll_logit(a$p, x, y)
#Probit model
ll_probit<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(pnorm(x[i,]%*%params)) + (1-y[i])*log(1-pnorm(x[i,]%*%params))
  }
  return(sum(f))
}
#Loglog model
ll_loglog<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(exp(-exp(x[i,]%*%params))) + (1-y[i])*log(1-exp(-exp(x[i,]%*%params)))
  }
  return(sum(f))
}
#Cloglog model
ll_cloglog<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(1-exp(-exp(x[i,]%*%params))) + (1-y[i])*log(exp(-exp(x[i,]%*%params)))
  }
  return(sum(f))
}
#Cauchit model
ll_cauchit<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log((atan(x[i,]%*%params)+pi/2)/pi) + (1-y[i])*log(1-(atan(x[i,]%*%params)+pi/2)/pi)
  }
  return(sum(f))
}



p <- ncol(x)
n <- nrow(x)
info_criteria <- function(k){
  if(a$link == "logit")
  return(-2*ll_logit() + k * p)
}
if(criteria == "AIC"){
  return(info_criteria(2))
}
if(criteria == "BIC"){
  return(info_criteria(log(n)))
}

##%in% 用法 a %in% table  a值是否包含于table中，为真输出TURE，否者输出FALSE

##forward stepwise
p <- ncol(x)
forward <- function(){
  for(i in 3:p){
    result <- list()
    criteria <- vector()
    result[[i]] <- frm(y, x[,c(2,i)], linkfrac="logit")
    criteria[i] <- ll_logit(result[[i]]$p, x[i], y)
  }
  return(criteria)
}
for(i in 3:p){
  result <- list()
  criteria <- vector()
  result[[i]] <- frm(y, x[,c(2,i)], linkfrac="logit")
  criteria[i] <- ll_logit(result[[i]]$p, x[i], y)
  criteria
}

a1 <- frm(y, x[,c(2,3)], linkfrac="logit")
criteria1 <- ll_logit(a1$p, x[,1:3], y)
criteria1

a2 <- frm(y, x[,c(2,4)], linkfrac="logit")
criteria2 <- ll_logit(a1$p, x[,c(1:2,4)], y)
criteria2

a <- list(a1, a2)
a
