#case：
library(frm)
setwd("E:/Master/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
#colnames(x)
#head(x)
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
f <- frm(y, x, linkfrac = "logit", table = FALSE)
summary(f)
f
x <- as.data.frame(x, row.names = TRUE)
head(x)
#head(x)
#Full Sample
#OLS
ols <- lm(y ~ x$mrate + x$lemp + x$lemp2 + x$age + x$age2 +x$sole)
ols
library(MASS)
stepAIC(ols, direction = "both")
#Subsample
#OLS
subdat <- dat[x[,1]<=1,]
head(subdat)
nrow(subdat)
suby <- subdat[,1]
subx <- x[x[,1]<=1,]
class(subx)
subols <- lm(suby ~ subx$mrate + subx$lemp + subx$lemp2 + subx$age + subx$age2 + subx$sole)
subols
stepAIC(subols, direction = "backward")
stepAIC(subols, scope = list(lower = ~1, upper = ~ subx$mrate * subx$lemp + subx$lemp2 + subx$age + subx$age2 + subx$sole),direction = "forward")
#So without sole the AIC will be smallest!


#Full Sample
#FRM
xnew <- cbind(dat[,2], dat[,2]^2, log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
colnames(xnew) <- c("mrate", "mrate2", "lemp", "lemp2","age", "age2","sole")
frmselect(xnew, y, criterion = "AIC",linkfrac = "logit", method = "both")
frmselect(xnew, y, criterion = "BIC",linkfrac = "logit", method = "allsubsets")
#By model selecting the first five variables are selected in the model which is the best one.

#forward stepwise
forward <- function(x,y, criterion = c("AIC", "BIC")){
  x <- cbind(rep(1,nrow(x)), x)
  p <- ncol(x)
  n <- nrow(x)
  if(criterion == "AIC"){
    k <- 2
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
  forward.1 <- function(x,y){
    result <- list()
    criteria <- vector()
    for(i in 2:p){
      result[[i]] <- frm(y, x[,i,drop = FALSE], linkfrac="logit", table = FALSE)
      criteria[i] <- -2*ll_logit(result[[i]]$p, x[,c(1,i)], y)
    }
    min_criteria <- min(na.omit(criteria))
    index <- which.min(criteria)-1
    variable <- paste0("X", index)
    return( list(min = min_criteria, index = index, variable = variable))
  }
  first_step <- forward.1(x,y)
  ind <- first_step$index + 1
  criteria_old <- first_step$min + k * 2
  m <- 2
  repeat{
    result <- list()
    criteria <- vector()
    for(i in 2:p){
      if(!(i %in% ind)){
        result[[i]] <- frm(y, x[, sort(c(ind,i)), drop = FALSE], linkfrac = "logit", table = FALSE)
        criteria[i] <- -2*ll_logit(result[[i]]$p, x[,sort(c(1,ind,i))], y) + k * (m+1)
      }
    }
    criteria_new <- min(na.omit(criteria))
    if(criteria_new >= criteria_old){
      break
    }
    ind <- c(ind, which.min(criteria))
    criteria_old <- criteria_new
    m <- m + 1
  }
  min_criteria <- criteria_old
  index <- ind - 1
  variable <- colnames(x)[ind]
  return(list(min = min_criteria, index = index, variable = variable))
}

forward(x, y, criterion = "AIC")
forward(x, y, criterion = "BIC")

#backward stepwise
backward <- function(x,y, criterion = c("AIC", "BIC")){
  n <- nrow(x)
  if(criterion == "AIC"){
    k <- 2
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
  first_step <- frm(y, x, linkfrac="logit", table = FALSE)
  x <- cbind(rep(1,nrow(x)), x)
  p <- ncol(x)
  criteria_old <- -2*ll_logit(first_step$p, x, y) + k*p
  ind <- NULL
  m <- 1
  repeat{
    result <- list()
    criteria <- vector()
    for(i in 2:p){
      if(!(i %in% ind)){
        result[[i]] <- frm(y, x[,-sort(c(1,ind,i)), drop = FALSE], linkfrac="logit", table = FALSE)
        criteria[i] <- -2*ll_logit(result[[i]]$p, x[,-sort(c(ind,i))], y) + k*(p-m)
      }
    }
    criteria_new <- min(na.omit(criteria))
    if(criteria_new >= criteria_old){
      break
    }
    ind <- c(ind, which.min(criteria))
    criteria_old <- criteria_new
    m <- m + 1
  }
  min_criteria <- criteria_old
  if(is.null(ind)){
    index <- ind
    variable <- colnames(x)
  }else{
    ind <- ind - 1
    index <- seq(1, p-1)[-ind]
    variable <- colnames(x)[index+1]
  }
  return(list(min = min_criteria, index = index, variable = variable))
}
backward(x, y, criterion = "AIC")
backward(x, y, criterion = "BIC")

#allsubsets
allsubsets <- function(x, y, criterion = c("AIC", "BIC")){
  x <- cbind(rep(1,nrow(x)), x)
  p <- ncol(x)
  n <- nrow(x)
  if(p-1 > 8){
    stop("Too many variables, using all-subsets method is time-consuming!")
  }
  if(criterion == "AIC"){
    k <- 2
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
  criteria_min_part <- NULL
  ind <- list()
  for(i in 1:(p-1)){
    count <- choose(p-1, i)
    com <- combn(p-1, i)
    result <- list()
    criteria <- vector()
    for(j in 1:count){
      result[[j]] <- frm(y, x[, com[,j]+1, drop = FALSE], linkfrac="logit", table = FALSE)
      criteria[j] <- -2*ll_logit(result[[j]]$p, x[,(c(1,com[,j]+1))], y)+ k*(i+1)
    }
    criteria_min_part <- c(criteria_min_part, min(criteria))
    ind <- c(ind, list(com[, which.min(criteria)]))
  }
  min_criteria <- min(criteria_min_part)
  index <- ind[[which.min(criteria_min_part)]]
  variable <- colnames(x)[index+1]
  return(list(min = min_criteria, index = index, variable = variable))
}
allsubsets(x, y, criterion = "AIC")
allsubsets(x, y, criterion = "BIC")

#both
both <- function(x,y, criterion = c("AIC", "BIC")){
  x <- cbind(rep(1,nrow(x)), x)
  p <- ncol(x)
  n <- nrow(x)
  if(criterion == "AIC"){
    k <- 2
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
  both.1 <- function(x,y){
    result <- list()
    criteria <- vector()
    for(i in 2:p){
      result[[i]] <- frm(y, x[,i,drop = FALSE], linkfrac="logit", table = FALSE)
      criteria[i] <- -2*ll_logit(result[[i]]$p, x[,c(1,i)], y)
    }
    min_criteria <- min(na.omit(criteria))
    index <- which.min(criteria)-1
    return(list(min = min_criteria, index = index))
  }
  first_step <- both.1(x,y)
  ind <- first_step$index + 1
  criteria_old <- first_step$min + k * 2
  m <- 2
  repeat{
    result_for <- list()
    criteria_for <- vector()
    result_back <- list()
    criteria_back <- vector()
    criteria_back_min <- vector()
    ind_for <- list()
    ind_back <- list()
    for(i in 2:p){
      if(!(i %in% ind)){
        result_for[[i]] <- frm(y, x[, sort(c(ind,i)), drop = FALSE], linkfrac = "logit", table = FALSE)
        criteria_for[i] <- -2*ll_logit(result_for[[i]]$p, x[,sort(c(1,ind,i))], y) + k*(m+1)
        ind_for[[i]] <- c(ind,i)
        l <- length(ind)
        if(l > 1){
          com <- combn(ind,l-1)
          for(j in 1:l){
            result_back[[j]] <- frm(y, x[, sort(c(com[,j],i)), drop = FALSE], linkfrac="logit", table = FALSE)
            criteria_back[j] <- -2*ll_logit(result_back[[j]]$p, x[,sort(c(1,com[,j],i))], y)+ k*m
          }
          criteria_back_min[i] <- min(na.omit(criteria_back))
          ind_back[[i]] <- c(list(com[, which.min(criteria_back)]),i)
          if(criteria_back_min[i] < criteria_for[i]){
            criteria_for[i] <- criteria_back_min[i]
            ind_for[[i]] <- ind_back[[i]]
          }
        }
      }
    }
    criteria_new <- min(na.omit(criteria_for))
    if(criteria_new >= criteria_old){
      break
    }
    ind <- ind_for[[which.min(criteria_for)]]
    criteria_old <- criteria_new
    m <- m + 1
  }
  min_criteria <- criteria_old
  index <- ind - 1
  variable <- colnames(x)[ind]
  return(list(min = min_criteria, index = index, variable = variable))
}
both(x, y, criterion = "AIC")
both(x, y, criterion = "BIC")


