setwd("E:/Master/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
newdat <- data.frame(y, x)
subdat <- newdat[y<1,]
suby <- subdat[,1]
subx <- subdat[,-1]
subx <- as.matrix(subdat[,-1])

betaselect <- function(x, y, criterion = c("AIC", "BIC", "HQ"),
                       link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                       method = c("forward", "backward", "both", "allsubsets")){
  if(any(missing(x) || missing(y))){
    stop("Error: Missing data")
  }
  if(is.null(colnames(x))){
    stop("Error: Please name the columns of x")
  }
  if(class(x) == "data.frame"){
    x <- as.matrix(x)
  }
  if(class(x)!= "matrix" ){
    stop("Error: x should be a matrix or data frame")
  }
  if(length(y)!= nrow(x)){
    stop("Error: The length of x and y are different")
  }
  if(missing(link)){
    link <- "logit"
  }
  p <- ncol(x)
  n <- nrow(x)
  if(any(criterion == "AIC" || missing(criterion))){
    k <- 2
    criterion <- "AIC"
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
  if(criterion == "HQ"){
    k <- log(log(n))
  }
  #forward stepwise
  forward <- function(x,y){
    forward.1 <- function(x,y){
      result <- list()
      criteria <- vector()
      for(i in 1:p){
        result[[i]] <- betareg::betareg(y ~ x[,i], link = link)
        criteria[i] <- -2 * result[[i]]$loglik + k*(1 + 2)
      }
      min_criteria <- min(criteria)
      index <- which.min(criteria)
      variable <- colnames(x)[index]
      return( list(criteria = criteria, min = min_criteria, index = index, variable = variable))
    }
    first_step <- forward.1(x,y)
    ind <- first_step$index
    criteria_old <- first_step$min
    m <- 3
    repeat{
      result <- list()
      criteria <- vector()
      for(i in 1:p){
        if(!(i %in% ind)){
          result[[i]] <- betareg::betareg(y ~ x[,sort(c(ind,i))], link = link)
          criteria[i] <- -2* result[[i]]$loglik + k * (m + 1)
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
    variable <- colnames(x)[ind]
    index <- sort(ind)
    est <- betareg::betareg(y ~ x[,index], link = link)
    coefficient <- coef(est)
    return(list(min_criteria = min_criteria, index = ind, variable = variable, coefficient = coefficient))
  }
  #backward stepwise
  backward <- function(x,y){
    first_step <- betareg::betareg(y ~ x, link = link)
    criteria_old <- -2 * first_step$loglik + k * (p + 2)
    ind <- NULL
    m <- 1
    repeat{
      result <- list()
      criteria <- vector()
      for(i in 1:p){
        if(!(i %in% ind)){
          result[[i]] <- betareg::betareg(y ~ x[,-sort(c(ind,i))], link = link)
          criteria[i] <- -2*result[[i]]$loglik + k * (p + 2 - m)
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
      est <- betareg::betareg(y ~ x, link = link)
      coefficient <- coef(est)
    }else{
      index <- seq(1, p)[-ind]
      variable <- colnames(x)[index]
      est <- betareg::betareg(y ~ x[,index], link = link)
      coefficient <- coef(est)
    }
    return(list(min_criteria = min_criteria, index = index, variable = variable, coefficient = coefficient))
  }
  #both
  both <- function(x,y){
    both.1 <- function(x,y){
      result <- list()
      criteria <- vector()
      for(i in 1:p){
        result[[i]] <- betareg::betareg(y ~ x[,i], link = link)
        criteria[i] <- -2 * result[[i]]$loglik + k*(1 + 2)
      }
      min_criteria <- min(na.omit(criteria))
      index <- which.min(criteria)
      return(list(min_criteria = min_criteria, index = index))
    }
    first_step <- both.1(x,y)
    ind <- first_step$index
    criteria_old <- first_step$min
    m <- 3
    repeat{
      result_for <- list()
      criteria_for <- vector()
      result_back <- list()
      criteria_back <- vector()
      criteria_back_min <- vector()
      ind_for <- list()
      ind_back <- list()
      for(i in 1:p){
        if(!(i %in% ind)){
          result_for[[i]] <- betareg::betareg(y ~ x[,sort(c(ind,i))], link = link)
          criteria_for[i] <- -2* result_for[[i]]$loglik + k * (m + 1)
          ind_for[[i]] <- c(ind,i)
          l <- length(ind)
          if(l > 1){
            com <- combn(ind,l-1)
            for(j in 1:l){
              result_back[[j]] <- betareg::betareg(y ~ x[, sort(c(com[,j],i))], link = link)
              criteria_back[j] <- -2*result_back[[j]]$loglik + k * m
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
    variable <- colnames(x)[ind]
    index <- sort(ind)
    est <- betareg::betareg(y ~ x[,index], link = link)
    coefficient <- coef(est)
    return(list(min_criteria = min_criteria, index = ind, variable = variable, coefficient = coefficient))
  }
  #allsubsets
  allsubsets <- function(x, y){
    if(p-1 > 8){
      stop("Too many variables, using all-subsets method is time-consuming!")
    }
    criteria_min_part <- NULL
    ind <- list()
    for(i in 1:p){
      count <- choose(p, i)
      com <- combn(p, i)
      result <- list()
      criteria <- vector()
      for(j in 1:count){
        result[[j]] <- betareg::betareg(y ~ x[, com[,j]], link = link)
        criteria[j] <- -2 * result[[j]]$loglik + k * (i + 2)
      }
      criteria_min_part <- c(criteria_min_part, min(criteria))
      ind <- c(ind, list(com[, which.min(criteria)]))
    }
    min_criteria <- min(criteria_min_part)
    index <- ind[[which.min(criteria_min_part)]]
    variable <- colnames(x)[index]
    est <- betareg::betareg(y ~ x[,index], link = link)
    coefficient <- coef(est)
    return(list(min_criteria = min_criteria, index = index, variable = variable, coefficient = coefficient))
  }
  if(any(method == "forward" || missing(method))){
    result <- forward(x, y)
    method <- "forward"
  }
  if(method == "backward"){
    result <- backward(x, y)
  }
  if(method == "allsubsets"){
    result <- allsubsets(x, y)
  }
  if(method == "both"){
    result <- both(x, y)
  }
  return(c(criterion = criterion, link = link, method = method, result))
}
betaselect(subx, suby)
betaselect(subx, suby, criterion = "HQ", method = "backward", link = "probit")
betaselect(subx, suby, criterion = "BIC", method = "both")
is.null(colnames(subx))
colnames(subx)

p <- ncol(subx)
n <-nrow(subx)
criteria <- function(x, y){
  if(p > 8){
    stop("Too many variables, find the minimal criterions is time-consuming!")
  }
  criterion <- NULL
  ind <- list()
  for(i in 1:p){
    count <- choose(p, i)
    com <- combn(p, i)
    result <- list()
    criteria <- vector()
    for(j in 1:count){
      result[[j]] <- betareg::betareg(y ~ x[, com[,j]])
      criteria[j] <- -2 * result[[j]]$loglik
    }
    criterion[i] <- min(criteria)
  }
  return(min_criteria = criterion)
}
num <- 1:p
cri <- criteria(subx, suby)
criteria_aic <- cri + 2 * (num + 2)
criteria_bic <- cri + log(n) * (num + 2)
criteria_hq <- cri + log(log(n)) * (num + 2)
plot(num, criteria_aic, type = "b", lwd = 2, main = "Three Criterions of Numbers of Variables ",
     xlab = "the number of variables", ylab = "criterion",
     ylim = c(min(criteria_aic), max(criteria_bic)))
lines(num, criteria_bic, type = "b", lwd = 2, col = "red")
lines(num, criteria_hq, type = "b", lwd = 2, col = "blue")
legend("topright", legend = c("AIC", "BIC", "HQ"), col = c("black", "red", "blue"), lwd = 2)


