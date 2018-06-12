library(frm)
set.seed(123)
N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2","X3")

ym <- exp(X[,1]+X[,2]+X[,3]+u)/(1+exp(X[,1]+X[,2]+X[,3]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1
#frm estimation of a logit fractional regression model
a<-frm(y, X, linkfrac="logit")

#Logit model
ll_logit<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(exp(x[i,]%*%params)/(1+exp(x[i,]%*%params))) + (1-y[i])*log(1-exp(x[i,]%*%params)/(1+exp(x[i,]%*%params)))
  }
  return(sum(f))
}

##backwand stepwise
x <- cbind(rep(1,nrow(X)), X)
p <- ncol(x)
n <- nrow(x)

backward.1 <- function(x, y){
  first_step <- frm(y, X, linkfrac="logit", table = FALSE)
  criteria <- -2*ll_logit(first_step$p, x, y) + 2*p
  return(criteria)
}
backward.1(x,y)
##########
head(x)
head(x[,-c(1,2)])
##########
backward.2 <- function(x, y){
  first_step <- backward.1(x, y)
  result <- list()
  criteria <- vector()
  for(i in 2:p){
    result[[i]] <- frm(y, x[,-c(1,i)], linkfrac="logit", table = FALSE)
    criteria[i] <- -2*ll_logit(result[[i]]$p, x[,-i], y)
  }
  min_criteria <- min(na.omit(criteria))
  index <- which.min(criteria)-1
  ind <- seq(1,p-1)[-index]
  variable <- paste0("X", ind)
  return(list(min = min_criteria, index = index, variable = variable))
}
backward.2(x, y)

################
index <- 2
c(1:4)[-index]
#################
a1 <- frm(y, x[,c(2,3)], linkfrac="logit", table = FALSE)
criteria1 <- -2*ll_logit(a1$p, x[,1:3], y)
criteria1

a2 <- frm(y, x[,c(2,4)], linkfrac="logit", table = FALSE)
criteria2 <- -2*ll_logit(a1$p, x[,c(1:2,4)], y)
criteria2

a3 <- frm(y, x[,c(3,4)], linkfrac="logit", table = FALSE)
criteria3 <- -2*ll_logit(a1$p, x[,c(1,3:4)], y)
criteria3

-2*ll_logit(a$p, x, y) + 2*p
#######################
x[,sort(c(3,1))]
x[,-sort(c(3,1))]
#######################
backward.3 <- function(x, y){
  first_step <- backward.2(x, y)
  ind <- first_step$index + 1
  result <- list()
  criteria <- vector()
  for(i in 2:p){
    if(!(i %in% ind)){
      result[[i]] <- frm(y, x[,-sort(c(1,ind,i)), drop = FALSE], linkfrac="logit", table = FALSE)
      criteria[i] <- -2*ll_logit(result[[i]]$p, x[,-sort(c(ind,i))], y)
    }
  }
  min_criteria <- min(na.omit(criteria))
  index.new <- which.min(criteria)-1
  index <- c(first_step$index, index.new)
  ind <- seq(1,p-1)[-index]
  variable <- paste0("X", ind)
  return(list(min = min_criteria, index = index, variable = variable))
}
backward.3(x, y)


backward <- function(x,y, criterion = c("AIC", "BIC")){
  k <- NULL
  if(criterion == "AIC"){
    k <- 2
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
  first_step <- frm(y, X, linkfrac="logit", table = FALSE)
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
    variable <- paste0("X", 1:p)
  }else{
    index <- ind - 1
    ind <- seq(1, p-1)[-index]
    variable <- paste0("X", ind)
  }
  return(list(min = min_criteria, index = index, variable = variable))
}
backward(x, y, criterion = "BIC")
#####################
seq(1, p-1)[--1]
ind <- seq(1, p-1)[-NULL]
paste0("X", 1:p)
#######################
