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
##forward stepwise
x <- cbind(rep(1,nrow(X)), X)
p <- ncol(x)
n <- nrow(x)
##k=1  logit model
forward.1 <- function(x,y){
  result <- list()
  criteria <- vector()
  for(i in 2:p){
    result[[i]] <- frm(y, x[,i,drop = FALSE], linkfrac="logit", table = FALSE)
    criteria[i] <- -ll_logit(result[[i]]$p, x[,c(1,i)], y)
  }
  min_criteria <- min(na.omit(criteria))
  index <- which.min(criteria)-1
  variable <- paste0("X", index)
  return( list(min = min_criteria, index = index, variable = variable))
}

aa <- forward.1(x,y)
aa

#####
head(x)
colnames(x)
x1<-x[,2,drop = FALSE]
colnames(x1)
head(x1)
#####
forward.2 <- function(x,y){
  first_step <- forward.1(x,y)
  ind <- first_step$index + 1
  result <- list()
  criteria <- vector()
  for(i in 2:p){
    if( i != ind){
      result[[i]] <- frm(y, x[,sort(c(ind,i))], linkfrac="logit", table = FALSE)
      criteria[i] <- -2*ll_logit(result[[i]]$p, x[,sort(c(1,ind,i))], y)
    }
  }
  min_criteria <- min(na.omit(criteria))
  index.new <- which.min(criteria)-1
  index <- c(first_step$index, index.new)
  variable <- paste0("X", index)
  #return(result$p)
  return(list(min = min_criteria, index = index, variable = variable))
}
forward.2(x,y)

ind <- 3
aaa <- x[,sort(c(1,ind,2))]
head(aaa)


forward.3 <- function(x,y){
  first_step <- forward.2(x,y)
  ind <- first_step$index + 1
  result <- list()
  criteria <- vector()
  for(i in 2:p){
    if(!(i %in% ind)){
      result[[i]] <- frm(y, x[,sort(c(ind,i))], linkfrac="logit", table = FALSE)
      criteria[i] <- -2*ll_logit(result[[i]]$p, x[,sort(c(1,ind,i))], y)
    }
  }
  min_criteria <- min(na.omit(criteria))
  index.new <- which.min(criteria)-1
  index <- c(first_step$index, index.new)
  variable <- paste0("X", index)
  return(list(min = min_criteria, index = index, variable = variable))
}
forward.3(x,y)

#########
num <- 3
num %in% ind
!(num %in% ind)
###########
forward <- function(x,y, criterion = c("AIC", "BIC")){
  first_step <- forward.1(x,y)
  ind <- first_step$index + 1
  criteria_old <- first_step$min
  repeat{
    result <- list()
    criteria <- vector()
    for(i in 2:p){
      if(!(i %in% ind)){
        result[[i]] <- frm(y, x[,sort(c(ind,i))], linkfrac="logit", table = FALSE)
        criteria[i] <- -2*ll_logit(result[[i]]$p, x[,sort(c(1,ind,i))], y)
      }
    }
    criteria_new <- min(na.omit(criteria))
    if(criterion == "AIC"){
      if(criteria_new - criteria_old >= 2){
        break
      }
    }
    if(criterion == "BIC"){
      if(criteria_new - criteria_old >= log(n)){
        break
      }
    }
    ind <- c(ind, which.min(criteria))
    criteria_old <- criteria_new
  }
  min_criteria <- criteria_old
  index <- ind - 1
  variable <- paste0("X", index)
  return(list(min = min_criteria, index = index, variable = variable))
}
forward(x, y, criterion = "AIC")

log(n)

setwd("E:/研究生/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), (log(dat[,3]))^2, dat[,4], (dat[,4])^2, dat[,5])
head(x)
colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
head(x)
a<-frm(y, X, linkfrac="logit")
