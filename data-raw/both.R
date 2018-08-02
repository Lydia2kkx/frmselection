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
#Logit Model
ll_logit<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(exp(x[i,]%*%params)/(1+exp(x[i,]%*%params))) + (1-y[i])*log(1-exp(x[i,]%*%params)/(1+exp(x[i,]%*%params)))
  }
  return(sum(f))
}
##both
x <- cbind(rep(1,nrow(X)), X)
p <- ncol(x)
n <- nrow(x)
#In this procedure, you start with an empty model and build up sequentially just like in forward selection.
#The only caveat is that every time you add a new variable, Xnew,
#you have to check to see if any of the other variables that are already in the model should be dropped after Xnew is included.
#In this approach, you can end up searching "nonlinearly" through all the different models.

both <- function(x,y, criterion = c("AIC", "BIC")){
  if(criterion == "AIC"){
    k <- 2
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
  for(i in 2:p){
    m <- 1
    result <- list()
    criteria <- vector()
    result[[i]] <- frm(y, x[, i, drop = FALSE], linkfrac = "logit", table = FALSE)
    criteria[i] <- -2*ll_logit(result[[i]]$p, x[,c(1,i)], y) + k * m
  }
  criteria_old <- min(criteria)
  ind <- which.min(criteria)
  repeat{
    result_for <- list()
    criteria_for <- vector()
    result_back <- list()
    criteria_back <- vector()
    for(i in 2:p){
      if(!(i %in% ind)){
        result_for[[i]] <- frm(y, x[, sort(c(ind,i)), drop = FALSE], linkfrac = "logit", table = FALSE)
        criteria_for[i] <- -2*ll_logit(result_for[[i]]$p, x[,sort(c(1,ind,i))], y) + k * m
        for(j in c(ind, i)){
          if(is.na(j)){
            next
          }else{
            result_back[[j]] <- frm(y, x[,sort(c(ind,i)[-j]), drop = FALSE], linkfrac="logit", table = FALSE)
            criteria_back[j] <- -2*ll_logit(result_back[[j]]$p, x[,sort(c(ind,i)[-j])], y) + k*(p-m)
          }
        }
        criteria_back_min <- min(na.omit(criteria_back))
      }
    }
  }
}




x <- vector()
for(i in c(2,4)){
  x[i] <- i^2
}
x

j <- 2
ind <- c(6,3,7,2)
i <- 5
sort(c(ind,i)[-j])


x <- cbind(rep(1,nrow(x)), x)
p <- ncol(x)
n <- nrow(x)

both.1 <- function(x,y){
  x <- cbind(rep(1,nrow(x)), x)
  p <- ncol(x)
  result <- list()
  criteria <- vector()
  for(i in 2:p){
    result[[i]] <- frm(y, x[,i,drop = FALSE], linkfrac="logit", table = FALSE)
    criteria[i] <- -2*ll_logit(result[[i]]$p, x[,c(1,i)], y) + 4
  }
  min_criteria <- min(na.omit(criteria))
  index <- which.min(criteria)-1
  variable <- paste0("X", index)
  return( list(min = min_criteria, index = index, variable = variable))
}
both.1(x,y)

both.2 <- function(x,y){
  first_step <- both.1(x,y)
  ind <- first_step$index + 1
  x <- cbind(rep(1,nrow(x)), x)
  p <- ncol(x)
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
  return(list(min = min_criteria, index = index, variable = variable))
}
both.2(x,y)

#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
both.3 <- function(x,y){
  first_step <- both.2(x,y)
  ind <- first_step$index + 1
  x <- cbind(rep(1,nrow(x)), x)
  p <- ncol(x)
  result_for <- list()
  criteria_for <- vector()
  result_back <- list()
  criteria_back <- vector()
  criteria_back_min <- NULL
  ind_back_min <- vector()
  m <- 3
  k <- 2
  for(i in 2:p){
    if(!(i %in% ind)){
      result_for[[i]] <- frm(y, x[, sort(c(ind,i)), drop = FALSE], linkfrac = "logit", table = FALSE)
      criteria_for[i] <- -2*ll_logit(result_for[[i]]$p, x[,sort(c(1,ind,i))], y) + k * m
      for(j in sort(ind)){
        if(is.na(j)){
          next
        }else{
          ind_back <- matrix()
          result_back[[j]] <- frm(y, x[,sort(c(i,ind[-which(j==ind)])), drop = FALSE], linkfrac="logit", table = FALSE)
          criteria_back[j] <- -2*ll_logit(result_back[[j]]$p, x[,sort(c(1,ind[-which(j==ind)],i))], y) + k*(m-1)
          ind_back[j,] <- sort(c(i,ind[-which(j==ind)]))
        }
      }
      criteria_back_min[i] <- min(na.omit(criteria_back))
      ind_back_min[i,] <- ind_back[which(criteria_back_min[i]),]
    }
  }
  criteria_min <- min(na.omit(c(criteria_for, criteria_back_min)))
  return(criteria_min)
}
both.3(x,y)
#How to output the index of selected variables where we using both forward and backward method in the one loop!
traceback()


matrix(c(1,2,NA,NA), nrow = 2, byrow = TRUE)

head(x)
ind1 <- c(4,3,6)
com <- combn(ind1, 2)
length(ind1)
c(com[,2],5)
i <- 2
j <- 2
sort(c(i,ind1[-which(j==ind1)]))
sort(c(i,ind1[-which(j==ind1)]))
sort(c(1,ind1[-which(j==ind1)],i))
head(x[,sort(c(ind,i)[-j]), drop = FALSE])
sort(c(ind,i))
result_back[[2]] <- frm(y, x[,sort(c(ind,i))[-which(j==c(ind,i))], drop = FALSE], linkfrac="logit", table = FALSE)
result_back[[2]]$p
criteria_back[2] <- -2*ll_logit(result_back[[2]]$p, x[,sort(c(1,ind,i))[-which(j==c(ind,i))]], y) + 2
criteria_back[2]
both.1(x,y)
head(x[,sort(c(1,ind,i))[-j]])

v <- NULL
mm <- c(ind1, i)
for(j in mm){
  v[j] <- j^2
}
v

both.3 <- function(x,y){
  first_step <- both.2(x,y)
  ind <- first_step$index + 1
  result_for <- list()
  criteria_for <- vector()
  result_back <- list()
  criteria_back <- vector()
  criteria_back_min <- vector()
  ind_for <- list()
  ind_back <- list()
  m <- 3
  k <- 2
  for(i in 2:p){
    if(!(i %in% ind)){
      result_for[[i]] <- frm(y, x[, sort(c(ind,i)), drop = FALSE], linkfrac = "logit", table = FALSE)
      criteria_for[i] <- -2*ll_logit(result_for[[i]]$p, x[,sort(c(1,ind,i))], y) + k * m
      ind_for[[i]] <- c(ind,i)
      n <- length(ind)
      if(n <=1){
        next
      }else{
        com <- combn(ind,n-1)
        for(j in 1:n){
          result_back[[j]] <- frm(y, x[, sort(c(com[,j],i)), drop = FALSE], linkfrac="logit", table = FALSE)
          criteria_back[j] <- -2*ll_logit(result_back[[j]]$p, x[,sort(c(1,com[,j],i))], y)+ k*(m-1)
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
  criteria_min <- min(na.omit(criteria_for))
  ind <- ind_for[[which.min(criteria_for)]]-1
  variable <- paste0("X", ind)
  return(list(min = criteria_min, index = ind, variable = variable))
}

both.3(x,y)

i <- 2
ind <- c(4,1,3)
n <- length(ind)
com <- combn(ind,n-1)
j <- 3
sort(c(com[,j],i))
