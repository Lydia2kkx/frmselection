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
ll_logit<-function(params,x,y){
  f <- NA
  for(i in 1:nrow(x)){
    f[i] <- y[i]*log(exp(x[i,]%*%params)/(1+exp(x[i,]%*%params))) + (1-y[i])*log(1-exp(x[i,]%*%params)/(1+exp(x[i,]%*%params)))
  }
  return(sum(f))
}
###All-subsets selection

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
