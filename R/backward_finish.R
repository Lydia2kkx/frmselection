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
backward <- function(x,y, criterion = c("AIC", "BIC")){
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

