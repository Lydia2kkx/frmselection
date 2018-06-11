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

forward <- function(x,y, criterion = c("AIC", "BIC")){
  if(criterion == "AIC"){
    k <- 2
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
  ind <- NULL
  m <- 1
  result_old <- frm(y, x[,2,drop = FALSE], linkfrac="logit", table = FALSE)
  criteria_old <- -2*ll_logit(result_old$p, x[,sort(c(1:2))], y) + k * m
  repeat{
    result <- list()
    criteria <- vector()
    for(i in 2:p){
      if(!(i %in% ind)){
        result[[i]] <- frm(y, x[, sort(c(ind,i)), drop = FALSE], linkfrac = "logit", table = FALSE)
        criteria[i] <- -2*ll_logit(result[[i]]$p, x[,sort(c(1,ind,i))], y) + k * m
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
  variable <- paste0("X", index)
  return(list(min = min_criteria, index = index, variable = variable))
}

forward(x, y, criterion = "AIC")
forward(x, y, criterion = "BIC")

