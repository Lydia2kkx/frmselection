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
  variable <- paste0("X", index)
  return(list(min = min_criteria, index = index, variable = variable))
}
allsubsets(x, y, criterion = "AIC")
allsubsets(x, y, criterion = "BIC")

#####################
choose(p-1, 2)
mm <- combn(8,3)
com <- combn(p-1, 2)
com
head(x[,com[,1]+1])
a <- min(3,2,5)
list(com[,1], mm[,2])[[sort(a)]]
aa <- NA
aa <- list(unlist(aa), c(1,2,3))
aa
aaa <- NULL
aaa <-c(aaa, c(1,2,3))
aaa

se<-c(1,2,3,4)
nl<-length(se)
res<-lapply(1:nl,function(i) combn(4,i))
res

ind <- list()
com <- combn(4,1)
ind <- c(ind, list(com[, 3]))
ind
com <- combn(4,2)
ind <- c(ind, list(com[, 3]))
ind
com <- combn(4,3)
ind <- c(ind, list(com[, 2]))
ind

ind <- list(rep(1,3),8,5,6,3,4,7)
min_criteria <- min(c(8,5,6,3,4,7))
xx <- which.min(c(8,5,6,3,4,7))
index <- ind[[xx]]
index
com[, which.min(xx)]
variable <- paste0("X", com[, which.min(xx)])

com <- combn(p-1, 3)
com

#######################
allsubsets <- function(x, y, criterion = c("AIC", "BIC")){
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
  subsets <- lapply(1:(p-1),function(i) combn((p-1),i))
  for(i in 1:(p-1)){
    result <- list()
    criteria <- vector()
    for(j in 1:subsets[[i]]){
      result[[j]] <- frm(y, x[, com[,j]+1, drop = FALSE], linkfrac="logit", table = FALSE)
      criteria[j] <- -2*ll_logit(result[[j]]$p, x[,(c(1,com[,j]+1))], y)+ k*(i+1)
    }
    criteria_min_part <- c(criteria_min_part, min(criteria))
    ind <- list(unlist(ind), com[, which.min(criteria)])
  }
  min_criteria <- min(criteria_min_part)
  index <- ind[[sort(min_criteria)]]
  variable <- paste0("X", index)
  return(list(min = min_criteria, index = index, variable = variable))
}
allsubsets(x, y, criterion = "AIC")
allsubsets(x, y, criterion = "BIC")
