g <- betareg::betareg(suby ~ mrate + lemp + lemp2 + age + age2 + sole, data = subx)

setwd("E:/Master/Programming with R")
dat <- read.table("401kjae.txt")
y <- dat[,1]
x <- cbind(dat[,2], log(dat[,3]), log(dat[,3])^2, dat[,4], (dat[,4])^2, dat[,5])
#colnames(x) <- c("mrate", "lemp", "lemp2","age", "age2","sole")
#newdat <- cbind(y, x)
newdat <- as.data.frame(cbind(y, x))
subdat <- newdat[y<1,]
suby <- subdat[,1]
subx <- subdat[,-1]
subx <- as.matrix(subdat[,-1])
gg <- betareg::betareg(suby ~ subx)
gg$coefficients
coef(gg)
gg$loglik

g <- betareg::betareg(suby ~ subx[,c(1,3,5)])
betareg::betareg(suby ~ subx[,1], data = subx)  #这样会报错，但是去掉data=subx就不会报错
betareg::betareg(suby ~ subx[,1])
betareg::betareg(suby ~ subx[,-sort(c(3,1,4))])

forward <- function(x,y, link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                    criterion = c("AIC", "BIC")){
  n <- nrow(x)
  p <- ncol(x)
  if(criterion == "AIC"){
    k <- 2
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
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
  return(list(min = min_criteria, index = ind, variable = variable))
}
forward(subx, suby, criterion = "AIC")
forward(subx, suby, link = "logit", criterion = "BIC")
forward(subx, suby, link = "probit", criterion = "AIC")

ind <- c(1,2)
xi <- cbind(subx[,3], subx[,1])
class(xi)
betareg::betareg(suby ~ xi + subx[,4])

forward.1 <- function(x,y){
  p <- ncol(x)
  result <- list()
  criteria <- vector()
  for(i in 1:p){
    result[[i]] <- betareg::betareg(y ~ x[,i], data = x)
    criteria[i] <- -2* result[[i]]$loglik + 2*(result[[i]]$df.null - result[[i]]$df.residual + 2)
  }
  min_criteria <- min(criteria)
  index <- which.min(criteria)
  variable <- colnames(x)[index]
  return( list(criteria = criteria, min = min_criteria, index = index, variable = variable))
}
forward.1(subx1, suby1)

result[1] <- betareg::betareg(y ~ x[,1], data = x)
