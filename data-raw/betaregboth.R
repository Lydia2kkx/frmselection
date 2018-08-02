both <- function(x,y, link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                 criterion = c("AIC", "BIC")){
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
    for(i in 1:p){
      result[[i]] <- betareg::betareg(y ~ x[,i], link = link)
      criteria[i] <- -2 * result[[i]]$loglik + k*(1 + 2)
    }
    min_criteria <- min(na.omit(criteria))
    index <- which.min(criteria)
    return(list(min = min_criteria, index = index))
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
  index <- ind
  variable <- colnames(x)[index]
  return(list(min = min_criteria, index = index, variable = variable))
}
both(subx, suby, criterion = "AIC")
both(subx, suby, criterion = "BIC")
