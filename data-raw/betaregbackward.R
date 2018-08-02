backward <- function(x,y, link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                     criterion = c("AIC", "BIC")){
  p <- ncol(x)
  n <- nrow(x)
  if(criterion == "AIC"){
    k <- 2
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
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
  }else{
    index <- seq(1, p)[-ind]
    variable <- colnames(x)[index]
  }
  return(list(min = min_criteria, index = index, variable = variable))
}
backward(subx, suby, link = "logit", criterion = "AIC")
