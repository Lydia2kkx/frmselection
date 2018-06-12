#allsubsets
allsubsets <- function(x, y, link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                       criterion = c("AIC", "BIC")){
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
  for(i in 1:p){
    count <- choose(p, i)
    com <- combn(p, i)
    result <- list()
    criteria <- vector()
    for(j in 1:count){
      result[[j]] <- betareg::betareg(y ~ x[, com[,j]], link = link)
      criteria[j] <- -2 * result[[j]]$loglik + k * (i + 2)
    }
    criteria_min_part <- c(criteria_min_part, min(criteria))
    ind <- c(ind, list(com[, which.min(criteria)]))
  }
  min_criteria <- min(criteria_min_part)
  index <- ind[[which.min(criteria_min_part)]]
  variable <- colnames(x)[index]
  return(list(min = min_criteria, index = index, variable = variable))
}
allsubsets(subx, suby, criterion = "AIC")
allsubsets(subx, suby, criterion = "BIC")
