both <- function(x,y, criterion = c("AIC", "BIC")){
  x <- cbind(rep(1,nrow(x)), x)
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
    for(i in 2:p){
      result[[i]] <- frm(y, x[,i,drop = FALSE], linkfrac="logit", table = FALSE)
      criteria[i] <- -2*ll_logit(result[[i]]$p, x[,c(1,i)], y)
    }
    min_criteria <- min(na.omit(criteria))
    index <- which.min(criteria)-1
    return(list(min = min_criteria, index = index))
  }
  first_step <- both.1(x,y)
  ind <- first_step$index + 1
  criteria_old <- first_step$min + k * 2
  m <- 2
  repeat{
    result_for <- list()
    criteria_for <- vector()
    result_back <- list()
    criteria_back <- vector()
    criteria_back_min <- vector()
    ind_for <- list()
    ind_back <- list()
    for(i in 2:p){
      if(!(i %in% ind)){
        result_for[[i]] <- frm(y, x[, sort(c(ind,i)), drop = FALSE], linkfrac = "logit", table = FALSE)
        criteria_for[i] <- -2*ll_logit(result_for[[i]]$p, x[,sort(c(1,ind,i))], y) + k*(m+1)
        ind_for[[i]] <- c(ind,i)
        l <- length(ind)
        if(l > 1){
          com <- combn(ind,l-1)
          for(j in 1:l){
            result_back[[j]] <- frm(y, x[, sort(c(com[,j],i)), drop = FALSE], linkfrac="logit", table = FALSE)
            criteria_back[j] <- -2*ll_logit(result_back[[j]]$p, x[,sort(c(1,com[,j],i))], y)+ k*m
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
  index <- ind - 1
  variable <- colnames(x)[ind]
  return(list(min = min_criteria, index = index, variable = variable))
}
both(x, y, criterion = "AIC")
both(x, y, criterion = "BIC")
