#' fractional model selection function based on beta regression
#'
#' betaselect is used to select a formula-based model by information criterions of fractional regression models.
#'
#' @usage betaselect(x, y, criterion, link, method)
#'
#' @param x a numeric matrix or data frame, with column names, containing the values of the covariates.
#' @param y a numeric vector containing the values of the response variable. It should be between 0 and 1.
#' @param criterion model selection critetion. Available options: AIC, BIC, HQ.The default value is AIC.
#' @param link link function, Available options: logit, probit, loglog, cloglog, log, cauchit.The default value is logit.
#' @param method the mode of stepwise search and allsubsets. Available options: forward, backward, both, allsubsets.The default value is forward.
#'
#' @return A list contains information criterion, link function, model selection method,
#' minimal value of information criterion, the order of variables which are chosen in the model,
#' the names of corresponding variables and the estimated coefficients of them.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' library(bamlss)
#' d <- GAMart()
#' x <- data.frame(d$x1, d$x2, d$x3)
#' y <- d$bnum
#' betaselect(x, y)


betaselect <- function(x, y, criterion = "AIC",link = "logit", method = "forward"){
  #In the first part, check the conditions whether the given command is satisfied with the requirements.
  if(any(missing(x) || missing(y))){
    stop("Error: Missing data")
  }
  if(is.null(colnames(x))){
    stop("Error: Please name the columns of x")
  }
  if(class(x) == "data.frame"){
    x <- as.matrix(x)
  }
  if(class(x)!= "matrix" ){
    stop("Error: x should be a matrix or data frame")
  }
  if(length(y)!= nrow(x)){
    stop("Error: The length of x and y are different")
  }
  if(any(y>=1 | y<=0)){
    stop("Error: invalid dependent variable, all observations must be in (0, 1) ")
  }
  linkfun <- c("logit", "probit", "cloglog", "cauchit", "log", "loglog")
  if (!link %in% linkfun){
    stop("The supported link functions are 'logit', 'probit', 'cloglog', 'loglog', 'cauchit' and 'log'!")
  }
  criterions <- c("AIC", "BIC", "HQ")
  if (!criterion %in% criterions){
    stop("The supported criterions are 'AIC', 'BIC' and 'HQ'!")
  }
  sup_methods <- c("forward", "backward", "both", "allsubsets")
  if (!method %in% sup_methods){
    stop("The supported methods are 'forward', 'backward', 'both' and 'allsubsets'!")
  }

  #The second part, three information criterions.
  p <- ncol(x)
  n <- nrow(x)
  if(criterion == "AIC"){
    k <- 2
    criterion <- "AIC"
  }
  if(criterion == "BIC"){
    k <- log(n)
  }
  if(criterion == "HQ"){
    k <- log(log(n))
  }

  #The third part, four different methods.
  #forward stepwise
  forward <- function(x,y){
    forward.1 <- function(x,y){
      result <- list()
      criteria <- numeric()
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
      criteria <- numeric()
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
    index <- sort(ind)
    est <- betareg::betareg(y ~ x[,index], link = link)
    coefficient <- coef(est)
    return(list(min_criteria = min_criteria, index = ind, variable = variable, coefficient = coefficient))
  }

  #backward stepwise
  backward <- function(x,y){
    first_step <- betareg::betareg(y ~ x, link = link)
    criteria_old <- -2 * first_step$loglik + k * (p + 2)
    ind <- NULL
    m <- 1
    repeat{
      result <- list()
      criteria <- numeric()
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
      est <- betareg::betareg(y ~ x, link = link)
      coefficient <- coef(est)
    }else{
      index <- seq(1, p)[-ind]
      variable <- colnames(x)[index]
      est <- betareg::betareg(y ~ x[,index], link = link)
      coefficient <- coef(est)
    }
    return(list(min_criteria = min_criteria, index = index, variable = variable, coefficient = coefficient))
  }

  #both
  both <- function(x,y){
    both.1 <- function(x,y){
      result <- list()
      criteria <- numeric()
      for(i in 1:p){
        result[[i]] <- betareg::betareg(y ~ x[,i], link = link)
        criteria[i] <- -2 * result[[i]]$loglik + k*(1 + 2)
      }
      min_criteria <- min(na.omit(criteria))
      index <- which.min(criteria)
      return(list(min_criteria = min_criteria, index = index))
    }
    first_step <- both.1(x,y)
    ind <- first_step$index
    criteria_old <- first_step$min
    m <- 3
    repeat{
      result_for <- list()
      criteria_for <- numeric()
      result_back <- list()
      criteria_back <- numeric()
      criteria_back_min <- numeric()
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
    variable <- colnames(x)[ind]
    index <- sort(ind)
    est <- betareg::betareg(y ~ x[,index], link = link)
    coefficient <- coef(est)
    return(list(min_criteria = min_criteria, index = ind, variable = variable, coefficient = coefficient))
  }

  #allsubsets
  allsubsets <- function(x, y){
    if(p-1 > 8){
      stop("Too many variables, using all-subsets method is time-consuming!")
    }
    criteria_min_part <- NULL
    ind <- list()
    for(i in 1:p){
      count <- choose(p, i)
      com <- combn(p, i)
      result <- list()
      criteria <- numeric()
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
    est <- betareg::betareg(y ~ x[,index], link = link)
    coefficient <- coef(est)
    return(list(min_criteria = min_criteria, index = index, variable = variable, coefficient = coefficient))
  }

  #The default value of method is forward stepwise.
  if(method == "forward"){
    result <- forward(x, y)
    method <- "forward"
  }
  if(method == "backward"){
    result <- backward(x, y)
  }
  if(method == "allsubsets"){
    result <- allsubsets(x, y)
  }
  if(method == "both"){
    result <- both(x, y)
  }

  #Return a list of results.
  return(c(criterion = criterion, link = link, method = method, result))
}
