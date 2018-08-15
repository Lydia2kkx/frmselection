#' fractional model selection function based on quasi-likelihood
#'
#' frmselect is used to select a formula-based model by information criterions of fractional regression models.
#'
#' @usage frmselect(x, y, criterion="AIC",linkfrac="logit", method="forward")
#'
#' @param x a numeric matrix, with column names, containing the values of the covariates.
#' @param y a numeric vector containing the values of the response variable. It should be between 0 and 1.
#' @param criterion model selection critetion. Available options: AIC, BIC, HQ.The default value is AIC.
#' @param linkfrac link function, Available options: logit, probit, loglog, cloglog, cauchit.The default value is logit.
#' @param method the mode of stepwise search and allsubsets. Available options: forward, backward, both, allsubsets.The default value is forward.
#'
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
#' x <- as.matrix(data.frame(d$x1, d$x2, d$x3))
#' y <- d$bnum
#' frmselect(x, y)

frmselect <- function(x,y, criterion = "AIC", linkfrac = "logit", method = "forward"){
  #The first part, check the conditions whether the given command is satisfied with the requirements.
  if(any(missing(x) || missing(y))){
    stop("Error: Missing data")
  }
  if(is.null(colnames(x))){
    stop("Error: Please name the columns of x")
  }
  if(class(x)!="matrix"){
    stop("Error: x should be a matrix")
  }
  if(length(y)!= nrow(x)){
    stop("Error: The length of x and y are different")
  }
  if(any(y>1 | y<0)){
    stop("Error: invalid dependent variable, all observations must be in [0, 1] ")
  }
  linkfun <- c("logit", "probit", "cloglog", "cauchit", "loglog")
  if (!linkfrac %in% linkfun){
    stop("The supported link functions are 'logit', 'probit', 'cloglog', 'loglog' and 'cauchit'!")
  }
  criterions <- c("AIC", "BIC", "HQ")
  if (!criterion %in% criterions){
    stop("The supported criterions are 'AIC', 'BIC' and 'HQ'!")
  }
  sup_methods <- c("forward", "backward", "both", "allsubsets")
  if (!method %in% sup_methods){
    stop("The supported methods are 'forward', 'backward', 'both' and 'allsubsets'!")
  }

  #The second part is to calculate the quasi-log-likelihood of supported link functions
  loglik <- function(params, x, y, linkfrac = c("logit", "probit", "loglog", "cloglog", "cauchit")){
    f <- NA
    for(i in 1:nrow(x)){
      if(linkfrac == "logit"){
        f[i] <- y[i]*log(exp(x[i,]%*%params)/(1+exp(x[i,]%*%params))) + (1-y[i])*log(1-exp(x[i,]%*%params)/(1+exp(x[i,]%*%params)))
      }
      if(linkfrac == "probit"){
        f[i] <- y[i]*log(pnorm(x[i,]%*%params)) + (1-y[i])*log(1-pnorm(x[i,]%*%params))
      }
      if(linkfrac == "loglog"){
        f[i] <- y[i]*log(exp(-exp(x[i,]%*%params))) + (1-y[i])*log(1-exp(-exp(x[i,]%*%params)))
      }
      if(linkfrac == "cloglog"){
        f[i] <- y[i]*log(1-exp(-exp(x[i,]%*%params))) + (1-y[i])*log(exp(-exp(x[i,]%*%params)))
      }
      if(linkfrac == "cauchit"){
        f[i] <- y[i]*log((atan(x[i,]%*%params)+pi/2)/pi) + (1-y[i])*log(1-(atan(x[i,]%*%params)+pi/2)/pi)
      }
    }
    return(sum(f))
  }

  #The third part, three information criterions.
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

  #The fourth part, four different methods.
  #forward stepwise
  forward <- function(x,y){
    x <- cbind(rep(1,nrow(x)), x)
    p <- ncol(x)
    forward.1 <- function(x,y){
      result <- list()
      criteria <- vector()
      for(i in 2:p){
        result[[i]] <- frm::frm(y, x[,i,drop = FALSE], linkfrac= linkfrac, table = FALSE)
        criteria[i] <- -2*loglik(result[[i]]$p, x[,c(1,i)], y, linkfrac= linkfrac)
      }
      min_criteria <- min(na.omit(criteria))
      index <- which.min(criteria)-1
      return(list(min_criteria = min_criteria, index = index))
    }
    first_step <- forward.1(x,y)
    ind <- first_step$index + 1
    criteria_old <- first_step$min + k * 2
    m <- 2
    repeat{
      result <- list()
      criteria <- vector()
      for(i in 2:p){
        if(!(i %in% ind)){
          result[[i]] <- frm::frm(y, x[, sort(c(ind,i)), drop = FALSE], linkfrac = linkfrac, table = FALSE)
          criteria[i] <- -2*loglik(result[[i]]$p, x[,sort(c(1,ind,i))], y, linkfrac= linkfrac) + k * (m+1)
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
    variable <- colnames(x)[ind]
    est <- frm::frm(y, x[, sort(ind), drop = FALSE], linkfrac = linkfrac, table = FALSE)
    coefficient <- est$p
    return(list(min_criteria = min_criteria, index = index, variable = variable, coefficient = coefficient))
  }

  #backward stepwise
  backward <- function(x,y){
    first_step <- frm::frm(y, x, linkfrac= linkfrac, table = FALSE)
    x <- cbind(rep(1,nrow(x)), x)
    p <- ncol(x)
    criteria_old <- -2*loglik(first_step$p, x, y, linkfrac= linkfrac) + k*p
    ind <- NULL
    m <- 1
    repeat{
      result <- list()
      criteria <- vector()
      for(i in 2:p){
        if(!(i %in% ind)){
          result[[i]] <- frm::frm(y, x[,-sort(c(1,ind,i)), drop = FALSE], linkfrac= linkfrac, table = FALSE)
          criteria[i] <- -2*loglik(result[[i]]$p, x[,-sort(c(ind,i))], y, linkfrac= linkfrac) + k*(p-m)
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
      est <- frm::frm(y, x[, -1, drop = FALSE], linkfrac = linkfrac, table = FALSE)
      variable <- colnames(x)
    }else{
      ind <- ind - 1
      index <- seq(1, p-1)[-ind]
      variable <- colnames(x)[index+1]
      est <- frm::frm(y, x[, sort(index + 1), drop = FALSE], linkfrac = linkfrac, table = FALSE)
    }
    coefficient <- est$p
    return(list(min_criteria = min_criteria, index = index, variable = variable, coefficient = coefficient))
  }

  #allsubsets
  allsubsets <- function(x, y){
    x <- cbind(rep(1,nrow(x)), x)
    p <- ncol(x)
    if(p-1 > 8){
      stop("Too many variables, using all-subsets method is time-consuming!")
    }
    criteria_min_part <- NULL
    ind <- list()
    for(i in 1:(p-1)){
      count <- choose(p-1, i)
      com <- combn(p-1, i)
      result <- list()
      criteria <- vector()
      for(j in 1:count){
        result[[j]] <- frm::frm(y, x[, com[,j]+1, drop = FALSE], linkfrac= linkfrac, table = FALSE)
        criteria[j] <- -2*loglik(result[[j]]$p, x[,(c(1,com[,j]+1))], y, linkfrac= linkfrac)+ k*(i+1)
      }
      criteria_min_part <- c(criteria_min_part, min(criteria))
      ind <- c(ind, list(com[, which.min(criteria)]))
    }
    min_criteria <- min(criteria_min_part)
    index <- ind[[which.min(criteria_min_part)]]
    variable <- colnames(x)[index+1]
    est <- frm::frm(y, x[, sort(index+1), drop = FALSE], linkfrac = linkfrac, table = FALSE)
    coefficient <- est$p
    return(list(min_criteria = min_criteria, index = index, variable = variable, coefficient = coefficient))
  }

  #both
  both <- function(x,y){
    x <- cbind(rep(1,nrow(x)), x)
    p <- ncol(x)
    both.1 <- function(x,y){
      result <- list()
      criteria <- vector()
      for(i in 2:p){
        result[[i]] <- frm::frm(y, x[,i,drop = FALSE], linkfrac= linkfrac, table = FALSE)
        criteria[i] <- -2*loglik(result[[i]]$p, x[,c(1,i)], y, linkfrac= linkfrac)
      }
      min_criteria <- min(na.omit(criteria))
      index <- which.min(criteria)-1
      return(list(min_criteria = min_criteria, index = index))
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
          result_for[[i]] <- frm::frm(y, x[, sort(c(ind,i)), drop = FALSE], linkfrac = linkfrac, table = FALSE)
          criteria_for[i] <- -2*loglik(result_for[[i]]$p, x[,sort(c(1,ind,i))], y, linkfrac= linkfrac) + k*(m+1)
          ind_for[[i]] <- c(ind,i)
          l <- length(ind)
          if(l > 1){
            com <- combn(ind,l-1)
            for(j in 1:l){
              result_back[[j]] <- frm::frm(y, x[, sort(c(com[,j],i)), drop = FALSE], linkfrac= linkfrac, table = FALSE)
              criteria_back[j] <- -2*loglik(result_back[[j]]$p, x[,sort(c(1,com[,j],i))], y, linkfrac= linkfrac)+ k*m
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
    est <- frm::frm(y, x[, sort(ind), drop = FALSE], linkfrac = linkfrac, table = FALSE)
    coefficient <- est$p
    return(list(min_criteria = min_criteria, index = index, variable = variable, coefficient = coefficient))
  }

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
  return(c(criterion = criterion, linkfrac = linkfrac, method = method, result))
}



