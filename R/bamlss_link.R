#' Bamlss quasi-binomial family
#'
#'This is a quasi-binomial family which is suitable for bamlss() function in bamlss package.
#'This quasi-binomial family only contains two link functions, which are used mostly in the economics,
#'logit and probit function.
#'
#' @param link link function, Available options: logit, probit
#' @param ... Arguments to pass to bamlss
#'
#' @details It is a family of bamlss package, can be used in the bamlss() function. The data type is
#' the fractional data which is between 0 and 1. Especially for probit link, the response data should be
#' strictly smaller than 1, otherwise it might be a warning that the backfitting algorithm did not converge!
#'
#' @export
#'
#' @examples
#' library("bamlss")
#' set.seed(123)
#' d <- GAMart()
#' b <- bamlss(bnum ~ x1 + x2 + x3, data = d, family = frm_bamlss(link = "logit") ,sampler = FALSE, multiple = FALSE)



frm_bamlss <- function(link = "logit", ...)
{
  linkfun <- c("logit", "probit")
  if (!link %in% linkfun){
   stop("The supported link functions are 'logit' and 'probit'!")
  }
  #link probit
  binomial2_bamlss <- function(...)
  {
    rval <- list(
      "family" = "quasibinomial",
      "names" = "pi",
      "links" = c(pi = "probit"),
      "valid.response" = function(x) {
        if(!is.factor(x)) {
          if(any(x > 1 || x < 0))
            stop("response should between 0 and 1!", call. = FALSE)
        } else {
          stop("response should not be factor levels!", call. = FALSE)
        }
        TRUE
      },
      "bayesx" = list(
        "pi" = c("quasibinomial_probit", "mu")
      ),
      "bugs" = list(
        "dist" = "dbern",
        "eta" = BUGSeta,
        "model" = BUGSmodel
      ),
      "mu" = function(par, ...) {
        par$pi
      },
      "d" = function(y, par, log = FALSE) {
        dnorm(y)
      },
      "p" = function(y, par, ...) {
        pnorm(y)
      },
      "r" = function(n, par) {
        rnorm(n)
      },
      "score" = list(
        "mu" = function(y, par, ...) {
          y - par$pi
        }
      ),
      "hess" = list(
        "pi" = function(y, par, ...) {
          par$pi * (1 - par$pi)
        }
      )
    )
    class(rval) <- "family.bamlss"
    rval
  }
  if(link != "logit"){
      return(binomial2_bamlss(link, ...))
  }

  rval <- list(
    "family" = "quasibinomial",
    "names" = "pi",
    "links" = c(pi = "logit"),
    "valid.response" = function(x) {
      if(!is.factor(x)) {
        if(any(x > 1 || x < 0))
          stop("response should between 0 and 1!", call. = FALSE)
      } else {
        stop("response should not be factor levels!", call. = FALSE)
      }
      TRUE
    },
    "bayesx" = list(
      "y" = c("binomial_logit", "mu")
    ),
    "bugs" = list(
      "dist" = "dbern",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "mu" = function(par, ...) {
      par$pi
    },
    "d" = function(y, par, log = FALSE) {
      choose(1, round(y))*par$pi^y*(1-par$pi)^(1-round(y))
    },
    "p" = function(y, par, ...) {
      pbinom(y, size = 1, prob = par$pi, ...)
    },
    "r" = function(n, par) {
      rbinom(n, size = 1, prob = par$pi)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        y - par$pi
      }
    ),
    "hess" = list(
      "pi" = function(y, par, ...) {
        par$pi * (1 - par$pi)
      }
    )
  )
  class(rval) <- "family.bamlss"
  rval
}




