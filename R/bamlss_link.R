#link logit
frm_bamlss <- function(link = c("logit", "probit"), ...)
{
  if(link != "logit")
    return(binomial2_bamlss(link, ...))

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


#link loglog, but it doesn't work
binomial2_bamlss <- function(...)
{
    rval <- list(
      "family" = "quasibinomial",
      "names" = "pi",
      "links" = c(pi = "loglog"),
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
        "pi" = c("quasibinomial_loglog", "mu")
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
        pevd(y)
      },
      "r" = function(n, par) {
        revd(n)
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

