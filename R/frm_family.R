process_factor_response <- function(x)
{
  if(!is.null(attr(x, "contrasts"))) {
    if(!is.null(dim(x)))
      x <- x[, ncol(x)]
  }
  if(is.factor(x))
    x <- as.integer(x) - 1L
  as.integer(x)
}

frm_bamlss <- function(link = "logit", ...)
{
  if(link != "logit")
    return(binomial2_bamlss(...))

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
      #y <- process_factor_response(y)
      choose(1, round(y))*par$pi^y*(1-par$pi)^(1-round(y))
    },
    "p" = function(y, par, ...) {
      #y <- process_factor_response(y)
      pbinom(y, size = 1, prob = par$pi, ...)
    },
    "r" = function(n, par) {
      rbinom(n, size = 1, prob = par$pi)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        #y <- process_factor_response(y)
        y - par$pi
      }
    ),
    "hess" = list(
      "pi" = function(y, par, ...) {
        par$pi * (1 - par$pi)
      }
    ),
    "initialize" = list("pi" = function(y, ...) {
      y <- process_factor_response(y)
      (y + 0.5) / 2
    })
  )
  class(rval) <- "family.bamlss"
  rval
}

frm_bamlss <- function(link = "logit", ...)
{
  if(link != "logit")
    return(binomial2_bamlss(...))

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
      "y" = c("quasibinomial_logit", "mu")
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

library("bamlss")
#library("stats")
set.seed(123)
d <- GAMart()
head(d)
#d$id
#f <- bnum ~ la(~x1+x2+x3)+la(id)
b <- bamlss(bnum ~ x1 + x2 + x3, data = d, family = "frm_bamlss"(link = "logit") ,sampler = FALSE,
            multiple = FALSE)
b <- bamlss(bnum ~ x1 + x2 + x3, data = d, family = "frm_bamlss"(link = "probit") ,sampler = FALSE,
             multiple = FALSE)
# higher derivatives of link "loglog" not available!
b <- bamlss(bnum ~ x1 + x2 + x3, data = d, family = "frm_bamlss"(link = "loglog") ,sampler = FALSE,
            multiple = FALSE)
summary(b)
coef(b)
#Check whether the function is right, comparing the estimated coefficients with frm()
x <- cbind(d$x1, d$x2, d$x3)
colnames(x) <- c("x1", "x2", "x3")
a <- frm::frm(d$bnum, x, linkfrac = "logit")
a1 <- frm::frm(d$bnum, x, linkfrac = "probit")
#If the factor covariate id is included, these two functions would treat id differently and coefficients would be different
b <- bamlss(bnum ~ x1 + x2 + x3 + id, data = d, family = "frm_bamlss", sampler = FALSE,
            multiple = FALSE)
coef(b)

x <- cbind(d$x1, d$x2, d$x3, d$id)
colnames(x) <- c("x1", "x2", "x3", "id")
a <- frm::frm(d$bnum, x, linkfrac = "logit")


x <- 2
attr(x, "contrasts")

any(c(x>3, x<5))
any(x>3, x<5)
any(x>3||x<5)
all(x>3,x<5)

x <- 3.3
round(x)
