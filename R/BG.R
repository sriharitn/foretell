#' Beta Geometric (BG) Model for Projecting Customer Retention.
#'
#' \code{BG} is a beta geometric model implemented based on \code{Fader and Hardie} probability based projection methedology. The survivor function for \code{BG} is \deqn{Beta(a,b+t)/Beta(a,b)}
#'
#' @param surv_value a numeric vector of historical customer retention percentage should start at 100 and non-starting values should be between 0 and less than 100
#' @param h forecasting horizon
#' @param lower lower limit used in \code{R} \code{optim} rotuine. Default is \code{c(1e-3,1e-3)}.
#'
#' @return
#' \item{fitted:}{Fitted Values based on historical data} \item{max.likelihood:}{Maximum Likelihood of Beta Geometric}
#' \item{params - a, b:}{Returns a and b paramters from maximum likelihood estimation for beta distribution}
#'
#' @examples
#' surv_value <- c(100,86.9,74.3,65.3,59.3)
#' h <- 6
#' BG(surv_value,h)
#'
#' @references {Fader  PS and Hardie BGS (2007), How to project customer retention. Volume 21, Issue 1. Journal of Interactive Marketing}
#' @export


BG <- function(surv_value,h,lower = c(1e-3,1e-3)){

  surv <- surv_value

  if(surv[1] != 100) stop("Starting Value should be 100")

  if(any(surv[-1] >= 100) | any(surv[-1] < 0)) stop("Starting Value should be 100 and non-starting value should be between 0 and less than 100")

  t <- length(surv)

  die <- diff(-surv)

  s <- rep(NA,length(surv))
  p <- rep(NA,length(surv))

  bg.log.lik<-function(params) {
    a<-params[1]
    b<-params[2]


    i = 0:(t-1)

    s <- beta(a,b+i)/beta(a,b)

    p <- diff(-s)

    ll_ <- (die[i])*log(p[i])

    ll <- sum(ll_)+(surv[t])*log(s[t])

    return(-ll)
  }


  max.lik.dbw <- tryCatch({
    stats::optim(c(1,2),fn=bg.log.lik,lower = lower,method="L-BFGS-B")
  }, error = function(error_condition) {
    message("stats::optim not working switching to nloptr::slsqp for maximum likelihood optimization")
    nloptr::slsqp(c(1,2),fn=bg.log.lik,lower = lower)
  })


  a <- max.lik.sgb$par[1]
  b <- max.lik.sgb$par[2]

  k <- 0:(t+h)

  sbg <- (beta(a, b+(k)) / beta(a, b))*100

  list(fitted = sbg[1:t],projected = sbg[(t+1):(t+h)],max.likelihood = max.lik.sgb$value, params = c(a = a,b = b))


}



