#' Beta Geometric (BG) Model for Projecting Customer Retention.
#'
#' \code{BG} is a beta geometric model implemented based on \code{Fader and Hardie} probability based projection methedology.
#'
#' @param surv_value a numeric vector of historical customer retention should start at 1 and values should be between 0 and 1
#' @param h forecasting horizon
#' @param lower lower limit used in \code{R} \code{optim} rotuine. Default is \code{c(1e-3,1e-3)}.
#'
#' @return
#' \item{fitted:}{Fitted Values based on historical data} \item{max.likelihood:}{Maximum Likelihood of Beta Geometric}
#' \item{params - alpha, beta:}{Returns alpha and beta paramters from maximum likelihood estimation for beta distribution}
#'
#' @examples
#' surv_value <- c(1,0.869,0.743,0.653,0.593)
#' h <- 6
#' BG(surv_value,h)
#'
#' @references {Fader  PS and Hardie BGS (2007), How to project customer retention. Volume 21, Issue 1. Journal of Interactive Marketing}
#' @export


BG <- function(surv_value,h,lower = c(1e-3,1e-3)){

  surv <- surv_value[-1]


  ## sBG likelihood
  t <- length(surv)

  i = 1:t

  survx <- c(1,surv)

  die <- survx[i]-survx[i+1]


  s <- rep(NA,length(surv))
  p <- rep(NA,length(surv))

  bg.log.lik<-function(params) {
    a<-params[1]
    b<-params[2]


    i = 1:t

    s <- beta(a,b+i)/beta(a,b)

    st <- c(1,s)

    p <- st[i]-st[i+1]

    ll_ <- (die[i]*1000)*log(p[i])

    ll <- sum(ll_)+(surv[t]*1000)*log(s[t])

    return(-ll)
  }


  max.lik.sgb <- stats::optim(c(1,2),fn=bg.log.lik,lower = lower,method="L-BFGS-B")


  a <- max.lik.sgb$par[1]
  b <- max.lik.sgb$par[2]

  k <- 1:(t+h)

  sbg <- beta(a, b+(k)) / beta(a, b)

  list(fitted = c(1,sbg[1:t]),projected = sbg[(t+1):(t+h)],max.likelihood = max.lik.sgb$value, params = c(alpha = a,beta = b))


}



