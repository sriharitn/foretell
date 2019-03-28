#' Beta discrete Weibull (BdW) Model for Projecting Customer Retention.
#'
#' \code{BdW} is a beta discrete weibull model implemented based on \code{Fader and Hardie} probability based projection methedology
#'
#' @param surv_value a numeric vector of historical customer retention should start at 1 and values should be between 0 and 1
#' @param h forecasting horizon
#' @param lower lower limit used in \code{R} \code{optim} rotuine. Default is \code{c(1e-3,1e-3)}.
#' @param upper upper limit used in \code{R} \code{optim} rotuine. Default is \code{c(10000,10000,10000)}.
#'
#' @return
#' \item{fitted values:}{Fitted Values based on historical data} \item{max.likelihood:}{Maximum Likelihood of Beta Geometric}
#' \item{params - alpha, beta c:}{Returns alpha and beta paramters from maximum likelihood estimation for beta distribution and c}
#'
#' @examples
#' surv_value <- c(1.0,0.869,0.743,0.653,0.593)
#' h <- 6
#' BdW(surv_value,h)
#'
#' @references {Fader  PS and Hardie BGS (2007), How to project customer retention. Volume 21, Issue 1. Journal of Interactive Marketing.}
#' @references {Fader  PS and Hardie BGS et al. (2018), How to Project Customer Retention Revisited: The Role of Duration Dependence. Volume 43, Journal of Interactive Marketing.}
#'
#' @export



BdW <- function(surv_value,h, lower = c(0.001,0.001,0.001),upper = c(10000,10000,10000)){

  surv <- surv_value[-1]


  ## sBG likelihood
  t <- length(surv)

  i = 1:t

  survx <- c(1,surv)

  die <- survx[i]-survx[i+1]


  s <- rep(NA,length(surv))
  p <- rep(NA,length(surv))

  dbw.log.lik<-function(params) {
    a<-params[1]
    b<-params[2]
    c<-params[3]


    i = 1:t

    s <- beta(a,b+i^c)/beta(a,b)


    st <- c(1,s)

    p <- st[i]-st[i+1]

    ll_ <- (die[i]*1000)*log(p[i])

    ll <- sum(ll_)+(surv[t]*1000)*log(s[t])

    return(-ll)
  }


  max.lik.dbw <- stats::optim(c(1,1,1),fn=dbw.log.lik,lower =lower, upper = upper,method="L-BFGS-B")


  a <- max.lik.dbw$par[1]
  b <- max.lik.dbw$par[2]
  c <- max.lik.dbw$par[3]


  k <- 1:(t+h)

  dbw <- beta(a, b+(k^c)) / beta(a, b)

  list(fitted = c(1,dbw[1:t]),projected = dbw[(t+1):(t+h)],max.likelihood = max.lik.dbw$value, params = c(alpha = a,beta = b,c = c))

}



