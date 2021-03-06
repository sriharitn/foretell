#' Beta discrete Weibull (BdW) Model for Projecting Customer Retention.
#'
#' \code{BdW} is a beta discrete weibull model implemented based on \code{Fader and Hardie} probability based projection methedology. The survivor function for \code{BdW} is \deqn{Beta(a,b+t^c)/Beta(a,b)}
#'
#' @param surv_value a numeric vector of historical customer retention percentage should start at 100 and non-starting values should be between 0 and less than 100
#' @param h forecasting horizon
#' @param lower lower limit used in \code{R} \code{optim} rotuine. Default is \code{c(1e-3,1e-3)}.
#' @param upper upper limit used in \code{R} \code{optim} rotuine. Default is \code{c(10000,10000,10000)}.
#'
#' @return
#' \item{fitted:}{Fitted values based on historical data}
#' \item{projected:}{Projected \code{h} values based on historical data}
#' \item{max.likelihood:}{Maximum Likelihood of Beta discrete Weibull}
#' \item{params - a, b and c:}{Returns a and b paramters from maximum likelihood estimation for beta distribution and c}
#'
#' @examples
#' surv_value <- c(100,86.9,74.3,65.3,59.3)
#' h <- 6
#' BdW(surv_value,h)
#'
#' @references {Fader P, Hardie B. How to project customer retention. Journal of Interactive Marketing. 2007;21(1):76-90.}
#' @references {Fader P, Hardie B, Liu Y, Davin J, Steenburgh T. "How to Project Customer Retention" Revisited: The Role of Duration Dependence. Journal of Interactive Marketing. 2018;43:1-16.}
#' @export



BdW <- function(surv_value,h, lower = c(0.001,0.001,0.001),upper = c(10000,10000,10000)){

  surv <- surv_value

  if(surv[1] != 100) stop("Starting Value should be 100")

  if(any(surv[-1] >= 100) | any(surv[-1] < 0)) stop("Starting Value should be 100 and non-starting value should be between 0 and less than 100")

  t <- length(surv)

  die <- diff(-surv)


  s <- rep(NA,length(surv))
  p <- rep(NA,length(surv))

  dbw.log.lik<-function(params) {
    a<-params[1]
    b<-params[2]
    c<-params[3]


    i = 0:(t-1)

    s <- beta(a,b+i^c)/beta(a,b)

    p <- diff(-s)

    ll_ <- (die[i])*log(p[i])

    ll <- sum(ll_)+(surv[t])*log(s[t])

    return(-ll)
  }

  max.lik.dbw <- tryCatch({
    stats::optim(c(1,1,1),fn=dbw.log.lik,lower =lower, upper = upper,method="L-BFGS-B")
  }, error = function(error_condition) {
    message("Note: stats::optim not working switching to nloptr::slsqp for maximum likelihood optimization")
    nloptr::slsqp(c(1,1,1),fn=dbw.log.lik,lower =lower, upper = upper)
  })


  a <- max.lik.dbw$par[1]
  b <- max.lik.dbw$par[2]
  c <- max.lik.dbw$par[3]


  k <- 0:(t+h)

  dbw <- (beta(a, b+(k^c)) / beta(a, b))*100

  list(fitted = dbw[1:t],projected = dbw[(t+1):(t+h)],max.likelihood = max.lik.dbw$value, params = c(a = a,b = b,c = c))

}



