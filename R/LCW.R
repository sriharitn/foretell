#' Latent Class Weibull (LCW) Model for Projecting Customer Retention
#'
#' \code{LCW} is a latent class weibull model implementation based on \code{Fader and Hardie} probability based projection methedology. The survivor function for \code{LCW} is \deqn{wS(t|t1,c1)+(1-w)S(t|t2,c2), 0<w<1}
#'
#' @param surv_value a numeric vector of historical customer retention percentage should start at 100 and non-starting values should be between 0 and less than 100
#' @param h forecasting horizon
#' @param lower lower limit used in \code{R} \code{optim} rotuine. Default is \code{c(0.001,0.001,0.001,0.001,0.001)}.
#' @param upper upper limit used in \code{R} \code{optim} rotuine. Default is \code{c(0.99999,10000,0.999999,10000,0.99999)}.
#' @param subjects Total number of customers or subject default 1000
#'
#' @return
#' \item{fitted:}{Fitted Values based on historical data}
#' \item{projected:}{Projected \code{h} values based on historical data}
#' \item{max.likelihood:}{Maximum Likelihood of LCW}
#' \item{params - t1,t2,c1,c2,w:}{Returns t1,c1,t2,c2,w paramters from maximum likelihood estimation}
#'
#' @examples
#' surv_value <- c(100,86.9,74.3,65.3,59.3,55.1,51.7,49.1,46.8,44.5,42.7,40.9,39.4)
#' h <- 6
#' LCW(surv_value,h)
#'
#' @references {Fader P, Hardie B. How to project customer retention. Journal of Interactive Marketing. 2007;21(1):76-90.}
#' @references {Fader P, Hardie B, Liu Y, Davin J, Steenburgh T. "How to Project Customer Retention" Revisited: The Role of Duration Dependence. Journal of Interactive Marketing. 2018;43:1-16.}
#' @export



LCW <- function(surv_value,h, lower = c(0.001,0.001,0.001,0.001,0.001),upper = c(0.99999,10000,0.999999,10000,0.99999),
                subjects = 1000){

  surv <- surv_value

  if(surv[1] != 100) stop("Starting Value should be 100")

  if(any(surv[-1] >= 100) | any(surv[-1] < 0)) stop("Starting Value should be 100 and non-starting value should be between 0 and less than 100")

  surv <- (surv/100)*subjects

  t <- length(surv)

  die <- diff(-surv)

  s <- rep(NA,length(surv))
  p <- rep(NA,length(surv))

  cw.lik<-function(params) {

    t1<-params[1]
    c1<-params[2]
    t2<-params[3]
    c2<-params[4]
    w <-params[5]

    i = 0:(t-1)

    s <- w*(1-t1)^(i^c1)+(1-w)*(1-t2)^(i^c2)

    p <- diff(-s)

    ll_ <- (die[i])*log(p[i])

    ll <- sum(ll_)+(surv[t])*log(s[t])

    return(-ll)
  }

  max.lik.cw <- tryCatch({
    stats::optim(c(0.5,2,0.5,1,0.6),fn=cw.lik,lower = lower,
                 upper = upper,method="L-BFGS-B")
  }, error = function(error_condition) {
    message("Note: stats::optim not working switching to nloptr::slsqp for maximum likelihood optimization")
    nloptr::slsqp(c(0.5,2,0.5,1,0.6),fn=cw.lik,lower = lower,
                  upper = upper)
  })


  t1 <- max.lik.cw$par[1]
  c1 <- max.lik.cw$par[2]
  t2 <- max.lik.cw$par[3]
  c2 <- max.lik.cw$par[4]
  w  <- max.lik.cw$par[5]


  k <- 0:(t+h)

  cw <- (w*(1-t1)^(k^c1)+(1-w)*(1-t2)^(k^c2))*100

  projected <- if (h > 0) {
    projected <- cw[(t + 1):(t + h)]
  } else {
    message("Forecast horizon h is: ",h,", No Forecast generated.")
    projected <- numeric(0) # Return an empty numeric vector if h = 0
  }

  list(fitted = cw[1:t],projected = projected,max.likelihood = max.lik.cw$value, params = c(t1=t1,c1=c1,t2=t2,c2=c2,w=w))

}



