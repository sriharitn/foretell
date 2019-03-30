#' Latent Class Weibull (LCW) Model for Projecting Customer Retention
#'
#' \code{LCW} is a latent class weibull model implementation based on \code{Fader and Hardie} probability based projection methedology. The survivor function for \code{LCW} is \deqn{wS(t|t1,c1)+(1-w)S(t|t2,c2), 0<w<1}
#'
#' @param surv_value a numeric vector of historical customer retention percentage should start at 100 and values should be between 0 and 100
#' @param h forecasting horizon
#' @param lower lower limit used in \code{R} \code{optim} rotuine. Default is \code{c(0.001,0.001,0.001,0.001,0.001)}.
#' @param upper upper limit used in \code{R} \code{optim} rotuine. Default is \code{c(0.99999,10000,0.999999,10000,0.99999)}.
#'
#' @return
#' \item{fitted:}{Fitted Values based on historical data} \item{max.likelihood:}{Maximum Likelihood of Beta Geometric}
#' \item{params - t1,t2,c1,c2,w:}{Returns t1,c1,t2,c2,w paramters from maximum likelihood estimation}
#'
#' @examples
#' surv_value <- c(100,86.9,74.3,65.3,59.3,55.1,51.7,49.1,46.8,44.5,42.7,40.9,39.4)
#' h <- 6
#' LCW(surv_value,h)
#'
#' @references {Fader  PS and Hardie BGS (2007), How to project customer retention. Volume 21, Issue 1. Journal of Interactive Marketing}
#' @references {Fader  PS and Hardie BGS et al. (2018), How to Project Customer Retention Revisited: The Role of Duration Dependence. Volume 43, Journal of Interactive Marketing}
#'
#' @export



LCW <- function(surv_value,h, lower = c(0.001,0.001,0.001,0.001,0.001),upper = c(0.99999,10000,0.999999,10000,0.99999)){

  surv <- surv_value

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


  max.lik.cw  <- stats::optim(c(0.5,2,0.5,1,0.6),fn=cw.lik,lower = lower,
                         upper = upper,method="L-BFGS-B")


  t1 <- max.lik.cw$par[1]
  c1 <- max.lik.cw$par[2]
  t2 <- max.lik.cw$par[3]
  c2 <- max.lik.cw$par[4]
  w  <- max.lik.cw$par[5]


  k <- 0:(t+h)

  cw <- w*(1-t1)^(k^c1)+(1-w)*(1-t2)^(k^c2)

  list(fitted = 1,cw[1:t],projected = cw[(t+1):(t+h)],max.likelihood = max.lik.cw$value, params = c(t1=t1,c1=c1,t2=t2,c2=c2,w=w))

}



