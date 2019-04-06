#' Excel based trendlines for projecting customer retention.
#'
#' \code{exltrend} generates Microsoft(r) Excel(r) based linear, logarthmic, exponential, polynomial of order 2, power trends.
#'
#' @param surv_value a numeric vector of historical customer retention percentage should start at 100 and non-starting values should be between 0 and less than 100
#' @param h forecasting horizon
#'
#' @return
#' \item{fitted:}{A data frame of fitted Values based on historical data for linear (lin.p), exponential (exp.p), logarthmic (log.p), polynomial (poly.p) of order 2 and power (pow.p) trends.}
#' \item{projected:}{A data frame of projected \code{h} values based on historical data for linear (lin.p), exponential (exp.p), logarthmic (log.p), polynomial (poly.p) of order 2 and power (pow.p) trends.}
#'
#' @examples
#' surv_value <- c(100,86.9,74.3,65.3,59.3)
#' h <- 6
#' exltrend(surv_value,h)
#'
#' @export



exltrend <- function(surv_value,h){

  surv <- surv_value

  if(surv[1] != 100) stop("Starting Value should be 100")

  if(any(surv[-1] >= 100) | any(surv[-1] < 0)) stop("Starting Value should be 100 and non-starting value should be between 0 and less than 100")

  t <- length(surv)

  ## Input dataset
  x <- 1:t
  new <- data.frame(x = 1:(t+h))

  ##Linear Trend
  fit.lm <- stats::glm(surv~x)

  ##Exponential
  f.exp <- function(x,a,b) {a * exp(b * x)}
  fit.exp <- stats::nls(surv ~ f.exp(x,a,b), start = c(a=100, b=-1))

  ##Log
  f.log <- function(x,a,b) {a * log(x) + b}
  fit.log <- stats::nls(surv ~ f.log(x,a,b), start = c(a=-25, b=100))

  # polynomial
  f.poly <- function(x,a,b,d) {(a*x^2) + (b*x) + d}
  fit.poly <- stats::nls(surv ~ f.poly(x,a,b,d), start = c(a=1, b=-15, d=100))

  ##Power
  f.pow <- function(x,a,b){a*x^b}
  fit.pow <- stats::nls(surv ~ f.pow(x,a,b), start = c(a=100, b=-1))

  ##Prediction

  pred.lm   <- stats::predict(fit.lm,newdata=new)
  pred.exp  <- stats::predict(fit.exp,newdata=new)
  pred.log  <- stats::predict(fit.log,newdata=new)
  pred.poly <- stats::predict(fit.poly,newdata=new)
  pred.pow  <- stats::predict(fit.pow,newdata=new)

  pred <- data.frame(lin.p = pred.lm, exp.p=pred.exp, log.p=pred.log,poly.p=pred.poly,pow.p=pred.pow)

  row.names(pred) <- c()


  return(list(fitted = pred[1:t,],projected = pred[(t+1):(t+h),]))

}



