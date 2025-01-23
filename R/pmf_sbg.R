#' 	Shifted Beta-geometric (sbg) Distribution Family Function
#'
#' Provides functions for the probability mass function (PMF), cumulative distribution function (CDF), quantile function, and random variate generation for the SBG distribution.
#'
#' @usage
#' dsbg(x, shape1, shape2, log = FALSE)
#' psbg(x, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
#' qsbg(p, shape1, shape2)
#' rsbg(n, shape1, shape2)
#'
#' @param x Vector of non-negative integers for `dsbg` and `psbg`.
#' @param p Vector of probabilities (0 <= p <= 1) for `qsbg`.
#' @param n Number of random variates to generate for `rsbg`.
#' @param shape1 First shape parameter (must be > 0).
#' @param shape2 Second shape parameter (must be > 0).
#' @param log Logical; if TRUE, probabilities are returned on the log scale (for `dsbg` and `psbg`).
#' @param lower.tail Logical; if TRUE (default), probabilities are P(X <= x), otherwise P(X > x) (for `psbg`).
#' @param log.p Logical; if TRUE, probabilities are returned on the log scale (for `psbg` and `qsbg`).
#' @return
#' - `dsbg`: A numeric vector of PMF values.
#' - `psbg`: A numeric vector of CDF values.
#' - `qsbg`: A numeric vector of quantile values.
#' - `rsbg`: A numeric vector of random variates.
#'
#' @examples
#' # PMF example
#' dsbg(1:5, shape1 = 2, shape2 = 3)
#'
#' # CDF example
#' psbg(1:5, shape1 = 2, shape2 = 3)
#'
#' # Quantile example
#' qsbg(c(0.1, 0.5, 0.9), shape1 = 2, shape2 = 3)
#'
#' # Random variates
#' rsbg(10, shape1 = 2, shape2 = 3)
#'
#' @name sbg
NULL

#' @export
dsbg <- function(x, shape1, shape2, log = FALSE) {
  non_integers <- x[x != as.integer(x)]
  if (length(non_integers) > 0) {
    warning(paste0("The following values of x are non-integers: ", paste(non_integers, collapse = ", ")))
  }
  if (any(x < 0)) {
    stop(paste0("All values of x must be non-negative integers. Found: ", paste(x[x < 0], collapse = ", ")))
  }
  if (shape1 <= 0 || shape2 <= 0) {
    stop("Both shape1 and shape2 must be greater than 0.")
  }
  cdf_value1 <- 1 - (beta(shape1, shape2 + x) / beta(shape1, shape2))
  cdf_value0 <- 1 - (beta(shape1, shape2 + (x - 1)) / beta(shape1, shape2))
  pdf_value <- cdf_value1 - cdf_value0
  pdf_value <- ifelse(x == 0, 0, pdf_value)
  if (log) {
    pdf_value <- log(pdf_value)
  }
  return(pdf_value)
}

#' @export

psbg <- function(x, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
  non_integers <- x[x != as.integer(x)]
  if (length(non_integers) > 0) {
    warning(paste0("The following values of x are non-integers: ", paste(non_integers, collapse = ", ")))
  }
  if (any(x < 0)) {
    stop(paste0("All values of x must be non-negative integers. Found: ", paste(x[x < 0], collapse = ", ")))
  }
  if (shape1 <= 0 || shape2 <= 0) {
    stop("Both shape1 and shape2 must be greater than 0.")
  }
  cdf_value <- 1 - (beta(shape1, shape2 + x) / beta(shape1, shape2))
  if (!lower.tail) {
    cdf_value <- 1 - cdf_value
  }
  if (log.p) {
    cdf_value <- log(cdf_value)
  }
  return(cdf_value)
}

#' @export

qsbg <- function(p, shape1, shape2) {
  if (any(p < 0 | p > 1)) {
    stop("p must be between 0 and 1")
  }
  if (shape1 <= 0 || shape2 <= 0) {
    stop("Both shape1 and shape2 must be greater than 0.")
  }

  compute_quantile_scalar <- function(p) {
    if (p == 0) {
      return(0)
    } else {
      result <- pracma::newton(
        function(t) beta(shape1, shape2 + t) / beta(shape1, shape2) - p,
        x0 = shape1 / (shape1 + shape2)
      )
      return(round(result$root))
    }
  }

  vectorized_compute <- Vectorize(compute_quantile_scalar)
  return(vectorized_compute(p))
}

#' @export

rsbg <- function(n, shape1, shape2) {
  if (!is.numeric(n) || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  if (shape1 <= 0 || shape2 <= 0) {
    stop("Both shape1 and shape2 must be greater than 0.")
  }
  theta <- rbeta(n, shape1, shape2)
  theta <- pmin(pmax(theta, 0.0000001), 0.9999999)
  return(rgeom(n, prob = theta) + 1)
}
