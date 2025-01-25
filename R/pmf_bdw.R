#' 	Beta discrete Weibull (bdw) Distribution Family Function
#'
#' Provides functions for the probability mass function (PMF), cumulative distribution function (CDF), quantile function, and random variate generation for the BdW distribution.
#'
#' @usage
#' dbdw(x, shape1, shape2, shape3, log = FALSE)
#' pbdw(x, shape1, shape2, shape3, lower.tail = TRUE, log.p = FALSE)
#' qbdw(p, shape1, shape2, shape3, lower.tail = TRUE, log.p = FALSE)
#' rbdw(n, shape1, shape2, shape3)
#'
#' @param x Vector of non-negative integers for `dbdw` and `pbdw`.
#' @param p Vector of probabilities (0 <= p <= 1) for `qbdw`.
#' @param n Number of random variates to generate for `rbdw`.
#' @param shape1 First shape parameter "a" (must be > 0).
#' @param shape2 Second shape parameter "b" (must be > 0).
#' @param shape3 Second shape parameter "c" (must be > 0).
#' @param log Logical; if TRUE, probabilities are returned on the log scale (for `dsbg` and `psbg`).
#' @param lower.tail Logical; if TRUE (default), probabilities are P(X <= x), otherwise P(X > x) (for `psbg`).
#' @param log.p Logical; if TRUE, probabilities are returned on the log scale (for `psbg` and `qsbg`).
#' @return
#' * `dbdw`: A numeric vector of PMF values.
#' * `pbdw`: A numeric vector of CDF values.
#' * `qbdw`: A numeric vector of quantile values.
#' * `rbdw`: A numeric vector of random variates.
#' @md
#'
#' @references {Fader P, Hardie B. How to project customer retention. Journal of Interactive Marketing. 2007;21(1):76-90.}
#' @references {Fader P, Hardie B, Liu Y, Davin J, Steenburgh T. "How to Project Customer Retention" Revisited: The Role of Duration Dependence. Journal of Interactive Marketing. 2018;43:1-16.}

#'
#' @examples
#' # PMF example
#' dbdw(1:5, shape1 = 2, shape2 = 3,shape3 = 0.5)
#'
#' # CDF example
#' pbdw(1:5, shape1 = 2, shape2 = 3, shape3 = 0.5)
#'
#' # Quantile example
#' qbdw(c(0.1, 0.5, 0.9), shape1 = 2, shape2 = 3 , shape3 = 0.5)
#'
#' # Random variates
#' rbdw(10, shape1 = 2, shape2 = 3,  shape3 = 0.5)
#'
#' @name BetadiscreteWeibull
NULL

#' @export
dbdw <- function(x, shape1, shape2, shape3, log = FALSE) {
  non_integers <- x[x != as.integer(x)]
  if (length(non_integers) > 0) {
    warning(paste0(
      "The following values of x are non-integers: ",
      paste(non_integers, collapse = ", ")
    ))
  }
  if (any(x < 0)) {
    stop(paste0(
      "All values of x must be non-negative integers. Found: ",
      paste(x[x < 0], collapse = ", ")
    ))
  }
  if (shape1 <= 0 || shape2 <= 0 || shape3 <= 0) {
    stop("shape1, shape2 and shape3 must be greater than 0.")
  }
  cdf_value1 <-
    1 - (beta(shape1, shape2 + x ^ shape3) / beta(shape1, shape2))
  cdf_value0 <-
    1 - (beta(shape1, shape2 + (x - 1) ^ shape3) / beta(shape1, shape2))
  pdf_value <- cdf_value1 - cdf_value0
  pdf_value <- ifelse(x == 0, 0, pdf_value)
  if (log) {
    pdf_value <- log(pdf_value)
  }
  return(pdf_value)
}

#' @export

pbdw <-
  function(x,
           shape1,
           shape2,
           shape3,
           lower.tail = TRUE,
           log.p = FALSE) {
    non_integers <- x[x != as.integer(x)]
    if (length(non_integers) > 0) {
      warning(paste0(
        "The following values of x are non-integers: ",
        paste(non_integers, collapse = ", ")
      ))
    }
    if (any(x < 0)) {
      stop(paste0(
        "All values of x must be non-negative integers. Found: ",
        paste(x[x < 0], collapse = ", ")
      ))
    }
    if (shape1 <= 0 || shape2 <= 0 || shape3 <= 0) {
      stop("shape1, shape2 and shape3 must be greater than 0.")
    }
    cdf_value <-
      1 - (beta(shape1, shape2 + x ^ shape3) / beta(shape1, shape2))
    if (!lower.tail) {
      cdf_value <- 1 - cdf_value
    }
    if (log.p) {
      cdf_value <- log(cdf_value)
    }
    return(cdf_value)
  }

#' @export

qbdw <-
  function(p,
           shape1,
           shape2,
           shape3,
           lower.tail = TRUE,
           log.p = FALSE) {
    if (log.p) {
      p <- exp(p)
    }

    if (!lower.tail) {
      p <- 1 - p
    }

    if (any(p < 0 | p > 1)) {
      stop("p must be between 0 and 1")
    }
    if (shape1 <= 0 || shape2 <= 0 || shape3 <= 0) {
      stop("shape1, shape2 and shape3 must be greater than 0.")
    }



    compute_quantile_scalar <- function(p) {
      if (p == 0) {
        return(0)
      } else if (p == 1) {
        return(Inf)
      } else {
        # Define the target function
        target_function <- function(t) {
          beta_val <- beta(shape1, shape2 + t ^ shape3)
          base_beta <- beta(shape1, shape2)
          return(beta_val / base_beta - p)
        }

        result <-
          uniroot(
            f = target_function,
            lower = 1e-8,
            # Small positive value to avoid invalid domain
            upper = 1e200,
            # Adjust upper limit based on problem scale
            extendInt = "yes"
          )$root

        return(round(result))
      }
    }

    vectorized_compute <- Vectorize(compute_quantile_scalar)
    return(vectorized_compute(p))
  }


#' @export

rbdw <- function(n, shape1, shape2, shape3) {
  if (!is.numeric(n) || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  if (shape1 <= 0 || shape2 <= 0 || shape3 <= 0) {
    stop("shape1, shape2 and shape3 must be greater than 0.")
  }
  theta <- rbeta(n, shape1, shape2)
  theta <- pmin(pmax(theta, 0.0000001), 0.9999999)

  theta <- rbeta(1, shape1, shape2)

  # Generate a uniform random number for inversion
  u <- runif(n, 0, 1)

  # Invert the CDF to solve for t
  t <- ceiling(((log(1 - u) / log(1 - theta)) ^ (1 / shape3)))+1

  return(t)



}
