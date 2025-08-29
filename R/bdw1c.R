## =========================== Beta Discrete Weibull ===========================
## -------- helpers --------
#' @keywords internal
.bdw_sigmoid <- function(x) 1 / (1 + exp(-x))
.bdw_logit   <- function(p) log(p / (1 - p))
.bdw_clip01  <- function(p, eps = 1e-8) pmin(1 - eps, pmax(eps, p))

# survival probabilities S(0..n_last) for sgm with a,b, + one-and-done & cure
#' @keywords internal
.bdw_surv_prob_from_pars <- function(a, b, c, d = 0, cure = 0, n_last) {
  stopifnot(is.finite(a), is.finite(b), is.finite(c), a > 0, b > 0, c > 0)
  k <- 0:n_last
  S_base <- exp(lbeta(a, b + k^c) - lbeta(a, b))
  if (n_last >= 1) {
    # apply one-and-done + cure modifiers for t>=1
    S_base[-1] <- (1 - d) * (cure + (1 - cure) * S_base[-1])
  }
  S_base
}

# scale conversion + label (kept identical behavior)
#' @keywords internal
.bdw_convert_scale <- function(S, object, scale = c("input","percent","count","prob")) {
  scale <- match.arg(scale)
  if (scale == "prob")    return(S)
  if (scale == "percent") return(S * 100)
  if (scale == "count")   return(S * object$N0)
  if (object$input == "percent") S * 100 else S * object$N0
}
.bdw_label_for_scale <- function(scale, object) {
  scale <- match.arg(scale, c("input","percent","count","prob"))
  if (scale == "input") {
    if (object$input == "percent") "% surviving" else "# surviving"
  } else c(percent = "% surviving", count = "# surviving", prob = "Survival probability")[scale]
}


# PD repair for covariance matrices (nearPD -> eigen floor)
#' @keywords internal
.bdw_make_pd <- function(S, eps_floor = 1e-8) {
  S <- as.matrix((S + t(S)) / 2)
  ch <- try(chol(S), silent = TRUE)
  if (!inherits(ch, "try-error")) return(S)
  if (requireNamespace("Matrix", quietly = TRUE)) {
    S2 <- as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)
    ch2 <- try(chol(S2), silent = TRUE)
    if (!inherits(ch2, "try-error")) return(S2)
    S <- S2
  }
  eg  <- eigen(S, symmetric = TRUE)
  eps <- max(1e-10, eps_floor * max(eg$values, 0))
  vals <- pmax(eg$values, eps)
  Sfix <- eg$vectors %*% diag(vals, length(vals)) %*% t(eg$vectors)
  ch3 <- try(chol(Sfix), silent = TRUE)
  if (inherits(ch3, "try-error")) stop("Could not repair Sigma to PD.")
  Sfix
}


# Jacobians (theta -> params / reparam) for vcov mapping (uses numDeriv)
# theta = (th_p[, th_d][, th_cure])
#' @keywords internal
.bdw_J_params <- function(theta, est_d, est_cure) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Package 'numDeriv' is required for Jacobians.")
  f <- function(th) {
    i <- 1L
    m      <- .bdw_clip01(.bdw_sigmoid(th[i])); i <- i + 1L
    q      <- .bdw_clip01(.bdw_sigmoid(th[i])); i <- i + 1L
    c_par  <- exp(th[i]);                          i <- i + 1L  # c > 0 via log-transform
    d      <- if (est_d)    { v <- .bdw_clip01(.bdw_sigmoid(th[i])); i <- i + 1L; v } else 0
    cure   <- if (est_cure)   .bdw_clip01(.bdw_sigmoid(th[i])) else 0
    ab_sum <- (1/q) - 1
    a <- m * ab_sum
    b <- (1 - m) * ab_sum
    c(a, b, c_par, if (est_d) d, if (est_cure) cure)
  }
  numDeriv::jacobian(f, theta)
}
.bdw_J_repar <- function(theta, est_d, est_cure) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Package 'numDeriv' is required for Jacobians.")
  f <- function(th) {
    i <- 1L
    m      <- .bdw_clip01(.bdw_sigmoid(th[i])); i <- i + 1L
    q      <- .bdw_clip01(.bdw_sigmoid(th[i])); i <- i + 1L
    c_par  <- exp(th[i]);                          i <- i + 1L
    d      <- if (est_d)    { v <- .bdw_clip01(.bdw_sigmoid(th[i])); i <- i + 1L; v } else 0
    cure   <- if (est_cure)   .bdw_clip01(.bdw_sigmoid(th[i])) else 0
    c(m = m, q = q, c = c_par, if (est_d) d = d, if (est_cure) cure = cure)
  }
  numDeriv::jacobian(f, theta)
}

# -----------------------------------------------------------------------
# Optimizer wrapper: BFGS -> (BFGS again as fallback) -> nloptr::slsqp
# -----------------------------------------------------------------------
.bdw_optimize <- function(theta0, nll) {
  out <- try(stats::optim(theta0, fn = nll, control = list(reltol = 1e-10)), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  out <- try(stats::optim(theta0, fn = nll, method = "BFGS"), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  if (!requireNamespace("nloptr", quietly = TRUE))
    stop("Need 'nloptr' for SLSQP fallback.")
  res <- nloptr::slsqp(x0 = theta0, fn = nll)
  list(par = res$par %||% res$solution, value = res$value, convergence = res$convergence,
       counts = NA, message = res$message)
}
`%||%` <- function(a,b) if (!is.null(a)) a else b



#' Beta-Discrete-Weibull (BDW) survival: One-and-Done and Cure
#'
#' Fits the Beta-Discrete-Weibull (BDW) survival model to a monotonically
#' nonincreasing survival series using unconstrained maximum likelihood
#' estimation. Internal reparameterization ensures \eqn{a>0}, \eqn{b>0}, and
#' \eqn{c>0} (e.g., via \eqn{\log c}), with optional one-and-done mass \eqn{d}
#' and cure fraction \eqn{\pi}. When \eqn{c=1}, the BDW reduces to the shifted
#' Beta-Geometric (sBG).
#'
#' @param surv_value Numeric survival series (percent starting at 100, or counts).
#' @param h Integer forecast horizon (>= 0).
#' @param N0 Cohort size if `input = "percent"`.
#' @param one_and_done Logical; include one-and-done parameter \eqn{d}.
#' @param cure Logical; include cure fraction parameter (fraction never churning).
#' @param starts_m,starts_q Starting values for \eqn{m,q} on the original scale in \eqn{(0,1)}.
#' @param starts_c Starting values for \eqn{c} on the original scale.
#' @param starts_d,starts_cure Starting values for \eqn{d} and \eqn{\pi} on the original scale in \eqn{[0,1]}.
#' @param compute_se Logical; compute robust SEs via Hessian + delta method.
#' @param input One of `"auto"`, `"percent"`, `"count"`. If `"auto"`, infer from `surv_value` and `N0`.
#' @param percent_tol,boundary_tol Tolerances for monotonicity checks and boundary SE skipping.
#'
#' @details
#' The baseline BDW survival (no one-and-done, no cure) at discrete times \eqn{t=0,1,\dots} is
#' \deqn{S(t) = \frac{B(a, b+t^{c})}{B(a,b)}.}
#' With one-and-done mass \eqn{d} and cure fraction \eqn{\pi}, survival for \eqn{t \ge 1} is
#' \deqn{S(t) = (1-d)\left[\pi + (1-\pi)\frac{B(a, b+t^{c})}{B(a,b)}\right], \qquad S(0)=1.}
#' Parameters \eqn{(a,b,c)} (and optionally \eqn{d,\pi}) are estimated by maximizing the
#' likelihood implied by the observed survival series using an internal unconstrained
#' parameterization. If \code{compute_se = TRUE}, standard errors on the original scale
#' are obtained via the delta method.
#'
#' \strong{Beta reparameterization by mean and polarization index:}
#' The model uses a numerically stable and interpretable reparameterization of
#' the Beta distribution \eqn{Beta(a,b)} in terms of a **mean**
#' \eqn{m \in (0,1)} and a **polarization (concentration) index**
#' \eqn{q \in (0,1)}.
#'
#' Standard Beta density:
#' \deqn{f(x \mid a,b) = \frac{x^{a-1}(1-x)^{b-1}}{B(a,b)}, \quad x\in(0,1),\; a,b>0.}
#'
#' Mapping to mean and polarization:
#' \deqn{m = \frac{a}{a+b}, \qquad q = \frac{1}{1+a+b}.}
#' Here, \eqn{m} is the expected probability, and \eqn{q} summarizes concentration:
#' small \eqn{q} means large \eqn{a+b} (high concentration); large \eqn{q}
#' means small \eqn{a+b} (diffuse).
#'
#' Inverse mapping (used internally for estimation):
#' \deqn{a = m \left(\frac{1}{q}-1\right), \qquad
#'       b = (1-m)\left(\frac{1}{q}-1\right).}
#'
#' This transformation is bijective for \eqn{m \in (0,1)}, \eqn{q \in (0,1)},
#' and guarantees \eqn{a>0, b>0}.
#'
#' \strong{Starting values:}
#' The arguments \code{starts_m} and \code{starts_q} provide starting values for
#' \eqn{m} and \eqn{q}, respectively. They are converted to \eqn{(a,b)} via the
#' inverse mapping above:
#' \deqn{a_{\text{start}} = \text{starts\_m}\left(\frac{1}{\text{starts\_q}}-1\right), \quad
#'       b_{\text{start}} = \left(1-\text{starts\_m}\right)\left(\frac{1}{\text{starts\_q}}-1\right).}
#' This parameterization typically improves optimization stability and makes
#' starting values more interpretable.
#'
#' @return
#' An object of class \code{bdw_fit} with the elements shown below:
#' \describe{
#'   \item{\code{y}}{Numeric vector of observed survival values on the input scale.}
#'   \item{\code{input}}{Character scalar, one of \code{"input"}, \code{"percent"},
#'     \code{"count"}, or \code{"prob"}; echoes how \code{y} was interpreted.}
#'   \item{\code{N0}}{Scalar. Reference cohort size used for count/percent scaling.}
#'   \item{\code{flags}}{Named logical vector with model switches, e.g.
#'     \code{c(one_and_done = TRUE/FALSE, cure = TRUE/FALSE)}.}
#'   \item{\code{coef_params}}{Named numeric vector of parameters on the natural scale.
#'     For sBG: \code{c(a,b)} and optionally \code{d} (one-and-done) and \code{cure}.
#'     For BDW: \code{c(a,b,c)} and optionally \code{d}, \code{cure}.}
#'   \item{\code{coef_repar}}{Named numeric vector with the mean-polarization reparameterization,
#'     typically \code{c(m,q)} (and \code{c} for BDW if modeled on the log/exp scale).}
#'   \item{\code{logits}}{Named numeric vector of unconstrained optimization variables
#'     (e.g., \code{theta_m}, \code{theta_q}, \code{theta_d}, \code{theta_cure}, \code{theta_c}).}
#'   \item{\code{fitted}}{Numeric vector \eqn{S(0{:}T_{\mathrm{obs}}-1)} on the probability scale; the in-sample fit.}
#'   \item{\code{projected}}{Numeric vector \eqn{S(0{:}T_{\mathrm{last}})} (fit + horizon) on the probability scale.}
#'   \item{\code{logLik}}{Maximized log-likelihood.}
#'   \item{\code{convergence}}{Optimizer return code (0 indicates successful convergence).}
#'   \item{\code{optim_out}}{Raw optimizer list (as returned by \code{optim()}); useful for debugging.}
#'   \item{\code{vcov_theta}}{Variance-covariance matrix on the unconstrained scale (logits).}
#'   \item{\code{vcov_params}}{Variance-covariance matrix for \code{coef_params} (delta-method mapped from \code{vcov_theta}).}
#'   \item{\code{vcov_repar}}{Variance-covariance matrix for \code{coef_repar} (e.g., \code{m,q}).}
#'   \item{\code{se_params}}{Named vector of standard errors for \code{coef_params}.}
#'   \item{\code{se_repar}}{Named vector of standard errors for \code{coef_repar}.}
#'   \item{\code{se_note}}{Character note if SEs are approximate/unstable (e.g., near-PD fix).}
#' }
#'
#'
#' @importFrom stats optim optimHess pnorm logLik
#' @importFrom graphics plot lines points legend abline
#' @importFrom grDevices adjustcolor
#'
#' @examples
#' \dontrun{
#' N0 <- 500; S <- c(100, 86.9, 74.3, 65.3, 59.3, 55.1, 51.7, 49.1)
#' fit <- sbg1c(S, h=6, N0=N0, input="percent")
#' summary(fit)
#' plot(fit, scale="percent")
#' }
#'
#' @references {Fader P, Hardie B. How to project customer retention. Journal of Interactive Marketing. 2007;21(1):76-90.}
#' @references {Fader P, Hardie B, Liu Y, Davin J, Steenburgh T. "How to Project Customer Retention" Revisited: The Role of Duration Dependence. Journal of Interactive Marketing. 2018;43:1-16.}


#' @export
#' @name bdw1c
bdw1c <- function(
    surv_value, h,
    N0 = NULL,                          # required if input="percent"
    one_and_done = FALSE,
    cure        = FALSE,                # include cure fraction parameter?
    starts_m  = 0.6,
    starts_q  = 0.1,
    starts_c  = 1.0,                    # NEW: starting value for c (>0)
    starts_d  = 0.05,
    starts_cure = 0.05,
    compute_se = TRUE,
    input = c("auto","percent","count"),
    percent_tol = 1e-8,
    boundary_tol = 1e-3                 # skip SEs if d/cure near 0/1
) {
  input <- match.arg(input)
  y_raw <- as.numeric(surv_value)
  if (length(y_raw) < 2L) stop("Need at least two survival points.")
  if (any(!is.finite(y_raw))) stop("surv_value must be finite.")
  if (any(diff(y_raw) > percent_tol)) stop("Survival must be non-increasing.")

  # Detect & convert to counts
  if (input == "auto") {
    if (abs(y_raw[1] - 100) <= percent_tol) input <- "percent"
    else if (y_raw[1] > 100 + percent_tol)  input <- "count"
    else stop("Can't auto-detect input: first value not 100 (percent) or >100 (count).")
  }
  if (input == "percent") {
    if (abs(y_raw[1] - 100) > percent_tol) stop("Percent input must start at 100.")
    if (is.null(N0) || N0 <= 0)            stop("When input='percent', supply positive N0.")
    y_counts <- y_raw / 100 * N0
  } else {
    if (y_raw[1] <= 0) stop("Count input must start with positive cohort size.")
    N0 <- y_raw[1]; y_counts <- y_raw
  }

  Tobs   <- length(y_counts)         # t=0..Tobs-1
  die    <- -diff(y_counts)          # length Tobs-1
  surv_T <-  y_counts[Tobs]

  est_d    <- isTRUE(one_and_done)
  est_cure <- isTRUE(cure)

  # logits + log(c)
  theta0 <- c(
    .bdw_logit(starts_m),
    .bdw_logit(starts_q),
    log(pmax(starts_c, 1e-6)),
    if (est_d)    .bdw_logit(starts_d)    else NULL,
    if (est_cure) .bdw_logit(starts_cure) else NULL
  )

  nll <- function(theta) {
    i <- 1L
    th_m     <- theta[i]; i <- i + 1L
    th_q     <- theta[i]; i <- i + 1L
    th_logc  <- theta[i]; i <- i + 1L
    th_d     <- if (est_d)    { v <- theta[i]; i <- i + 1L; v } else -Inf
    th_cur   <- if (est_cure)   theta[i] else -Inf

    m    <- .bdw_sigmoid(th_m);    q <- .bdw_sigmoid(th_q)
    cpar <- exp(th_logc)
    d    <- if (est_d)    .bdw_sigmoid(th_d)   else 0
    cure <- if (est_cure) .bdw_sigmoid(th_cur) else 0

    ab_sum <- (1/q) - 1
    a <- max(m * ab_sum, 1e-12)
    b <- max((1 - m) * ab_sum, 1e-12)

    S  <- .bdw_surv_prob_from_pars(a, b, cpar, d, cure, n_last = Tobs - 1)
    p  <- pmax(-diff(S), 1e-12)
    ST <- max(S[Tobs], 1e-12)

    -(sum(die * log(p)) + surv_T * log(ST))
  }

  # optimize
  opt <- .bdw_optimize(theta0, nll)

  # decode
  i <- 1L
  th_m     <- opt$par[i]; i <- i + 1L
  th_q     <- opt$par[i]; i <- i + 1L
  th_logc  <- opt$par[i]; i <- i + 1L
  th_d     <- if (est_d)    { v <- opt$par[i]; i <- i + 1L; v } else -Inf
  th_cur   <- if (est_cure)   opt$par[i] else -Inf

  m_hat    <- .bdw_sigmoid(th_m)
  q_hat    <- .bdw_sigmoid(th_q)
  c_hat    <- exp(th_logc)
  d_hat    <- if (est_d)    .bdw_sigmoid(th_d)   else 0
  cure_hat <- if (est_cure) .bdw_sigmoid(th_cur) else 0

  ab_sum_hat <- (1 / q_hat) - 1
  a_hat <- m_hat * ab_sum_hat
  b_hat <- (1 - m_hat) * ab_sum_hat

  # fitted & projected
  S_fit <- .bdw_surv_prob_from_pars(a_hat, b_hat, c_hat, d_hat, cure_hat, n_last = Tobs - 1)
  fitted <- if (input == "percent") S_fit * 100 else S_fit * N0

  if (h > 0) {
    S_all  <- .bdw_surv_prob_from_pars(a_hat, b_hat, c_hat, d_hat, cure_hat, n_last = Tobs - 1 + h)
    S_proj <- S_all[(Tobs + 1):(Tobs + h)]
    projected <- if (input == "percent") S_proj * 100 else S_proj * N0
  } else projected <- numeric(0)

  # SEs via Hessian in theta-space
  vcov_theta <- vcov_params <- vcov_repar <- NULL
  se_params <- se_repar <- NULL
  se_note <- NULL
  do_se <- isTRUE(compute_se)

  if (do_se) {
    probs <- c()
    if (est_d)    probs <- c(probs, d_hat)
    if (est_cure) probs <- c(probs, cure_hat)
    if (length(probs)) {
      if (any(!is.finite(probs))) {
        do_se  <- FALSE
        se_note <- "SEs disabled: d/cure not estimable (NA/Inf)."
      } else if (any(probs < boundary_tol | probs > (1 - boundary_tol))) {
        do_se  <- FALSE
        se_note <- sprintf("SEs skipped (d/cure within %.3f of 0 or 1).", boundary_tol)
      }
    }
  }

  if (do_se) {
    H <- try(stats::optimHess(opt$par, nll), silent = TRUE)
    if (!inherits(H, "try-error")) {
      H <- (H + t(H)) / 2
      inv <- try(solve(H), silent = TRUE)
      pseudo_used <- FALSE
      if (inherits(inv, "try-error")) {
        if (requireNamespace("Matrix", quietly = TRUE)) {
          Hpd <- as.matrix(Matrix::nearPD(H, corr = FALSE)$mat)
          inv <- try(solve(Hpd), silent = TRUE)
          if (inherits(inv, "try-error")) {
            eg <- eigen(Hpd, symmetric = TRUE)
            tol <- max(1e-10, 1e-8 * max(eg$values, 0))
            inv_vals <- ifelse(eg$values > tol, 1/eg$values, 0)
            inv <- eg$vectors %*% diag(inv_vals, length(inv_vals)) %*% t(eg$vectors)
          }
        } else {
          eg <- eigen(H, symmetric = TRUE)
          tol <- max(1e-10, 1e-8 * max(eg$values, 0))
          inv_vals <- ifelse(eg$values > tol, 1/eg$values, 0)
          inv <- eg$vectors %*% diag(inv_vals, length(inv_vals)) %*% t(eg$vectors)
        }
        pseudo_used <- TRUE
      }
      vcov_theta <- inv
      attr(vcov_theta, "pseudo_inverse") <- pseudo_used

      if (!requireNamespace("numDeriv", quietly = TRUE)) {
        warning("numDeriv not available; skipping SEs.")
      } else {
        theta_hat <- c(th_m, th_q, th_logc, if (est_d) th_d, if (est_cure) th_cur)
        Jp <- .bdw_J_params(theta_hat, est_d, est_cure)
        Jr <- .bdw_J_repar (theta_hat, est_d, est_cure)
        vcov_params <- Jp %*% vcov_theta %*% t(Jp)
        vcov_repar  <- Jr %*% vcov_theta %*% t(Jr)
        se_params   <- sqrt(pmax(diag(vcov_params), 0))
        se_repar    <- sqrt(pmax(diag(vcov_repar),  0))
        names(se_params) <- c("a","b","c", if (est_d) "d", if (est_cure) "cure")
        names(se_repar)  <- c("m","q","c", if (est_d) "d", if (est_cure) "cure")
      }
    }
  }

  obj <- list(
    call  = match.call(),
    y     = y_raw,
    input = input, N0 = N0, h = h,
    flags = c(one_and_done = est_d, cure = est_cure),
    coef_params = setNames(c(a_hat, b_hat, c_hat, if (est_d) d_hat, if (est_cure) cure_hat),
                           c("a","b","c", if (est_d) "d", if (est_cure) "cure")),
    coef_repar  = setNames(c(m_hat, q_hat, c_hat, if (est_d) d_hat, if (est_cure) cure_hat),
                           c("m","q","c", if (est_d) "d", if (est_cure) "cure")),
    logits      = setNames(c(th_m, th_q, th_logc, if (est_d) th_d, if (est_cure) th_cur),
                           c("theta_m","theta_q","theta_logc", if (est_d) "theta_d", if (est_cure) "theta_cure")),
    fitted      = fitted,
    projected   = projected,
    logLik      = -opt$value,
    convergence = opt$convergence,
    optim_out   = opt,
    # SE stuff
    vcov_theta  = vcov_theta,
    vcov_params = vcov_params,
    vcov_repar  = vcov_repar,
    se_params   = se_params,
    se_repar    = se_repar,
    se_note     = se_note
  )
  class(obj) <- "bdw_fit"
  obj
}

#' @export
#' @method print bdw_fit
print.bdw_fit <- function(x, digits = 4, ...) {
  cat("Beta-Discrete-Weibull (bDW)\n")
  cat("One-and-done:", x$flags["one_and_done"], " | Cure:", x$flags["cure"], "\n")
  cat("LogLik:", round(x$logLik, digits), "  Convergence:", x$convergence, "\n")
  if (!is.null(x$se_note)) cat("Note:", x$se_note, "\n")
  cat("\nCoefficients (original scale):\n")
  print(round(x$coef_params, digits))
  cat("\nReparameterization (m, q, c",
      if (x$flags["one_and_done"]) ", d" else "",
      if (x$flags["cure"])        ", cure" else "", "):\n", sep = "")
  print(round(x$coef_repar, digits))
  invisible(x)
}


#' @export
#' @method summary bdw_fit
summary.bdw_fit <- function(object, digits = 4, ...) {
  make_tab <- function(est, se) {
    if (is.null(est)) return(NULL)
    if (is.null(se)) {
      data.frame(Estimate = est, check.names = FALSE)
    } else {
      data.frame(Estimate = est, `Std. Error` = se,check.names = FALSE)
    }
  }
  tab_params <- make_tab(object$coef_params, object$se_params)
  tab_repar  <- make_tab(object$coef_repar,  object$se_repar)

  out <- list(
    call        = object$call,
    logLik      = object$logLik,
    convergence = object$convergence,
    note        = object$se_note,
    params      = tab_params,
    repar       = tab_repar,
    vcov_theta  = object$vcov_theta,
    vcov_params = object$vcov_params,
    vcov_repar  = object$vcov_repar
  )
  class(out) <- "summary.bdw_fit"
  out
}


#' @export
#' @method print summary.bdw_fit
print.summary.bdw_fit <- function(x, digits = 4, ...) {
  cat("Beta-Discrete-Weibull (bDW) summary\n")
  cat("LogLik:", round(x$logLik, digits), "  Convergence:", x$convergence, "\n")
  if (!is.null(x$note)) cat("Note:", x$note, "\n")
  cat("\nOriginal parameters (a, b, c[, d][, cure]):\n")
  if (!is.null(x$params)) print(round(x$params, digits)) else cat("(not available)\n")
  cat("\nReparameterization (m, q, c[, d][, cure]):\n")
  if (!is.null(x$repar))  print(round(x$repar, digits))  else cat("(not available)\n")
  invisible(x)
}

#' @export
#' @method coef bdw_fit
coef.bdw_fit <- function(object, reparam = FALSE, ...) {
  if (isTRUE(reparam)) object$coef_repar else object$coef_params
}

#' @export
#' @method vcov bdw_fit
vcov.bdw_fit <- function(object, what = c("params","repar","theta"), ...) {
  what <- match.arg(what)
  if (what == "params") return(object$vcov_params)
  if (what == "repar")  return(object$vcov_repar)
  object$vcov_theta
}

#' @export
#' @method logLik bdw_fit
logLik.bdw_fit <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- length(object$logits)
  class(val) <- "logLik"
  val
}

#' @export
#' @method fitted bdw_fit
fitted.bdw_fit <- function(object, scale = c("input","percent","count","prob"), ...) {
  scale <- match.arg(scale)
  Tobs <- length(object$y)
  p <- object$coef_params
  a <- p["a"]; b <- p["b"]; cpar <- p["c"]
  d    <- if ("d"    %in% names(p)) p["d"]    else 0
  cure <- if ("cure" %in% names(p)) p["cure"] else 0
  S <- .bdw_surv_prob_from_pars(a, b, cpar, d, cure, n_last = Tobs - 1)
  .bdw_convert_scale(S[1:Tobs], object, scale)
}

#' @export
#' @method predict bdw_fit
predict.bdw_fit <- function(object, n_ahead = NULL, scale = c("input","percent","count","prob"), ...) {
  scale <- match.arg(scale)
  Tobs <- length(object$y)
  h <- if (is.null(n_ahead)) object$h else as.integer(n_ahead)
  if (is.null(h) || h <= 0) return(numeric(0))
  p <- object$coef_params
  a <- p["a"]; b <- p["b"]; cpar <- p["c"]
  d    <- if ("d"    %in% names(p)) p["d"]    else 0
  cure <- if ("cure" %in% names(p)) p["cure"] else 0
  S <- .bdw_surv_prob_from_pars(a, b, cpar, d, cure, n_last = Tobs - 1 + h)
  .bdw_convert_scale(S[(Tobs + 1):(Tobs + h)], object, scale)
}

#' @export
#' @method plot bdw_fit
plot.bdw_fit <- function(x,
                         scale = c("input","percent","count","prob"),
                         bands = NULL,  # list: $fundamental, $estimation, $both (each: time, lower, upper, S_hat?, part)
                         show_observed_points = TRUE, show_observed_line = TRUE,
                         show_model_line = TRUE, show_projection_line = TRUE,
                         show_boundary = TRUE, connect_boundary = TRUE,
                         # lines
                         col_observed = "gray40",
                         col_model    = "steelblue",
                         col_proj     = "tomato",
                         lwd_observed = 1.2, lwd_model = 2, lwd_proj = 2,
                         # polygon fills + borders
                         col_fundamental_fill = "forestgreen",
                         col_estimation_fill  = "darkorange",
                         col_both_fill        = "purple",
                         alpha_fundamental = 0.20,
                         alpha_estimation  = 0.20,
                         alpha_both        = 0.20,
                         col_fundamental_border = "forestgreen",
                         col_estimation_border  = "darkorange",
                         col_both_border        = "purple",
                         lwd_band_border = 1.2,
                         draw_band_center = TRUE,
                         lty_band_center  = 3,
                         lwd_band_center  = 1.5,
                         main = NULL, xlab = "Period", ylab = NULL, xlim = NULL, ylim = NULL, ...) {
  scale <- match.arg(scale)
  y <- x$y; Tobs <- length(y)
  h <- x$h; if (is.null(h)) h <- length(x$projected)

  # observed -> requested scale
  if (scale == "prob") {
    y_obs <- if (x$input == "percent") (y / 100) else (y / x$N0)
  } else if (scale == "percent") {
    y_obs <- if (x$input == "percent") y else (y / x$N0) * 100
  } else if (scale == "count") {
    y_obs <- if (x$input == "percent") (y / 100) * x$N0 else y
  } else y_obs <- y

  # model params
  p <- x$coef_params
  a <- p["a"]; b <- p["b"]; cpar <- p["c"]
  d    <- if ("d"    %in% names(p)) p["d"]    else 0
  cure <- if ("cure" %in% names(p)) p["cure"] else 0

  # fitted/proj
  S_fit <- .bdw_surv_prob_from_pars(a, b, cpar, d, cure, n_last = Tobs - 1)
  y_fit <- .bdw_convert_scale(S_fit[1:Tobs], x, scale)

  y_proj <- numeric(0); t_prj <- integer(0)
  if (!is.null(h) && h > 0) {
    S_all  <- .bdw_surv_prob_from_pars(a, b, cpar, d, cure, n_last = Tobs - 1 + h)
    S_proj <- S_all[(Tobs + 1):(Tobs + h)]
    y_proj <- .bdw_convert_scale(S_proj, x, scale)
    t_prj  <- Tobs:(Tobs + h - 1)
  }

  t_obs <- 0:(Tobs - 1)

  if (is.null(ylab)) ylab <- .bdw_label_for_scale(scale, x)
  if (is.null(main)) {
    od <- if (isTRUE(x$flags["one_and_done"])) " + one-and-done" else ""
    cu <- if (isTRUE(x$flags["cure"]))        " + cure"         else ""
    main <- paste0("Beta-Discrete-Weibull", od, cu)
  }

  # axis limits: include polygons if present
  if (is.null(ylim)) {
    y_all <- c(y_obs, y_fit, y_proj)
    for (nm in c("fundamental","estimation","both")) {
      if (!is.null(bands[[nm]])) y_all <- c(y_all, bands[[nm]]$lower, bands[[nm]]$upper)
    }
    ylim <- range(y_all, na.rm = TRUE)
    if (scale %in% c("percent","prob") || (scale == "input" && x$input == "percent")) {
      ylim[1] <- max(0, ylim[1])
      cap <- if (scale == "percent" || (scale == "input" && x$input == "percent")) 100 else 1
      ylim[2] <- min(cap, max(ylim[2], cap))
    }
  }
  if (is.null(xlim)) xlim <- c(0, if (length(t_prj)) max(t_prj) else max(t_obs))

  # canvas
  plot(t_obs, y_obs, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)

  draw_poly_band <- function(df, col_fill, col_border, alpha) {
    if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
    polygon(c(df$time, rev(df$time)),
            c(df$lower, rev(df$upper)),
            col    = grDevices::adjustcolor(col_fill, alpha.f = alpha),
            border = NA)
    lines(df$time, df$lower, col = col_border, lwd = lwd_band_border)
    lines(df$time, df$upper, col = col_border, lwd = lwd_band_border)
    if (draw_band_center && !is.null(df$S_hat)) {
      lines(df$time, df$S_hat, col = col_border, lwd = lwd_band_center, lty = lty_band_center)
    }
    invisible(NULL)
  }
  draw_poly_band(bands$both,        col_both_fill,        col_both_border,        alpha_both)
  draw_poly_band(bands$estimation,  col_estimation_fill,  col_estimation_border,  alpha_estimation)
  draw_poly_band(bands$fundamental, col_fundamental_fill, col_fundamental_border, alpha_fundamental)

  if (show_observed_line)   lines(t_obs, y_obs, col = col_observed, lwd = lwd_observed)
  if (show_observed_points) points(t_obs, y_obs, col = col_observed, pch = 16)
  if (show_model_line)      lines(t_obs, y_fit, col = col_model, lwd = lwd_model)

  if (length(y_proj) && show_projection_line) {
    if (connect_boundary) lines(c(Tobs - 1, Tobs), c(y_fit[Tobs], y_proj[1]), col = col_proj, lwd = lwd_proj, lty = 2)
    lines(t_prj, y_proj, col = col_proj, lwd = lwd_proj, lty = 2)
  }
  if (show_boundary && length(y_proj)) abline(v = Tobs - 1, lty = 3, col = "darkgray")

  # legend
  leg_items <- c(); leg_cols <- c(); leg_lty <- c(); leg_lwd <- c(); leg_pch <- c()

  # Observed
  if (show_observed_points || show_observed_line) {
    leg_items <- c(leg_items, "Observed")
    leg_cols  <- c(leg_cols, col_observed)
    leg_lty   <- c(leg_lty, if (show_observed_line) 1 else NA)
    leg_lwd   <- c(leg_lwd, if (show_observed_line) lwd_observed else NA)
    leg_pch   <- c(leg_pch, if (show_observed_points) 16 else NA)
  }

  # Fitted (no points)
  if (show_model_line) {
    leg_items <- c(leg_items, "Fitted")
    leg_cols  <- c(leg_cols, col_model)
    leg_lty   <- c(leg_lty, 1)
    leg_lwd   <- c(leg_lwd, lwd_model)
    leg_pch   <- c(leg_pch, NA)
  }

  # Projected (no points in legend)
  if (length(y_proj) && show_projection_line) {
    leg_items <- c(leg_items, "Projected")
    leg_cols  <- c(leg_cols, col_proj)
    leg_lty   <- c(leg_lty, 2)
    leg_lwd   <- c(leg_lwd, lwd_proj)
    leg_pch   <- c(leg_pch, NA)
  }

  # Bands (no points)
  if (!is.null(bands$fundamental)) {
    leg_items <- c(leg_items, "Fundamental band")
    leg_cols  <- c(leg_cols, col_fundamental_border)
    leg_lty   <- c(leg_lty, 1)
    leg_lwd   <- c(leg_lwd, lwd_band_border)
    leg_pch   <- c(leg_pch, NA)
  }
  if (!is.null(bands$estimation)) {
    leg_items <- c(leg_items, "Estimation band")
    leg_cols  <- c(leg_cols, col_estimation_border)
    leg_lty   <- c(leg_lty, 1)
    leg_lwd   <- c(leg_lwd, lwd_band_border)
    leg_pch   <- c(leg_pch, NA)
  }
  if (!is.null(bands$both)) {
    leg_items <- c(leg_items, "Both band")
    leg_cols  <- c(leg_cols, col_both_border)
    leg_lty   <- c(leg_lty, 1)
    leg_lwd   <- c(leg_lwd, lwd_band_border)
    leg_pch   <- c(leg_pch, NA)
  }

  if (length(leg_items) > 0) {
    legend("topright",
           legend = leg_items,
           col    = leg_cols,
           lty    = leg_lty,
           lwd    = leg_lwd,
           pch    = leg_pch,
           merge  = TRUE,      # show point+line on "Observed" if both present
           pt.cex = 1.2,
           bty    = "n")
  }

  invisible(bands)
}

