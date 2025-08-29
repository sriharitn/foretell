## =========================== Geometric (Geom) ===========================
## -------- helpers --------
#' @keywords internal
.geom_sigmoid <- function(x) 1 / (1 + exp(-x))
.geom_logit   <- function(p) log(p / (1 - p))
.geom_clip01  <- function(p, eps = 1e-8) pmin(1 - eps, pmax(eps, p))

# survival probabilities S(0..n_last) for Geom with p, + one-and-done & cure
#' @keywords internal
.geom_surv_prob_from_pars <- function(p, d = 0, cure = 0, n_last) {
  stopifnot(is.finite(p), p > 0, p < 1)
  k <- 0:n_last
  S_base <- (1 - p)^k
  if (n_last >= 1) {
    # apply one-and-done + cure modifiers for t>=1 (same pattern as bDW/sBG)
    S_base[-1] <- (1 - d) * (cure + (1 - cure) * S_base[-1])
  }
  S_base
}

# scale conversion + label (kept identical behavior)
#' @keywords internal
.geom_convert_scale <- function(S, object, scale = c("input","percent","count","prob")) {
  scale <- match.arg(scale)
  if (scale == "prob")    return(S)
  if (scale == "percent") return(S * 100)
  if (scale == "count")   return(S * object$N0)
  if (object$input == "percent") S * 100 else S * object$N0
}
.geom_label_for_scale <- function(scale, object) {
  scale <- match.arg(scale, c("input","percent","count","prob"))
  if (scale == "input") {
    if (object$input == "percent") "% surviving" else "# surviving"
  } else c(percent = "% surviving", count = "# surviving", prob = "Survival probability")[scale]
}

# PD repair for covariance matrices (nearPD -> eigen floor)
#' @keywords internal
.geom_make_pd <- function(S, eps_floor = 1e-8) {
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
.geom_J_params <- function(theta, est_d, est_cure) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Package 'numDeriv' is required for Jacobians.")
  f <- function(th) {
    i <- 1L
    p    <- .geom_clip01(.geom_sigmoid(th[i])); i <- i + 1L
    d    <- if (est_d)    { v <- .geom_clip01(.geom_sigmoid(th[i])); i <- i + 1L; v } else 0
    cure <- if (est_cure)   .geom_clip01(.geom_sigmoid(th[i])) else 0
    c(p, if (est_d) d, if (est_cure) cure)
  }
  numDeriv::jacobian(f, theta)
}
.geom_J_repar <- function(theta, est_d, est_cure) {
  # For geometric we take the reparameterization to be identical (p[,d][,cure])
  .geom_J_params(theta, est_d, est_cure)
}

# -----------------------------------------------------------------------
# Optimizer wrapper: BFGS -> (BFGS again as fallback) -> nloptr::slsqp
# -----------------------------------------------------------------------
#' @keywords internal
.geom_optimize <- function(theta0, nll) {
  if (length(theta0) == 1L) {
    # 1D problem: use Brent search on logit(p) (-Inf, Inf)
    out <- try(stats::optim(theta0, fn = nll, method = "Brent",
                            lower = -20, upper = 20), silent = TRUE)
    return(out)
  }
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
#' @keywords internal
`%||%` <- function(a,b) if (!is.null(a)) a else b



#' Geometric survival: One and Done and Cure Models
#'
#' Fit a geometric survival model (optionally with one-and-done and cure),
#' using logit parameterization, robust SEs from the Hessian,
#' and S3 methods for printing, summarizing, predicting, plotting.
#'
#' @param surv_value Numeric survival series (percent starting at 100, or counts).
#' @param h Integer forecast horizon (>= 0).
#' @param N0 Cohort size if `input="percent"`.
#' @param one_and_done Logical; include one-and-done parameter \eqn{d}.
#' @param cure Logical; include cure fraction parameter.
#' @param starts_p Starting values for \eqn{p} on the original scale in \eqn{(0,1)}.
#' @param starts_d,starts_cure Starting values for \eqn{d} and \eqn{\pi} on the original scale in \eqn{[0,1]}.
#' @param compute_se Logical; compute robust SEs via Hessian + delta method.
#' @param input One of `"auto"`, `"percent"`, `"count"`.
#' @param percent_tol,boundary_tol Tolerances for monotonicity and boundary SE skipping.
#'
#' @details
#' Let \eqn{p \in (0,1)} be the per-period churn probability in a geometric model.
#' The baseline geometric survival (no one-and-done, no cure) at discrete times
#' \eqn{t = 0,1,\dots} is
#' \deqn{S(t) = (1-p)^t, \qquad S(0)=1.}
#' With one-and-done mass \eqn{d} and cure fraction \eqn{\pi}, survival for \eqn{t \ge 1} is
#' \deqn{S(t) = (1-d)\left[\pi + (1-\pi)\,(1-p)^{\,t-1}\right], \qquad S(0)=1.}
#'
#' Parameters are estimated by maximizing the likelihood implied by the observed
#' survival series using an internal unconstrained parameterization for numerical
#' stability (e.g., a logit transform for \eqn{p,d,\pi}). If \code{compute_se = TRUE},
#' standard errors on the original scale are obtained via the delta method from
#' the Hessian on the working scale.
#'
#' @return An object of class `geom_fit`.
#'
#' @importFrom stats optim optimHess rgeom rbinom pnorm logLik
#' @importFrom graphics plot lines points legend abline
#' @importFrom grDevices adjustcolor
#' @examples
#' \dontrun{
#' N0 <- 500; S <- c(100,94,90,87,85,84,83,82,81,80,79)
#' fit <- geom1c(S, h=6, N0=N0, input="percent")
#' summary(fit)
#' plot(fit, scale="percent")
#' }

#' @export
#' @name geom1c
geom1c <- function(
    surv_value, h,
    N0 = NULL,
    one_and_done = FALSE,
    cure        = FALSE,
    starts_p    = 0.15,
    starts_d    = 0.05,
    starts_cure = 0.05,
    compute_se  = TRUE,
    input = c("auto","percent","count"),
    percent_tol = 1e-8,
    boundary_tol = 1e-3
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

  # logits
  theta0 <- c(
    .geom_logit(starts_p),
    if (est_d)    .geom_logit(starts_d)    else NULL,
    if (est_cure) .geom_logit(starts_cure) else NULL
  )

  nll <- function(theta) {
    i <- 1L
    th_p   <- theta[i]; i <- i + 1L
    th_d   <- if (est_d)    { v <- theta[i]; i <- i + 1L; v } else -Inf
    th_cur <- if (est_cure)   theta[i] else -Inf

    p    <- .geom_sigmoid(th_p)
    d    <- if (est_d)    .geom_sigmoid(th_d)   else 0
    cure <- if (est_cure) .geom_sigmoid(th_cur) else 0

    S  <- .geom_surv_prob_from_pars(p, d, cure, n_last = Tobs - 1)
    f  <- pmax(-diff(S), 1e-12)
    ST <- max(S[Tobs], 1e-12)
    -(sum(die * log(f)) + surv_T * log(ST))
  }

  # optimize
  opt <- .geom_optimize(theta0, nll)

  # decode
  i <- 1L
  th_p   <- opt$par[i]; i <- i + 1L
  th_d   <- if (est_d)    { v <- opt$par[i]; i <- i + 1L; v } else -Inf
  th_cur <- if (est_cure)   opt$par[i] else -Inf

  p_hat    <- .geom_sigmoid(th_p)
  d_hat    <- if (est_d)    .geom_sigmoid(th_d)   else 0
  cure_hat <- if (est_cure) .geom_sigmoid(th_cur) else 0

  # fitted & projected
  S_fit <- .geom_surv_prob_from_pars(p_hat, d_hat, cure_hat, n_last = Tobs - 1)
  fitted <- if (input == "percent") S_fit * 100 else S_fit * N0

  if (h > 0) {
    S_all  <- .geom_surv_prob_from_pars(p_hat, d_hat, cure_hat, n_last = Tobs - 1 + h)
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
        theta_hat <- c(th_p, if (est_d) th_d, if (est_cure) th_cur)
        Jp <- .geom_J_params(theta_hat, est_d, est_cure)
        Jr <- .geom_J_repar (theta_hat, est_d, est_cure)
        vcov_params <- Jp %*% vcov_theta %*% t(Jp)
        vcov_repar  <- Jr %*% vcov_theta %*% t(Jr)
        se_params   <- sqrt(pmax(diag(vcov_params), 0))
        se_repar    <- sqrt(pmax(diag(vcov_repar),  0))
        names(se_params) <- c("p", if (est_d) "d", if (est_cure) "cure")
        names(se_repar)  <- c("p", if (est_d) "d", if (est_cure) "cure")
      }
    }
  }

  obj <- list(
    call  = match.call(),
    y     = y_raw,
    input = input, N0 = N0, h = h,
    flags = c(one_and_done = est_d, cure = est_cure),
    coef_params = setNames(c(p_hat, if (est_d) d_hat, if (est_cure) cure_hat),
                           c("p", if (est_d) "d", if (est_cure) "cure")),
    coef_repar  = setNames(c(p_hat, if (est_d) d_hat, if (est_cure) cure_hat),
                           c("p", if (est_d) "d", if (est_cure) "cure")),
    logits      = setNames(c(th_p, if (est_d) th_d, if (est_cure) th_cur),
                           c("theta_p", if (est_d) "theta_d", if (est_cure) "theta_cure")),
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
  class(obj) <- "geom_fit"
  obj
}

#' @export
#' @method print geom_fit
print.geom_fit <- function(x, digits = 4, ...) {
  cat("Geometric (Geom) unconstrained (logits)\n")
  cat("One-and-done:", x$flags["one_and_done"], " | Cure:", x$flags["cure"], "\n")
  cat("LogLik:", round(x$logLik, digits), "  Convergence:", x$convergence, "\n")
  if (!is.null(x$se_note)) cat("Note:", x$se_note, "\n")
  cat("\nCoefficients (original scale):\n")
  print(round(x$coef_params, digits))
  cat("\nReparameterization (p",
      if (x$flags["one_and_done"]) ", d" else "",
      if (x$flags["cure"])        ", cure" else "", "):\n", sep = "")
  print(round(x$coef_repar, digits))
  invisible(x)
}


#' @export
#' @method summary geom_fit
summary.geom_fit <- function(object, digits = 4, ...) {
  make_tab <- function(est, se) {
    if (is.null(est)) return(NULL)
    if (is.null(se)) {
      data.frame(Estimate = est, check.names = FALSE)
    } else {
      z <- est / se
      p <- 2 * stats::pnorm(-abs(z))
      data.frame(
        Estimate   = est,
        `Std. Error` = se,
         check.names = FALSE
      )
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
  class(out) <- "summary.geom_fit"
  out
}

#' @export
#' @method print summary.geom_fit
print.summary.geom_fit <- function(x, digits = 4, ...) {
  cat("Geometric summary\n")
  cat("LogLik:", round(x$logLik, digits), "  Convergence:", x$convergence, "\n")
  if (!is.null(x$note)) cat("Note:", x$note, "\n")
  cat("\nOriginal parameters (p, [, d][, cure]):\n")
  if (!is.null(x$params)) print(round(x$params, digits)) else cat("(not available)\n")
  cat("\nReparameterization (p[, d][, cure]):\n")
  if (!is.null(x$repar))  print(round(x$repar, digits))  else cat("(not available)\n")
  invisible(x)
}

#' @export
#' @method coef geom_fit
coef.geom_fit <- function(object, reparam = FALSE, ...) {
  if (isTRUE(reparam)) object$coef_repar else object$coef_params
}

#' @export
#' @method vcov geom_fit
vcov.geom_fit <- function(object, what = c("params","repar","theta"), ...) {
  what <- match.arg(what)
  if (what == "params") return(object$vcov_params)
  if (what == "repar")  return(object$vcov_repar)
  object$vcov_theta
}

#' @export
#' @method logLik geom_fit
logLik.geom_fit <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- length(object$logits)
  class(val) <- "logLik"
  val
}

#' @export
#' @method fitted geom_fit
fitted.geom_fit <- function(object, scale = c("input","percent","count","prob"), ...) {
  scale <- match.arg(scale)
  Tobs <- length(object$y)
  p <- object$coef_params
  ph <- p["p"]
  d    <- if ("d"    %in% names(p)) p["d"]    else 0
  cure <- if ("cure" %in% names(p)) p["cure"] else 0
  S <- .geom_surv_prob_from_pars(ph, d, cure, n_last = Tobs - 1)
  .geom_convert_scale(S[1:Tobs], object, scale)
}

#' @export
#' @method predict geom_fit
predict.geom_fit <- function(object, n_ahead = NULL, scale = c("input","percent","count","prob"), ...) {
  scale <- match.arg(scale)
  Tobs <- length(object$y)
  h <- if (is.null(n_ahead)) object$h else as.integer(n_ahead)
  if (is.null(h) || h <= 0) return(numeric(0))
  p <- object$coef_params
  ph <- p["p"]
  d    <- if ("d"    %in% names(p)) p["d"]    else 0
  cure <- if ("cure" %in% names(p)) p["cure"] else 0
  S <- .geom_surv_prob_from_pars(ph, d, cure, n_last = Tobs - 1 + h)
  .geom_convert_scale(S[(Tobs + 1):(Tobs + h)], object, scale)
}

#' @export
#' @method plot geom_fit
plot.geom_fit <- function(x,
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
  ph <- p["p"]
  d    <- if ("d"    %in% names(p)) p["d"]    else 0
  cure <- if ("cure" %in% names(p)) p["cure"] else 0

  # fitted/proj
  S_fit <- .geom_surv_prob_from_pars(ph, d, cure, n_last = Tobs - 1)
  y_fit <- .geom_convert_scale(S_fit[1:Tobs], x, scale)

  y_proj <- numeric(0); t_prj <- integer(0)
  if (!is.null(h) && h > 0) {
    S_all  <- .geom_surv_prob_from_pars(ph, d, cure, n_last = Tobs - 1 + h)
    S_proj <- S_all[(Tobs + 1):(Tobs + h)]
    y_proj <- .geom_convert_scale(S_proj, x, scale)
    t_prj  <- Tobs:(Tobs + h - 1)
  }

  t_obs <- 0:(Tobs - 1)

  if (is.null(ylab)) ylab <- .geom_label_for_scale(scale, x)
  if (is.null(main)) {
    od <- if (isTRUE(x$flags["one_and_done"])) " + one-and-done" else ""
    cu <- if (isTRUE(x$flags["cure"]))        " + cure"         else ""
    main <- paste0("Geometric", od, cu)
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


