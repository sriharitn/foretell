# sim.sbg1c.R
# =======================================================================
# Shifted Beta Geometric: simulator + nested bootstrap bands (package-ready)
# =======================================================================

# ----- INTERNAL HELPERS (not exported) ---------------------------------

# lifetimes measured in discrete periods: 1,2,3,..., Inf (cured)
#' @keywords internal
rsbg_mix <- function(n, a, b, d = 0, cure = 0) {
  # coerce to scalars and sanitize
  n    <- as.integer(n)[1]; a <- as.numeric(a)[1]; b <- as.numeric(b)[1]
  d    <- as.numeric(d)[1];  cure <- as.numeric(cure)[1]
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) stop("Invalid (a,b).")
  if (!is.finite(n) || n <= 0) stop("Invalid n.")
  eps <- 1e-12
  if (!is.finite(d)) d <- 0; if (!is.finite(cure)) cure <- 0
  d    <- min(max(d,    0), 1 - eps)
  cure <- min(max(cure, 0), 1 - eps)

  drop1 <- rbinom(n, 1, d) == 1L
  ten <- integer(n); ten[drop1] <- 1L
  stay <- !drop1
  if (any(stay)) {
    cured <- rbinom(sum(stay), 1, cure) == 1L
    idx_cured <- which(stay)[cured]; idx_nc <- which(stay)[!cured]
    if (length(idx_cured)) ten[idx_cured] <- Inf
    if (length(idx_nc)) {
      theta <- rbeta(length(idx_nc), a, b)
      theta <- pmin(pmax(theta, 1e-7), 1 - 1e-7)
      ten[idx_nc] <- rgeom(length(idx_nc), prob = theta) + 1L
    }
  }
  ten
}

# Survival curve from individual tenures
#' @keywords internal
.sbg_surv_from_tenures <- function(ten, Tlast) {
  # S(t) = mean(tenure > t) for t = 0..Tlast
  sapply(0:Tlast, function(tt) mean(ten > tt))
}

# Make covariance PD and draw from MVN in theta-space (logits)
# Depends on .sbg_make_pd() defined elsewhere in your package.
#' @keywords internal
.sbg_draw_thetas <- function(B, mu, Sigma) {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required for nested bootstrap.")
  Sig <- .sbg_make_pd(Sigma)
  MASS::mvrnorm(B, mu = mu, Sigma = Sig)
}


# ----- EXPORTED USER FUNCTION ------------------------------------------

#' Parametric bootstrap uncertainty bands for a shifted geometric survival fit
#'
#' Computes uncertainty bands for a fitted \code{sbg_fit} object, returning
#' up to three band types:
#' \itemize{
#'   \item \strong{fundamental} - process (sampling) variability at fixed parameters
#'   \item \strong{estimation}  - parameter uncertainty mapped to the mean survival curve
#'   \item \strong{both}        - predictive bands (parameter draw + one process draw)
#' }
#'
#' Set \code{B2 = 0} to skip process simulation inside the estimation bands (fast,
#' pure estimation-only) to average process noise in estimation bands.
#'
#' @param object A fitted \code{sbg_fit} from your geometric model constructor
#'   (e.g., \code{geom_extended()}), containing fields \code{y}, \code{N0}, \code{h},
#'   \code{coef_params}, \code{logits}, \code{vcov_theta}, and \code{flags}.
#' @param B1 Integer, number of parameter draws for estimation and both bands.
#' @param B2 Integer, Average process noise in estimation bands (set \code{B2 = 0}
#'   for pure estimation bands).
#' @param level Confidence level for bands (default \code{0.95}).
#' @param scale Output scale: \code{"prob"}, \code{"percent"}, or \code{"count"}.
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical; print minimal progress messages.
#'
#' @return A list with data frames named \code{fundamental}, \code{estimation},
#'   and \code{both} (whichever were computed). Each data frame has columns:
#'   \code{time}, \code{S_hat}, \code{lower}, \code{upper}, \code{part}.
#'
#' @examples
#' \dontrun{
#' N0 <- 400
#' S  <- c(100, 86.9, 74.3, 65.3, 59.3, 55.1, 51.7, 49.1)
#' fit <- geom1c(S, h = 6, N0 = N0, input = "percent")
#' b   <- sim.geom1c(fit, B1 = 200, B2 = 50, level = 0.90, scale = "percent", seed = 123)
#' plot(fit, scale = "percent", bands = b)
#' }
#' @references {King G. , Tomz M. and Wittenberg J. Making the most of statistical analyses: Improving interpretation and presentation.
#' American journal of political science, 2000;347-61.}
#'
#' @importFrom stats rgeom rbinom quantile optimHess pnorm setNames
#' @importFrom grDevices adjustcolor
#' @importFrom graphics plot lines points legend abline
#' @export
sim.sbg1c <- function(
    object,
    B1 = 1000, B2 = 0, level = 0.95,
    scale = c("percent","count","prob"),
    seed = NULL, verbose = FALSE
) {
  stopifnot(inherits(object, "sbg_fit"))
  scale <- match.arg(scale)
  if (!is.null(seed)) set.seed(seed)

  Tobs  <- length(object$y)
  Tlast <- (Tobs - 1) + max(0L, as.integer(object$h))
  N0    <- object$N0
  p     <- object$coef_params
  a0    <- p["a"]; b0 <- p["b"]
  d0    <- if ("d"    %in% names(p)) p["d"]    else 0
  c0    <- if ("cure" %in% names(p)) p["cure"] else 0

  to_scale <- function(Sprob) switch(scale,
                                     prob = Sprob, percent = Sprob * 100, count = Sprob * N0
  )
  alpha <- 1 - level
  res <- list()

  ## Fundamental only
  if (isTRUE(verbose)) message("Computing fundamental bands (B1=", B1, ")...")
  SmatF <- matrix(NA_real_, nrow = B1, ncol = Tlast + 1L)
  for (b in seq_len(B1)) {
    ten <- rsbg_mix(N0, a0, b0, d0, c0)
    SmatF[b, ] <- .sbg_surv_from_tenures(ten, Tlast)
  }
  S_hat <- colMeans(SmatF)
  res$fundamental <- data.frame(
    time  = 0:Tlast,
    S_hat = to_scale(S_hat),
    lower = to_scale(apply(SmatF, 2, stats::quantile, probs = alpha/2)),
    upper = to_scale(apply(SmatF, 2, stats::quantile, probs = 1 - alpha/2)),
    part  = ifelse(0:Tlast <= (Tobs - 1), "fitted", "projected")
  )

  ## Estimation only (avg over B1 process draws)
  if (is.null(object$vcov_theta)) stop("Need vcov_theta from fit for estimation uncertainty.")
  if (isTRUE(verbose)) message("Computing estimation bands (B1=", B1,
                               if (B2 > 0) paste0(", averaging over B2=", B2) else ", closed-form", ")...")
  TH <- .sbg_draw_thetas(B1, mu = unname(object$logits), Sigma = object$vcov_theta)

  SmatE <- matrix(NA_real_, nrow = B1, ncol = Tlast + 1L)
  for (bix in seq_len(B1)) {
    th <- TH[bix, ]
    i <- 1L
    m  <- .sbg_clip01(.sbg_sigmoid(th[i])); i <- i + 1L
    q  <- .sbg_clip01(.sbg_sigmoid(th[i])); i <- i + 1L
    dd <- 0
    if (isTRUE(object$flags["one_and_done"])) {
      dd <- .sbg_clip01(.sbg_sigmoid(th[i])); i <- i + 1L
    }
    cu <- 0
    if (isTRUE(object$flags["cure"])) {
      cu <- .sbg_clip01(.sbg_sigmoid(th[i])); i <- i + 1L
    }
    ab_sum <- (1/q) - 1
    aa <- max(m * ab_sum, 1e-12)
    bb <- max((1 - m) * ab_sum, 1e-12)


    if (B2 <= 0) {
      # Pure estimation: no process simulation, compute S directly
      SmatE[bix, ] <- .sbg_surv_prob_from_pars(aa, bb, dd, cu, n_last = Tlast)
    } else {
      # Average over B2 simulated paths for each theta draw
      S_acc <- replicate(B2, .sbg_surv_from_tenures(rsbg_mix(N0, aa, bb, dd, cu), Tlast))
      SmatE[bix, ] <- rowMeans(S_acc)
    }


  }
  S_hat <- colMeans(SmatE)
  res$estimation <- data.frame(
    time  = 0:Tlast,
    S_hat = to_scale(S_hat),
    lower = to_scale(apply(SmatE, 2, stats::quantile, probs = alpha/2)),
    upper = to_scale(apply(SmatE, 2, stats::quantile, probs = 1 - alpha/2)),
    part  = ifelse(0:Tlast <= (Tobs - 1), "fitted", "projected")
  )

  ## Both (param draw + single process draw)
  SmatB <- matrix(NA_real_, nrow = B1, ncol = Tlast + 1L)
  if (isTRUE(verbose)) message("Computing predictive (both) bands (B1=", B1, ", single process draw each)...")
  for (bix in seq_len(B1)) {
    th <- TH[bix, ]
    i <- 1L
    m  <- .sbg_clip01(.sbg_sigmoid(th[i])); i <- i + 1L
    q  <- .sbg_clip01(.sbg_sigmoid(th[i])); i <- i + 1L
    dd <- 0
    if (isTRUE(object$flags["one_and_done"])) { dd <- .sbg_clip01(.sbg_sigmoid(th[i])); i <- i + 1L }
    cu <- 0
    if (isTRUE(object$flags["cure"]))        { cu <- .sbg_clip01(.sbg_sigmoid(th[i])); i <- i + 1L }
    ab_sum <- (1/q) - 1
    aa <- max(m * ab_sum, 1e-12)
    bb <- max((1 - m) * ab_sum, 1e-12)
    SmatB[bix, ] <- .sbg_surv_from_tenures(rsbg_mix(N0, aa, bb, dd, cu), Tlast)
  }
  S_hat <- colMeans(SmatB)
  res$both <- data.frame(
    time  = 0:Tlast,
    S_hat = to_scale(S_hat),
    lower = to_scale(apply(SmatB, 2, stats::quantile, probs = alpha/2)),
    upper = to_scale(apply(SmatB, 2, stats::quantile, probs = 1 - alpha/2)),
    part  = ifelse(0:Tlast <= (Tobs - 1), "fitted", "projected")
  )

  res
}

