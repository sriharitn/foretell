% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmf_bdw.R
\name{BetadiscreteWeibull}
\alias{BetadiscreteWeibull}
\title{Beta discrete Weibull (bdw) Distribution Family Function}
\usage{
dbdw(x, shape1, shape2, shape3, log = FALSE)
pbdw(x, shape1, shape2, shape3, lower.tail = TRUE, log.p = FALSE)
qbdw(p, shape1, shape2, shape3, lower.tail = TRUE, log.p = FALSE)
rbdw(n, shape1, shape2, shape3)
}
\arguments{
\item{x}{Vector of non-negative integers for \code{dbdw} and \code{pbdw}.}

\item{p}{Vector of probabilities (0 <= p <= 1) for \code{qbdw}.}

\item{n}{Number of random variates to generate for \code{rbdw}.}

\item{shape1}{First shape parameter "a" (must be > 0).}

\item{shape2}{Second shape parameter "b" (must be > 0).}

\item{shape3}{Second shape parameter "c" (must be > 0).}

\item{log}{Logical; if TRUE, probabilities are returned on the log scale (for \code{dsbg} and \code{psbg}).}

\item{lower.tail}{Logical; if TRUE (default), probabilities are P(X <= x), otherwise P(X > x) (for \code{psbg}).}

\item{log.p}{Logical; if TRUE, probabilities are returned on the log scale (for \code{psbg} and \code{qsbg}).}
}
\value{
\itemize{
\item \code{dbdw}: A numeric vector of PMF values.
\item \code{pbdw}: A numeric vector of CDF values.
\item \code{qbdw}: A numeric vector of quantile values.
\item \code{rbdw}: A numeric vector of random variates.
}
}
\description{
Provides functions for the probability mass function (PMF), cumulative distribution function (CDF), quantile function, and random variate generation for the BdW distribution.
}
\examples{
# PMF example
dbdw(1:5, shape1 = 2, shape2 = 3,shape3 = 0.5)

# CDF example
pbdw(1:5, shape1 = 2, shape2 = 3, shape3 = 0.5)

# Quantile example
qbdw(c(0.1, 0.5, 0.9), shape1 = 2, shape2 = 3 , shape3 = 0.5)

# Random variates
rbdw(10, shape1 = 2, shape2 = 3,  shape3 = 0.5)

}
\references{
{Fader P, Hardie B. How to project customer retention. Journal of Interactive Marketing. 2007;21(1):76-90.}

{Fader P, Hardie B, Liu Y, Davin J, Steenburgh T. "How to Project Customer Retention" Revisited: The Role of Duration Dependence. Journal of Interactive Marketing. 2018;43:1-16.}
}
