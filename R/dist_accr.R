#' Accrual
#'
#' Density, distribution function, quantile function,
#' and random generation for the accrual distribution.
#'
#' @param x,q             vector of quantiles.
#' @param p               vector of probabilities.
#' @param n               number of observations.
#' @param arm             object of class 'arm'.
#' @param lower.tail      logical; if TRUE, probabilities are \eqn{P(X \le x)};
#'   otherwise, \eqn{P(X > x)}.
#' @return \code{daccr} gives the density, \code{paccr} gives the distribution
#'   function, \code{qaccr} gives the quantile function, and \code{raccr} generates
#'   random deviates.
#' @seealso \code{\link{create_arm}} and \code{\link{create_arm_lachin}}
#'   for creating an object of class 'arm'.
#' @export
daccr <- function(x, arm) {

  if (arm$accr_dist=="pieceuni" & length(arm$accr_param)==1) { # Uniform
    stats::dunif(x, min=0, max=arm$accr_time)
  } else if (arm$accr_dist=="pieceuni") { # Piecewise uniform
    (x <= arm$accr_time) *
      c(0, arm$accr_param / diff(arm$accr_interval))[pmin(rowSums(interval_atrisk(x, arm$accr_interval)),
                                                          length(arm$accr_interval)-1)+1]
  } else  { # Truncated-exponential
    (x >= 0) *
      (x <= arm$accr_time) *
      arm$accr_param *
      exp(-arm$accr_param * x) /
      (1 - exp(-arm$accr_param * arm$accr_time))
  }

}

#' @rdname daccr
#' @export
paccr <-function(q, arm, lower.tail=T) {

  q <- pmax(pmin(q, arm$accr_time), 0)

  if (arm$accr_dist=="pieceuni" & length(arm$accr_param)==1) { # Uniform
    cdf <- stats::punif(q=q, min=0, max=arm$accr_time)
  } else if (arm$accr_dist=="pieceuni") { # Piecewise-uniform
    cdf <- ifelse(q < max(arm$accr_interval),
                  rowSums(t(t(interval_ptime(q, arm$accr_interval)) /
                              diff(arm$accr_interval) *
                              arm$accr_param)),
                  1)
  } else { # Truncated-exponential
    cdf <- (1 - exp(-arm$accr_param * q)) / (1 - exp(-arm$accr_param * arm$accr_time))
  }

  lower.tail * cdf + (1 - lower.tail) * (1 - cdf)

}

#' @rdname daccr
#' @export
qaccr <- function(p, arm) {

  if (arm$accr_dist=="pieceuni" & length(arm$accr_param)==1) { # Uniform
    p * arm$accr_time
  } else if (arm$accr_dist=="pieceuni") { # Piecewise uniform
    interval_probrisk <- c(1, utils::head(1-cumsum(arm$accr_param), -1)) # prob of being at risk (for accrual) at the beginning of each accrual interval
    interval_occur    <- rowSums(t(matrix(rep(interval_probrisk, length(p)), ncol=length(p))) >= (1-p)) # interval in which accrual occured
    margin <- (p - 1 + interval_probrisk[interval_occur]) / # (margin time width) = (margin prob width) /
      - diff(c(interval_probrisk,0))[interval_occur] * # (interval probability width) *
      diff(arm$accr_interval)[interval_occur] # (interval time width)
    arm$accr_interval[interval_occur] + margin
  } else { # Truncated exponential
    -log( 1 - p * (1 - exp(-arm$accr_param * arm$accr_time)) ) / arm$accr_param
  }

}

#' @rdname daccr
#' @export
raccr <- function(n=1, arm) {
  qaccr(stats::runif(n), arm)
}
