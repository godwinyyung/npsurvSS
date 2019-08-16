#' Survival
#'
#' Density, distribution function, hazard function, quantile function,
#' and random generation for the survival distribution.
#'
#' @param x,q             vector of quantiles.
#' @param p               vector of probabilities.
#' @param n               number of observations.
#' @param arm             object of class arm.
#' @param include_cured   logical; if TRUE, mixture distribution of cured
#'   and uncured patients is considered; otherwise, only the distribution
#'   for uncured patients is considered.
#' @param lower.tail      logical; if TRUE, probabilities are \eqn{P(X \le x)};
#'   otherwise, \eqn{P(X > x)}.
#' @seealso \code{\link{create_arm}} and \code{\link{create_arm_lachin}}
#'   for creating an object of class 'arm'.
#' @export
dsurv <- function(x, arm, include_cured=T) {

  # pdf for uncured patients
  if (length(arm$surv_interval)==2  ) { # Weibull
    dens <- stats::dweibull(x,
                            shape=arm$surv_shape,
                            scale=1/arm$surv_scale)
  } else { # Piecewise-exponential
    dens <- hsurv(x, arm, include_cured=F) * psurv(x, arm, lower.tail=F, include_cured=F)
  }

  # pdf for overall population
  if (include_cured) {
    dens <- (1 - arm$surv_cure) * dens
  }

  return(dens)
}

#' @rdname dsurv
#' @export
psurv <- function(q, arm, include_cured=T, lower.tail=T) {

  # cdf for uncured patients
  if ( length(arm$surv_interval)==2 ) { # Weibull
    surv <- stats::pweibull(q,
                            shape=arm$surv_shape,
                            scale=1/arm$surv_scale,
                            lower.tail=F)
  } else { # Piecewise-exponential
    surv <- as.vector(exp( -interval_ptime(q, arm$surv_interval) %*% arm$surv_scale ))

  }

  # cdf for overall population
  if (include_cured) {
    surv <- arm$surv_cure + (1 - arm$surv_cure) * surv
  }

  return(lower.tail * (1-surv) + (1-lower.tail) * surv) # cdf or 1-cdf
}

#' @rdname dsurv
#' @export
hsurv <- function(x, arm, include_cured=T) {

  # hazard for uncured patients
  if (length(arm$surv_interval)==2  ) { # Weibull
    hazard <- pmax(0,
                   (arm$surv_scale * arm$surv_shape) * (arm$surv_scale * x)^(arm$surv_shape - 1))
  } else { # piecewise-exponential
    hazard <- c(0, arm$surv_scale)[rowSums(interval_atrisk(x, arm$surv_interval))+1]
  }

  # hazard function for overall population
  if (include_cured) {
    surv_uncured <- psurv(x, arm, lower.tail=F, include_cured=F)
    hazard <- hazard * (1-arm$surv_cure) * surv_uncured /
      (arm$surv_cure + (1-arm$surv_cure) * surv_uncured)
  }

  return(hazard)
}

#' @rdname dsurv
#' @export
qsurv <- function(p, arm, include_cured=T, lower.tail=T) {

  if (!include_cured) {
    arm$surv_cure <- 0
  }
  if (!lower.tail) {
    p <- 1 - p
  }

  if ( length(arm$surv_interval)==2 ) { # Weibull
    p <- ifelse(p > (1 - arm$surv_cure), 1, 1 - (1 - arm$surv_cure - p) / (1 - arm$surv_cure))
    q <- stats::qweibull(p, shape=arm$surv_shape, scale=1/arm$surv_scale)
  } else { # Piecewise-exponential
    p2      <- 1 - (1 - arm$surv_cure - p) / (1 - arm$surv_cure)
    k       <- length(p2)
    interval_probrisk <- c(1, cumprod( exp(-arm$surv_scale * diff(arm$surv_interval)) )) # probability of being at risk at the beginning of each survival interval
    interval_occur    <- rowSums(t(matrix(rep(interval_probrisk, k), ncol=k)) >= (1-p2))
    margin  <- suppressWarnings(-log((1-p2) / interval_probrisk[interval_occur]) / arm$surv_scale[interval_occur])
    q       <- arm$surv_interval[interval_occur] + margin
    q       <- ifelse(p==0, -Inf,
                      ifelse(p < 1 - arm$surv_cure, q, Inf))
  }

  return(q)

}

#' @rdname dsurv
#' @export
rsurv <- function(n=1, arm, include_cured=T) {
  # inverse transform sampling
  qsurv(p=stats::runif(n),
        arm=arm,
        include_cured=include_cured)
}
