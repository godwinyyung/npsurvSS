#' Maximum observed time
#'
#' Density, distribution function, and expected value for the maximum observed
#' time in a single arm of patients.
#'
#' @param x,q             vector of quantiles.
#' @param arm             object of class 'arm'.
#' @param include_cens    logical; if TRUE, include time-to-censoring as potential
#'   observed time; otherwise, observed time equals time-to-event.
#' @param lower.tail      logical; if TRUE, probabilities are \eqn{P(X \le x)};
#'   otherwise, \eqn{P(X > x)}.
#' @return \code{dmaxU} gives the density, \code{pmaxU} gives the distribution
#'   function, and \code{emaxU} gives the expected value.
#' @details Given a patient's time-to-event \eqn{T_i} and time-to-censoring
#' \eqn{C_i}, \eqn{U_i=\min(T_i, C_i)} defines the patient's observed time. The
#' maximum observed time over patients of a single arm is then \eqn{\max_i U_i}.
#' @seealso \code{\link{create_arm}} and \code{\link{create_arm_lachin}}
#'   for creating an object of class 'arm'.
#' @export
dmaxU <- function(x, arm, include_cens=T) {
  if (include_cens) {
    arm$size * (1 - prob_risk(arm, x)) ^ (arm$size - 1) *
      ( dsurv(x, arm) * ploss(x, arm, lower.tail=F) * paccr(pmin(arm$accr_time, arm$total_time-x), arm) +
          psurv(x, arm, lower.tail=F) * dloss(x, arm) * paccr(pmin(arm$accr_time, arm$total_time-x), arm) +
          psurv(x, arm, lower.tail=F) * ploss(x, arm, lower.tail=F) * daccr(pmin(arm$accr_time, arm$total_time-x), arm) * (x > arm$follow_time)
      )
  } else {
    arm$size * (1 - (prob_event(arm) - prob_event(arm, tmax=x))) ^ (arm$size - 1) *
      dens_event(arm, x)
  }
}

#' @rdname dmaxU
#' @export
pmaxU <- function(q, arm, include_cens=T, lower.tail=T) {
  if (include_cens) {
    cdf <- (1 - prob_risk(arm, q)) ^ arm$size
  } else {
    cdf <- (1 - (prob_event(arm) - prob_event(arm, tmax=q))) ^ arm$size
  }
  lower.tail * cdf + (1 - lower.tail) * (1-cdf)
}

#' @rdname dmaxU
#' @export
emaxU <- function(arm, include_cens=T) {
  stats::integrate(function(x) x * dmaxU(x, arm, include_cens=include_cens),
                   lower=0,
                   upper=arm$total_time)$value
}

#' Minimax observed time
#'
#' Density, distribution function, quantile function, and expected value for the
#' minimum of the maximum observed time over two treatment arms.
#'
#' @param x,q             vector of quantiles.
#' @param arm0            object of class 'arm'.
#' @param arm1            object of class 'arm'.
#' @param include_cens    logical; if TRUE, include time-to-censoring as potential
#'   observed time; otherwise, observed time equals time-to-event.
#' @param lower.tail      logical; if TRUE, probabilities are \eqn{P(X \le x)};
#'   otherwise, \eqn{P(X > x)}.
#' @param p               vector of probabilities.
#' @param margin          margin of accuracy.
#' @return \code{dminimaxU} gives the density, \code{pminimaxU} gives the distribution
#'   function, \code{qminimaxU} gives the quantile function, and \code{eminimaxU}
#'   gives the expected value.
#' @details Given a patient in arm \eqn{X_i=j} with time-to-event \eqn{T_i} and time-to-censoring
#' \eqn{C_i}, \eqn{U_i=\min(T_i, C_i)} defines the patient's observed time. The
#' maximum observed time over patients of arm \eqn{j} is then \eqn{\max_{i:X_i=j} U_i},
#' and the minimax observed time over two arms is \eqn{\min_j (\max_{i:X_i=j} U_i)}.
#' @seealso \code{\link{create_arm}} and \code{\link{create_arm_lachin}}
#'   for creating an object of class 'arm'.
#' @export
dminimaxU <- function(x, arm0, arm1, include_cens=T) {
  dmaxU(x, arm0, include_cens=include_cens) *
    pmaxU(x, arm1, include_cens=include_cens, lower.tail=F) +
    dmaxU(x, arm1, include_cens=include_cens) *
    pmaxU(x, arm0, include_cens=include_cens, lower.tail=F)
}

#' @rdname dminimaxU
#' @export
pminimaxU <- function(q, arm0, arm1, include_cens=T, lower.tail=T) {
  surv <- pmaxU(q, arm0, include_cens=include_cens, lower.tail=F) *
    pmaxU(q, arm1, include_cens=include_cens, lower.tail=F)
  lower.tail * (1 - surv) + (1 - lower.tail) * surv
}

#' @rdname dminimaxU
#' @export
qminimaxU <- function(p, arm0, arm1, include_cens=T, margin=0.01) {
  time.vec <- seq(0, arm0$total_time, margin)
  if (length(p) == 1) {
    time.vec[which.min(pminimaxU(q=time.vec, arm0, arm1, include_cens=include_cens) < p)-1]
  } else {
    sapply(p, function(p2) qminimaxU(p2, arm0, arm1, include_cens=include_cens))
  }
}

#' @rdname dminimaxU
#' @export
eminimaxU <- function(arm0, arm1, include_cens=T) {
  stats::integrate(function(x) x * dminimaxU(x, arm0, arm1, include_cens=include_cens),
                   lower=0,
                   upper=arm0$total_time)$value
}
