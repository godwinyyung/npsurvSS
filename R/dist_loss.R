#' Loss to follow-up
#'
#' Density, distribution function, hazard function, quantile function,
#' and random generation for the loss to follow-up distribution.
#'
#' @param x,q             vector of quantiles.
#' @param p               vector of probabilities.
#' @param n               number of observations.
#' @param arm             object of class 'arm'.
#' @param lower.tail      logical; if TRUE, probabilities are \eqn{P(X \le x)};
#'   otherwise, \eqn{P(X > x)}.
#' @return \code{dloss} gives the density, \code{ploss} gives the distribution
#'   function, \code{hloss} gives the hazard function, \code{qloss} gives the
#'   quantile function, and \code{rloss} generates random deviates.
#' @seealso \code{\link{create_arm}} and \code{\link{create_arm_lachin}}
#'   for creating an object of class 'arm'.
#' @export
dloss <- function(x, arm) {
  stats::dweibull(x, shape=arm$loss_shape, scale=1/arm$loss_scale)
}

#' @rdname dloss
#' @export
ploss <- function(q, arm, lower.tail=T) {
  stats::pweibull(q, shape=arm$loss_shape, scale=1/arm$loss_scale, lower.tail=lower.tail)
}

#' @rdname dloss
#' @export
hloss <- function(x, arm) {
  dloss(x, arm) / ploss(x, arm, lower.tail=F)
}

#' @rdname dloss
#' @export
qloss <- function(p, arm, lower.tail=T) {
  stats::qweibull(p, shape=arm$loss_shape, scale=1/arm$loss_scale, lower.tail=lower.tail)
}

#' @rdname dloss
#' @export
rloss <- function(n=1, arm) {
  stats::rweibull(n, shape=arm$loss_shape, scale=1/arm$loss_scale)
}
