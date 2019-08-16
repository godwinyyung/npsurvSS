#' Expected number of events
#'
#' Given one or two treatment arms, calculate the expected number
#' of events \eqn{d} at time \eqn{\tau}.
#'
#' @param arm0  object of class 'arm'.
#' @param arm1  object of class 'arm'.
#' @param tau   vector of times. Defaults to total study duration.
#' @examples
#' arm0 <- create_arm(size=120, accr_time=6, surv_scale=0.05, loss_scale=0.005, follow_time=12)
#' arm1 <- create_arm(size=120, accr_time=6, surv_scale=0.03, loss_scale=0.005, follow_time=12)
#' exp_events(arm0)
#' exp_events(arm0, arm1)
#' exp_events(arm0, tau=c(10,NA))
#' exp_events(arm0, arm1, tau=c(10,NA))
#' @seealso \code{\link{exp_duration}} for calculating time to achieve expected
#'   events d, \code{\link{create_arm}} and \code{\link{create_arm_lachin}}
#'   for creating an object of class 'arm'.
#' @export
exp_events <- function(arm0=NA, arm1=NA, tau=NA) {
  if (length(tau) == 1) { # only one calculation is desired
    d <- 0
    if (!is.na(arm0[1])) {
      if (is.na(tau)) {
        tau <- arm0$total_time
      }
      arm0$total_time   <- tau
      arm0$follow_time  <- arm0$total_time - arm0$accr_time
      d <- d + arm0$size * prob_event(arm0)
    }
    if (!is.na(arm1[1])) {
      if (is.na(tau)) {
        tau <- arm1$total_time
      }
      arm1$total_time   <- tau
      arm1$follow_time  <- arm1$total_time - arm1$accr_time
      d <- d + arm1$size * prob_event(arm1)
    }
    return(d)
  } else {
    sapply(tau, function(x) exp_events(arm0, arm1, x))
  }
}

#' Expected duration
#'
#' Given one or two treatment arms, calculate the time \eqn{\tau} at which
#' the expected number of events equals \eqn{d}.
#'
#' @param arm0  object of class 'arm'.
#' @param arm1  object of class 'arm'.
#' @param d     vector of number of events.
#' @param search_start  value at which the search for duration tau starts.
#' @param search_prec   value controlling the desired precision before
#'   terminating the search.
#' @param max_duration  maximum \eqn{\tau} for consideration.
#' @examples
#' arm0 <- create_arm(size=120, accr_time=6, surv_scale=0.05, loss_scale=0.005, follow_time=12)
#' arm1 <- create_arm(size=120, accr_time=6, surv_scale=0.03, loss_scale=0.005, follow_time=12)
#' exp_duration(arm0, d=61)
#' exp_duration(arm0, arm1, d=103)
#' exp_duration(arm0, d=c(35,61))
#' exp_duration(arm0, arm1, d=c(57,103))
#' @seealso \code{\link{exp_events}} for calculating expected events d at time tau,
#'   \code{\link{create_arm}} and \code{\link{create_arm_lachin}}
#'   for creating an object of class 'arm'.
#' @export
exp_duration <- function(arm0=NA, arm1=NA, d, search_start=10, search_prec=1e-2, max_duration=1000) {

  if (length(d) == 1) { # only one calculation is desired

    # test if desired number of events is attainable
    if (exp_events(arm0, arm1, tau=max_duration) < d) {
      stop(paste("Desired number of events is not attainable by time", max_duration), call.=F)
    }

    # continue, since desired number of events is attainable
    tau   <- search_start  # total study time
    delta <- search_start  # incremental increase or decrease to tau
    less  <- T           # expected num of events is less than d
    cont  <- T           # continue algorithm
    while(cont) {

      # calculate expected number of events
      temp <- exp_events(arm0, arm1, tau)

      # stop or update necessary parameters
      if (abs(delta) <= search_prec+1e-10) { # delta within desired precision
        if ((less & temp >= d) | (!less & temp <= d)) {
          cont <- F
        } else {
          tau <- tau + delta
        }
      } else { # delta not within desired precision
        if ((less & temp >= d) | (!less & temp <= d)) {
          delta <- -0.1 * delta
          less  <- !less
        }
        tau <- tau + delta
      }
    }

    # return calculated total study time
    return(tau)

  } else { # calculate multiple total study times
    sapply(d, function(x) exp_duration(arm0, arm1, x, search_start, search_prec, max_duration))
  }

}
