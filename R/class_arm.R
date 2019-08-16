#' Create an 'arm' object
#'
#' Create an object of class 'arm' by specifying the trial parameters for
#' a single arm, including the sample size, accrual distribution,
#' survival distribution, loss to follow-up distribution, and study duration.
#'
#' @param size sample size. If total sample size is unknown, provide
#'   the integer sample size relative to the opposing arm, e.g. 1 for
#'   1:2 randomization ratio or 2 for 2:3.
#' @param accr_time accrual duration.
#' @param accr_dist accrual distribution. Default is piecewise uniform.
#'   Alternatively, 'truncexp' allows for a truncated exponential distribution
#'   as proposed by Lachin and Foulkes (1986). Depending on the value of
#'   \code{accr_param}, this distribution can be either convex or concave.
#' @param accr_interval accrual intervals. Defaults to the single interval
#'   spanning from 0 to \code{accr_time}. If a piecewise uniform accrual
#'   with more than one interval is desired, specify \code{accr_interval}
#'   as the vector of increasing changepoints (knots) starting from 0 and
#'   ending with \code{accr_time}, e.g. c(0, 2, 4) defines a piecewise
#'   uniform distribution with two intervals, [0, 2) and [2, 4].
#' @param accr_param additional accrual parameter(s). For a piecewise uniform accrual
#'   with more than one interval, specify \code{accr_param} as the vector of probabilities
#'   a patient is enrolled in each interval. The probabilities should naturally sum to 1.
#'   For \code{accr_dist}='truncexp', specify \code{accr_param} as a single number.
#'   >0 results in a convex distribution and <0 results in a concave distribution.
#' @param surv_cure proportion of patients that are cured.
#' @param surv_interval survival intervals. Defaults to the single interval
#'   spanning from 0 to infinity. If a piecewise exponential survival is
#'   desired for uncured patients, specify \code{surv_interval} as the vector
#'   of increasing changepoints (knots) starting from 0 and ending with infinity,
#'   e.g. c(0, 6, 10, Inf).
#' @param surv_shape Weibull shape parameter for the survival distribution of
#'   uncured patients.
#' @param surv_scale Weibull scale parameter for the survival distrubition of
#'   uncured patients. Piecewise exponential survival may be defined by specifying
#'   \code{surv_shape}=1 and \code{surv_scale} as the vector of piecewise hazard
#'   rates.
#' @param loss_shape Weibull shape parameter for the loss to follow-up distribution.
#' @param loss_scale Weibull scale parameter for the loss to follow-up distribution.
#' @param follow_time follow-up duration.
#' @param total_time total study duration. Only 1 of the 2 parameters, \code{follow_time}
#'   or \code{total_time}, need to be defined. If neither is defined, \code{total_time}
#'   is defaulted to max value 1e6.
#' @examples
#' # Example 1
#' example <- create_arm(size=120,
#'   accr_time=6,                   # uniform accrual
#'   surv_scale=0.05,               # exponential survival
#'   loss_scale=0.005,              # exponential loss to follow-up
#'   follow_time=12)
#' class(example)                   # this example also satisfies properties of subclass 'lachin'
#'
#' # Example 2
#' create_arm(size=120,
#'   accr_time=6,                   # truncated exponential accrual
#'   accr_dist="truncexp",
#'   accr_param=0.1,
#'   surv_shape=2,                  # weibull survival
#'   surv_scale=0.05,
#'   loss_shape=1.5,                # weilbull loss to follow-up
#'   loss_scale=0.005,
#'   total_time=18)
#'
#' # Example 3
#' create_arm(size=120,
#'   accr_time=6,
#'   accr_interval=c(0,2,4,6),      # piecewise uniform accrual
#'   accr_param=c(0.2,0.3,0.5),
#'   surv_cure=0.1,                 # 10% cure fraction
#'   surv_interval=c(0,6,10,Inf),   # piecewise exponential survival for uncured patients
#'   surv_scale=c(0.05,0.04,0.03),
#'   loss_shape=0.7,                # weibull loss to follow-up
#'   loss_scale=0.005,
#'   total_time=18)
#' @seealso \code{\link{create_arm_lachin}} for creating an object of subclass 'lachin'.
#' @references
#' Lachin, J. M. and Foulkes, M. A. (1986) Evaluation of sample size and power for analyses of
#'   survival with allowance for nonuniform patient entry, losses to follow-up,
#'   noncompliance, and stratification. \emph{Biometrics}, \strong{42}, 507-519.
#' @export
create_arm <- function(size,
                       accr_time,
                       accr_dist = "pieceuni",
                       accr_interval = c(0, accr_time),
                       accr_param = NA,
                       surv_cure = 0,
                       surv_interval = c(0, Inf),
                       surv_shape=1,
                       surv_scale,
                       loss_shape=1,
                       loss_scale,
                       follow_time = Inf,
                       total_time = Inf) {

  # Check accr_dist
  if (! accr_dist %in% c("pieceuni", "truncexp")) {
    stop("Please specify a valid accrual distribution.", call.=F)
  }

  # Ensure accr_interval is properly defined
  accr_interval = sort(unique(c(0, accr_interval, accr_time)))
  if (min(accr_interval) < 0 | max(accr_interval) > accr_time) {
    stop("accr_interval is out of range.", call.=F)
  }

  # Check accr_param
  if (accr_dist == "pieceuni") {
    if (length(accr_param) != length(accr_interval) - 1) {
      stop("Number of accrual intervals (accr_interval) does not match number of
         accrual parameters (accr_param).", call.=F)
    }
    if (length(accr_interval) > 2 & sum(accr_param) != 1) {
      stop("accr_param must sum to 1.", call.=F)
    }
  } else if (is.na(accr_param) | length(accr_param) > 1) {
    stop("Truncated exponential is a one-parameter family distribution.", call.=F)
  }

  # Ensure surv_interval is properly defined
  surv_interval = sort(unique(c(0, surv_interval, Inf)))
  if (min(surv_interval) < 0) {
    stop("surv_interval is out of range.", call.=F)
  }

  # Check surv_shape
  if (surv_shape != 1 & length(surv_scale) > 1) {
    surv_shape = 1
    warning("Piecewise weibull is not supported. surv_shape defaulted to 1.", call.=F)
  }

  # Check surv_scale
  if (length(surv_scale) != length(surv_interval) - 1) {
    stop("Number of survival intervals (surv_interval) does not match number of
         piecewise hazards (surv_scale).", call.=F)
  }

  # Check loss to follow-up distribution
  if (length(loss_shape) > 1 | length(loss_scale) > 1) {
    loss_shape = loss_shape[1]
    loss_scale = loss_scale[1]
    warning("Only Weibull loss to follow-up is supported. First number in loss_shape
            and loss_scale are considered. The rest are ignored.", call.=F)
  }

  # Check total study time
  if (is.infinite(follow_time) & is.infinite(total_time)) {
    total_time = 1e6
    follow_time = total_time - accr_time
    warning("Neither follow_time nor total_time were defined. Therefore, total_time is
            defaulted to max value.", call.=F)
  } else if (!is.infinite(follow_time) & !is.infinite(total_time) & accr_time+follow_time != total_time) {
    total_time = accr_time + follow_time
    warning("follow_time and total_time were inconsistently defined.
            total_time will be ignored.", call.=F)
  } else if (is.infinite(follow_time)) {
    follow_time = total_time - accr_time
  } else {
    total_time = accr_time + follow_time
  }

  # prepare object for return
  arm <- list(size = size,
              accr_time = accr_time,
              accr_dist = accr_dist,
              accr_interval = accr_interval,
              accr_param = accr_param,
              surv_cure = surv_cure,
              surv_interval = surv_interval,
              surv_shape = surv_shape,
              surv_scale = surv_scale,
              loss_shape = loss_shape,
              loss_scale = loss_scale,
              follow_time = follow_time,
              total_time = total_time)

  # check if object satisfies properties of class 'lachin'
  if (length(accr_param)==1 & # uniform or truncated exponential accrual
      length(surv_interval)==2 & surv_shape==1 & # exponential survival
      loss_shape==1) { # exponential loss to follow-up
    class(arm) <- append(class(arm), "lachin")
  }

  # attach class 'arm' to object
  class(arm) <- append(class(arm), "arm")

  return(arm)

}

#' Create a 'lachin' object
#'
#' Create an object of class 'lachin' by specifying the trial parameters for
#' a single arm, including the sample size, accrual distribution,
#' survival distribution, loss to follow-up distribution, and study duration.
#' 'Lachin' objects are also 'arm' objects, but with accrual limited
#' to the uniform and truncated exponential distributions, and survival and
#' loss to follow-up limited to the exponential distribution. 'Lachin' objects
#' have the advantage that expectations for certain counting processes have
#' closed form equations and can therefore be calculated more efficiently
#' (Lachin and Foulkes, 1986).
#'
#' @param size sample size. If total sample size is unknown, provide
#'   the integer sample size relative to the opposing arm, e.g. 1 for
#'   1:2 randomization ratio or 2 for 2:3.
#' @param accr_time accrual duration.
#' @param accr_dist accrual distribution. Default is uniform (piecewise uniform
#'   with one interval). Alternatively, 'truncexp' allows for a truncated
#'   exponential distribution as proposed by Lachin and Foulkes (1986).
#'   Depending on the value of \code{accr_param}, this distribution can be
#'   either convex or concave.
#' @param accr_param additional accrual parameter for \code{accr_dist}='truncexp'.
#'\code{accr_param}>0 specifies a convex distribution and \code{accr_param}<0
#'specifies a concave distribution.
#' @param surv_median median survival.
#' @param surv_exphazard exponential hazard rate for the survival distribution.
#' @param surv_milestone a tuple c(milestone, probability) that uniquely defines the
#'   exponential survival distribution, e.g. c(12, 0.8) corresponds to the exponential
#'   distribution with 80\% survival rate at 12 months.
#' @param loss_median median loss to follow-up.
#' @param loss_exphazard exponential hazard rate for the loss to follow-up distribution.
#' @param loss_milestone a tuple c(milestone, probability) that uniquely defines the
#'   exponential loss to follow-up distribution, e.g. c(12, 0.99) corresponds to the
#'   exponential distribution with 1\% loss to follow-up at 12 months.
#' @param follow_time Follow-up duration. Either \code{follow_time} or
#'   \code{total_time} (below) should be specified.
#' @param total_time Total study duration. Either \code{follow_time} (above) or
#'   \code{total_time} should be specified.
#' @seealso \code{\link{create_arm}} for creating an object of class 'arm'.
#' @references Lachin, J. M. and Foulkes, M. A. (1986) Evaluation of sample size and power for analyses of
#'   survival with allowance for nonuniform patient entry, losses to follow-up,
#'   noncompliance, and stratification. \emph{Biometrics}, \strong{42}, 507-519.
#' @examples
#' # 3 arms with similar survival and loss to follow-up
#' create_arm_lachin(size=120, accr_time=6,
#'   surv_median=14,
#'   loss_median=140,
#'   follow_time=12)
#' create_arm_lachin(size=120, accr_time=6,
#'   surv_exphazard=0.05,
#'   loss_exphazard=0.005,
#'   follow_time=12)
#' create_arm_lachin(size=120, accr_time=6,
#'   accr_dist="truncexp",
#'   accr_param=0.1,
#'   surv_milestone=c(14, 0.5),
#'   loss_milestone=c(140, 0.5),
#'   total_time=18)
#' @export
create_arm_lachin <- function(size,
                              accr_time,
                              accr_dist = "pieceuni",
                              accr_param = NA,
                              surv_median = NA,
                              surv_exphazard = NA,
                              surv_milestone = NA,
                              loss_median = NA,
                              loss_exphazard = NA,
                              loss_milestone = NA,
                              follow_time = Inf,
                              total_time = Inf) {

  # Check accr_param
  if (accr_dist == "pieceuni" & !is.na(accr_param)) {
    accr_param = NA
    warning("accr_param is ignored.", call.=F)
  }

  # Survival
  if (sum(!is.na(c(surv_median, surv_exphazard, surv_milestone[1]))) > 1) {
    stop("Please specify just one of surv_median, surv_exphazard, or surv_milestone.", call.=F)
  } else if (!is.na(surv_median)) {
    surv_scale = per2haz(surv_median)
  } else if (!is.na(surv_exphazard)) {
    surv_scale = surv_exphazard
  } else {
    surv_scale = per2haz(x=surv_milestone[1], per=1-surv_milestone[2])
  }

  # Loss to follow-up
  if (sum(!is.na(c(loss_median, loss_exphazard, loss_milestone[1]))) > 1) {
    stop("Please specify just one of loss_median, loss_exphazard, or loss_milestone.", call.=F)
  } else if (!is.na(loss_median)) {
    loss_scale = per2haz(loss_median)
  } else if (!is.na(loss_exphazard)) {
    loss_scale = loss_exphazard
  } else {
    loss_scale = per2haz(x=loss_milestone[1], per=1-loss_milestone[2])
  }

  arm <- create_arm(size = size,
                    accr_time = accr_time,
                    accr_dist = accr_dist,
                    accr_param = accr_param,
                    surv_scale = surv_scale,
                    loss_scale = loss_scale,
                    follow_time = follow_time,
                    total_time = total_time)

  return(arm)

}

#' Convert exponential parameters
#'
#' Convert exponential survival percentile or hazard rate to the other.
#' @param x     survival percentile or exponential hazard rate
#' @param per   (per)th percentile
#' @details \deqn{y=-log(1-per)/x}
#' @examples
#' per2haz(14)          # hazard rate for exponential with 14-month median
#' per2haz(0.05)        # median survival for exponential with hazard rate 0.05
#' per2haz(14, 0.8)     # hazard rate for exponential with 80th percentile survival at 14 months
#' per2haz(0.27, 0.8)   # 80th percentile survival for exponential with hazard rate 0.27
#' @export
per2haz <- function(x, per=0.5) {
  -log(1-per)/x
}
