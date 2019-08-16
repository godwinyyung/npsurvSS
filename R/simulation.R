#' Simulate complete data for a single arm
#'
#' Simulate the complete data for a single arm, including time to accrual,
#' event, and loss of follow-up. No cutoff (by number of events or time) is
#' applied. Hence, no patients are administratively censored.
#'
#' @param arm     object of class 'arm'.
#' @param label   numeric label for the simulated arm, e.g. 0 for control,
#'   1 for treatment
#' @return
#'   \item{arm}{label}
#'   \item{time.accr}{time to accrual}
#'   \item{time.obs}{time to observation from accrual}
#'   \item{time.total}{time to observation from start of study}
#'   \item{censor}{0=censor, 1=event}
#'   \item{reason}{event description ('[experience ]event', '[loss to ]followup', 'administration[ censoring]')}
#'   \item{time.surv}{time to event}
#'   \item{time.loss}{time to loss of follow-up}
#' @seealso \code{\link{create_arm}} for creating an object of class 'arm'.
#' @examples
#' arm0 <- create_arm(size=120, accr_time=6, surv_scale=0.05, loss_scale=0.005, follow_time=12)
#' simulate_arm(arm0, label=0)
#' @export
simulate_arm <- function(arm, label=1) {

  n         <- ceiling(arm$size)

  # simulate time of accrual, time to event, and time to loss of follow-up
  time.accr <- raccr(n, arm)
  time.surv <- rsurv(n, arm)
  time.loss <- rloss(n, arm)

  censor    <- 1 * (time.loss >= time.surv) # censoring indicator
  reason    <- ifelse(censor==1, "event", "followup")
  time.obs  <- pmin(time.loss, time.surv) # observed time (from enrollment)
  time.total<- time.accr + time.obs # observed time (from start of study)

  # create and return simulated data
  data      <- data.frame(arm = rep(label, n),
                          time.accr,
                          time.obs,
                          time.total,
                          censor,
                          reason,
                          time.surv,
                          time.loss,
                          stringsAsFactors = F)
  return(data)
}

#' Simulate a clinical trial
#'
#' Simulate a single- or two-arm clinical trial, where end of study (EOS)
#' is triggered after a number of events has been observed or a certain time
#' has elapsed. Whereas \code{simulate_arm} provides complete data for patients,
#' including time to event and loss of follow-up, \code{simulate_trial} mimicks
#' an actual survival study by providing only the observed time (minimum of time
#' to event or censoring) and censoring indicator.
#'
#' @param arm0      object of class 'arm'.
#' @param arm1      object of class 'arm'.
#' @param events    number of required events to trigger end of study; overrides
#'   study duration defined within \code{arm0} and \code{arm1}.
#' @param duration  time from first-patient-in to trigger end of study; overrides
#'   study duration defined within \code{arm0} and \code{arm1}. If both \code{events}
#'   and \code{duration} are specified, end of study is triggered by either criteria,
#'   whichever occurs first.
#' @return
#'   \item{arm}{0=\code{arm0}, 1=\code{arm1}}
#'   \item{time.accr}{time to accrual}
#'   \item{time.obs}{time to observation from accrual}
#'   \item{time.total}{time to observation from start of study}
#'   \item{censor}{0=censor, 1=event}
#'   \item{reason}{event description ('[experience ]event', '[loss to ]followup', 'administration[ censoring]')}
#' @seealso \code{\link{simulate_arm}} for simulating complete data for a single
#' arm, \code{\link{create_arm}} for creating an object of class 'arm'.
#' @examples
#' arm0 <- create_arm(size=120, accr_time=6, surv_scale=0.05, loss_scale=0.005, follow_time=12)
#' arm1 <- create_arm(size=120, accr_time=6, surv_scale=0.03, loss_scale=0.005, follow_time=12)
#' simulate_trial(arm0, duration=10)
#' simulate_trial(arm0, arm1, events=50)
#' @export
simulate_trial <- function(arm0=NA,
                           arm1=NA,
                           events=NA,
                           duration=Inf) {

  # Simulate arms
  data <- c()
  if (!is.na(arm0[1])) {
    data <- simulate_arm(arm0, label=0)
  }
  if (!is.na(arm1[1])) {
    data <- rbind(data, simulate_arm(arm1, label=1))
  }
  data      <- data[order(data$time.total), # order by first to last observation
                    c("arm", "time.accr", "time.obs", "time.total", "censor", "reason")] # keep only observable data

  # Calculate time of analysis
  if (is.na(events)) { # required events not defined
    duration <- min(duration, max(data$time.total)) # trigger = duration
  } else {
    if (events <= sum(data$censor)) { # events achievable with infinite follow-up
      index    <- which.max(cumsum(data$censor) >= events)
      duration <- min(duration, data$time.total[index]) # trigger = duration or events, whichever first
    } else { # events not achievable with infinite follow-up
      duration <- max(data$time.total) # trigger = duration
    }
  }

  # Update survival data
  censor.sfu              <- data$time.total > duration   # subjects being followed at time of censoring
  data$censor[censor.sfu] <- 0                # censor observations that go beyond study period
  data$reason[censor.sfu] <- "administration"
  events             <- sum(data$censor)
  data$time.total         <- data$time.total * (data$time.total<=duration) + duration * (data$time.total > duration)
  data$time.obs           <- data$time.total - data$time.accr

  return(data)
}
