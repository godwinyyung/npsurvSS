###
# CLASS ARM
###

# For a subject in the provided arm, calculate the probability he or
# she is observed to be at risk at time=teval after enrollment.
prob_risk <- function(arm, teval) {
  psurv(teval, arm, lower.tail=F) *
    ploss(teval, arm, lower.tail=F) *
    paccr(pmin(arm$accr_time, arm$total_time-teval), arm)
}

# For a subject in the provided arm, calculate the density of event
# at time=teval after enrollment.
dens_event <- function(arm, teval) {
  dsurv(teval, arm) *
    ploss(teval, arm, lower.tail=F) *
    paccr(pmin(arm$accr_time, arm$total_time-teval), arm)
}

# For a subject in the provided arm, calculate the probability he or
# she is observed to have experienced an event by time=teval after enrollment.
prob_event <- function(arm, tmin=0, tmax=arm$total_time) {
  UseMethod("prob_event", arm)
}

# prob_event for arm of class "arm"
prob_event.arm <- function(arm, tmin=0, tmax=arm$total_time) {
  l = length(tmax)
  if (l==1) {
    return(stats::integrate(function(x) dens_event(arm, x), lower=tmin, upper=tmax)$value)
  } else {
    if (length(tmin)==1) {
      tmin = rep(tmin, l)
    }
    return(sapply(seq(l), function(i) prob_event(arm, tmin[i], tmax[i])))
  }
}

# For a subject in the provided arm, calculate the probability he or
# she is observed to be censored by time=teval after enrollment due to
# loss of follow-up, administration, or either.
prob_censor <- function(arm, tmin=0, tmax=arm$total_time, reason="either") {
  UseMethod("prob_censor", arm)
}

# prob_censor for arm of class "arm"
prob_censor.arm <- function(arm, tmin=0, tmax=arm$total_time, reason="either") {
  l = length(tmax)
  if (l == 1) {
    prob <- 0
    if (reason %in% c("followup", "either")) {
      prob <- prob + stats::integrate(function(x) dloss(x, arm) *
                                        psurv(x, arm, lower.tail=F) *
                                        paccr(pmin(arm$accr_time, arm$total_time-x), arm),
                                      lower=tmin,
                                      upper=tmax)$value
    }
    if (reason %in% c("administration", "either")) {
      prob <- prob + stats::integrate(function(x) daccr(arm$total_time-x, arm) *
                                        psurv(x, arm, lower.tail=F) *
                                        ploss(x, arm, lower.tail=F),
                                      lower=max(tmin, arm$follow_time),
                                      upper=tmax)$value
    }
    return(prob)
  } else {
    if (length(tmin)==1) {
      tmin = rep(tmin, l)
    }
    return(sapply(seq(l), function(i) prob_censor(arm, tmin[i], tmax[i], reason)))
  }
}

###
# CLASS LACHIN
###

# prob_event for arm of class "lachin"
prob_event.lachin <- function(arm, tmin=0, tmax=arm$total_time) {
  tmid <- pmax(pmin(tmax, arm$follow_time), tmin)
  rate <- arm$surv_scale + arm$loss_scale

  if (arm$accr_dist=="pieceuni") { # uniform
    arm$surv_scale / rate *
      (exp(-rate * tmin) -
         exp(-rate * tmid) *
         (1 - (arm$total_time - tmid) / arm$accr_time + 1 / rate / arm$accr_time) -
         exp(-rate * tmax) *
         ((arm$total_time - tmax) / arm$accr_time - 1 / rate / arm$accr_time)
      )
  } else { # truncated exponential
    rate2 <- rate - arm$accr_param
    arm$surv_scale / rate * (exp(-rate * tmin) - exp(-rate * tmid)) +
      arm$surv_scale / rate / (1 - exp(-arm$accr_param * arm$accr_time)) * (exp(-rate * tmid) - exp(-rate * tmax)) +
      arm$surv_scale / rate2 / (1 - exp(-arm$accr_param * arm$accr_time)) * (exp(-rate2 * tmax - arm$accr_param * arm$total_time) - exp(-rate2 * tmid - arm$accr_param * arm$total_time))
  }
}

# prob_censor for arm of class "lachin"
prob_censor.lachin <- function(arm, tmin=0, tmax=arm$total_time, reason="either") {
  prob <- 0
  prob_event <- prob_event(arm, tmin, tmax)
  if (reason %in% c("followup", "either")) {
    prob <- prob + arm$loss_scale / arm$surv_scale * prob_event
  }
  if (reason %in% c("administration", "either")) {
    prob <- prob + 1 - prob_event * (1 + arm$loss_scale / arm$surv_scale)
  }
  return(prob)
}
