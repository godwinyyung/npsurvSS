###
# Kaplan-Meier based tests
###

# Cumulative hazard
sigma2j_cumh <- function(arm, teval) {
  if (teval < 0 | teval >= arm$total_time) {
    stop(paste("Time of evaluation must be in the interval [0, ", arm$total_time, ").", sep=""))
  }

  stats::integrate(function(x) hsurv(x, arm) / prob_risk(arm, x),
                   lower=0,
                   upper=teval)$value
}

# Milestone survival
sigma2j_surv <- function(arm, teval) {
  psurv(teval, arm, lower.tail=F)^2 * sigma2j_cumh(arm, teval)
}

# Percentile survival
sigma2j_perc <- function(arm, perc) {
  teval <- qsurv(perc, arm)
  sigma2j_cumh(arm, teval) / hsurv(teval, arm)^2
}

# Restricted mean survival time
deltaj_rmst <- function(x, arm) {
  stats::integrate(function(y) psurv(y, arm, lower.tail=F),
                   lower=0,
                   upper=x)$value
}

sigma2j_rmst <- function(arm, teval) {
  inner <- function(x) {
    sapply(x, function(x1) stats::integrate(function(x2) psurv(x2, arm, lower.tail=F),
                                            lower=x1,
                                            upper=teval)$value
    )
  }
  stats::integrate(function(x) inner(x)^2 * hsurv(x, arm) / prob_risk(arm, x),
                   lower=0,
                   upper=teval)$value
}

###
# Weighted log-rank test
###

weight_wlr <- function(x, arm0, arm1, weight="1") { # weight function

  n   <- arm0$size + arm1$size
  p1 <- arm1$size / n
  p0 <- 1 - p1

  if (weight=="1") {
    1
  } else if (weight=="n") {
    n * (p0 * prob_risk(arm0, x) + p1 * prob_risk(arm1, x))
  } else if (weight=="sqrtN") {
    sqrt( n * (p0 * prob_risk(arm0, x) + p1 * prob_risk(arm1, x)) )
  } else if (grepl("FH", weight)) {
    weight <- strsplit(weight, "_")[[1]]
    p <- as.numeric(substring(weight[2], 2))
    q <- as.numeric(substring(weight[3], 2))
    esurv <- p0 * psurv(x, arm0) + p1 * psurv(x, arm1)
    esurv^p * (1-esurv)^q
  } else {
    stop("Please specify valid weight function.")
  }
}

delta_wlr <- function(arm0, arm1, weight="1", approx="asymptotic") {

  p1 <- arm1$size / (arm0$size + arm1$size)
  p0 <- 1 - p1

  if (approx == "event driven") {
    if (sum(arm0$surv_shape != arm1$surv_shape) > 0 |
        length( unique(arm1$surv_scale / arm0$surv_scale) ) > 1) {
      stop("Hazard is not proportional over time.", call.=F)
    } else if (weight != "1") {
      stop("Weight must equal `1`.", call.=F)
    }
    theta <- c(arm0$surv_shape * log( arm1$surv_scale / arm0$surv_scale ))[1] # log hazard ratio
    nu    <- p0 * prob_event(arm0) + p1 * prob_event(arm1) # probability of event
    delta <- theta * p0 * p1 * nu
  } else if (approx == "asymptotic") {
    delta <- stats::integrate(function(x) weight_wlr(x, arm0, arm1, weight) *
                                (1 / p0 / prob_risk(arm0, x) + 1 / p1 / prob_risk(arm1, x)) ^ (-1) *
                                ( hsurv(x, arm1) - hsurv(x, arm0) ),
                              lower=0,
                              upper=arm0$total_time)$value
  } else if (approx == "generalized schoenfeld") {
    delta <- stats::integrate(function(x) weight_wlr(x, arm0, arm1, weight) *
                                log( hsurv(x, arm1) / hsurv(x, arm0) ) *
                                p0 * prob_risk(arm0, x) * p1 * prob_risk(arm1, x) /
                                ( p0 * prob_risk(arm0, x) + p1 * prob_risk(arm1, x) )^2 *
                                ( p0 * dens_event(arm0, x) + p1 * dens_event(arm1, x)),
                              lower=0,
                              upper=arm0$total_time)$value
  } else {
    stop("Please specify a valid approximation for the mean.", call.=F)
  }

  return(delta)

}

deltaj_wlr <- function(j, arm0, arm1, weight="1") {

  p1 <- arm1$size / (arm0$size + arm1$size)
  p0 <- 1 - p1
  if (j==0) {
    arm <- arm0
  } else {
    arm <- arm1
  }
  stats::integrate(function(x) weight_wlr(x, arm0, arm1, weight) *
                     (j - p1 * prob_risk(arm1, x) /
                        ( p0 * prob_risk(arm0, x) + p1 * prob_risk(arm1, x) )) *
                     (dens_event(arm, x) -
                        prob_risk(arm, x) *
                        ( p0 * dens_event(arm0, x) + p1 * dens_event(arm1, x) ) /
                        ( p0 * prob_risk(arm0, x) + p1 * prob_risk(arm1, x) )),
                   lower=0,
                   upper=arm0$total_time)$value

}

sigma2_wlr <- function(arm0, arm1, weight="1", approx="asymptotic") {

  p1 <- arm1$size / (arm0$size + arm1$size)
  p0 <- 1 - p1

  if (approx == "event driven") {
    nu      <- p0 * prob_event(arm0) + p1 * prob_event(arm1)
    sigma2  <- p0 * p1 * nu
  } else if (approx %in% c("asymptotic", "generalized schoenfeld")) {
    sigma2  <- stats::integrate(function(x) weight_wlr(x, arm0, arm1, weight)^2 *
                                  p0 * prob_risk(arm0, x) * p1 * prob_risk(arm1, x) /
                                  ( p0 * prob_risk(arm0, x) + p1 * prob_risk(arm1, x) )^2 *
                                  ( p0 * dens_event(arm0, x) + p1 * dens_event(arm1, x)),
                                lower=0,
                                upper=arm0$total_time)$value
  } else {
    stop("Please specify a valid approximation for the mean.", call.=F)
  }

  return(sigma2)

}

sigma2j_wlr <- function(j, arm0, arm1, weight="1") {

  p1 <- arm1$size / (arm0$size + arm1$size)
  p0 <- 1 - p1
  if (j==0) {
    arm <- arm0
  } else {
    arm <- arm1
  }

  inner <- function(x) {
    sapply(x, function(x1)
      stats::integrate(function(x2) weight_wlr(x2, arm0, arm1, weight) *
                         (j - p1 * prob_risk(arm1, x2) / ( p0 * prob_risk(arm0, x2) + p1 * prob_risk(arm1, x2) )) *
                         ( p0 * dens_event(arm0, x2) + p1 * dens_event(arm1, x2) ) /
                         ( p0 * prob_risk(arm0, x2) + p1 * prob_risk(arm1, x2) ),
                       lower=0,
                       upper=x1)$value
    )
  }
  sigma2jA <- stats::integrate(function(x) weight_wlr(x, arm0, arm1, weight)^2 *
                                 (j - p1 * prob_risk(arm1, x) / ( p0 * prob_risk(arm0, x) + p1 * prob_risk(arm1, x) ))^2 *
                                 dens_event(arm, x),
                               lower=0,
                               upper=arm0$total_time)$value
  sigma2jB <- -1 * deltaj_wlr(j, arm0, arm1, weight)^2
  sigma2jC <- 2 * stats::integrate(function(x) inner(x) *
                                     weight_wlr(x, arm0, arm1, weight) *
                                     (j - p1 * prob_risk(arm1, x) / ( p0 * prob_risk(arm0, x) + p1 * prob_risk(arm1, x) )) *
                                     (( p0 * dens_event(arm0, x) + p1 * dens_event(arm1, x) ) /
                                        ( p0 * prob_risk(arm0, x) + p1 * prob_risk(arm1, x) ) *
                                        prob_risk(arm, x) -
                                        dens_event(arm, x)),
                                   lower=0,
                                   upper=arm0$total_time)$value
  sigma2jA + sigma2jB + sigma2jC

}

tsigma2_wlr <- function(arm0, arm1, weight="1", approx="block") {

  p1 <- arm1$size / (arm0$size + arm1$size)
  p0 <- 1 - p1
  if (approx %in% c("block", "simple")) {
    tsigma2 <- p0 * sigma2j_wlr(0, arm0, arm1, weight) +
      p1 * sigma2j_wlr(1, arm0, arm1, weight)
    if (approx == "simple") {
      tsigma2 <- p0 *
        p1 *
        (deltaj_wlr(0, arm0, arm1, weight) - deltaj_wlr(1, arm0, arm1, weight))^2 +
        tsigma2
    }
  } else {
    stop("Please specify a valid approximation for the variance.", call.=F)
  }

  return(tsigma2)

}

# Cox log hazard ratio
sigma2_clhr <- function(arm0, arm1) {
  if (sum(arm0$surv_shape!=arm1$surv_shape)>0 |
      length(unique(arm1$surv_scale/arm0$surv_scale))>1) {
    warning("Hazard is not proportional over time.")
  }

  p1 <- arm1$size / (arm0$size + arm1$size)
  1 / stats::integrate(function(x) ( 1 / (1 - p1) / dens_event(arm0, x) +
                                       1/ p1 / dens_event(arm1, x) ) ^ (-1),
                       lower=0,
                       upper=arm0$total_time)$value
}
