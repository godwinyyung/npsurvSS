# Calculate design parameters delta, sigma2, and tsigma2
calc_design <- function(arm0, arm1, test) {

  p1  <- arm1$size / (arm0$size + arm1$size)
  p0  <- 1 - p1
  tsigma2 <- NULL

  # survival difference and ratio
  if (grepl("survival", test$test)) {

    if (! "milestone" %in% names(test)) {
      stop(paste("Please provide milestone for ",
                 test$test,
                 ".",
                 sep=""),
           call.=F)
    }
    milestone <- test$milestone
    surv0 <- psurv(milestone, arm0, lower.tail=F)
    surv1 <- psurv(milestone, arm1, lower.tail=F)

    if (test$test == "survival difference") { # survival difference
      delta   <- surv0 - surv1
      sigma2  <- sigma2j_surv(arm0, milestone) / p0 + sigma2j_surv(arm1, milestone) / p1
    } else if (test$test == "survival ratio") { # survival ratio
      delta   <- log(surv0) - log(surv1)
      sigma2  <- sigma2j_cumh(arm0, milestone) / p0 + sigma2j_cumh(arm1, milestone) / p1
    } else {
      stop("Please specify valid survival contrast.", call.=F)
    }

  } else if (grepl("rmst", test$test)) { # rmst difference and ratio

    if (! "milestone" %in% names(test)) {
      stop(paste("Please provide milestone for ",
                 test$test,
                 ".",
                 sep=""),
           call.=F)
    }
    milestone <- test$milestone
    rmst0 <- deltaj_rmst(milestone, arm0)
    rmst1 <- deltaj_rmst(milestone, arm1)

    if (test$test == "rmst difference") { # rmst difference
      delta   <- rmst0 - rmst1
      sigma2  <- sigma2j_rmst(arm0, milestone) / p0 + sigma2j_rmst(arm1, milestone) / p1
    } else if (test$test == "rmst ratio") { # rmst ratio
      delta   <- log(rmst0) - log(rmst1)
      sigma2  <- sigma2j_rmst(arm0, milestone) / p0 / rmst0^2 +
        sigma2j_rmst(arm1, milestone) / p1 / rmst1^2
    } else {
      stop("Please specify a valid rmst contrast.", call.=F)
    }

  } else if (grepl("percentile", test$test)) { # percentile difference or ratio

    if (! "percentile" %in% names(test)) {
      stop(paste("Please provide percentile for ",
                 test$test,
                 ".",
                 sep=""),
           call.=F)
    }
    percentile <- test$percentile
    perc0 <- qsurv(percentile, arm0)
    perc1 <- qsurv(percentile, arm1)

    if (test$test == "percentile difference") { # percentile difference
      delta   <- perc0 - perc1
      sigma2  <- sigma2j_perc(arm0, percentile) / p0 + sigma2j_perc(arm1, percentile) / p1
    } else if (test$test == "percentile ratio") { # percentile ratio
      delta   <- log(perc0) - log(perc1)
      sigma2  <- sigma2j_perc(arm0, percentile) / p0 / perc0^2 +
        sigma2j_perc(arm1, percentile) / p1 / perc1^2
    } else {
      stop("Please specify a valid percentile contrast.", call.=F)
    }
  } else if (test$test == "hazard ratio") {
    delta     <- c(arm0$surv_shape * log( arm1$surv_scale / arm0$surv_scale ))[1] # log hazard-ratio
    sigma2    <- sigma2_clhr(arm0, arm1)
  } else if (test$test == "weighted logrank") {
    # weight
    if (! "weight" %in% names(test)) {
      test$weight <- "1"
    }
    # delta and sigma2
    if (! "mean.approx" %in% names(test)) {
      test$mean.approx <- "asymptotic"
    }
    delta     <- delta_wlr(arm0, arm1, test$weight, test$mean.approx)
    sigma2    <- sigma2_wlr(arm0, arm1, test$weight, test$mean.approx)
    # tsigma2
    if (! "var.approx" %in% names(test)) {
      test$var.approx <- "1"
    }
    if (test$var.approx == "1") {
      tsigma2 <- sigma2
    } else {
      tsigma2 <- tsigma2_wlr(arm0, arm1, test$weight, test$var.approx)
    }
  } else {
    stop("Please specify a valid test.", call.=F)
  }
  if(is.null(tsigma2)) {
    tsigma2 <- sigma2
  }

  return(list(delta=delta, sigma2=sigma2, tsigma2=tsigma2))

}

#' Power
#'
#' Calculate power for a two-arm survival study.
#' @param arm0  object of class 'arm'.
#' @param arm1  object of class 'arm'.
#' @param test  list or list of lists. Each list must contain at minimum
#'   the key 'test' describing the type of statistical test. Default test
#'   is the "weighted logrank". Kaplan-Meier based tests ("survival difference",
#'   "survival ratio", "rmst difference", "rmst ratio", "percentile difference",
#'   and "percentile ratio") require the user to define an additional key,
#'   either the desired 'milestone' or 'percentile'. The weighted log-rank test
#'   does not require additional keys. However, user may choose which weight function
#'   ("1"=unweighted, "n"=Gehan-Breslow, "sqrtN"=Tarone-Ware, "FH_p[a]_q[b]"=
#'   Fleming-Harrington with p=a and q=b) and which approximation for the
#'   large-sample mean ("asymptotic", "generalized schoenfeld", "event driven",
#'   "freedman", "rubinstein") and variance ("1", "block[ randomization]", "simple[ randomization]") 
#'   they wish to use. Default choice is 'weight'="1", 'mean.approx'="asymptotic", and 'var.approx'="1".
#'   For more details regarding the different mean and variance approximations
#'   for the weight log-rank test, please see Yung and Liu (2020). If there are multiple lists, 
#'   then users may provide a 'label' for each list to be displayed in the output.
#' @param alpha type 1 error rate
#' @param sides 1=1-sided test, 2=2-sided test
#' @return power.
#' @seealso \code{\link{create_arm}} for creating an object of class 'arm'.
#' @references Yung, G and Liu, Y. (2020). Sample size and power for the weighted
#' log-rank test and Kaplan-Meier based tests with allowance for non-proportional
#' hazards. \emph{Biometrics} 76(3):939-950.
#' @examples
#' arm0 <- create_arm(size=120, accr_time=6, surv_scale=0.05, loss_scale=0.005, follow_time=12)
#' arm1 <- create_arm(size=120, accr_time=6, surv_scale=0.03, loss_scale=0.005, follow_time=12)
#' power_two_arm(arm0, arm1)
#' power_two_arm(arm0, arm1, list(test="weighted logrank",
#'   weight="n",
#'   mean.approx="generalized schoenfeld",
#'   var.approx="block"))
#' power_two_arm(arm0, arm1, list(test="survival difference", milestone=12))
#' power_two_arm(arm0, arm1, list(test="rmst ratio", milestone=12))
#' power_two_arm(arm0, arm1, list(test="percentile difference", percentile=0.25))
#' power_two_arm(arm0, arm1, list(
#'   list(test="weighted logrank", label="Logrank"),
#'   list(test="survival difference", milestone=12, label="12-month survival difference")))
#' @export
power_two_arm <- function(arm0,
                          arm1,
                          test=list(test="weighted logrank"),
                          alpha=0.025,
                          sides=1) {

  if (! inherits(test[[1]], "list")) { # one test to perform

    n       <- arm0$size + arm1$size
    design  <- calc_design(arm0, arm1, test)
    out     <- stats::pnorm((sqrt(design$sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(n) * design$delta) /
                              sqrt(design$tsigma2),
                            lower.tail=F) +
      (sides==2) * stats::pnorm((sqrt(design$sigma2) * stats::qnorm(1 - alpha / sides) - sqrt(n) * design$delta) /
                           sqrt(design$tsigma2),
                         lower.tail=F)
    return(out)

  } else { # multiple tests to perform

    out <- c()
    for (i in 1:length(test)) {
      label = ifelse("label" %in% names(test[[i]]), test[[i]]$label, i)
      out <- rbind(out, c(label, power_two_arm(arm0, arm1, test[[i]], alpha, sides)))
    }
    out <- data.frame(out)
    names(out) <- c("test", "power")
    return(out)

  }

}

#' Sample size
#'
#' Calculate required sample size and expected number of events for a
#' two-arm survival study.
#' @param arm0  object of class 'arm'.
#' @param arm1  object of class 'arm'.
#' @param test  list or list of lists. Each list must contain at minimum
#'   the key 'test' describing the type of statistical test. Default test
#'   is the "weighted logrank". Kaplan-Meier based tests ("survival difference",
#'   "survival ratio", "rmst difference", "rmst ratio", "percentile difference",
#'   and "percentile ratio") require the user to define an additional key,
#'   either the desired 'milestone' or 'percentile'. The weighted log-rank test
#'   does not require additional keys. However, user may choose which weight function
#'   ("1"=unweighted, "n"=Gehan-Breslow, "sqrtN"=Tarone-Ware, "FH_[a]_[b]"=
#'   Fleming-Harrington with p=a and q=b) and which approximation for the
#'   large-sample mean ("asymptotic", "generalized schoenfeld", "event driven",
#'   "freedman", "rubinstein") and variance ("1", "block[ randomization]", "simple[ randomization]") 
#'   they wish to use. Default choice is 'weight'="1", 'mean.approx'="asymptotic", and 'var.approx'="1".
#'   For more details regarding the different mean and variance approximations
#'   for the weight log-rank test, please see Yung and Liu (2020). If there are multiple lists, 
#'   then users may provide a 'label' for each list to be displayed in the output.
#' @param power 1 - type 2 error rate
#' @param alpha type 1 error rate
#' @param sides 1=1-sided test, 2=2-sided test
#' @return
#'   \item{n0}{sample size for \code{arm0}}
#'   \item{n1}{sample size for \code{arm1}}
#'   \item{n}{total sample size}
#'   \item{d0}{expected number of events for \code{arm0}}
#'   \item{d1}{expected number of events for \code{arm1}}
#'   \item{d}{total expected number of events; can be used to convert a time-driven
#'   trial to an event-driven trial.}
#' @seealso \code{\link{create_arm}} for creating an object of class 'arm'.
#' @references Yung, G and Liu, Y. (2020). Sample size and power for the weighted
#' log-rank test and Kaplan-Meier based tests with allowance for non-proportional
#' hazards. \emph{Biometrics} 76(3):939-950.
#' @examples
#' arm0 <- create_arm(size=120, accr_time=6, surv_scale=0.05, loss_scale=0.005, follow_time=12)
#' arm1 <- create_arm(size=120, accr_time=6, surv_scale=0.03, loss_scale=0.005, follow_time=12)
#' size_two_arm(arm0, arm1)
#' size_two_arm(arm0, arm1, list(test="weighted logrank",
#'   weight="n",
#'   mean.approx="generalized schoenfeld",
#'   var.approx="block"))
#' size_two_arm(arm0, arm1, list(test="survival difference", milestone=12))
#' size_two_arm(arm0, arm1, list(test="rmst ratio", milestone=12))
#' size_two_arm(arm0, arm1, list(test="percentile difference", percentile=0.25))
#' size_two_arm(arm0, arm1, list(
#'   list(test="weighted logrank", label="Logrank"),
#'   list(test="survival difference", milestone=12, label="12-month survival difference")))

#' @export
size_two_arm <- function(arm0,
                         arm1,
                         test=list(test="weighted logrank"),
                         power=0.8,
                         alpha=0.025,
                         sides=1) {

  if (! inherits(test[[1]], "list")) { # one test to perform

    # sample size for 1-sided test
    p1      <- arm1$size / (arm0$size + arm1$size)
    p0      <- 1 - p1
    design  <- calc_design(arm0, arm1, test)
    out     <- ( sqrt(design$sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(design$tsigma2) * stats::qnorm(power) )^2 /
      design$delta^2 *
      c(p0, p1)
    # out     <- ceiling(out) # rounding

    # refine sample size for 2-sided test
    if (sides==2) {
      i         <- 1
      arm0$size <- out[1]
      arm1$size <- out[2]
      cont      <- T
      while (cont) {
        temp <- power_two_arm(arm0, arm1, test, alpha, sides)
        if (temp > power) {
          i         <- i + 1
          arm0$size <- out[1] - i * p0
          arm1$size <- out[2] - i * p1
        } else {
          out   <- out - (i-1) * c(p0, p1)
          cont  <- F
        }
      }
    }

    out <- c(out, # n0, n1
             sum(out), # n
             out[1] * prob_event(arm0), # d0
             out[2]  *prob_event(arm1), # d1
             out[1] * prob_event(arm0) + out[2] * prob_event(arm1)) # d
    # if (out[4] %% 1 < out[5] %% 1) { # rounding
    #   out[4] <- floor(out[4])
    #   out[5] <- ceiling(out[5])
    # } else {
    #   out[4] <- ceiling(out[4])
    #   out[5] <- floor(out[5])
    # }
    names(out) <- c("n0", "n1", "n", "d0", "d1", "d")
    return(out)

  } else { # multiple tests to perform

    out <- c()
    for (i in 1:length(test)) {
      label = ifelse("label" %in% names(test[[i]]), test[[i]]$label, i)
      out <- rbind(out, c(label, size_two_arm(arm0, arm1, test[[i]], power, alpha, sides)))
    }
    out <- data.frame(out)
    names(out)[1] <- "test"
    return(out)

  }

}
