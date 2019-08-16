# Indicates whether a person (row) was at risk at the beginning of
# the survival interval (column)
interval_atrisk <- function(x, surv_interval = c(0, Inf)) {

  k <- length(x)
  t(matrix(rep(surv_interval, k), ncol=k)) <= x

}

# Observed person-time for each person (row) across each
# survival interval (column)
interval_ptime <- function(x, surv_interval = c(0, Inf)) {

  atrisk <- interval_atrisk(x, surv_interval)
  t(diff(
    t(atrisk) * c(utils::head(surv_interval, -1), 0) +
      t((1-atrisk) * x)
  ))

}
