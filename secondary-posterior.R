# function for fitting and extracting posteriors
secondary_posterior <- function(obs, target_date, start_date, window = 14,
                                prior, prior_scale = 1, args = NULL) {
  if (missing(prior)) {
    prior <- NULL
  }
  if (missing(start_date)) {
    start_date <- NULL
  }

  snapshot <- copy(obs)[date <= target_date]
  if (!is.null(start_date)) {
    snapshot <- snapshot[date >= start_date]
  }

  burn_in <- max(14, as.integer(nrow(snapshot) - window))

  if (!is.null(prior)) {
     prior <- prior[, sd := sd * prior_scale]
     args <- update_secondary_args(args, posterior = prior)
  }

  model <- do.call(estimate_secondary, c(
    list(reports = snapshot, burn_in = burn_in),
    args
  ))

  posterior <- extract_stan_param(
    model$fit,
    CrIs = c(seq(0.1, 0.9, 0.1), 0.95)
  )

  out <- model
  out$posterior <- posterior
  return(out)
}

secondary_snapshot <- function(obs, dates, priors, iterate_prior = FALSE,
                               obs_weeks = 8, obs_fit = 2, ...) {
  fits <- vector("list", length(dates))
  names(fits) <- dates
  for (i in seq_along(dates)) {
    # set target, start date for fitting, and window size
    target_date <- dates[i]
    window <- obs_fit * 7
    start_date <- target_date - weeks(obs_weeks)

    # set up prior
    prior <- NULL
    if (!missing(priors)) {
      prior <- priors[[as.character(target_date)]]$posterior
    }
    if (iterate_prior) {
      prior <- fits[[as.character(target_date - weeks(1))]]$posterior
    }
    # restrict to current data of interest
    message("Estimating for the: ", as.character(target_date))
    fits[[as.character(target_date)]] <- secondary_posterior(
        obs,
        target_date,
        window = obs_fit * 7,
        start_date = start_date,
        prior = prior,
        ...
    )
  }
  return(fits)
}