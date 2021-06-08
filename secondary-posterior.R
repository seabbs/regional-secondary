# function for fitting and extracting posteriors
secondary_posterior <- function(obs, target_date, window = 14, prior,
                                prior_scale = 1, args = NULL) {
  if (missing(prior)) {
    prior <- NULL
  }

  snapshot <- copy(obs)[date <= target_date]
  burn_in <- as.integer(nrow(snapshot) - window)

  if (!is.null(prior)) {
     prior <- prior[, sd := sd * prior_scale]
     args <- update_secondary_args(args, posterior = prior)
  }

  model <- do.call(estimate_secondary, c(
    list(reports = obs, burn_in = burn_in),
    args
  ))

  posterior <- extract_stan_param(
    model$fit,
    CrIs = c(seq(0.1, 0.9, 0.1), 0.95)
  )
  return(posterior)
}