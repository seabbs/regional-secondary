# load required packages
library(EpiNow2)
library(future.apply)
library(purrr)
library(data.table)

# inner function for forecasting a single region
forecast_region <- function(target_region, reports, case_forecast, verbose = TRUE, 
                            return_fit = TRUE, return_plots = TRUE, ...) {
  if (verbose) {
    message("Forecasting for: ", target_region)
  }
  # filter for target region
  target_obs <- reports[region == target_region][, region := NULL]
  pred_cases <- case_forecast[region == target_region]
  pred_cases <- pred_cases[, .(date, sample, value = cases)]
  pred_cases <- pred_cases[date > max(target_obs$date)]
  
  # estimate relationship fitting to just the last month of data
  cases_to_deaths <- estimate_secondary(target_obs, verbose = FALSE, ...)
  out <- list()
  if (return_plots) {
      out$plots$fit <- plot(cases_to_deaths)
  }

  deaths_forecast <- forecast_secondary(cases_to_deaths, pred_cases, samples = 1000)
  if (return_plots) {
    out$plots$forecast <- plot(deaths_forecast, from = max(target_obs$date) - 7)
  }
  # link in previous observations to forecast
  obs_as_samples <- target_obs[, .(date, value = secondary, sample = list(unique(deaths_forecast$samples$sample)))]
  obs_as_samples <- obs_as_samples[, .(sample = as.numeric(unlist(sample))), by = c("date", "value")]
  deaths_forecast$samples <- rbindlist(list(
    obs_as_samples,
    deaths_forecast$samples
  ), use.names = TRUE)
  
  # return samples, summary + estimated fit
  out$samples <- deaths_forecast$samples
  out$summarised <- deaths_forecast$predictions
  out$summarised_posterior <- extract_stan_param(cases_to_deaths$fit, CrIs = c(seq(0.1, 0.9, 0.1), 0.95))
  if (return_fit) {
    out$estimate_secondary <- cases_to_deaths
  }
  if (verbose) {
    message("Completed forecast for: ", target_region)
  }
  return(out)
}

# wrapper for forecasting across regions
# additional arguments are passed to estimate_secondary
regional_secondary <- function(reports, case_forecast, verbose = interactive(), 
                               return_fit = TRUE, return_plots = TRUE,  ...) {
  
  # Convert to data.table
  reports <- as.data.table(reports)
  case_forecast <- as.data.table(case_forecast)
  
  # run the forecast safely in case of failure
  safe_forecast_region <- safely(forecast_region)
  
  # forecast all regions
  forecasts <- future_lapply(unique(reports$region), 
                             safe_forecast_region,
                             reports = reports, 
                             case_forecast = case_forecast,
                             verbose = verbose,
                             return_fit = return_fit,
                             return_plots = return_plots,
                             future.seed = TRUE, 
                             future.scheduling = Inf,
                             ...)
  # pick out error messages
  errors <- map(forecasts, ~ .[[2]])
  names(errors) <- unique(reports$region)
  # pick out results and name
  forecasts <- map(forecasts, ~ .[[1]])
  names(forecasts) <- unique(reports$region)
  
  # format output
  out <- list()
  out$region <- forecasts
  out$samples <- rbindlist(map(forecasts, ~ .$samples), idcol = "region")
  out$summarised <- rbindlist(map(forecasts, ~ .$summarised), idcol = "region")
  out$errors <- errors
  return(out)
}