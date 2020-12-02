# load required packages
library(EpiNow2)
library(future.apply)
library(purrr)
library(data.table)

# inner function for forecasing a single region
forecast_region <- function(target_region, reports, case_forecast, verbose = TRUE, ...) {
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
  out$plots$fit <- plot(cases_to_deaths)
  
  deaths_forecast <- forecast_secondary(cases_to_deaths, pred_cases, samples = 1000)
  out$plots$forecast <- plot(deaths_forecast, from = max(target_obs$date) - 7)
  
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
  out$estimate_secondary <- cases_to_deaths
  if (verbose) {
    message("Completed forecast for: ", target_region)
  }
  return(out)
}

# wrapper for forecasting across regions
# additional arguments are passed to estimate_secondary
regional_secondary <- function(reports, case_forecast, verbose = interactive(), ...) {
  # run the forecast safely in case of failure
  safe_forecast_region <- safely(forecast_region)
  
  # forecast all regions
  forecasts <- future_lapply(unique(reports$region), 
                             safe_forecast_region,
                             reports = reports, 
                             case_forecast = case_forecast,
                             verbose = verbose,
                             future.seed = TRUE, 
                             future.scheduling = Inf,
                             ...)
  # pick out results and name
  forecasts <- purrr::map(forecasts, ~ .[[1]])
  names(forecasts) <- unique(reports$region)
  
  # format output
  out <- list()
  out$region <- forecasts
  out$samples <- rbindlist(map(forecasts, ~ .$samples), idcol = "region")
  out$summarised <- rbindlist(map(forecasts, ~ .$summarised), idcol = "region")
  return(out)
}