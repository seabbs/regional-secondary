library(data.table)
library(purrr)
library(EpiNow2)
library(rstan)
library(ggplot2)

# simulate data according to a convolution model
weight_cmf <- function(x, meanlog, sdlog) {
  cmf <- cumsum(dlnorm(1:length(x), meanlog, sdlog)) -
    cumsum(dlnorm(0:(length(x) - 1), meanlog, sdlog))
  cmf <- cmf / plnorm(length(x), meanlog, sdlog)
  conv <- sum(x * rev(cmf), na.rm = TRUE)
  return(conv)
}

simulate_secondary <- function(data, type = "incidence",
                               obs_model = "poisson", delay_max = 30, ...) {
  type <- match.arg(type, choices = c("incidence", "prevalence"))
  obs_model <- match.arg(obs_model, choices = c("none", "poisson", "negbin"))
  data <- as.data.table(data)
  data <- copy(data)
  data <- data[, index := 1:.N]
  # apply scaling
  data <- data[, scaled := scaling * primary]
  # add convolution
  data <- data[,
    conv := pmap_dbl(list(i = index, m = meanlog, s = sdlog),
     function(i, m, s) {
       weight_cmf(scaled[max(1, i - delay_max):i],
                  meanlog = m, sdlog = s)
     })]
  # build model
  if (type == "incidence") {
    data <- data[, secondary := conv]
  }else if (type == "prevalence") {
    data <- data[1, secondary := scaled]
    for (i in 2:nrow(data)) {
      index <-
        data[c(i - 1, i)][, secondary := shift(secondary, 1) - conv]
      index <- index[secondary < 0, secondary := 0]
      data[i, ] <- index[2][, secondary := secondary + scaled]
    }
  }
  # check secondary is greater that zero
  data <- data[secondary < 0, secondary := 0]
  data <- data[!is.na(secondary)]
  # apply observation model
  if (obs_model == "poisson") {
    data <- data[, secondary := map_dbl(secondary, ~ rpois(1, .))]
  }else if (obs_model == "negbin") {
    data <- data[, secondary := map_dbl(secondary, ~ rnbinom(1, mu = .), ...)]
  }
  data <- data[, secondary := as.integer(secondary)]
  return(data)
}

# summarise simulated scenarios
summarise_scenario <- function(hosp, window = 1, dates) {
  summarised_scenarios <- copy(hosp)[, .(date, scaling, meanlog, sdlog)]
  summarised_scenarios <- melt(summarised_scenarios, id.vars = "date")
  cris <- function(index, window, x) {
    x <- data.table(value = x[max(1, index - window + 1):index], type = "temp")
    cris <-
      EpiNow2::calc_summary_measures(x, CrIs = c(seq(0.1, 0.9, 0.1), 0.95))
    return(cris)
  }
  summarised_scenarios <- summarised_scenarios[,
    .(date = date, summary = map(1:.N, ~ cris(index = ., window, x = value))),
      by = c("variable")]
  summarised_scenarios <- summarised_scenarios[,
    rbindlist(summary), by = c("date", "variable")][, type := NULL]
  if (!missing(dates)) {
    summarised_scenarios <-
    summarised_scenarios[as.character(date) %in% as.character(dates)]
  }

  setnames(summarised_scenarios, "date", "target_date")
  return(summarised_scenarios)
}

# join multiple simulations together
join_simulations <- function(simulations, labels, to_week = FALSE) {
  simulations <- map2(simulations, labels, ~ .x[, target := .y])
  simulations <- rbindlist(simulations, fill = TRUE, use.names = TRUE)

  if (to_week) {
    simulations <- simulations[,
      date := floor_date(date, "week", week_start = 1)]
    simulations <- simulations[, lapply(.SD, sum),
                                by = c("date", "target"),
                                .SDcols = c("secondary"), ]
  }
  simulations <- simulations[, target : as.factor(target)]
  return(simulations)
}


# plot parameter posteriors using output from summarise_parameter_posteriors
plot_trace <- function(draws, data = NULL,
                       samples = 100,
                       alpha = 0.01, obs_alpha = 0.8,
                       scale = "continuous", x_axis = TRUE,
                       scale_label = "scaling") {
  scale <- match.arg(scale, choice = c("continuous", "log", "percent"))
  draws <- as.data.table(draws)[, Source := "Model"]
  draws <- draws[sample <= samples]
  if (!is.null(data)) {
    data <- data[, Source := "Simulation"]
  }

  plot <- ggplot(draws) +
    aes(x = date, y = value, group = sample) +
    geom_line(size = 1.1, alpha = alpha) +
    theme_minimal() +
    labs(x = "Date", y = scale_label) +
    guides(size = NULL) +
    theme_cowplot() +
    scale_x_date(date_breaks = "2 week", date_labels = "%b %d") + 
    theme(axis.text.x = ggplot2::element_text(angle = 90)) + 
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~target, ncol = 1, scales = "free_y")

  if (!is.null(data)) {
    plot <- plot +
      geom_point(aes(group = NULL), alpha = obs_alpha, size = 1.2,
                 data = data, colour =  "black")
  }

  if (scale %in% "log") {
    plot <- plot + 
      scale_y_continuous(labels = comma, trans = log_trans()) +
      labs(x = "Date", y = paste0(scale_label, "(log scale)"))
  }else if (scale %in% "continuous") {
    plot <- plot +
      scale_y_continuous(labels = comma) +
      labs(x = "Date", y = scale_label)
  }else if (scale %in%  "percent") {
        plot <- plot +
      scale_y_continuous(labels = percent) +
      labs(x = "Date", y = scale_label)
  }

  if (!x_axis) {
    plot <- plot +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
  }

  return(plot)
}

# add the difference between variables
add_diff_variable <- function(dt, variable, label, by, fill = 0,
                              shift_col = "value",
                              exact_cols = c("median", "mean", "secondary",
                                         "value"),
                              partial_cols = c("lower_", "upper_")) {
  dt <- copy(dt)
  dt_alt <- dt[target %in% variable]

  cols <- colnames(dt_alt)
  target_cols <- intersect(cols, exact_cols)
  target_cols <- c(target_cols, grep(
    paste(partial_cols, collapse = "|"),
    cols, fixed = FALSE, value = TRUE
  ))

  if (missing(by)) {
    by <- "across"
    dt_alt[, across := 1]
  }

  dt_alt <- dt_alt[,
    (target_cols) := map(.SD, ~ . - shift(get(shift_col), fill = fill)),
    .SDcols = target_cols, by = by]
  dt_alt <- dt_alt[, target := label]
  dt_alt <- suppressWarnings(dt_alt[, across := NULL])
  dt <- rbind(dt, dt_alt)
  return(dt)
}


# plot posterior predictions
plot_predictions <- function (simulations, predictions,
                              variable, log = TRUE, keep_x = TRUE) {
  simulations <- copy(simulations)
  predictions <- copy(predictions)

  if (!missing(variable)) {
    simulations <- simulations[target %in% variable]
    predictions <- predictions[target %in% variable]
  }

  plot <- ggplot(simulations) +
    aes(x = date, y = secondary) +
    geom_ribbon(data = predictions,
                aes(ymin = lower_30, ymax = upper_30, group = target_date),
                alpha = 0.15, fill = "#1B9E77") +
    geom_ribbon(data = predictions,
                aes(ymin = lower_60, ymax = upper_60, group = target_date),
                alpha = 0.15, fill = "#1B9E77") +
    geom_ribbon(data = predictions,
                aes(ymin = lower_90, ymax = upper_90, group = target_date),
                alpha = 0.15, fill = "#1B9E77") +
    geom_point(alpha = 0.2, size = 0.8) +
    theme_cowplot() +
    scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position = "none") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")

  if (log) {
    plot <- plot + 
      scale_y_continuous(labels = comma, trans = log_trans()) +
      labs(x = "Date", y = "Notifications (log scale)")
  }else{
    plot <- plot +
      scale_y_continuous(labels = comma) +
      labs(x = "Date", y = "Notifications")
  }

  if (!keep_x) {
    plot <- plot +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
  }
  return(plot)
}
