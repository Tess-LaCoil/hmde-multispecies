diagnostic_list_single_species <- list(
  #Global
  plot1=list(
    pars=c("global_error_sigma"),
    name="canham_sigma"
  ),

  plot2 <- list(
    pars = "ind_max_growth",
    name = "ind_max_growth"
  ),
  plot3 <- list(
    pars = "ind_size_at_max_growth",
    name = "ind_size_at_max_growth"
  ),
  plot4 <- list(
    pars = "ind_k",
    name = "ind_k"
  ),

  #Species
  plot5=list(
    pars=c("pop_max_growth_mean",
           "pop_max_growth_sd",
           "pop_size_at_max_growth_mean",
           "pop_size_at_max_growth_sd",
           "pop_k_mean",
           "pop_k_sd"),
    name="canham_pars"
  )
)

extract_fit_estimates <- function(sp_codes,
                                  fit_file_path,
                                  rstan_file_path,
                                  plot_diagnostics,
                                  warmup = 1500){
  for(i in 1:length(sp_codes)){
    print(paste0("Extracting: ", sp_codes[i]))
    fit <- readRDS(paste0(fit_file_path,
                          sp_codes[i],
                          "_SRSWR_fit.rds"))
    fit@sim$warmup <- warmup
    fit@model_name <- "canham_multi_ind"

    rstan_data <- readRDS(paste0(rstan_file_path,
                                 sp_codes[i],
                                 "_SRSWR_stan_data.rds"))
    obs_data <- tibble(
      ind_id = rstan_data$ind_id,
      time = rstan_data$time,
      y_obs = rstan_data$time,
      obs_index = rstan_data$obs_index
    )

    if(plot_diagnostics){
      plot_all_diagnostics(fit, sp_code = sp_codes[i],
                           plot_list = diagnostic_list_single_species)
    }

    ests <- hmde_extract_estimates(fit = fit,
                                   input_measurement_data = obs_data)

    saveRDS(ests, paste0("output/data/", sp_codes[i], "_ests.rds"))
  }
}

plot_all_diagnostics <- function(fit, sp_code, plot_list, inc_warmup = FALSE){
  for(i in 1:length(plot_list)){
    plot_diagnostic_trace(fit,
                          pars=plot_list[[i]]$pars,
                          name=plot_list[[i]]$name,
                          inc_warmup,
                          sp_code)
  }

  plot_diagnostic_rhat_hist(fit, sp_code, name = "RHat_hist")
}

# Build and export histogram of RHats
plot_diagnostic_rhat_hist <- function(fit, sp_code, name){
  #Histogram of R_hat
  fit_summary <- summary(fit)
  Rhat <- fit_summary$summary[,10]
  n_NaN <- length(which(is.na(Rhat)))
  n <- length(Rhat)

  filename <- paste0("output/figures/diagnostic/", paste(sp_code,  name, sep="_"), "_RHat_Hist.png")
  png(filename, width=500, height=300)
  hist(Rhat, main=paste(n_NaN,"NaN values from", n, "Rhats"))
  dev.off()
}

#Save diagnostic plot to file
plot_diagnostic_trace <- function(fit, pars, name, inc_warmup, sp_code){
  if(grepl("ind", name)){
    width <- 4000
    height <- 4000
  } else {
    width <- 800
    height <- 800
  }

  filename <- paste0("output/figures/diagnostic/", paste(sp_code,  name, sep="_"), ".png")
  png(filename, width=width, height=height)
  print(traceplot(fit, pars=pars, inc_warmup=inc_warmup))
  dev.off()
}
