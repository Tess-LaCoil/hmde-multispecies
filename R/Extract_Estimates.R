extract_fit_estimates <- function(sp_codes,
                                  fit_file_path,
                                  rstan_file_path){
  for(i in 1:length(sp_codes)){
    print(paste0("Extracting: ", sp_codes[i]))
    fit <- readRDS(paste0(fit_file_path,
                          sp_codes[i],
                          "_SRSWR_fit.rds"))
    rstan_data <- readRDS(paste0(rstan_file_path,
                                 sp_codes[i],
                                 "_SRSWR_stan_data.rds"))
    obs_data <- tibble(
      ind_id = rstan_data$ind_id,
      time = rstan_data$time,
      y_obs = rstan_data$time,
      obs_index = rstan_data$obs_index
    )

    fit@model_name <- "canham_multi_ind"

    ests <- hmde_extract_estimates(fit = fit,
                                   input_measurement_data = obs_data)

    saveRDS(ests, paste0("output/data/", sp_codes[i], "_ests.rds"))
  }
}
