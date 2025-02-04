fit_hierarchical_models <- function(sp_codes,
                                    rstan_file_path){
  for(i in 1:length(sp_codes)){
    data_filename <- paste0(rstan_file_path, sp_code, "_SRSWR_stan_data.rds")
    rstan_data <- readRDS(data_filename)

    fit <- hmde_model("canham_multi_ind") %>%
      hmde_assign_data(n_obs = rstan_data$n_obs,
                       n_ind = rstan_data$n_ind,
                       y_obs = rstan_data$y_obs,
                       obs_index = rstan_data$obs_index,
                       time = rstan_data$time,
                       ind_id = rstan_data$ind_id) %>%
      hmde_run(iter=3000,
               warmup = 1500,
               chains = 4,
               cores = 4,
               control = list(max_treedepth = 14))

    saveRDS(fit, file = paste0("output/data/", sp_codes, "_SRSWR_fit.rds"))
  }
}

fit_species_level_models <- function(sp_codes, rstan_file_path){
  for(i in 1:length(sp_codes)){
    data_filename <- paste0(rstan_file_path, sp_codes[i], "_SRSWR_stan_data.rds")
    input_data <- readRDS(data_filename)

    model_data_tibble <-
      tibble(
        ind_id = input_data$ind_id,
        time = input_data$time,
        y_obs_init = input_data$y_obs
      ) %>%
      group_by(ind_id) %>%
      mutate(
        interval = lead(time) - time,
        delta_obs = lead(y_obs_init) - y_obs_init,
        y_obs = (lead(y_obs_init) + y_obs_init)/2
      ) %>%
      ungroup() %>%
      na.omit()

    rstan_data <- list(
      n_obs = nrow(model_data_tibble),
      y_obs = model_data_tibble$y_obs,
      delta_obs = model_data_tibble$delta_obs,
      interval = model_data_tibble$interval
    )

    model <- stan_model(file="stan/canham_species_level.stan")

    print(paste(sp_codes[i], " model fit start at:", Sys.time()))
    fit <- sampling(model,
                    data=rstan_data,
                    iter=3000,
                    warmup = 1500,
                    chains = 4,
                    cores = 4,
                    control = list(max_treedepth = 14))

    print(paste("Model fit end at:", Sys.time()))

    save_filename <- paste0("output/", sp_codes[i], "_SRSWR_species_level_fit.rds")
    saveRDS(fit, file=save_filename)
  }
}
