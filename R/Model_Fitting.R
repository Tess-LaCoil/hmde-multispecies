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
               warmup = 1000,
               chains = 4,
               cores = 4,
               control = list(max_treedepth = 14))
  }
  saveRDS(fit, file = paste0("output/data/", sp_codes, "_SRSWR_fit.rds"))
}



   data_filename <-
   rstan_data <-

   model <- stan_model(file="stan/canham_multi_ind.stan")

   print(paste("Model fit start at:", Sys.time()))
   fit <- sampling(model,
                   data=rstan_data,
                   iter=3000,
                   warmup = 1000,
                   chains = 4,
                   cores = 4,
                   control = list(max_treedepth = 14))

   print(paste("Model fit end at:", Sys.time()))

   save_filename <- paste0("output/", sp_code, "_SRSWR_fit.rds")
   saveRDS(fit, file=save_filename)
 }
