run_analysis <- function(sp_codes,
                         col_vec,
                         est_file_path){
  #Construct tibbles containing all species data
  full_data <- construct_full_data_tibbles(sp_codes, est_file_path)
  saveRDS(full_data, file = "output/data/full_est_data.rds")

  #Analysis of single ind
  individual_growth_analysis(ind_data = full_data$ind_data_full)

  #Plots plots plots
  fig_1_plots(full_data, sp_code = "gar2in")
  plot_3d_scatter(ind_data = full_data$ind_data_full, col_vec)

}


# Construct data tibbles for measurement, individual, and population-level for all species
construct_full_data_tibbles <- function(sp_codes,
                                        est_file_path){
  measurement_data_full <- tibble(
    ind_id = c(),
    BCI_ind_id = c(),
    species = c(),
    sp_code = c(),
    y_obs = c(),
    y_hat = c(),
    time = c()
  )

  ind_data_full <- tibble(
    ind_id = c(),
    BCI_ind_id = c(),
    ind_max_growth = c(),
    ind_size_at_max_growth = c(),
    ind_k = c(),
    S_initial = c(),
    S_final = c(),
    sp_code = c(),
    species = c()
  )

  sp_data_full <- tibble(
    sp_code=c(),
    species=c(),
    par_name = c(),
    mean=c(),
    median=c(),
    CI_lower=c(),
    CI_upper=c()
  )

  for(i in 1:length(sp_codes)){
    ests <- readRDS(paste0(est_file_path, sp_codes[i],"_ests.rds"))
    sample_data <- readRDS(paste0("input/",sp_codes[i],"_SRSWR_sample.rds"))

    ind_id_tibble <- sample_data %>%
      select(treeid, treeid_factor) %>%
      distinct() %>%
      arrange(treeid_factor)

    #Measurement data
    ests$measurement_data <- ests$measurement_data %>%
      group_by(ind_id) %>%
      arrange(obs_index) %>%
      mutate(BCI_ind_id = ind_id_tibble$treeid[ind_id],
             sp_code = sp_codes[i],
             species = sample_data$species[1],
             S_initial = first(y_hat),
             S_final = last(y_hat)) %>%
      arrange(ind_id, obs_index) %>%
      ungroup()

    measurement_data_temp <- tibble(
      ind_id = ests$measurement_data$ind_id,
      BCI_ind_id = ests$measurement_data$BCI_ind_id,
      species = ests$measurement_data$species,
      sp_code =  ests$measurement_data$sp_code,
      y_obs = ests$measurement_data$y_obs,
      y_hat = ests$measurement_data$y_hat,
      time = ests$measurement_data$time
    )

    measurement_data_full <- rbind(measurement_data_full,
                                   measurement_data_temp)

    #Individual-level data
    ind_initial_final <- ests$measurement_data %>%
      select(ind_id, S_initial, S_final, sp_code, species) %>%
      mutate(ind_id = as.integer(ind_id)) %>%
      distinct()

    ests$individual_data <- ests$individual_data %>%
      group_by(ind_id) %>%
      mutate(BCI_ind_id = ind_id_tibble$treeid[ind_id],
             sp_code = sp_codes[i],
             species = sample_data$species[1]) %>%
      ungroup()

    #Median rather than mean due to sample skewness
    ind_data_temp <- tibble(
      ind_id = ests$individual_data$ind_id,
      BCI_ind_id = ests$individual_data$BCI_ind_id,
      ind_max_growth = ests$individual_data$ind_max_growth_median,
      ind_size_at_max_growth = ests$individual_data$ind_size_at_max_growth_median,
      ind_k = ests$individual_data$ind_k_median,
        ) %>%
      left_join(ind_initial_final, by = "ind_id")

    ind_data_full <- rbind(ind_data_full, ind_data_temp)

    #Population-level data
    sp_data_temp <- tibble(
      sp_code = sp_codes[i],
      species = sample_data$species[1],
      par_name = ests$population_data$par_name,
      mean = ests$population_data$mean,
      median = ests$population_data$median,
      CI_lower = ests$population_data$CI_lower,
      CI_upper = ests$population_data$CI_upper
    )

    sp_data_full <- rbind(sp_data_full, sp_data_temp)
  }

  return_data <- list(
    measurement_data_full = measurement_data_full,
    ind_data_full = ind_data_full,
    sp_data_full = sp_data_full
  )
}

#Analysis of single ind
individual_growth_analysis <- function(ind_data){

}

fig_1_plots <- function(sp_code)

#3D scatter plot of individuals coloured by species
plot_3d_scatter(ind_data, col_vec){
  # 3D scatter plot of parameters
  fig <- plot_ly(ind_data,
                 x = ~ind_max_growth,
                 y = ~ind_size_at_max_growth,
                 z = ~ind_k,
                 color = ~species,
                 colors = col_vec,
                 size = 0.1,
                 opacity = 0.5) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'g_max',
                                     type = "log"),
                        yaxis = list(title = 'S_max',
                                     type = "log"),
                        zaxis = list(title = 'k',
                                     type = "log")))

  htmlwidgets::saveWidget(
    widget = fig, #the plotly object
    file = "output/figures/3dParScatter.html", #the path & file name
    selfcontained = TRUE #creates a single html file
  )
}

