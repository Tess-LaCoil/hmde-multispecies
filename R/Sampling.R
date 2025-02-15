sample_data <- function(sample_size,
                        data_path,
                        out_path,
                        sp_codes,
                        col_vec,
                        plot_size_histograms = TRUE){
  data_trees <- read.csv(data_path)

  data_trees_transformed <-
    data_trees %>%
    arrange(treeid, julian_date) %>%
    group_by(treeid) %>%
    mutate(
      #Dates
      census_interval = (julian_date - lag(julian_date)) / 365.25,
      census = rank(julian_date),
      zeroed_time = (julian_date-min(julian_date)) / 365.25,

      #Individual flags
      max_census = max(census),
      same_stem = max(stemid) - min(stemid)
    ) %>%
    ungroup()

  data_trees_transformed$census_interval[is.na(data_trees_transformed$census_interval)] <- 0 #Replace NAs with 0s for first observation

  data <- data_trees_transformed %>%
    filter(max_census == 6) %>% #Trees that has been around the entire time
    filter(same_stem == 0) %>% #Removes records that may have different stems
    filter(sp %in% sp_codes)


  sample_size <- 300
  for(i in 1:length(sp_codes)){
    species_data <- data %>%
      filter(sp == sp_codes[i])

    #If there are enough individuals take a sample otherwise take the lot.
    if(length(unique(species_data$treeid)) > sample_size){
      sample_IDs <- sample(unique(species_data$treeid),
                           size = sample_size)

      sample_data <- species_data %>%
        filter(treeid %in% sample_IDs) %>%
        mutate(y_obs = dbh,
               time = zeroed_time,
               treeid_factor = as.numeric(as.factor(treeid))) %>%
        select(treeid, treeid_factor, sp, species, y_obs, census, time)

    } else {
      sample_data <- species_data %>%
        mutate(y_obs = dbh,
               time = zeroed_time,
               treeid_factor = as.numeric(as.factor(treeid))) %>%
        select(treeid, treeid_factor, sp, species, y_obs, census, time)
    }

    #Plot histograms of size distributions
    if(plot_size_histograms){
      plot_data <- data_trees_transformed %>%
        filter(same_stem == 0) %>% #Removes records that may have different stems
        filter(sp == sp_codes[i])

      full_data_size <- ggplot(data = plot_data, aes(x=dbh)) +
        geom_density(fill = col_vec[i],
                       colour = "black") +
        labs(x = "Observed DBH",
             title = paste0(sample_data$species[1], ": all trees")) +
        theme_classic()

      sample_data_size <- ggplot(data = sample_data, aes(x=y_obs)) +
        geom_density(fill = col_vec[i],
                       colour = "black") +
        labs(x = "Observed DBH",
             title = "Sampled trees") +
        xlim(min(plot_data$dbh, na.rm = TRUE), max(plot_data$dbh, na.rm = TRUE)) +
        theme_classic()

      dens_plot <- plot_grid(
        full_data_size,
        sample_data_size,
        nrow = 2,
        labels = paste0("Size density for ", sample_data$species[1])
      )

      file_name <- paste0("output/figures/", sp_codes[i], "_SizeDens.svg")
      ggsave(file_name, plot=dens_plot, width=100, height=120, units="mm", device = "svg")
    }

    saveRDS(sample_data, file = paste0(out_path,
                                       sp_codes[i],
                                       "_SRSWR_sample.rds"))

    stan_data <- list(
      n_obs = nrow(sample_data),
      n_ind = length(unique(sample_data$treeid)),
      y_obs = sample_data$y_obs,
      obs_index = sample_data$census,
      time = sample_data$time,
      ind_id = sample_data$treeid_factor
    )
    saveRDS(stan_data, file = paste0(out_path,
                                     sp_codes[i],
                                     "_SRSWR_stan_data.rds"))
  }
}

