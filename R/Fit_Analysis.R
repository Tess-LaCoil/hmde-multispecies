#Preliminaries
scatterplot_list <- list(
  gmax_plot = list(
    par="Log g_max",
    mean="pop_max_growth_mean",
    sd="pop_max_growth_sd"
  ),

  Smax_plot = list(
    par="Log S_max",
    mean="pop_size_at_max_growth_mean",
    sd="pop_size_at_max_growth_sd"
  ),

  k_plot = list(
    par="Log k",
    mean="pop_k_mean",
    sd="pop_k_sd"
  )
)

#Control function
run_analysis <- function(sp_codes,
                         col_vec,
                         est_file_path){
  #Construct tibbles containing all species data
  full_data <- construct_full_data_tibbles(sp_codes, est_file_path)
  saveRDS(full_data, file = "output/data/full_data.rds")

  #Tables
  species_summary_table <- build_sp_tables(full_data)
  write.csv(species_summary_table, file = "output/data/species_summary_table.csv")

  #Trait stuff
  trait_analysis(species_summary_table)

  #Plots plots plots
  plot_obs_est_size_scatter(full_data,
                            col_vec = c("#f8766d", "#e7861b", "#cf9400",
                                        "#afa100", "#00b81f", "#00bf7d",
                                        "#00c1aa", "#00bdd0", "#00b4ef",
                                        "#00a5ff", "#77a5ff", "#9590ff",
                                        "#cf78ff", "#f066ea", "#ff62bc",
                                        "#fc717f"),
                            sp_codes = c("alsebl", "beilpe", "cordbi", "faraoc",
                                         "hirttr", "jac1co", "pri2co", "protpa",
                                         "protte", "quaras", "swars1", "swars2",
                                         "simaam", "tachve", "tet2pa", "tri2tu"))

  plot_3d_scatter(ind_data = full_data$ind_data_full, col_vec)

  plot_ridges(ind_par_data = full_data$ind_data_full,
              growth_par_names = c("ind_max_growth",
                                   "ind_size_at_max_growth",
                                   "ind_k"),
              plot_par_names = c("Max growth rate (g_max) cm/yr",
                                 "Diameter at max growth (S_max) cm",
                                 "Spread parameter (k)"),
              measurement_data = full_data$measurement_data,
              col_vec = col_vec)

  plot_ind_par_pairs(full_data$ind_data_full,
                     c("ind_max_growth",
                       "ind_size_at_max_growth",
                       "ind_k"),
                     c("Max growth rate (g_max) cm/yr",
                       "Diameter at max growth (S_max) cm",
                       "Spread parameter (k)"))

  plot_growth_pieces(individual_data = full_data$ind_data_full,
                     measurement_data = full_data$measurement_data_full,
                     species_level_data = full_data$sp_model_data_full,
                     growth_par_names = c("ind_max_growth",
                                          "ind_size_at_max_growth",
                                          "ind_k"),
                     growth_function = hmde_model_des("canham_single_ind"),
                     col_vec,
                     exclude_vec = c(21260, 83062)) #Exclude extremes from G. recondita, H. Triandra

  six_species_focus(measurement_data = full_data$measurement_data_full,
                    individual_data = full_data$ind_data_full,
                    plot_initial_sizes = FALSE,
                    exclude_vec = c(21260, 83062), #Excluding two extreme individuals
                    focus_ind_list <- list(
                      faraoc = c(6, 37, 14, 38, 26), #F. occidentalis - faraoc done
                      gar2in = c(3, 12, 72, 126, 160), #G. recondita - gar2in done
                      hirttr = c(57, 58, 28, 8, 2), #H. triandra - hirttr done
                      jac1co = c(58, 50, 49, 1, 7), #J. copaia - jac1co done
                      simaam =  c(21, 23, 11, 163, 5), #S. amara  - simaam done
                      tachve = c(35, 31, 4, 8, 34) #T. panamensis - tachve done
                    ),
                    sp_codes <- c("faraoc", "gar2in", "hirttr",
                                  "jac1co", "simaam", "tachve"),
                    sp_names <- c("Faramea occidentalis", "Garcinia recondita",
                                  "Hirtella triandra", "Jacaranda copaia",
                                  "Simarouba amara", "Tachigali panamensis"),
                    col_vec <- c("#afa100", "#72b000", "#00b81f",
                                 "#00bf7d", "#cf78ff", "#f066ea"))

  fig_1_plots(chosen_sp_code = "gar2in",
              colour = "#72b000",
              focus_ind_vec = c(3, 12, 72, 126, 160),
              full_data,
              species_summary_table,
              exclude_vec = c(21260, 83062)
              )#Excluding two extreme individuals

  canham_amenability(
    chosen_sp_codes = c("faraoc", "gar2in", "hirttr",
                        "jac1co", "simaam", "tachve"),
    full_data,
    exclude_vec = c(21260, 83062),
    colours = c("#afa100", "#72b000", "#00b81f",
                "#00bf7d", "#cf78ff", "#f066ea")
  )
}

#-----------------------------------------------------------------------------#
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

  sp_model_data_full <- tibble(
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
      mutate(BCI_ind_id = treeid,
             ind_id = treeid_factor) %>%
      select(BCI_ind_id, ind_id) %>%
      distinct() %>%
      arrange(ind_id)

    #Measurement data
    ests$measurement_data <- ests$measurement_data %>%
      group_by(ind_id) %>%
      arrange(obs_index) %>%
      mutate(BCI_ind_id = ind_id_tibble$BCI_ind_id[ind_id],
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
      mutate(BCI_ind_id = ind_id_tibble$BCI_ind_id[ind_id],
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

    sp_model_data_full <- rbind(sp_model_data_full,
                                readRDS(paste0(est_file_path, sp_codes[i],
                                               "_sp_level_model_ests.rds")))
  }

  return_data <- list(
    measurement_data_full = measurement_data_full,
    ind_data_full = ind_data_full,
    sp_data_full = sp_data_full,
    sp_model_data_full = sp_model_data_full
  )
}


#-----------------------------------------------------------------------------#
#Construct tables for species-level data
build_sp_tables <- function(full_data){
  #Species names
  sp_names <- full_data$sp_data_full %>%
    select(sp_code, species) %>%
    distinct()

  #Convert to wide format
  hierarchical_summary_table <- full_data$sp_data_full %>%
    mutate(par_name_factor = as.factor(par_name)) %>%
    select(sp_code, par_name_factor, median) %>%
    spread(par_name_factor, median)

  #Renaming parameters to reflect log-normal distribution structure
  names(hierarchical_summary_table) <- c(
    "sp_code","pop_log_k_mean","pop_log_k_sd", "pop_log_max_growth_mean",
    "pop_log_max_growth_sd", "pop_log_size_at_max_growth_mean",
    "pop_log_size_at_max_growth_sd"
  )

  #Exponentiate means
  hierarchical_summary_table <- hierarchical_summary_table %>%
    mutate(pop_k_mean = exp(pop_log_k_mean),
           pop_max_growth_mean = exp(pop_log_max_growth_mean),
           pop_size_at_max_growth_mean = exp(pop_log_size_at_max_growth_mean))

  #Sample sizes
  ind_summary_size <- full_data$ind_data_full %>%
    group_by(sp_code) %>%
    summarise(samp = length(unique(BCI_ind_id)),
              median_ind_max_growth = median(ind_max_growth),
              median_ind_size_at_max_growth = median(ind_size_at_max_growth),
              median_ind_k = median(ind_k),
              gmax_95 = as.numeric(quantile(ind_max_growth, 0.95)),
              k_95 = as.numeric(quantile(ind_k, 0.95)),
              max_est_size = max(S_final),
              sd_ind_log_size_at_max_growth = sd(log(ind_size_at_max_growth)),
              sd_ind_log_k = sd(log(ind_k)))

  obs_growth_data <- full_data$measurement_data_full %>%
    group_by(BCI_ind_id) %>%
    mutate(interval = time - lag(time),
           delta_obs = (y_obs - lag(y_obs))/interval) %>%
    ungroup() %>%
    filter(!is.na(delta_obs)) %>%
    group_by(sp_code) %>%
    summarise(gobs_95 = as.numeric(quantile(delta_obs, 0.95)),
              Sobs_95 = as.numeric(quantile(y_obs, 0.95)),
              Sobs_max = max(y_obs)
    )

  #Get species-level model data
  sp_model_wide <- full_data$sp_model_data_full %>%
    mutate(par_name_factor = as.factor(par_name)) %>%
    select(sp_code, par_name_factor, median_est) %>%
    spread(par_name_factor, median_est) %>%
    rename(pop_level_model_k = pop_k,
           pop_level_model_max_growth = pop_max_growth,
           pop_level_model_size_at_max_growth = pop_size_at_max_growth)

  #Load in walltime data
  species_walltime <- read.csv("input/species_walltime.csv")

  #Trait data
  trait_table <-  read.csv("input/TraitTable.csv")

  species_summary_table <- left_join(sp_names, hierarchical_summary_table, by = "sp_code") %>%
    left_join(ind_summary_size, by = "sp_code") %>%
    left_join(obs_growth_data, by = "sp_code") %>%
    left_join(species_walltime, by = "sp_code") %>%
    left_join(sp_model_wide, by = "sp_code") %>%
    left_join(trait_table, by = "species") %>%
    mutate(
      log_pop_model_g_max = log(pop_level_model_max_growth),
      log_ind_g_max_95 = log(gmax_95),
      log_gobs_95 = log(gobs_95),
      log_pop_mean_max_growth = log(pop_max_growth_mean)
    )

  table_1 <- species_summary_table %>%
    select(species, sp_code, samp, walltime)
  write.csv(table_1, "output/data/Paper_Table_1.csv", row.names = FALSE)

  table_2 <- species_summary_table %>%
    select(species, pop_max_growth_mean, pop_log_max_growth_sd,
           pop_size_at_max_growth_mean, pop_log_size_at_max_growth_sd,
           pop_k_mean, pop_log_k_sd)
  table_2[,2:8] <- signif(table_2[,2:8], digits = 4)
  write.csv(table_2, "output/data/Paper_Table_2.csv", row.names = FALSE)

  table_3 <- species_summary_table %>%
    select(species,
           median_ind_size_at_max_growth,
           pop_size_at_max_growth_mean,
           pop_level_model_size_at_max_growth,
           max_est_size,
           Sobs_max)
  table_3[,2:5] <- round(table_3[,2:5], digits = 1)
  write.csv(table_3, "output/data/Paper_Table_3.csv", row.names = FALSE)

  table_4 <- species_summary_table %>%
    select(
      species,
      pop_max_growth_mean,
      gmax_95,
      gobs_95,
      pop_level_model_max_growth,
      median_ind_size_at_max_growth,
      pop_level_model_size_at_max_growth,
      median_ind_k,
      sd_ind_log_k,
      pop_level_model_k,
      Sobs_95, WD, Hmax
    )
  table_4_rounded <- table_4
  table_4_rounded[,2:11] <- signif(table_4_rounded[,2:11], digits = 4)
  write.csv(table_4_rounded, "output/data/Trait_Table.csv", row.names = FALSE)

  return(species_summary_table)
}

#-----------------------------------------------------------------------------#
trait_analysis <- function(species_summary_table){
  #All pairs
  pairs_comp_data <- species_summary_table %>%
    select(species,
           pop_max_growth_mean, gmax_95, gobs_95,
           pop_log_size_at_max_growth_sd, pop_level_model_max_growth,
           sd_ind_log_size_at_max_growth,
           WD, Hmax, mean_b_coeff
    )

  pairs_plots <- list()
  for(i in 8:10){
    for(j in 2:7){
      if(i != j){
        plot <- pairs_scatter_spearmans_plot(dataset = pairs_comp_data,
                                             trait = names(pairs_comp_data)[i],
                                             pair = names(pairs_comp_data)[j],
                                             x_lab = names(pairs_comp_data)[i],
                                             y_lab = names(pairs_comp_data)[j])
        pairs_plots <- append(pairs_plots, list(plot))
      }
    }
  }
  pairs_grid <- plot_grid(plotlist = pairs_plots,
                    nrow = 6,
                    byrow = FALSE,
                    align = "hv")
  file_name <- "output/figures/AllPairs_TraitScatter.svg"
  ggsave(file_name, plot=pairs_grid, width=12, height=18)

  #Specific comparisons for density, max size, and light response
  growth_comp_data <- species_summary_table %>%
    select(species,
           pop_level_model_max_growth,
           pop_max_growth_mean, gmax_95, gobs_95,
           WD, Hmax
    )

  #Plots and test output
  growth_scatterplot_set <- trait_analysis_stats(trait_data = growth_comp_data,
                       trait_name = c("WD", "Hmax"),
                       comp_pair_names = c("pop_level_model_max_growth",
                                           "gobs_95",
                                           "pop_max_growth_mean",
                                           "gmax_95"
                                           ),
                       plot_x_labs = c("Wood density", "Max height"),
                       plot_y_labs = c("Sp. avg. model g_max",
                                       "Obs. 95% growth",
                                       "Sp. mean g_max",
                                       "Ind. 95% g_max"
                                       ))
  growth_grid <- plot_grid(plotlist = growth_scatterplot_set,
            ncol = 2,
            byrow = FALSE,
            align = "hv")
  file_name <- "output/figures/SizeDens_TraitScatter.svg"
  ggsave(file_name, plot=growth_grid, width=7, height=14)

  #Table 4 of rank coeffs
  trait_corr_pairs <- list(
    c("pop_max_growth_mean","pop_level_model_max_growth"),
    c("pop_max_growth_mean","gmax_95"),
    c("pop_max_growth_mean","gobs_95"),
    c("gmax_95","gobs_95"),
    c("pop_log_size_at_max_growth_sd",
      "sd_ind_log_size_at_max_growth")
  )
  trait_corr_table <- tibble(par_names = c(), r_s = c(), p_val = c())

  for(i in 1:length(trait_corr_pairs)){
    test <- cor.test(species_summary_table[[trait_corr_pairs[[i]][1]]],
                     species_summary_table[[trait_corr_pairs[[i]][2]]],
                     method = "spearman")
    corr_temp <- tibble(
      par_names = paste0(trait_corr_pairs[[i]][1], ", ",
                         trait_corr_pairs[[i]][2]),
      r_s = test$estimate,
      p_val = test$p.value
    )
    trait_corr_table <- rbind(trait_corr_table, corr_temp)
  }
  write.csv(trait_corr_table, "output/data/trait_corr_table.csv")

  trait_corr_table_rounded <- trait_corr_table
  trait_corr_table_rounded[,2:3] <- signif(trait_corr_table[,2:3], digits = 4)
  write.csv(trait_corr_table_rounded, "output/data/trait_corr_table_rounded.csv")

  #Herault comparisons
  #Table of linear cors for g_max
  herault_growth_pars <- c(
    "log_pop_model_g_max",
    "log_ind_g_max_95"
  )

  #Plots and test output
  gmax_herault_scatterplot_set <- trait_analysis_stats(trait_data = species_summary_table,
                                                 trait_name = herault_growth_pars,
                                                 comp_pair_names = c("WD", "Hmax", "Sobs_95"),
                                                 plot_x_labs = c("Log Sp. avg. model g_max",
                                                                 "Log Ind. 95% g_max"),
                                                 plot_y_labs = c("Wood density",
                                                                 "Max height",
                                                                 "95% Obs. size")
                                                 )
  growth_grid <- plot_grid(plotlist = gmax_herault_scatterplot_set,
                           ncol = 2,
                           byrow = FALSE,
                           align = "hv")
  file_name <- "output/figures/log_gmax_Scatter.svg"
  ggsave(file_name, plot=growth_grid, width=7, height=10.5)

  h_max_cor <- species_summary_table %>%
    summarise(
      trait = "Hmax",
      r_log_g_max_95 = cor.test(Hmax,
                                log_ind_g_max_95,
                                method = "pearson")$estimate,
      r_log_g_max_95_p = cor.test(Hmax,
                                log_ind_g_max_95,
                                method = "pearson")$p.value,
      r_log_pop_model_g_max = cor.test(Hmax,
                                log_pop_model_g_max,
                                method = "pearson")$estimate,
      r_log_pop_model_g_max_p = cor.test(Hmax,
                                       log_pop_model_g_max,
                                       method = "pearson")$p.value
    )
  wd_cor <- species_summary_table %>%
    summarise(
      trait = "WD",
      r_log_g_max_95 = cor.test(WD,
                                log_ind_g_max_95,
                                method = "pearson")$estimate,
      r_log_g_max_95_p = cor.test(WD,
                                  log_ind_g_max_95,
                                  method = "pearson")$p.value,
      r_log_pop_model_g_max = cor.test(WD,
                                       log_pop_model_g_max,
                                       method = "pearson")$estimate,
      r_log_pop_model_g_max_p = cor.test(WD,
                                         log_pop_model_g_max,
                                         method = "pearson")$p.value
    )
  sobs_max_cor <- species_summary_table %>%
    summarise(
      trait = "sobs_95",
      r_log_g_max_95 = cor.test(Sobs_95,
                                log_ind_g_max_95,
                                method = "pearson")$estimate,
      r_log_g_max_95_p = cor.test(Sobs_95,
                                  log_ind_g_max_95,
                                  method = "pearson")$p.value,
      r_log_pop_model_g_max = cor.test(Sobs_95,
                                       log_pop_model_g_max,
                                       method = "pearson")$estimate,
      r_log_pop_model_g_max_p = cor.test(Sobs_95,
                                         log_pop_model_g_max,
                                         method = "pearson")$p.value
    )

  herault_growth_corr_table <- rbind(h_max_cor, wd_cor, sobs_max_cor)

  herault_growth_corr_table_rounded <- herault_growth_corr_table
  herault_growth_corr_table_rounded[,2:5] <- signif(herault_growth_corr_table[,2:5], digits = 4)
  write.csv(herault_growth_corr_table_rounded, "output/data/herault_growth_corr_table_rounded.csv")

  #Herault S_max comparisons
  S_max_pars <- c("median_ind_size_at_max_growth",
                  "pop_level_model_size_at_max_growth")
  #Plots and test output
  Smax_herault_scatterplot_set <- trait_analysis_stats(trait_data = species_summary_table,
                                                       trait_name = c("Hmax", "Sobs_max"),
                                                       comp_pair_names = S_max_pars,
                                                       plot_x_labs = c("Max height", "95% Obs. Size"),
                                                       plot_y_labs = c("Ind. 50% S_max",
                                                                       "Sp. Avg. model S_max"
                                                       ))
  smax_grid <- plot_grid(plotlist = Smax_herault_scatterplot_set,
                           ncol = 2,
                           byrow = FALSE,
                           align = "hv")
  file_name <- "output/figures/Smax_Scatter.svg"
  ggsave(file_name, plot=smax_grid, width=7, height=7)

  Hmax_cors <- species_summary_table %>%
    summarise(
      trait = "Hmax",
      r_median_ind_size_at_max_growth = cor.test(Hmax,
                              median_ind_size_at_max_growth,
                              method = "pearson")$estimate,
      r_median_ind_size_at_max_growth_p = cor.test(Hmax,
                                median_ind_size_at_max_growth,
                                method = "pearson")$p.value,
      r_pop_level_model_size_at_max_growth = cor.test(Hmax,
                                pop_level_model_size_at_max_growth,
                                method = "pearson")$estimate,
      r_pop_level_model_size_at_max_growth_p = cor.test(Hmax,
                                  pop_level_model_size_at_max_growth,
                                  method = "pearson")$p.value
    )
  s_ij_max_cors <- species_summary_table %>%
    summarise(
      trait = "size_max",
      r_median_ind_size_at_max_growth = cor.test(Sobs_max,
                                                 median_ind_size_at_max_growth,
                                                 method = "pearson")$estimate,
      r_median_ind_size_at_max_growth_p = cor.test(Sobs_max,
                                                   median_ind_size_at_max_growth,
                                                   method = "pearson")$p.value,
      r_pop_level_model_size_at_max_growth = cor.test(Sobs_max,
                                                      pop_level_model_size_at_max_growth,
                                                      method = "pearson")$estimate,
      r_pop_level_model_size_at_max_growth_p = cor.test(Sobs_max,
                                                        pop_level_model_size_at_max_growth,
                                                        method = "pearson")$p.value
    )
  s_ij_95_cors <- species_summary_table %>%
    summarise(
      trait = "size_95",
      r_median_ind_size_at_max_growth = cor.test(Sobs_95,
                                                 median_ind_size_at_max_growth,
                                                 method = "pearson")$estimate,
      r_median_ind_size_at_max_growth_p = cor.test(Sobs_95,
                                                   median_ind_size_at_max_growth,
                                                   method = "pearson")$p.value,
      r_pop_level_model_size_at_max_growth = cor.test(Sobs_95,
                                                      pop_level_model_size_at_max_growth,
                                                      method = "pearson")$estimate,
      r_pop_level_model_size_at_max_growth_p = cor.test(Sobs_95,
                                                        pop_level_model_size_at_max_growth,
                                                        method = "pearson")$p.value
    )
  Smax_cor <- rbind(Hmax_cors, s_ij_max_cors, s_ij_95_cors)
  herault_S_max_cor_table_rounded <- Smax_cor
  herault_S_max_cor_table_rounded[,2:5] <- signif(Smax_cor[,2:5], digits = 3)
  write.csv(herault_S_max_cor_table_rounded, "output/data/herault_S_max_cor_table_rounded.csv")


  #Herault k comparisons
  k_pars <- c("sd_ind_log_k", "pop_level_model_k")
  #Plots and test output
  k_herault_scatterplot_set <- trait_analysis_stats(trait_data = species_summary_table,
                                                    trait_name = c("WD"),
                                                    comp_pair_names = k_pars,
                                                    plot_x_labs = c("Wood density"),
                                                    plot_y_labs = c("SD Ind. k",
                                                                    "Sp. Avg. model k"
                                                       ))
  k_grid <- plot_grid(plotlist = k_herault_scatterplot_set,
                         ncol = 1,
                         byrow = FALSE,
                         align = "hv")
  file_name <- "output/figures/k_Scatter.svg"
  ggsave(file_name, plot=k_grid, width=3.5, height=7)

  wd_cor <- species_summary_table %>%
    summarise(
      trait = "WD",
      r_sd_ind_log_k = cor.test(WD,
                                sd_ind_log_k,
                                method = "pearson")$estimate,
      r_sd_ind_log_k_p = cor.test(WD,
                                  sd_ind_log_k,
                                  method = "pearson")$p.value,
      r_pop_level_model_k = cor.test(WD,
                               pop_level_model_k,
                               method = "pearson")$estimate,
      r_pop_level_model_k_p = cor.test(WD,
                                 pop_level_model_k,
                                 method = "pearson")$p.value,
    )
  herault_k_cor_table_rounded <- wd_cor
  herault_k_cor_table_rounded[,2:7] <- signif(wd_cor[,2:7], digits = 3)
  write.csv(herault_k_cor_table_rounded, "output/data/herault_k_cor_table_rounded.csv")
}

trait_analysis_stats <- function(trait_data,
                                 trait_name,
                                 comp_pair_names,
                                 plot_x_labs,
                                 plot_y_labs){
  scatterplot_set <- list()

  for(i in 1:length(trait_name)){
    for(j in 1:length(comp_pair_names)){
      print(paste0("Testing ", trait_name[i], " and ", comp_pair_names[j]))

      plot <- pairs_scatter_spearmans_plot(dataset = trait_data,
                                   trait = trait_name[i],
                                   pair = comp_pair_names[j],
                                   x_lab = plot_x_labs[i],
                                   y_lab = plot_y_labs[j])

      scatterplot_set <- append(scatterplot_set, list(plot))
    }
  }

  return(scatterplot_set)
}

pairs_scatter_spearmans_plot <- function(dataset,
                                         trait,
                                         pair,
                                         x_lab,
                                         y_lab,
                                         print_test = FALSE){
  test <- cor.test(dataset[[trait]],
                   dataset[[pair]], method = "spearman")

  if(print_test){
    print(test)
  }

  plot_data <- dataset[c(trait, pair)]
  names(plot_data) <- c("x", "y")

  if(test$p.value < 3^{-16}){
    plot_label <- paste(
      paste0("r_s = ", signif(test$estimate, digits = 3)),
      "p < 2.2e-16",
      sep = "\n"
    )
  } else {
    plot_label <- paste(
      paste0("r_s = ", signif(test$estimate, digits = 3)),
      paste0("p = ", formatC(test$p.value, format = "e", digits = 3)),
      sep = "\n"
    )
  }

  x_pos <- max(plot_data$x, na.rm = TRUE) -
    0.2* (max(plot_data$x, na.rm = TRUE) -
            min(plot_data$x, na.rm = TRUE))
  y_pos <- max(plot_data$y, na.rm = TRUE) -
    0.1* (max(plot_data$y, na.rm = TRUE) -
            min(plot_data$y, na.rm = TRUE))

  plot <- ggplot(data = plot_data, aes(x = x, y = y)) +
    geom_point(colour = "green4", size = 1) +
    labs(x = x_lab, y = y_lab) +
    annotate("text", x=x_pos, y=y_pos, label= plot_label) +
    theme_classic()

  return(plot)
}

#-----------------------------------------------------------------------------#
#Scatter plots of estimated and observed sizes for each species
plot_obs_est_size_scatter <- function(full_data, col_vec, sp_codes){
  plot_list <- list()
  for(i in 1:length(sp_codes)){
    plot_data <- full_data$measurement_data_full %>%
      filter(sp_code == sp_codes[i])

    plot_list[[i]] <- ggplot_obs_est_scatter(plot_data,
                                             colour = col_vec[i],
                                             title = plot_data$species[1])
  }

  obs_est_scatter <- plot_grid(
    plotlist = plot_list,
    nrow = 4,
    align = "hv"
  )

  file_name <- "output/figures/ObsEstScatter.svg"
  ggsave(file_name, plot=obs_est_scatter, width=1200, height=1200,
         units="px", device = "svg", dpi = 96)
}

ggplot_obs_est_scatter <- function(plot_data, colour, title){
  r_sq_est <- cor(plot_data$y_obs,
                  plot_data$y_hat)^2
  r_sq <- paste0("R^2 = ",
                 signif(r_sq_est,
                        digits = 3))

  plot <-
    ggplot(data=plot_data, aes(x=y_obs, y=y_hat)) +
    geom_point(colour = colour, shape = 19, alpha = 0.4, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                alpha = 1, colour = "black") +
    annotate("text",
             x = 0.3*max(plot_data$y_obs),
             y = 0.9*max(plot_data$y_hat),
             label = r_sq) +
    labs(y = "Estimated DBH cm",
         x="Observed DBH cm", title = title) +
    theme_classic()

  return(plot)
}

#-----------------------------------------------------------------------------#
# Ridge plots
plot_ridges <- function (ind_par_data, growth_par_names, plot_par_names,
                         measurement_data, col_vec){
  for(i in 1:length(growth_par_names)){
    plot_data <- ind_par_data %>%
      select(growth_par_names[i], BCI_ind_id, species, sp_code) %>%
      rename(val = growth_par_names[i]) %>%
      group_by(species) %>%
      mutate(median = median(val)) %>%
      ungroup()
    file_name <- paste("output/figures/", growth_par_names[i],".png", sep="")
    plot <- gg_plot_ridges(plot_data, col_vec, log10 = TRUE,
                           quantile = FALSE, par_name = plot_par_names[i])
    ggsave(file_name, plot=plot, width=130, height=100, units="mm")
  }

  plot_data <- measurement_data %>%
    mutate(val = (lead(y_hat) - y_hat)/(lead(time) - time)) %>%
    group_by(species) %>%
    mutate(median = median(val)) %>%
    ungroup()

  file_name <- "output/figures/GrowthRate_Density.png"
  plot <- gg_plot_ridges(plot_data, col_vec, quantile = FALSE,
                         par_name = "Annualised growth rates", log10=TRUE)
  ggsave(file_name, plot=plot, width=130, height=100, units="mm")
}

gg_plot_ridges <- function(plot_data, col_vec, log10 = FALSE, quantile = FALSE,
                           par_name){
  plot <- ggplot(plot_data, aes(x=val, y=reorder(x=species, X=median),
                                fill=species)) +
    geom_density_ridges() +
    xlab(par_name) +
    ylab("") +
    scale_fill_discrete(col_vec) +
    theme_classic() +
    theme(legend.position="none")

  if(log10){
    plot <- plot +
      scale_x_log10()
  }

  if(quantile){
    quant_data <- plot_data %>%
      group_by(species) %>%
      summarize(quant90 = quantile(val, 0.9),
                median = quantile(val, 0.5)) %>%
      arrange(median)

    for(j in 1:length(quant_data)){
      plot <- plot +
        geom_segment(data = quant_data,
                     aes(y = rank(median),
                         yend = rank(median) + 1,
                         x = quant90,
                         xend = quant90),
                     linewidth = 1)
    }
  }

  return(plot)
}

#-----------------------------------------------------------------------------#
#Pairs plots Ind Pars
plot_ind_par_pairs <- function(individual_data,
                               pairs_par_names,
                               pairs_plot_names){
  print("Producing pair plots.")
  #Iterate through all pairs of parameters
  for(i in 1:length(pairs_par_names)){
    for(j in  1:length(pairs_par_names)){
      if(i < j){
        #Build data frame
        plot_data <- tibble(species = individual_data$species,
                            sp_code = individual_data$sp_code)
        plot_data$x <- individual_data[[pairs_par_names[i]]]
        plot_data$y <- individual_data[[pairs_par_names[j]]]

        #Produce and export plot to file
        plot <- gg_pairplot(plot_data, pairs_plot_names[c(i,j)])
        file_name <- paste("output/figures/pair_plot_", i, j , ".png", sep="")
        ggsave(file_name, plot=plot, width=145, height=145, units="mm")

        #Produce grid of plots with facet_wrap()
        plot <- gg_pairplot(plot_data, pairs_plot_names[c(i,j)], grid=TRUE)
        file_name <- paste("output/figures/pair_plot_grid", i, j , ".png", sep="")
        ggsave(file_name, plot=plot, width=210, height=200, units="mm")
      }
    }
  }
}

gg_pairplot <- function(plot_data, axis_names, grid=FALSE){
  plot <- ggplot(data=plot_data, aes(x=x, y=y, col=species)) +
    geom_point(size=0.8, alpha=0.3) +
    labs(x = axis_names[1], y = axis_names[2]) +
    theme_classic()

  if(grid){
    plot <- plot +
      geom_point(data=plot_data[2:4], aes(x=x, y=y, group=sp_code),
                 colour="grey", alpha=0.5, size=0.8, shape=16) +
      geom_point(size=0.8, alpha=0.5, shape=10) +
      facet_wrap(~species, nrow = 4) +
      theme(legend.position="none")
  }

  plot <- plot + scale_x_log10() + scale_y_log10()

  return(plot)
}

#---------------------------------------------------------------------------#
# Growth function piece plot
plot_growth_pieces <- function(individual_data,
                               measurement_data,
                               species_level_data,
                               growth_par_names,
                               growth_function,
                               col_vec,
                               exclude_vec = c()){
  print("Producing plots of growth histories.")
  species_code <- unique(individual_data$sp_code)
  for(i in 1:length(species_code)){
    #Pull out data for  species data
    temp_ind_data <- individual_data %>%
      filter(sp_code == species_code[i],
             !BCI_ind_id %in% exclude_vec)
    n_ind <- nrow(temp_ind_data)

    #Get data frame of the growth parameter estimates
    growth_par_ests <- tibble(g_max = temp_ind_data$ind_max_growth,
                              s_max = temp_ind_data$ind_size_at_max_growth,
                              k = temp_ind_data$ind_k)

    sp_measurement_data <- measurement_data %>%
      filter(sp_code == species_code[i])

    sp_level_pars <- species_level_data %>%
      filter(sp_code == species_code[i])
    sp_par_ests <- c(
      g_max = sp_level_pars$median_est[which(sp_level_pars$par_name == "pop_max_growth")],
      s_max = sp_level_pars$median_est[which(sp_level_pars$par_name == "pop_size_at_max_growth")],
      k = sp_level_pars$median_est[which(sp_level_pars$par_name == "pop_k")]
    )

    #Plot growth history
    plot_ind_sp_level_models(growth_par_ests,
                             sp_measurement_data,
                             sp_par_ests,
                             growth_function,
                             n_ind,
                             S_0 = temp_ind_data$S_initial,
                             S_final = temp_ind_data$S_final,
                             colour = col_vec[i],
                             name = species_code[i])
  }
}

#Plot species-level model and growth functions pieces
plot_ind_sp_level_models <- function(growth_par_ests,
                                     sp_measurement_data,
                                     sp_par_ests,
                                     growth_function,
                                     n_ind,
                                     S_0,
                                     S_final,
                                     colour,
                                     name){
  #Growth function pieces
  print(paste("Plotting growth histories:", name, sep=" "))
  file_name <- paste("output/figures/", name,
                     "_Growth_History.png", sep="")
  function_plot <- ggplot_sample_growth_trajectories(post_pars = growth_par_ests,
                                                     growth_function,
                                                     max_growth_size = max(S_final),
                                                     min_growth_size = 0,
                                                     S_0,
                                                     S_final,
                                                     colour,
                                                     species = sp_measurement_data$species[1])
  ggsave(file_name, plot=function_plot, width=90, height=80, units="mm", device = "png")

  args_list <- list(pars = sp_par_ests)

  #Plot species average function over individual growth pieces
  function_plot_with_sp_avg <- function_plot +
    geom_function(fun=growth_function, args=args_list, alpha=1,
                  linetype = "dashed",
                  color="black", linewidth=0.8, xlim=c(1, max(S_final)))
  file_name <- paste("output/figures/", name, "_Growth_History_Sp_fn.png", sep="")
  ggsave(file_name, plot=function_plot_with_sp_avg, width=90, height=80, units="mm", device = "png")
}

ggplot_sample_growth_trajectories <- function(post_pars, growth_function,
                                              max_growth_size, min_growth_size,
                                              S_0, S_final, colour, species){
  plot <- ggplot() +
    xlim(min_growth_size, max_growth_size) +
    labs(x = "Size (DBH) cm", y="Growth rate cm/yr", title = species) +
    theme_classic()

  for(i in 1:nrow(post_pars)){
    args_list <- list(pars=post_pars[i,])
    plot <- plot +
      geom_function(fun=growth_function, args=args_list, alpha=0.2,
                    color=colour, linewidth=1, xlim=c(S_0[i], S_final[i]))
  }

  plot <- plot +
    geom_hline(yintercept=0, colour="black")

  return(plot)
}

#-----------------------------------------------------------------------------#
#Analysis of 6 focus species
six_species_focus <- function(measurement_data,
                              individual_data,
                              plot_initial_sizes = FALSE,
                              exclude_vec = c(),
                              focus_ind_list,
                              sp_codes,
                              sp_names,
                              col_vec){
  if(plot_initial_sizes){
    for(i in c(4)){ #1:length(sp_codes)){
      temp_measurement_data <- measurement_data %>%
        filter(sp_code == sp_codes[i])
      temp_ind_data <- individual_data %>%
        filter(sp_code == sp_codes[i])

      plot_obs_and_est_life_history(temp_measurement_data,
                                    nrow(temp_ind_data),
                                    sp_codes[i], colour = col_vec[i])
    }
  }

  for(i in 1:length(sp_codes)){
    focus_ind_vec <- focus_ind_list[[sp_codes[i]]]
    temp_measurement_data <- measurement_data %>%
      filter(sp_code == sp_codes[i],
             !BCI_ind_id %in% exclude_vec)
    temp_ind_data <- individual_data %>%
      filter(sp_code == sp_codes[i],
             !BCI_ind_id %in% exclude_vec)

    plot_focus_ind_figs(focus_inds = focus_ind_vec,
                        measurement_data = temp_measurement_data,
                        ind_data = temp_ind_data,
                        growth_function = hmde_model_des("canham_multi_ind"),
                        species = temp_ind_data$species[1],
                        colour = col_vec[i])
  }
}

#----------------------------------------------------------------------------#
# Observed and estimated life histories
#Plots a sample of individuals to show how the model estimates behave against the observed size
plot_obs_and_est_life_history <- function(plotting_data, n_ind, species, colour){
  for(j in 1:n_ind){
    if((j-1)%%100 == 0){print(paste0("Plotting sizes over time ", j))}

    sample <- plotting_data %>%
      filter(ind_id == j) %>%
      select(y_obs, time, BCI_ind_id, y_hat, ind_id) %>%
      mutate(time = 1990 + round(time, digits=0))

    #Build data frame
    data <- data.frame(size=sample$y_obs,
                       time=sample$time,
                       cond=rep("Obs.", times=length(sample$time)))
    data <- rbind (data, data.frame(size=sample$y_hat,
                                    time=sample$time,
                                    cond=rep("Est.", times=length(sample$time))))

    #Produce plot
    file_name <- paste("output/figures/sampled/", species,
                       "_Sizes_",
                       j,".png", sep="")

    plot <- ggplot_obs_and_est_life_history(data, title = sample$BCI_ind_id[1], colour)

    ggsave(file_name, plot=plot, width=130, height=100, units="mm")
  }
}

#Produces plot for plot_obs_and_est_life_history()
ggplot_obs_and_est_life_history<- function(data, title, colour){
  plot <- ggplot(data=data, aes(x, y)) +
    geom_line(aes(x=time, y=size, color=as.factor(cond),
                  group=cond, linetype=as.factor(cond)), linewidth=1.1) +
    geom_point(aes(x=time, y=size, color=as.factor(cond),
                   shape=as.factor(cond),
                   group=cond), size=2.5) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_color_manual(values = c("black", colour)) +
    xlab("Year") +
    ylab("Size (cm)") +
    ggtitle(title) +
    labs(color = NULL) +
    theme_classic()

  return(plot)
}

#--------------------------------------------------------------------------#
#Multipe ind. focus plots
plot_focus_ind_figs <- function(focus_inds, measurement_data, ind_data, growth_function, species, colour){
  #Plot sizes over time
  focus_ind_sizes <- measurement_data %>%
    filter(ind_id %in% focus_inds)
  plot_multi_ind_life_history(focus_ind_sizes, species, colour)

  #Plot growth function pieces
  plot_growth_focus(focus_inds, ind_data, growth_function, species, colour)
}

#Plot growth functions with focused individuals highlighted
plot_growth_focus <- function(focus_inds, ind_data, growth_function, species, colour){
  focus_ind_pars <- ind_data %>%
    filter(ind_id %in% focus_inds)

  post_pars <- data.frame(g_max = ind_data$ind_max_growth,
                          s_max = ind_data$ind_size_at_max_growth,
                          k = ind_data$ind_k)
  focus_post_pars <- data.frame(g_max = focus_ind_pars$ind_max_growth,
                                s_max = focus_ind_pars$ind_size_at_max_growth,
                                k = focus_ind_pars$ind_k)

  plot <- ggplot_sample_growth_trajectories(post_pars,
                                            growth_function,
                                            max_growth_size = max(ind_data$S_final),
                                            min_growth_size = 1,
                                            S_0 = ind_data$S_initial,
                                            S_final = ind_data$S_final,
                                            colour = colour,
                                            species = ind_data$species[1])
  for(i in 1:nrow(focus_post_pars)){
    args_list <- list(pars=focus_post_pars[i,])
    plot <- plot +
      geom_function(fun=growth_function, args=args_list, alpha=1,
                    color="#333333", linewidth=1, xlim=c(focus_ind_pars$S_initial[i],
                                                         focus_ind_pars$S_final[i]))
  }
  ggsave(paste0("output/figures/", species, "_GrowthFocus.svg"),
         plot=plot, width=99, height=72, units="mm", device = "svg")
}

#Plot multiple individuals life history in same figure
plot_multi_ind_life_history<- function(plotting_data, species, colour){
  #Build data frame
  data <- data.frame(size=plotting_data$y_obs,
                     time=plotting_data$time,
                     cond=rep("Observed", times=length(plotting_data$time)),
                     id_num = paste(plotting_data$BCI_ind_id, "Obs") )
  data <- rbind(data,
                data.frame(size=plotting_data$y_hat,
                           time=plotting_data$time,
                           cond=rep("Estimated", times=length(plotting_data$time)),
                           id_num = paste(plotting_data$BCI_ind_id, "Est")
                )
  )
  data <- data %>%
    mutate(time = 1990 + round(time, digits=0))

  #Produce plot
  file_name <- paste("output/figures/", species,
                     "_Sizes.svg", sep="")

  plot <- ggplot_multi_ind_life_history(data, species, colour) +
    theme(legend.position="none")

  ggsave(file_name, plot=plot, width=99, height=72, units="mm", device = "svg")

}

#Produces plot_multi_ind_life_history()
ggplot_multi_ind_life_history <- function(data, species, colour){
  plot <- ggplot(data=data, aes(x, y)) +
    geom_line(aes(x=time, y=size,
                  color=as.factor(cond),
                  alpha = as.factor(cond),
                  group = id_num,
                  linetype=as.factor(cond)), linewidth=1.1) +
    geom_point(aes(x=time, y=size, color=as.factor(cond),
                   shape=as.factor(cond)), size=2.5) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_color_manual(values = c(colour, "black")) +
    scale_alpha_manual(values = c(1, 0.5)) +
    xlab("Year") +
    ylab("Size (DBH) cm") +
    ggtitle(species) +
    labs(color = NULL, shape=NULL, linetype = NULL, alpha=NULL) +
    theme_classic() +
    theme(legend.position = "bottom")

  return(plot)
}

#-----------------------------------------------------------------------------#
#3D scatter plot of individuals coloured by species
plot_3d_scatter <- function(ind_data, col_vec){
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

#-----------------------------------------------------------------------------#
#Fig. 1 plot
fig_1_plots <- function(chosen_sp_code, colour, focus_ind_vec,
                        full_data, species_summary_table,
                        growth_function = hmde_canham_de,
                        exclude_vec = NULL){
  species_summary_vec <- species_summary_table %>%
    filter(sp_code == chosen_sp_code)

  pars_combo <- c(
    g_max = species_summary_vec$pop_level_model_max_growth,
    y_max = species_summary_vec$pop_level_model_size_at_max_growth,
    k = species_summary_vec$pop_level_model_k
  )

  ind_data <- full_data$ind_data_full %>%
    filter(sp_code == chosen_sp_code,
           !BCI_ind_id %in% exclude_vec)
  measurement_data <- full_data$measurement_data_full %>%
    filter(sp_code == chosen_sp_code,
           !BCI_ind_id %in% exclude_vec) %>%
    group_by(ind_id) %>%
    mutate(obs_index = rank(time)) %>%
    ungroup()

  #Get size estimates from species-level model
  sp_level_measurement_data <- tibble()
  for(i in ind_data$ind_id){ #Iterate through individuals and project forward from s_0
    temp <- measurement_data %>%
      filter(ind_id == i)

    yini <- c(Y = temp$y_obs[1])
    size_est <- ode(yini,
        times = temp$time,
        func = Canham_DE,
        parms = pars_combo,
        method = "ode45")[,2]
    temp$y_hat <- size_est

    sp_level_measurement_data <- rbind(sp_level_measurement_data, temp)
  }

  #Data for plots of sizes over time
  plotting_data <- sp_level_measurement_data %>%
    filter(ind_id %in% focus_ind_vec)

  data_h <- rbind(data.frame(size=plotting_data$y_obs,
                             time=plotting_data$time,
                             cond=rep("Observed", times=length(plotting_data$time)),
                             id_num = paste(plotting_data$BCI_ind_id, "Obs") ),
                  data.frame(size=plotting_data$y_hat,
                             time=plotting_data$time,
                             cond=rep("Estimated", times=length(plotting_data$time)),
                             id_num = paste(plotting_data$BCI_ind_id, "Est")
                  )) %>%
    mutate(time = 1990 + round(time, digits=0))

  data_sp <- rbind(data.frame(size=plotting_data$y_obs,
                              time=plotting_data$time,
                              cond=rep("Observed", times=length(plotting_data$time)),
                              id_num = paste(plotting_data$BCI_ind_id, "Obs") ),
                   data.frame(size=plotting_data$y_hat_species,
                              time=plotting_data$time,
                              cond=rep("Estimated", times=length(plotting_data$time)),
                              id_num = paste(plotting_data$BCI_ind_id, "Est")
                   )) %>%
    mutate(time = 1990 + round(time, digits=0))

  #Produce size over time plots
  h_sizes_plot <- ggplot_multi_ind_life_history(data_h,
                                                species = "G. recondita: hierarchical fit",
                                                colour = colour) +
    theme(text = element_text(size = 10),
          legend.position = "inside",
          legend.position.inside = c(0.1, 0.9))
  sp_sizes_plot <- ggplot_multi_ind_life_history(data_sp,
                                                 species = "G. recondita: species-level fit",
                                                 colour = colour) +
    theme(text = element_text(size = 10),
          legend.position = "inside",
          legend.position.inside = c(0.1, 0.9))

  #Produce growth function plots
  focus_ind_pars <- ind_data %>%
    filter(ind_id %in% focus_ind_vec)

  post_pars <- data.frame(g_max = ind_data$ind_max_growth,
                          s_max = ind_data$ind_size_at_max_growth,
                          k = ind_data$ind_k)
  focus_post_pars <- data.frame(g_max = focus_ind_pars$ind_max_growth,
                                s_max = focus_ind_pars$ind_size_at_max_growth,
                                k = focus_ind_pars$ind_k)

  h_g_plot <- ggplot_sample_growth_trajectories(post_pars,
                                                hmde_canham_de,
                                                max_growth_size = max(ind_data$S_final),
                                                min_growth_size = 1,
                                                S_0 = ind_data$S_initial,
                                                S_final = ind_data$S_final,
                                                colour = "#72b000",
                                                species = "Growth functions fit to individuals") +
    theme(text = element_text(size = 10))

  for(i in 1:nrow(focus_post_pars)){
    args_list <- list(pars=focus_post_pars[i,])
    h_g_plot <- h_g_plot +
      geom_function(fun=hmde_canham_de, args=args_list, alpha=1,
                    color="#333333", linewidth=1, xlim=c(focus_ind_pars$S_initial[i],
                                                         focus_ind_pars$S_final[i]))
  }

  #Extract y limits from h_g_plot
  ylims <- c(0, layer_scales(h_g_plot)$y$range$range[2])
  args_list <- list(pars=as.numeric(species_summary_vec[1,c(22, 23, 21)]))
  sp_g_plot <- ggplot() +
    geom_function(fun=hmde_canham_de, args=args_list, alpha=1,
                  color="#333333", linewidth=1, xlim=c(1,
                                                       max(ind_data$S_final))) +
    ylim(ylims) +
    labs(x = "DBH cm", y="Growth rate cm/yr", title = "Species average model") +
    geom_hline(yintercept=0, colour="black") +
    theme_classic() +
    theme(text = element_text(size = 10))

  #Scatter plot of size
  size_scatter_data <- rbind(data.frame(y_obs=measurement_data$y_obs,
                                        y_hat = measurement_data$y_hat,
                                        cond="Hierarchical fit"),
                             data.frame(y_obs=sp_level_measurement_data$y_obs,
                                        y_hat = sp_level_measurement_data$y_hat_species,
                                        cond="Species-level fit"))

  #Separate size scatter plots for model fits
  ylims <- c(
    min(size_scatter_data$y_hat),
    max(size_scatter_data$y_hat)
  )
  h_size_scatter <- ggplot_obs_est_scatter(plot_data = measurement_data,
                                           colour = "#72b000",
                                           title = "Hierarchical fit size estimation") +
    ylim(ylims[1], ylims[2]) +
    theme(text = element_text(size = 10))

  sp_size_scatter <- ggplot_obs_est_scatter(plot_data = sp_level_measurement_data,
                                            colour = "#72b000",
                                            title = "Species average fit size estimation") +
    ylim(ylims[1], ylims[2]) +
    theme(text = element_text(size = 10))

  fig_1_grid <- plot_grid(h_sizes_plot,
                          sp_sizes_plot,
                          h_g_plot,
                          sp_g_plot,
                          h_size_scatter,
                          sp_size_scatter,
                          nrow = 3,
                          align = "hv",
                          axis = "b",
                          rel_heights = c(0.4, 0.3, 0.3),
                          labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))

  file_name <- "output/figures/Fig1Grid_SizeSep.svg"
  ggsave(file_name, plot=fig_1_grid, width=700, height=1150,
         units="px", device = "svg", dpi = 96)
}

# Numerical integration RK45 setup
#Create DE function for deSolve
Canham_DE <- function(Time, State, Pars) { #Pars: g_max, y_max, k
  with(as.list(c(State, Pars)), {
    dY <- g_max * exp(-0.5 * (log(Y / y_max) / k)^2)

    return(list(c(dY)))
  })
}

#-----------------------------------------------------------------------------#
#Canham amenability comparison
canham_amenability <- function(chosen_sp_codes,
                               full_data,
                               exclude_vec,
                               colours){
  plot_list <- list()
  for(i in 1:length(chosen_sp_codes)){
    ind_data <- full_data$ind_data_full %>%
      filter(sp_code == chosen_sp_codes[i],
             !BCI_ind_id %in% exclude_vec) %>%
      mutate(size_rank = rank(S_final))

    #Focus the three biggest individuals
    focus_ind <- ind_data %>%
      filter(size_rank %in% (nrow(ind_data)-3):nrow(ind_data))

    post_pars <- data.frame(g_max = ind_data$ind_max_growth,
                            s_max = ind_data$ind_size_at_max_growth,
                            k = ind_data$ind_k)

    #Produce plot of all growth functions
    plot <- ggplot_sample_growth_trajectories(post_pars,
                                              growth_function = hmde_canham_de,
                                              max_growth_size = max(ind_data$S_final),
                                              min_growth_size = 1,
                                              S_0 = ind_data$S_initial,
                                              S_final = ind_data$S_final,
                                              colour = colours[i],
                                              species = ind_data$species[1])

    #Overlay largest individuals
    for(j in 1:3){
      args_list <- list(pars=as.numeric(focus_ind[j,3:5]))
      plot <- plot +
        geom_function(fun=hmde_canham_de, args=args_list, alpha=1,
                      color="#111111", linewidth=1, xlim=c(focus_ind$S_initial[j],
                                                             focus_ind$S_final[j])) +
        geom_function(fun=hmde_canham_de, args=args_list, alpha=1, linetype = "longdash",
                      color="#444444", linewidth=0.5, xlim=c(1, focus_ind$S_initial[j]))
    }


    plot_list[[i]] <- plot
  }

  amenability_plot <- plot_grid(plotlist = plot_list,
            nrow = 3,
            align = "hv",
            byrow = FALSE,
            labels = c("(a)", "(b)",
                       "(c)", "(d)",
                       "(e)", "(f)"))

  file_name <- "output/figures/CanhamAmenability.svg"
  ggsave(file_name, plot=amenability_plot, width=900, height=1050,
         dpi = 96, units="px", device = "svg")
}
