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
  saveRDS(full_data, file = "output/data/full_est_data.rds")

  #Plots plots plots
  #fig_1_plots(full_data, sp_code = "gar2in")
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
                     col_vec)
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
                               growth_function, col_vec){
  print("Producing plots of growth histories.")
  species_code <- unique(individual_data$sp_code)
  for(i in 1:length(species_code)){
    #Pull out data for  species data
    temp_ind_data <- individual_data %>%
      filter(sp_code == species_code[i])
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

#-----------------------------------------------------------------------------#
#Fig. 1 plot
fig_1_plots <- function(){

}


#-----------------------------------------------------------------------------#
# Numerical integration
#Runge-Kutta 4th order
rk4_est <- function(S_0, growth, pars, step_size, N_step){
  runge_kutta_int <- c(S_0)
  for(i in 2:N_step){
    k1 <- growth(runge_kutta_int[i-1], pars)
    k2 <- growth((runge_kutta_int[i-1] + step_size*k1/2), pars)
    k3 <- growth((runge_kutta_int[i-1] + step_size*k2/2), pars)
    k4 <- growth((runge_kutta_int[i-1] + step_size*k3), pars)

    runge_kutta_int[i] <- runge_kutta_int[i-1] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*step_size
  }
  return(runge_kutta_int)
}

#RK4 algorithm with substeps
rk4_step <- function(y,  pars, interval, growth_function){
  k1 <- growth_function(y, pars)
  k2 <- growth_function(y+interval*k1/2.0, pars)
  k3 <- growth_function(y+interval*k2/2.0, pars)
  k4 <- growth_function(y+interval*k3, pars)

  y_hat <- y + (1.0/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4) * interval

  return(y_hat)
}

rk4 <- function(y, pars, interval, step_size, growth_function){
  duration <-  0
  y_hat <- y

  while(duration < interval){
    #Determine the relevant step size
    step_size_temp <- min(c(step_size, interval-duration))

    #Get next size estimate
    y_hat <- rk4_step(y_hat, pars, step_size_temp, growth_function)

    #Increment observed duration
    duration <- duration + step_size_temp
  }

  return(y_hat)
}
