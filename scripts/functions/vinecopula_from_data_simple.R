#Wrapper function around rvinecopulib (spline-based)
library(kde1d)
library(rvinecopulib)

#estimate splines using fitdistrplus package and actuar package for extra distribution options
estimate_spline_marginal <- function(covariate, xmin = NaN) {
  covariate <- covariate[!is.na(covariate)]
  param <- kde1d(covariate, xmin = xmin)
  marg <- list(pdf = function(u) qkde1d(u, param),
               pit = function(x) pkde1d(x, param),
               rdist = function(n) rkde1d(n, param),
               density = function(x) dkde1d(x, param),
               dist = param)
  return(marg)
}


#object for estimation and methods contour(), plot() and simulate()

estimate_vinecopula_from_data <- function(dat, var_types = rep(c, ncol(dat)), ...) {
  
  dat_out <- dat

  #estimate the marginal splines for each variable and create uniform data set
  marginals <- list()
  dat_unif <- dat_out
  for (ii in 1:ncol(dat_out)) {
    marginals[[ii]] <- estimate_spline_marginal(dat_out[, ii])
    dat_unif[, variable] <- marginals[[ii]]$pit(dat_unif[, variable])
  }

  #estimate copula
  vine_coefs <- vinecop(dat_unif, ...)
  
  #create output object
  vine_output <- list(vine_copula = vine_coefs, marginals = marginals, 
                      polynomial = polynomial, 
                      names = c(ID_name = ID_name, time_name = time_name), 
                      time_range = time_range, 
                      variables_of_interest = variables_of_interest)
  if (keep_data) {
    vine_output <- append(vine_output, list(original_data = dat_out, uniform_data = dat_unif), after = 3)
  }
  
  class(vine_output) <- "estVineCopula"
  
  return(vine_output)
}

plot.estVineCopula <- function(vine_output, ...) {
  plot(vine_output$vine_copula, ...)
}

contour.estVineCopula <- function(vine_output, ...) {
  contour(vine_output$vine_copula, ...)
}

#simulation from estimated vine copula
#value_only = TRUE for longitudinal predictions, or FALSE for list with 
#   parameters and longitudinal predictions
simulate.estVineCopula <- function(vine_output, n, value_only = TRUE) {
  
  #Simulation
  if (any(vine_output$vine_copula$var_types == "d")) {
    vine_distribution <- vinecop_dist(vine_output$vine_copula$pair_copulas, vine_output$vine_copula$structure, var_types = vine_output$vine_copula$var_types)
    dat_sim <- as.data.frame(rvinecop(n, vine_distribution))
    names(dat_sim) <- vine_output$variables_of_interest[1:ncol(dat_sim)]
  } else {
    dat_sim <- as.data.frame(rvinecop(n, vine_output$vine_copula))
  }
  
  for (variable in names(dat_sim)) {
    ind_na <- is.na(dat_sim[, variable])
    dat_sim[!ind_na, variable] <- vine_output$marginals[[variable]]$pdf(dat_sim[!ind_na, variable])
  }
  if (any(vine_output$vine_copula$var_types == "d")) {
    dat_sim[, vine_output$vine_copula$var_types == "d"] <- round(dat_sim[, vine_output$vine_copula$var_types == "d"])
  }
  
  #polynomials
  if (vine_output$polynomial) {
    gest_times <- seq(vine_output$time_range[1], vine_output$time_range[2], length.out = 100)
    time_data <- as.data.frame(expand_grid(rownames(dat_sim), gest_times))
    names(time_data) <- vine_output$names
    suppressMessages(df_sim <- dat_sim %>% 
      rownames_to_column(vine_output$names["ID_name"]) %>%
      right_join(time_data))
    
    for (variable in vine_output$variables_of_interest) {
      col_ind <- grep(paste0("_", variable), names(df_sim))
      df_sim[, variable] <- df_sim[, col_ind[1]] + df_sim[, col_ind[2]]*df_sim[, vine_output$names["time_name"]] + 
        df_sim[, col_ind[3]]*df_sim[, vine_output$names["time_name"]]^2
    }
    output <- df_sim %>% dplyr::select(all_of(c(as.character(vine_output$names), vine_output$variables_of_interest)))
    if (value_only) {
      return(output)
    } else {
      return(list(values = output, parameters = dat_sim))
    }
  }
  
  return(dat_sim)
}

