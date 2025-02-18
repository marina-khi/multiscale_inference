size_and_power_calculations <- function(grid_type_ = "normal", seed_ = 1, n_ts_ = 2,
                                        beta_ = c(1, 1, 1), a_x_vec_ = c(0.25, 0.25, 0.25), phi_ = 0.25,
                                        a_ = 0.25, sigma_ = 0.25,
                                        rho_ = 0.25,
                                        n_rep_ = 1000, sim_runs_ = 1000,
                                        different_T_ = c(100, 250, 500), different_alpha_ = c(0.01, 0.05, 0.1), different_b_ = c(0),
                                        q_ = 25, r_ = 10, numCores_ = 2,
                                        filename_ext_ = ""){
  source("functions/functions.R")
  size_and_power_array <- array(NA, dim = c(length(different_T_),
                                            length(different_b_),
                                            length(different_alpha_)),
                                dimnames = list(t = different_T_,
                                                b = different_b_,
                                                alpha = different_alpha_))
  #Constructing the set of pairwise comparisons
  ijset <- expand.grid(i = 1:n_ts_, j = 1:n_ts_)
  ijset <- ijset[ijset$i < ijset$j, ]
  
  for (t_len in different_T_){
    set.seed(seed_)
    k <- match(t_len, different_T_)
    #Constructing the grid
    if (grid_type_ == "normal"){
      u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
      h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
      h_grid <- h_grid[h_grid > log(t_len) / t_len]
      grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
    } else if (grid_type_ == "dense"){
      u_grid <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
      h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 2 / t_len)
      h_grid <- h_grid[h_grid > log(t_len) / t_len]
      grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
    } else if (grid_type_ == "dyadic"){
      #Constructing the very sparse grid
      h_min  <- ceiling(log(t_len))/t_len
      K_seq  <- seq(from = 0, to = t_len, by = 1)
      K_seq  <- K_seq[2^K_seq * h_min < 0.25]
      h_grid <- 2^K_seq * h_min
      
      grid           <- list()
      grid$grid_type <- "non-default"
      grid$gset      <- data.frame()
      
      for (h in h_grid){
        s_seq     <- seq(from = 0, to = floor((1/ h - 1)/2), by = 1) 
        u_grid    <- (2 * s_seq + 1) * h  
        gset      <- expand.grid(u = u_grid, h = h)
        grid$gset <- rbind(grid$gset, gset)
      }
      
      grid$gset_full <- grid$gset
      grid$pos_full  <- rep(TRUE, dim(grid$gset)[1])
      
      grid$gset <- grid$gset[grid$pos_full, ]
      grid$bws  <- unique(grid$gset[, 2])
      grid$lens <- rep(NA, length(grid$bws))
      for (i in seq_len(length(grid$bws)))
        grid$lens[i] <- sum(grid$gset[, 2] == grid$bws[i])
    } else {
      cat("The type of grid is not supported")
    }
    
    #Calculating the Gaussian quantiles in parallel
    tic()
    cl <- makePSOCKcluster(numCores_)
    registerDoParallel(cl)
    foreach (val = 1:sim_runs_, .combine = "cbind") %dopar% {
      source("functions/functions.R")
      repl(rep_ = val, n_ts_ = n_ts_, t_len_ = t_len, grid_ = grid,
           gaussian_sim = TRUE)
      # Loop one-by-one using foreach
    } -> simulated_pairwise_gaussian
    stopCluster(cl)
    toc()
    
    simulated_gaussian <- apply(simulated_pairwise_gaussian, 2, max)
    
    probs      <- seq(0.5, 0.995, by = 0.005)
    quantiles  <- as.vector(quantile(simulated_gaussian, probs = probs))
    quantiles  <- rbind(probs, quantiles)
    
    colnames(quantiles) <- NULL
    rownames(quantiles) <- NULL
    
    quants <- as.vector(quantiles[2, ])
    
    tic()
    cl <- makePSOCKcluster(numCores_)
    registerDoParallel(cl)
    foreach (val = 1:n_rep_, .combine = "cbind") %dopar% {
      source("functions/functions.R")
      repl(rep_ = val, n_ts_ = n_ts_, t_len_ = t_len, grid_ = grid,
           a_ = a_, sigma_ = sigma_,
           beta_ = beta_, a_x_vec_ = a_x_vec_, phi_ = phi_,
           rho_ = rho_, different_b_ = different_b_,
           q_ = q_, r_ = r_, gaussian_sim = FALSE)
      # Loop one-by-one using foreach
    } -> simulated_pairwise_statistics
    stopCluster(cl)
    toc()
    
    for (j in 1:length(different_b_)){
      simulated_statistic <- apply(simulated_pairwise_statistics[((j - 1) * n_ts_ * n_ts_ + 1):(j * n_ts_ * n_ts_), ], 2, max)
      
      size_and_power_vec <- c()
      for (alpha in different_alpha_){
        if (sum(probs == (1 - alpha)) == 0)
          pos <- which.min(abs(probs - (1 - alpha)))
        if (sum(probs == (1 - alpha)) != 0)
          pos <- which.max(probs == (1 - alpha))    
        quant <- quants[pos]
        
        num_of_rej         <- sum(simulated_statistic > quant)/n_rep_
        size_and_power_vec <- c(size_and_power_vec, num_of_rej) 
        
        cat("Ratio of rejection is ", num_of_rej, "with b = ", different_b_[j],
            ", alpha = ", alpha, "and T = ", t_len, "\n")
      }
      
      #Storing the results in a 3D array
      size_and_power_array[k, j, ] <- size_and_power_vec
    }
  }
  
  
  #Output of the results
  for (b in different_b_){
    l   <- match(b, different_b_)
    tmp <- as.matrix(size_and_power_array[, l, ])
    if (b == 0){
      filename = paste0("output/revision/", n_ts, "_ts_", phi_ * 100, "_",
                        rho_ * 100, "_size", filename_ext_, ".tex")
    } else {
      filename = paste0("output/revision/", n_ts, "_ts_", phi_ * 100, "_",
                        rho_ * 100, "_power_b_", b * 100, filename_ext_, ".tex")
    }
    output_matrix(tmp, filename, numcols_ = 4)
    line <- paste0("%This simulation was done for the seed ", seed_,
                   ", for the following values of the parameters: n_ts = ", n_ts_,
                   ", with ", n_rep_, " simulations for calculating size and power and ", sim_runs_,
                   " simulations to calculate the Gaussian quantiles. Furthermore, for the error process we have a = ",
                   a_, " and sigma = ", sigma_, 
                   ". For the covariate process a_1 = a_2 = a_3 = ", a_x_vec_[1], " and phi = ", phi_,
                   ". For the fixed effect, we have rho = ", rho_,
                   ". The grid is ", grid_type_)     
    write(line, file = filename, append = TRUE)
  }
}