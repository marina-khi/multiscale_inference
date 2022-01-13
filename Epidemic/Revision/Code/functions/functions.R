#Local linear estimator of the trend function 
#using the rectangular kernel. 
nadaraya_watson_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    result = sum((abs((grid_p - u) / bw) <= 1) * data_p)
    norm = sum((abs((grid_p - u) / bw) <= 1))
    return(result/norm)
  }
}


produce_plots <- function (results, l, data_i, data_j,
                           gov_resp_i, gov_resp_j,
                           country_i, country_j, filename){
  Tlen <- length(data_i)
  gset <- results$gset_with_values[[l]]
  
  pdf(filename, width=5.5, height=13, paper="special")
  layout(matrix(c(1, 2, 3, 4),ncol=1), widths=c(2.2, 2.2, 2.2, 2.2),
         heights=c(1.5, 1.5, 1.5, 1.8), TRUE)

  #Setting the layout of the graphs

  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2.5, 0)) #Margins for each plot
  par(oma = c(0.2, 0.5, 1, 0.2)) #Outer margins


  plot(data_i, ylim=c(min(data_i, data_j), max(data_i, data_j)), type="l",
      col="black", ylab="", xlab="", mgp=c(1, 0.5, 0), cex.axis = 0.8)
  lines(data_j, col="red")

  title(main = "(a) observed new cases per day", font.main = 1, line = 0.5, cex.main = 0.8)
  legend("topright", inset = 0.02, legend=c(country_i, country_j),
         col = c("black", "red"), lty = 1, cex = 0.75, ncol = 1, bg = 'white')

  par(mar = c(0.5, 0.5, 2.5, 0)) #Margins for each plot

  #Plotting the smoothed version of the time series that we have
  grid_points <- seq(from = 1 / Tlen, to = 1, length.out = Tlen) #grid points for estimating
  smoothed_i  <- mapply(nadaraya_watson_smoothing, grid_points,
                        MoreArgs = list(data_i, grid_points, bw = 3.5 / Tlen))
  smoothed_j  <- mapply(nadaraya_watson_smoothing, grid_points,
                        MoreArgs = list(data_j, grid_points, bw = 3.5 / Tlen))
  
  plot(smoothed_i, ylim=c(min(data_i, data_j), max(data_i, data_j)), type="l",
       col="black", ylab="", xlab = "", mgp=c(1,0.5,0), cex.axis = 0.8)
  title(main = "(b) smoothed curves from (a)", font.main = 1, line = 0.5, cex.main = 0.8)
  lines(smoothed_j, col="red")
   
  plot(gov_resp_i, ylim=c(0, 100), type="l",
        col="black", ylab="", xlab = "", mgp=c(1, 0.5, 0), cex.axis = 0.8)
  title(main = "(c) government response index", font.main = 1, line = 0.5, cex.main = 0.8)
  lines(gov_resp_j, col="red")

  par(mar = c(2.7, 0.5, 2.5, 0)) #Margins for each plot

  a_t_set <- subset(gset, test == TRUE, select = c(u, h))
  if (nrow(a_t_set) > 0){
    p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h) * Tlen + 0.5,
                          'endpoint' = (a_t_set$u + a_t_set$h) * Tlen - 0.5, 'values' = 0)
    p_t_set$values <- (1:nrow(p_t_set))/nrow(p_t_set)
    
    #Produce minimal intervals
    p_t_set2  <- compute_minimal_intervals(p_t_set)

    plot(NA, xlim=c(0, Tlen),  ylim = c(0, 1 + 1 / nrow(p_t_set)), xlab="", mgp=c(2, 0.5, 0), yaxt = "n", cex.axis = 0.8)
    title(main = "(d) (minimal) intervals produced by our test", font.main = 1, cex.main = 0.8, line = 0.5)
    title(xlab = "days since first Monday after the hundredth case", line = 1.7, cex.lab = 0.8)
    segments(p_t_set2$startpoint, p_t_set2$values, p_t_set2$endpoint, p_t_set2$values, lwd = 2)
    segments(p_t_set$startpoint, p_t_set$values, p_t_set$endpoint, p_t_set$values, col = "gray")
    print.xtable(xtable(p_t_set, digits = c(3), align = paste(replicate(4, "c"), collapse = "")),
                 type="latex", file=paste0("plots/", country_i, "_vs_", country_j, ".tex"),
                 include.colnames = FALSE)
    print.xtable(xtable(p_t_set2, digits = c(3), align = paste(replicate(4, "c"), collapse = "")),
                 type="latex", file=paste0("plots/", country_i, "_vs_", country_j, "_min_intervals.tex"),
                 include.colnames = FALSE)
  } else {
    #If there are no intervals where the test rejects, we produce empty plots

    plot(NA, xlim=c(0, Tlen),  ylim = c(0, 1), xlab="", ylab = "", mgp=c(2,0.5,0), yaxt = "n", cex.axis = 0.7)
    title(main = "(d) (minimal) intervals produced by our test", font.main = 1, cex.main = 0.8, line = 0.5)
    title(xlab = "days since first Monday after the hundredth case", line = 1.7, cex.lab = 0.8)
  }
  mtext(paste0("Comparison of ", country_i, " and ", country_j), side = 3, line = -1, outer = TRUE, font = 1, cex = 0.9)
  dev.off()
}

# functions for data simulations
lambda_fct <- function(u, c = 1000, height = 5000, position = 10) {
  return (height * exp(-(position * u - 3) ^ 2 / 2) + c)
}

r_doublepois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta - 1))
}

simulate_data <- function(n_ts, t_len, lambda_vec, sigma) {
  data <- matrix(0, ncol = n_ts, nrow = t_len)
  for(t in 1:t_len) {
    data[t, ] <- r_doublepois(n = n_ts, lambda_vec[t], sigma^2)
  }
  return(data)
}

calculate_size <- function(t_len, n_ts, alpha_vec, lambda_vec = lambda_vec,
                           sigma = sigma,
                           n_sim = 1000, sim_runs = 1000){
  
  #Constructing the set of intervals
  grid  <- construct_weekly_grid(t_len)

  #Constructing the set of the pairwise comparisons
  ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
  ijset <- ijset[ijset$i < ijset$j, ]

  # compute critical value
  quantiles <- compute_quantiles(t_len = t_len, n_ts = n_ts, 
                                 grid = grid, ijset = ijset,
                                 sim_runs = sim_runs)
  
  probs <- as.vector(quantiles$quant[1, ])
  quant <- as.vector(quantiles$quant[2, ])
  
  crit_val <- c()
  for (alpha in alpha_vec){
    if (sum(probs == (1 - alpha)) == 0)
      pos <- which.min(abs(probs - (1 - alpha)))
    if (sum(probs == (1 - alpha)) != 0)
      pos <- which.max(probs == (1 - alpha))
    crit_val <- c(crit_val, quant[pos])
  }
  
  
  # carry out multiscale test
  test_res <- matrix(NA, ncol = length(alpha_vec), nrow = n_sim)
  
  for(sim in 1:n_sim) {
    Y <- simulate_data(n_ts = n_ts, t_len = t_len, lambda_vec = lambda_vec, sigma = sigma)
    sigma_vec <- rep(0, n_ts)
    for (i in 1:n_ts){
      sigma_squared <- sum((Y[2:t_len, i] - Y[1:(t_len - 1), i]) ^ 2) / (2 * sum(Y[, i]))
      sigma_vec[i] <- sqrt(sigma_squared)
    }

    sigmahat <- sqrt(mean(sigma_vec * sigma_vec))
    
    result <- compute_statistics(data = Y, sigma = sigmahat, n_ts = n_ts,
                                 grid = grid, ijset = ijset)
    
    test_stat       <- max(result$stat, na.rm = TRUE)
    test_res[sim, ] <- as.numeric(test_stat > crit_val)
  }
  print(paste("Empirical size: ",  colSums(test_res) / n_sim, sep=""))
  return(colSums(test_res) / n_sim)
}


calculate_power <- function(t_len, n_ts, alpha_vec, lambda_vec_1 = lambda_vec_1,
                            lambda_vec_2 = lambda_vec_2,
                            sigma = sigma,
                            n_sim = 1000, sim_runs = 1000){
  
  #Constructing the set of intervals
  grid  <- construct_weekly_grid(t_len)

  #Constructing the set of the pairwise comparisons
  ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
  ijset <- ijset[ijset$i < ijset$j, ]

  # compute critical value
  quantiles <- compute_quantiles(t_len = t_len, n_ts = n_ts, 
                                 grid = grid, ijset = ijset,
                                 sim_runs = sim_runs)
  
  probs <- as.vector(quantiles$quant[1, ])
  quant <- as.vector(quantiles$quant[2, ])
  
  crit_val <- c()
  for (alpha in alpha_vec){
    if (sum(probs == (1 - alpha)) == 0)
      pos <- which.min(abs(probs - (1 - alpha)))
    if (sum(probs == (1 - alpha)) != 0)
      pos <- which.max(probs == (1 - alpha))
    crit_val <- c(crit_val, quant[pos])
  }
  

  # carry out multiscale test
  test_res       <- matrix(NA, ncol = length(alpha_vec), nrow = n_sim)

  for(sim in 1:n_sim) {
    #The first time series has the different mean function then the others
    Y1 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_1, sigma = sigma)
    Y2 <- simulate_data(n_ts = n_ts - 1, t_len = t_len, lambda_vec = lambda_vec_2, sigma = sigma)
    Y  <- cbind(Y1, Y2)
    
    sigma_vec <- rep(0, n_ts)
    for (i in 1:n_ts){
      sigma_squared <- sum((Y[2:t_len, i] - Y[1:(t_len - 1), i]) ^ 2) / (2 * sum(Y[, i]))
      sigma_vec[i] <- sqrt(sigma_squared)
    }
    
    sigmahat <- sqrt(mean(sigma_vec * sigma_vec))
    
    result <- compute_statistics(data = Y, sigma = sigmahat, n_ts = n_ts,
                                 grid = grid, ijset = ijset)
    
    test_stat_group1 <- max(result$stat_pairwise[1, ], na.rm = TRUE)
    test_stat_group2 <- max(result$stat_pairwise[-1, ], na.rm = TRUE)
    
    test_res[sim, ]       <- as.numeric((test_stat_group1 > crit_val) & (test_stat_group2 <= crit_val))
  }
  print(paste("Power: ",  colSums(test_res) / n_sim, sep=""))
  return(colSums(test_res) / n_sim)
}
