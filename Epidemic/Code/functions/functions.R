#Local Linear estimator using the epanechnikov kernel. 
nadaraya_watson_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    result = sum((abs((grid_points - u) / bw) <= 1) * data_p)
    norm = sum((abs((grid_points - u) / bw) <= 1))
    #for (i in 1:T_size){
    #  result = result + epanechnikov_kernel((grid_p[i] - u) / bw) * data_p[i]
    #  norm = norm + epanechnikov_kernel((grid_p[i] - u) / bw)
    #}
    return(result/norm)
  }
}


produce_plots <- function (results, l, data_i, data_j, smoothed_i, smoothed_j,
                           gov_resp_i, gov_resp_j, lagged_gov_resp_i, lagged_gov_resp_j,
                           country_i, country_j){
  Tlen <- length(data_i)
  gset <- results$gset_with_vals[[l]]

  layout(matrix(c(1, 2, 3, 4),ncol=1), widths=c(2.2, 2.2, 2.2, 2.2),
         heights=c(1.5, 1.5, 1.5, 1.8), TRUE)
  #Setting the layout of the graphs

  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 2, 0.2)) #Outer margins

  plot(data_i, ylim=c(min(data_i, data_j), max(data_i, data_j)), type="l",
      col="blue", ylab="", xlab="", mgp=c(1, 0.5, 0))
  lines(data_j, col="red")
  title(main = "(a) observed new cases per day", font.main = 1, line = 0.5)
  legend("topright", inset = 0.02, legend=c(country_i, country_j),
         col = c("blue", "red"), lty = 1, cex = 0.95, ncol = 1)

  par(mar = c(0.5, 0.5, 3, 0)) #Margins for each plot

  plot(smoothed_i, ylim=c(min(data_i, data_j), max(data_i, data_j)), type="l",
     col="blue", ylab="", xlab = "", mgp=c(1,0.5,0))
  title(main = "(b) smoothed curve from (a)", font.main = 1, line = 0.5)
  lines(smoothed_j, col="red")

  plot(gov_resp_i, ylim=c(0, 100), type="l",
       col="blue", ylab="", xlab = "", mgp=c(1, 0.5, 0))
  title(main = "(c) government response index", font.main = 1, line = 0.5)
  lines(gov_resp_j, col="red")
  #lines(lagged_gov_resp_i, col="blue", lty = "dashed", lwd = 3)
  #lines(lagged_gov_resp_j, col="red", lty = "dashed", lwd = 3)

  par(mar = c(2.7, 0.5, 3, 0)) #Margins for each plot

  a_t_set <- subset(gset, test == 1, select = c(u, h))
  if (nrow(a_t_set) > 0){
    p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h) * Tlen + 0.5,
                          'endpoint' = (a_t_set$u + a_t_set$h) * Tlen - 0.5, 'values' = 0)
    p_t_set$values <- (1:nrow(p_t_set))/nrow(p_t_set)
    
    #Produce minimal intervals
    p_t_set2  <- compute_minimal_intervals(p_t_set)

    plot(NA, xlim=c(0, Tlen),  ylim = c(0, 1 + 1 / nrow(p_t_set)), xlab="", mgp=c(2, 0.5, 0), yaxt = "n")
    title(main = "(d) minimal intervals produced by our test", font.main = 1, line = 0.5)
    title(xlab = "days since the hundredth case", line = 1.7, cex.lab = 0.9)
    segments(p_t_set2$startpoint, p_t_set2$values, p_t_set2$endpoint, p_t_set2$values, lwd = 2)
    segments(p_t_set$startpoint, p_t_set$values, p_t_set$endpoint, p_t_set$values, col = "gray")
  } else {
    #If there are no intervals where the test rejects, we produce empty plots

    plot(NA, xlim=c(0, Tlen),  ylim = c(0, 1), xlab="", ylab = "", mgp=c(2,0.5,0), yaxt = "n")
    title(main = "(d) minimal intervals produced by our test", font.main = 1, line = 0.5)
    title(xlab = "days since the hundredth case", line = 1.7, cex.lab = 0.9)
  }
  mtext(paste0("Comparison of ", country_i, " and ", country_j), side = 3, line = 0, outer = TRUE, font = 1, cex = 1.2)
  print.xtable(xtable(p_t_set2, digits = c(3), align = paste(replicate(4, "c"), collapse = "")),
               type="latex", file=paste0("plots/", country_i, "_vs_", country_j, ".tex"),
               include.colnames = FALSE)
  
}

# functions to simulate data

lambda_fct <- function(u, c) {
  return (5000 * exp(-(10 * u - 3) ^ 2 / 2) + c)
}

r_doublepois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta - 1))
}

simulate_data <- function(n_ts, t_len, lambda_vec, sigma) {
  data <- matrix(0, ncol = n_ts, nrow = t_len)
  #for(i in 1:n_ts) {
  #  data[, i]  <- rdoublepois(n = t_len, m = lambda_vec, s = rep(1 / sigma^2, t_len))
  #}
  for(t in 1:t_len) {
    data[t, ] <- r_doublepois(n = n_ts, lambda_vec[t], sigma^2)
  }
  return(data)
}

simulate_data_iid <- function(n_ts, t_len, lambda_vec, sigma) {
  data_ <- matrix(0, ncol = n_ts, nrow = t_len)
  for(t in 1:t_len) {
    data_[t, ]  <- lambda_vec[t] + sigma * sqrt(lambda_vec[t]) * rnorm(n_ts)
  }
  data_[data_ < 0] <- 0
  return(data_)
}


calculate_size <- function(t_len, n_ts, alpha_vec, lambda_vec = lambda_vec,
                           sigma = sigma,
                           n_sim = 1000, sim_runs = 1000, iid = FALSE){
  
  #Constructing the set of intervals
  grid                   <- construct_weekly_grid(t_len)
  gset                   <- grid$gset
  gset_cpp               <- as.matrix(gset)
  gset_cpp               <- as.vector(gset_cpp)
  storage.mode(gset_cpp) <- "double"
  
  #Constructing the set of the pairwise comparisons
  ijset                   <- expand.grid(i = 1:n_ts, j = 1:n_ts)
  ijset                   <- ijset[ijset$i < ijset$j, ]
  ijset_cpp               <- as.matrix(ijset)
  ijset_cpp               <- as.vector(ijset_cpp)
  storage.mode(ijset_cpp) <- "integer"
  
  # compute critical value
  quantiles <- compute_quantiles(t_len = t_len, grid = grid, n_ts = n_ts, 
                                 sim_runs = sim_runs, epidem = TRUE)
  
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
    if (iid){
      Y <- simulate_data_iid(n_ts = n_ts, t_len = t_len, lambda_vec = lambda_vec, sigma = sigma)
    } else {
      Y <- simulate_data(n_ts = n_ts, t_len = t_len, lambda_vec = lambda_vec, sigma = sigma)
    }
    sigma_vec <- rep(0, n_ts)
    for (i in 1:n_ts){
      sigma_squared <- sum((Y[2:t_len, i] - Y[1:(t_len - 1), i]) ^ 2) / (2 * sum(Y[, i]))
      sigma_vec[i] <- sqrt(sigma_squared)
    }
    
    result <- compute_multiple_statistics(t_len = t_len, n_ts = n_ts, data = Y,
                                          gset = gset_cpp, ijset = ijset_cpp,
                                          sigma_vec = sigma_vec, epidem = TRUE)
    
    test_stat <- max(result$stat, na.rm = TRUE)
    test_res[sim, ] <- as.numeric(test_stat > crit_val)
  }
  print(paste("Empirical size: ",  colSums(test_res) / n_sim, sep=""))
  return(colSums(test_res) / n_sim)
}
