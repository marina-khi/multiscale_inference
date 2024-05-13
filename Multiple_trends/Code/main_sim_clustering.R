rm(list=ls())

library(car)
library(ggplot2)
library(MSinference)
library(Matrix)
library(foreach)
library(parallel)
library(iterators)
library(doParallel)
library(tictoc)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("functions/functions.R")

##############################
#Defining necessary constants#
##############################

n_ts     <- 15 #number of different time series for simulation
n_rep    <- 5000 #number of simulations for calculating size and power
sim_runs <- 5000 #number of simulations to calculate the Gaussian quantiles

different_T     <- c(100, 250, 500) #Different lengths of time series
different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels

a     <- 0.25 
sigma <- 0.25
q     <- 25 #Parameters for the estimation of long-run-variance
r     <- 10

#For the covariate process
beta    <- c(1, 1, 1)
a_x_vec <- c(0.25, 0.25, 0.25) #VAR(1) coefficients
phi     <- 0.25                 #dependence between the innovations

#For the fixed effects
rho     <- 0.25 #covariance between the fixed effects

numCores = round(parallel::detectCores() * .80)

seed <- 135798642
set.seed(seed)

correct_specification <- c(rep(1, (floor(n_ts / 3))),
                           rep(2, (floor(2 * n_ts / 3) - floor(n_ts / 3))),
                           rep(3, n_ts - floor(2 * n_ts / 3)))


###################################
#This is how time series look like#
###################################

#Plotting the trend functions
t_len <- 250

m1 <- numeric(t_len)
m2 <- numeric(t_len)
m1 <- 0.35 * b_function((1:t_len)/t_len, 0.25, 0.25) - 0.35 * b_function((1:t_len)/t_len, 0.75, 0.25)
m2 <- b_function((1:t_len)/t_len, 0.75, 0.025) - b_function((1:t_len)/t_len, 0.25, 0.025)

pdf(paste0("output/revision/clustering_functions.pdf"),
    width = 12, height = 12, paper="special")
par(mfrow = c(3, 1))
par(mar = c(4, 3, 0.5, 0)) #Margins for each plot
par(oma = c(0.5, 0.5, 0.5, 0.2)) #Outer margins

errors <- arima.sim(model = list(ar = a),
                    innov = rnorm(t_len, 0, sigma),
                    n = t_len)
plot(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
     y = rep(0, t_len), ylim = c(-1.5, 1.5),
     xlab = "", ylab = "", main = NULL,
     type = 'l', cex = 0.8)
lines(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
      y = errors, type = "l", col = "red")
mtext(side = 1, text = expression(f[1](t/T) + epsilon[it]), line = 2.3, cex = 1)


errors <- arima.sim(model = list(ar = a),
                    innov = rnorm(t_len, 0, sigma),
                    n = t_len)
plot(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
     y = m1, ylim = c(-1.5, 1.5),
     xlab = "", ylab = "", main = NULL,
     type = 'l', cex = 0.8)
lines(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
      y = m1 + errors, type = "l", col = "red")
mtext(side = 1, text = expression(f[2](t/T) + epsilon[it]), line = 2.3, cex = 1)

errors <- arima.sim(model = list(ar = a),
                    innov = rnorm(t_len, 0, sigma),
                    n = t_len)
plot(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
     y = m2, ylim = c(-1.5, 1.5),
     xlab = "", ylab = "", main = NULL,
     type = 'l', cex = 0.8)
lines(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
      y = m2 + errors, type = "l", col = "red")
mtext(side = 1, text = expression(f[3](t/T) + epsilon[it]), line = 2.3, cex = 1)

dev.off()

###################################################
#Simulating the data and performing the clustering#
###################################################

for (t_len in different_T){
  #Constructing the grid
  u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
  h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
  
  m_matrix <- matrix(0, nrow = t_len, ncol = n_ts)
  for (i in (floor(n_ts / 3) + 1):(floor(2 * n_ts / 3))){
    m_matrix[, i] <- 0.35 * b_function((1:t_len)/t_len, 0.25, 0.25) - 0.35 * b_function((1:t_len)/t_len, 0.75, 0.25)
  }
  for (i in (floor(2 * n_ts / 3) + 1):n_ts){
    m_matrix[, i] <- b_function((1:t_len)/t_len, 0.75, 0.025) - b_function((1:t_len)/t_len, 0.25, 0.025)
  }
  
  cat("Calculating the Gaussian quantiles\n")
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:n_rep, .combine = "cbind") %dopar% { 
    repl_clustering(rep = val, t_len_ = t_len, n_ts_ = n_ts,
                    grid_ = grid, sigma_ = 1,
                    gaussian_sim = TRUE) #Loop one-by-one using foreach
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
  
  cat("Calculating the distance measures for T = ", t_len,"\n")
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:n_rep, .combine = "cbind") %dopar% { 
    repl_clustering(rep = val, t_len_ = t_len, n_ts_ = n_ts,
                    grid_ = grid,
                    m_matrix_ = m_matrix, 
                    a_ = a, sigma_ = sigma,
                    beta_ = beta, a_x_vec_ = a_x_vec,
                    phi_ = phi,  rho_ = rho,
                    q_ = q, r_ = r, 
                    gaussian_sim = FALSE, comparison = FALSE) #Loop one-by-one using foreach
  } -> simulated_statistic
  stopCluster(cl)
  toc()

  cat("Performing HAC\n")
  
  tic()
  for (alpha in different_alpha){
    if (sum(probs == (1 - alpha)) == 0)
      pos <- which.min(abs(probs - (1 - alpha)))
    if (sum(probs == (1 - alpha)) != 0)
      pos <- which.max(probs == (1 - alpha))    
    quant <- quants[pos]
    
    groups_mat <- matrix(NA, ncol = n_rep, nrow = n_ts)
    colnames(groups_mat) <- paste0("rep_", 1:n_rep)
    rownames(groups_mat) <- paste0("ts_", 1:n_ts)
    
    number_of_groups_vec <- c()
    for (i in 1:n_rep){
      #Multiscale method with unknown number of clusters
      statistic_vector <- simulated_statistic[, i]
      statistic_value  <- max(statistic_vector)
      if (statistic_value > quant) {
        statistic_matrix  <- matrix(statistic_vector, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
        statistic_matrix  <- forceSymmetric(statistic_matrix, uplo = "U")
        statistic_matrix  <- as.dist(statistic_matrix)
        clustering        <- hclust(statistic_matrix, method = "complete")
        groups            <- cutree(clustering, h = quant)
        number_of_groups  <- max(groups)
      } else {
        number_of_groups <- 1
        groups           <- rep(1, n_ts)
      }
      groups_mat[, i]      <- groups
      number_of_groups_vec <- c(number_of_groups_vec, number_of_groups)
    }
    clustering_results  <- rbind(number_of_groups_vec, groups_mat)

    filename = paste0("output/revision/misc/results_for_T_", t_len, "_and_alpha_", alpha * 100, ".RData")
    save(clustering_results, file = filename)
  }
  toc()
}


######################################################
#Analysis of the clustering results for our procedure#
######################################################

cat("Analysis of the results for the multiscale method with unknown number of clusters\n")

correct_groups_vec <- c()
correct_structure_vec <- c()

group_count <- list()
error_count <- list()

j <- 0
for (t_len in different_T){
  for (alpha in different_alpha){
    filename = paste0("output/revision/misc/results_for_T_", t_len, "_and_alpha_", alpha * 100, ".RData")
    load(file = filename)
    results <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                                results_matrix_ = clustering_results,
                                correct_specification_ = correct_specification)
    if (alpha == 0.05){
      j <- j + 1
      group_count[[j]] <- table(factor(clustering_results[1, ], levels = 1:5))
      error_count[[j]] <- table(factor(results$num_of_errors, levels = 0:8))
    }
    correct_groups_vec    <- c(correct_groups_vec, results$correct_number_of_groups/n_rep)
    correct_structure_vec <- c(correct_structure_vec, results$correctly_specified_groups/n_rep)
  }
}

correct_groups    <- matrix(correct_groups_vec, ncol = length(different_alpha), byrow = TRUE)
correct_structure <- matrix(correct_structure_vec, ncol = length(different_alpha), byrow = TRUE)


#######################
#Output of the results#
#######################

produce_hist_plots(file_extension_ = "", different_T_ = different_T,
                   n_rep_ = n_rep, group_count_ = group_count,
                   error_count_ = error_count)

filename = paste0("output/revision/", n_ts, "_ts_correct_group_number.tex")
rownames(correct_groups) <- different_T
output_matrix(matrix_ = correct_groups, filename_ = filename, numcols_ = 4)

filename2 = paste0("output/revision/", n_ts, "_ts_correct_group_structure.tex")
rownames(correct_structure) <- different_T
output_matrix(matrix_ = correct_structure, filename_ = filename2, numcols_ = 4)