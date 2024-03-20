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
n_rep    <- 1000 #number of simulations for calculating size and power
sim_runs <- 1000 #number of simulations to calculate the Gaussian quantiles

different_T <- c(100, 250) #Different lengths of time series
alpha       <- 0.05
#different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels

a_hat <- 0.25 
sigma <- 0.25
q     <- 25 #Parameters for the estimation of long-run-variance
r     <- 10

numCores  = round(parallel::detectCores() * .70)

seed <- 135792468
set.seed(seed)

correct_specification <- c(rep(1, (floor(n_ts / 3))),
                           rep(2, (floor(2 * n_ts / 3) - floor(n_ts / 3))),
                           rep(3, n_ts - floor(2 * n_ts / 3)))


###################################
#This is how time series look like#
###################################


#Plotting the trend functions
t_len <- 100

m1 <- numeric(t_len)
m2 <- numeric(t_len)
m1 <- b_function((1:t_len)/t_len, 0.25, 0.25) - b_function((1:t_len)/t_len, 0.75, 0.25)
m2 <- 2 * b_function((1:t_len)/t_len, 0.75, 0.025) - 2 * b_function((1:t_len)/t_len, 0.25, 0.025)

pdf(paste0("output/revision/clustering_functions.pdf"),
    width = 12, height = 8, paper="special")
par(mfrow = c(2, 1))
par(mar = c(4, 3, 0.5, 0)) #Margins for each plot
par(oma = c(0.5, 0.5, 0.5, 0.2)) #Outer margins

errors <- arima.sim(model = list(ar = a_hat),
                    innov = rnorm(t_len, 0, sigma),
                    n = t_len)
plot(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
     y = m1, ylim = c(-4.5, 4.5),
     xlab = "", ylab = "", main = NULL,
     type = 'l', cex = 0.8)
lines(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
      y = m1 + errors, type = "l", col = "red")
#mtext(side = 1, text = paste0("b = ", b), line = 2.3, cex = 1)

errors <- arima.sim(model = list(ar = a_hat),
                    innov = rnorm(t_len, 0, sigma),
                    n = t_len)
plot(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
     y = m2, ylim = c(-4.5, 4.5),
     xlab = "", ylab = "", main = NULL,
     type = 'l', cex = 0.8)
lines(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
      y = m2 + errors, type = "l", col = "red")
#mtext(side = 1, text = paste0("b = ", b), line = 2.3, cex = 1)

dev.off()

###################################################
#Simulating the data and performing the clustering#
###################################################

for (t_len in different_T){
  #Constructing the grid
  u_grid <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
  h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len)
  
  m1 <- numeric(t_len)
  m2 <- numeric(t_len)
  m1 <- 0.5 * b_function((1:t_len)/t_len, 0.25, 0.25) - 0.5 * b_function((1:t_len)/t_len, 0.75, 0.25)
  m2 <- 2 * b_function((1:t_len)/t_len, 0.75, 0.025) - 2 * b_function((1:t_len)/t_len, 0.25, 0.025)
  
  cat("Calculating the distance measures for T = ", t_len,"\n")
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:n_rep, .combine = "cbind") %dopar% { 
    repl_clustering_revision(rep = val, t_len_ = t_len, n_ts_ = n_ts,
                             grid_ = grid, m1_ = m1, m2_ = m2, a_hat_ = a_hat,
                             sigma_ = sigma, q_ = q, r_ = r,
                             h_ = h_grid[1]) #Loop one-by-one using foreach
  } -> simulated_statistic
  stopCluster(cl)
  toc()

  cat("Calculating the Gaussian quantiles\n")
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:n_rep, .combine = "cbind") %dopar% { 
    repl_clustering_revision(rep = val, t_len_ = t_len, n_ts_ = n_ts,
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
  
  cat("Performing HAC\n")
  
  tic()
  if (sum(probs == (1 - alpha)) == 0)
    pos <- which.min(abs(probs - (1 - alpha)))
  if (sum(probs == (1 - alpha)) != 0)
    pos <- which.max(probs == (1 - alpha))    
  quant <- quants[pos]
  
  groups_mat <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_mat) <- paste0("rep_", 1:n_rep)
  rownames(groups_mat) <- paste0("ts_", 1:n_ts)
  
  groups_mat2 <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_mat2) <- paste0("rep_", 1:n_rep)
  rownames(groups_mat2) <- paste0("ts_", 1:n_ts)
  
  groups_benchmark_mat <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_benchmark_mat) <- paste0("rep_", 1:n_rep)
  rownames(groups_benchmark_mat) <- paste0("ts_", 1:n_ts)

  groups_benchmark2_mat <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_benchmark2_mat) <- paste0("rep_", 1:n_rep)
  rownames(groups_benchmark2_mat) <- paste0("ts_", 1:n_ts)
    
  number_of_groups_vec <- c()
  for (i in 1:n_rep){
    #Multiscale method with unknown number of clusters
    statistic_vector <- simulated_statistic[1:(nrow(simulated_statistic)/3), i]
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
    
    #Multiscale method with fixed number of clusters
    number_of_groups2  <- 3
    statistic_matrix2  <- matrix(statistic_vector, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix2  <- forceSymmetric(statistic_matrix2, uplo = "U")
    statistic_matrix2  <- as.dist(statistic_matrix2)
    clustering2        <- hclust(statistic_matrix2, method = "complete")
    groups2            <- cutree(clustering2, k = 3)

    groups_mat2[, i]      <- groups2

    #Benchmark method (L2 distance)
    statistic_vector_benchmark <- simulated_statistic[(nrow(simulated_statistic)/3 + 1):(2 * nrow(simulated_statistic) / 3), i]
    statistic_matrix_benchmark <- matrix(statistic_vector_benchmark, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark <- forceSymmetric(statistic_matrix_benchmark, uplo = "U")
    statistic_matrix_benchmark <- as.dist(statistic_matrix_benchmark)
    clustering_benchmark       <- hclust(statistic_matrix_benchmark, method = "complete")
    groups_benchmark           <- cutree(clustering_benchmark, k = 3)
    groups_benchmark_mat[, i]  <- groups_benchmark
    
    #Benchmark method 2 (max distance)
    statistic_vector_benchmark2 <- simulated_statistic[(2 * nrow(simulated_statistic)/3 + 1):nrow(simulated_statistic), i]
    statistic_matrix_benchmark2 <- matrix(statistic_vector_benchmark2, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark2 <- forceSymmetric(statistic_matrix_benchmark2, uplo = "U")
    statistic_matrix_benchmark2 <- as.dist(statistic_matrix_benchmark2)
    clustering_benchmark2       <- hclust(statistic_matrix_benchmark2, method = "complete")
    groups_benchmark2           <- cutree(clustering_benchmark2, k = 3)
    groups_benchmark2_mat[, i]  <- groups_benchmark2
    
  }
  clustering_results            <- rbind(number_of_groups_vec, groups_mat)
  clustering_results2           <- rbind(rep(3, n_rep), groups_mat2)
  clustering_results_benchmark  <- rbind(rep(3, n_rep), groups_benchmark_mat)
  clustering_results_benchmark2 <- rbind(rep(3, n_rep), groups_benchmark2_mat)
  
  filename = paste0("output/revision/misc/results_for_T_", t_len, ".RData")
  save(clustering_results, file = filename)
  filename2 = paste0("output/revision/misc/results_for_T_", t_len, "_3clusters.RData")
  save(clustering_results2, file = filename2)
  filename_benchmark = paste0("output/revision/misc/results_for_T_", t_len, "_benchmark.RData")
  save(clustering_results_benchmark, file = filename_benchmark)
  filename_benchmark2 = paste0("output/revision/misc/results_for_T_", t_len, "_benchmark2.RData")
  save(clustering_results_benchmark2, file = filename_benchmark2)
  toc()
}


######################################################
#Analysis of the clustering results for our procedure#
######################################################


cat("Analysis of the results for the multiscale method with unknown number of clusters\n")
correct_groups   <- c()
correct_structure <- c()

group_count <- list()
error_count <- list()

j <- 0

for (t_len in different_T){
  filename = paste0("output/revision/misc/results_for_T_", t_len, ".RData")
  load(file = filename)
  results <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                              results_matrix_ = clustering_results)

  j <- j + 1
  group_count[[j]] <- table(factor(clustering_results[1, ], levels = 1:5))
  error_count[[j]] <- table(factor(results$num_of_errors, levels = 0:8))

  correct_groups    <- c(correct_groups, results$correct_number_of_groups/n_rep)
  correct_structure <- c(correct_structure, results$correctly_specified_groups/n_rep)
}

cat("Analysis of the results for the multiscale method with known number of clusters (N = 3)\n")
correct_groups2   <- c()
correct_structure2 <- c()

group_count2 <- list()
error_count2 <- list()

j <- 0

for (t_len in different_T){
  filename = paste0("output/revision/misc/results_for_T_", t_len, "_3clusters.RData")
  load(file = filename)
  results <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                              results_matrix_ = clustering_results2)
  
  j <- j + 1
  group_count2[[j]] <- table(factor(clustering_results2[1, ], levels = 1:5))
  error_count2[[j]] <- table(factor(results$num_of_errors, levels = 0:8))
  
  correct_groups2    <- c(correct_groups2, results$correct_number_of_groups/n_rep)
  correct_structure2 <- c(correct_structure2, results$correctly_specified_groups/n_rep)
}


#################################################################
#Analysis of the clustering results for the benchmark procedures#
#################################################################
cat("Analysis of the results for the benchmark method with L2 distance\n")
correct_groups_benchmark   <- c()
correct_structure_benchmark <- c()

group_count_benchmark <- list()
error_count_benchmark <- list()

j <- 0

for (t_len in different_T){
  filename = paste0("output/revision/misc/results_for_T_", t_len, "_benchmark.RData")
  load(file = filename)

  results <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                              results_matrix_ = clustering_results_benchmark)

  j <- j + 1
  group_count_benchmark[[j]] <- table(factor(clustering_results_benchmark[1, ], levels = 1:5))
  error_count_benchmark[[j]] <- table(factor(results$num_of_errors, levels = 0:8))

  correct_groups_benchmark    <- c(correct_groups_benchmark, results$correct_number_of_groups/n_rep)
  correct_structure_benchmark <- c(correct_structure_benchmark, results$correctly_specified_groups/n_rep)
}

cat("Analysis of the results for the benchmark method with max distance\n")
correct_groups_benchmark2   <- c()
correct_structure_benchmark2 <- c()

group_count_benchmark2 <- list()
error_count_benchmark2 <- list()

j <- 0

for (t_len in different_T){
  filename = paste0("output/revision/misc/results_for_T_", t_len, "_benchmark2.RData")
  load(file = filename)
  
  results <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                              results_matrix_ = clustering_results_benchmark2)
  
  j <- j + 1
  group_count_benchmark2[[j]] <- table(factor(clustering_results_benchmark2[1, ], levels = 1:5))
  error_count_benchmark2[[j]] <- table(factor(results$num_of_errors, levels = 0:8))
  
  correct_groups_benchmark2    <- c(correct_groups_benchmark2, results$correct_number_of_groups/n_rep)
  correct_structure_benchmark2 <- c(correct_structure_benchmark2, results$correctly_specified_groups/n_rep)
}


#######################
#Output of the results#
#######################

labels_ <- c("(a)", "(b)", "(c)")

pdf(paste0("output/revision/histograms_groups.pdf"), width=8, height=2.9, paper="special")
par(mfrow = c(1,2))
par(mar = c(3, 2, 0.5, 1)) #Margins for each plot
par(oma = c(1.4, 1.5, 0.5, 0.2)) #Outer margins

for (j in 1:length(different_T)){
  t_len <- different_T[j]
  bp1 <- barplot(group_count[[j]], ylim = c(0, 1.05 * n_rep), xlab = "",
                 main = "", ylab = "", xaxt = 'n', space = 0)
  text(x = bp1, y = group_count[[j]], label = group_count[[j]], cex = 0.8, pos = 3)
  axis(1, at = bp1, labels = 1:5, tick=FALSE, line=-0.5, cex.axis=1)
  mtext(side=1, text= paste0(labels_[j], " T = ", t_len), line=3.4)
  title(xlab="number of groups", mgp=c(1.5,1,0), cex.lab=1)
}
dev.off()

pdf(paste0("output/revision/histograms_errors.pdf"), width=8, height=2.9, paper="special")
par(mfrow = c(1,2))
par(mar = c(3, 2, 0.5, 1)) #Margins for each plot
par(oma = c(1.4, 1.5, 0.5, 0.2)) #Outer margins

for (j in 1:length(different_T)){
  t_len <- different_T[j]
  bp2 <- barplot(error_count[[j]], ylim = c(0, 1.05 *  n_rep), xlab = "",
                 main = "", ylab = "", xaxt = 'n', space = 0)
  text(x = bp2, y = error_count[[j]], label = error_count[[j]], cex = 0.8, pos = 3)
  mtext(side=1, text= paste0(labels_[j], " T = ", t_len), line=3.4)
  axis(1, at = bp2, labels = 0:8, tick=FALSE, line=-0.5, cex.axis=1)
  title(xlab="number of errors", mgp=c(1.5,1,0), cex.lab=1)
}
dev.off()


# filename = paste0("output/tables/", n_ts, "_ts_correct_group_number.tex")
# rownames(correct_groups) <- different_T
# output_matrix(correct_groups, filename)
# 
# filename2 = paste0("output/tables/", n_ts, "_ts_correct_group_structure.tex")
# rownames(correct_structure) <- different_T
# output_matrix(correct_structure, filename2)
 


pdf(paste0("output/revision/histograms_groups_benchmark.pdf"), width=8, height=2.9, paper="special")
par(mfrow = c(1,2))
par(mar = c(3, 2, 0.5, 1)) #Margins for each plot
par(oma = c(1.4, 1.5, 0.5, 0.2)) #Outer margins

for (j in 1:length(different_T)){
  t_len <- different_T[j]
  bp1 <- barplot(group_count_benchmark[[j]], ylim = c(0, 1.05 * n_rep), xlab = "",
                 main = "", ylab = "", xaxt = 'n', space = 0)
  text(x = bp1, y = group_count_benchmark[[j]], label = group_count_benchmark[[j]], cex = 0.8, pos = 3)
  axis(1, at = bp1, labels = 1:5, tick=FALSE, line=-0.5, cex.axis=1)
  mtext(side=1, text= paste0(labels_[j], " T = ", t_len), line=3.4)
  title(xlab="number of groups for the benchmark", mgp=c(1.5,1,0), cex.lab=1)
}
dev.off()

pdf(paste0("output/revision/histograms_errors_benchmark.pdf"), width=8, height=2.9, paper="special")
par(mfrow = c(1,2))
par(mar = c(3, 2, 0.5, 1)) #Margins for each plot
par(oma = c(1.4, 1.5, 0.5, 0.2)) #Outer margins

for (j in 1:length(different_T)){
  t_len <- different_T[j]
  bp2 <- barplot(error_count_benchmark[[j]], ylim = c(0, 1.05 *  n_rep), xlab = "",
                 main = "", ylab = "", xaxt = 'n', space = 0)
  text(x = bp2, y = error_count_benchmark[[j]], label = error_count_benchmark[[j]], cex = 0.8, pos = 3)
  mtext(side=1, text= paste0(labels_[j], " T = ", t_len), line=3.4)
  axis(1, at = bp2, labels = 0:8, tick=FALSE, line=-0.5, cex.axis=1)
  title(xlab="number of errors for the benchmark", mgp=c(1.5,1,0), cex.lab=1)
}
dev.off()

# filename = paste0("output/tables/", n_ts, "_ts_correct_group_number_benchmark.tex")
# rownames(correct_groups_benchmark) <- different_T
# output_matrix(correct_groups_benchmark, filename)
# 
# filename2 = paste0("output/tables/", n_ts, "_ts_correct_group_structure_benchmark.tex")
# rownames(correct_structure_benchmark) <- different_T
# output_matrix(correct_structure_benchmark, filename)