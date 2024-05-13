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

different_T <- c(100, 250, 500) #Different lengths of time series
alpha       <- 0.05 #Different confidence levels

a     <- 0.25 
sigma <- 0.25

numCores = round(parallel::detectCores() * .80)

seed <- 13577531

correct_specification <- c(rep(1, (floor(n_ts / 3))),
                           rep(2, (floor(2 * n_ts / 3) - floor(n_ts / 3))),
                           rep(3, n_ts - floor(2 * n_ts / 3)))

###################################################
#Simulating the data and performing the clustering#
###################################################

for (t_len in different_T){
  set.seed(seed)
  #Constructing the grid
  u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
  h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
  
  m_matrix <- matrix(0, nrow = t_len, ncol = n_ts)
  for (i in (floor(n_ts / 3) + 1):(floor(2 * n_ts / 3))){
    m_matrix[, i] <- 0.3 * b_function((1:t_len)/t_len, 0.25, 0.25) - 0.3 * b_function((1:t_len)/t_len, 0.75, 0.25)
  }
  for (i in (floor(2 * n_ts / 3) + 1):n_ts){
    m_matrix[, i] <- b_function((1:t_len)/t_len, 0.75, 0.025) - b_function((1:t_len)/t_len, 0.25, 0.025)
  }
  
  cat("Calculating the distance measures for T = ", t_len,"\n")
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:n_rep, .combine = "cbind") %dopar% { 
    repl_clustering(rep = val, t_len_ = t_len, n_ts_ = n_ts,
                    grid_ = grid, m_matrix_ = m_matrix,
                    a_ = a, sigma_ = sigma, beta_ = NULL, q_ = 25, r_ = 10, 
                    gaussian_sim = FALSE, comparison = TRUE) #Loop one-by-one using foreach
  } -> simulated_statistic
  stopCluster(cl)
  toc()
  
  cat("Performing HAC\n")
  
  tic()
  groups_mat <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_mat) <- paste0("rep_", 1:n_rep)
  rownames(groups_mat) <- paste0("ts_", 1:n_ts)
  
  groups_benchmark_mat <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_benchmark_mat) <- paste0("rep_", 1:n_rep)
  rownames(groups_benchmark_mat) <- paste0("ts_", 1:n_ts)
  
  groups_benchmark_mat2 <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_benchmark_mat2) <- paste0("rep_", 1:n_rep)
  rownames(groups_benchmark_mat2) <- paste0("ts_", 1:n_ts)

  groups_benchmark_mat3 <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_benchmark_mat3) <- paste0("rep_", 1:n_rep)
  rownames(groups_benchmark_mat3) <- paste0("ts_", 1:n_ts)
  
  groups_benchmark_mat4 <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_benchmark_mat4) <- paste0("rep_", 1:n_rep)
  rownames(groups_benchmark_mat4) <- paste0("ts_", 1:n_ts)

  groups_benchmark_mat5 <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_benchmark_mat5) <- paste0("rep_", 1:n_rep)
  rownames(groups_benchmark_mat5) <- paste0("ts_", 1:n_ts)
  
  groups_benchmark_mat6 <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_benchmark_mat6) <- paste0("rep_", 1:n_rep)
  rownames(groups_benchmark_mat6) <- paste0("ts_", 1:n_ts)
  
  for (i in 1:n_rep){
    #Multiscale method with estimated lrv
    statistic_vector <- simulated_statistic[1:(nrow(simulated_statistic)/7), i]
    statistic_matrix <- matrix(statistic_vector, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix <- forceSymmetric(statistic_matrix, uplo = "U")
    statistic_matrix <- as.dist(statistic_matrix)
    clustering       <- hclust(statistic_matrix, method = "complete")
    groups_mat[, i]  <- cutree(clustering, k = 3)

    #Benchmark method (h1)
    statistic_vector_benchmark <- simulated_statistic[(nrow(simulated_statistic)/7 + 1):(2 * nrow(simulated_statistic) / 7), i]
    statistic_matrix_benchmark <- matrix(statistic_vector_benchmark, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark <- forceSymmetric(statistic_matrix_benchmark, uplo = "U")
    statistic_matrix_benchmark <- as.dist(statistic_matrix_benchmark)
    clustering_benchmark       <- hclust(statistic_matrix_benchmark, method = "complete")
    groups_benchmark_mat[, i]  <- cutree(clustering_benchmark, k = 3)
    
    #Benchmark method (h2)
    statistic_vector_benchmark2 <- simulated_statistic[(2 * nrow(simulated_statistic)/7 + 1):(3 * nrow(simulated_statistic) / 7), i]
    statistic_matrix_benchmark2 <- matrix(statistic_vector_benchmark2, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark2 <- forceSymmetric(statistic_matrix_benchmark2, uplo = "U")
    statistic_matrix_benchmark2 <- as.dist(statistic_matrix_benchmark2)
    clustering_benchmark2       <- hclust(statistic_matrix_benchmark2, method = "complete")
    groups_benchmark_mat2[, i]  <- cutree(clustering_benchmark2, k = 3)

    #Benchmark method (h3)
    statistic_vector_benchmark3 <- simulated_statistic[(3 * nrow(simulated_statistic)/7 + 1):(4 * nrow(simulated_statistic) / 7), i]
    statistic_matrix_benchmark3 <- matrix(statistic_vector_benchmark3, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark3 <- forceSymmetric(statistic_matrix_benchmark3, uplo = "U")
    statistic_matrix_benchmark3 <- as.dist(statistic_matrix_benchmark3)
    clustering_benchmark3       <- hclust(statistic_matrix_benchmark3, method = "complete")
    groups_benchmark_mat3[, i]  <- cutree(clustering_benchmark3, k = 3)

    #Benchmark method (h4)
    statistic_vector_benchmark4 <- simulated_statistic[(4 * nrow(simulated_statistic)/7 + 1):(5 * nrow(simulated_statistic) / 7), i]
    statistic_matrix_benchmark4 <- matrix(statistic_vector_benchmark4, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark4 <- forceSymmetric(statistic_matrix_benchmark4, uplo = "U")
    statistic_matrix_benchmark4 <- as.dist(statistic_matrix_benchmark4)
    clustering_benchmark4       <- hclust(statistic_matrix_benchmark4, method = "complete")
    groups_benchmark_mat4[, i]  <- cutree(clustering_benchmark4, k = 3)

    #Benchmark method (h5)
    statistic_vector_benchmark5 <- simulated_statistic[(5 * nrow(simulated_statistic)/7 + 1):(6 * nrow(simulated_statistic) / 7), i]
    statistic_matrix_benchmark5 <- matrix(statistic_vector_benchmark5, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark5 <- forceSymmetric(statistic_matrix_benchmark5, uplo = "U")
    statistic_matrix_benchmark5 <- as.dist(statistic_matrix_benchmark5)
    clustering_benchmark5       <- hclust(statistic_matrix_benchmark5, method = "complete")
    groups_benchmark_mat5[, i]  <- cutree(clustering_benchmark5, k = 3)
    
    #Benchmark method (h6)
    statistic_vector_benchmark6 <- simulated_statistic[(6 * nrow(simulated_statistic)/7 + 1):nrow(simulated_statistic), i]
    statistic_matrix_benchmark6 <- matrix(statistic_vector_benchmark6, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark6 <- forceSymmetric(statistic_matrix_benchmark6, uplo = "U")
    statistic_matrix_benchmark6 <- as.dist(statistic_matrix_benchmark6)
    clustering_benchmark6       <- hclust(statistic_matrix_benchmark6, method = "complete")
    groups_benchmark_mat6[, i]  <- cutree(clustering_benchmark6, k = 3)
  }
  
  clustering_results            <- rbind(rep(3, n_rep), groups_mat)
  clustering_results_benchmark  <- rbind(rep(3, n_rep), groups_benchmark_mat)
  clustering_results_benchmark2 <- rbind(rep(3, n_rep), groups_benchmark_mat2)
  clustering_results_benchmark3 <- rbind(rep(3, n_rep), groups_benchmark_mat3)
  clustering_results_benchmark4 <- rbind(rep(3, n_rep), groups_benchmark_mat4)
  clustering_results_benchmark5 <- rbind(rep(3, n_rep), groups_benchmark_mat5)
  clustering_results_benchmark6 <- rbind(rep(3, n_rep), groups_benchmark_mat6)
  
  filename = paste0("output/revision/misc/results_for_T_", t_len, "_comparison.RData")
  save(clustering_results, clustering_results_benchmark,
       clustering_results_benchmark2, clustering_results_benchmark3,
       clustering_results_benchmark4, clustering_results_benchmark5,
       clustering_results_benchmark6, file = filename)
}


######################################################
#Analysis of the clustering results for our procedure#
######################################################

cat("Analysis of the results for the multiscale method with unknown number of clusters\n")
correct_structure <- c()
error_count       <- list()

correct_structure_benchmark <- c()
error_count_benchmark       <- list()

correct_structure_benchmark2 <- c()
error_count_benchmark2       <- list()

correct_structure_benchmark3 <- c()
error_count_benchmark3       <- list()

correct_structure_benchmark4 <- c()
error_count_benchmark4       <- list()

correct_structure_benchmark5 <- c()
error_count_benchmark5       <- list()

correct_structure_benchmark6 <- c()
error_count_benchmark6       <- list()

j <- 1

for (t_len in different_T){
  filename = paste0("output/revision/misc/results_for_T_", t_len, "_comparison.RData")
  load(file = filename)
  cat("Analysis of the results for the multiscale method with the estimated lrv\n")
  results <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                              results_matrix_ = clustering_results,
                              correct_specification_ = correct_specification)
  
  error_count[[j]]  <- table(factor(results$num_of_errors, levels = 0:8))
  correct_structure <- c(correct_structure, results$correctly_specified_groups/n_rep)

  cat("Analysis of the results for the benchmark method with bandwidth h1\n")
  results_benchmark <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                                        results_matrix_ = clustering_results_benchmark,
                                        correct_specification_ = correct_specification)
   
  error_count_benchmark[[j]]  <- table(factor(results_benchmark$num_of_errors, levels = 0:8))
  correct_structure_benchmark <- c(correct_structure_benchmark, results_benchmark$correctly_specified_groups/n_rep)

  cat("Analysis of the results for the benchmark method with bandwidth h2\n")
  results_benchmark2 <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                                         results_matrix_ = clustering_results_benchmark2,
                                         correct_specification_ = correct_specification)
  
  error_count_benchmark2[[j]]  <- table(factor(results_benchmark2$num_of_errors, levels = 0:8))
  correct_structure_benchmark2 <- c(correct_structure_benchmark2, results_benchmark2$correctly_specified_groups/n_rep)

  cat("Analysis of the results for the benchmark method with bandwidth h3\n")
  results_benchmark3 <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                                         results_matrix_ = clustering_results_benchmark3,
                                         correct_specification_ = correct_specification)
  
  error_count_benchmark3[[j]]  <- table(factor(results_benchmark3$num_of_errors, levels = 0:8))
  correct_structure_benchmark3 <- c(correct_structure_benchmark3, results_benchmark3$correctly_specified_groups/n_rep)
  
  cat("Analysis of the results for the benchmark method with bandwidth h4\n")
  results_benchmark4 <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                                         results_matrix_ = clustering_results_benchmark4,
                                         correct_specification_ = correct_specification)
  
  error_count_benchmark4[[j]]  <- table(factor(results_benchmark4$num_of_errors, levels = 0:8))
  correct_structure_benchmark4 <- c(correct_structure_benchmark4, results_benchmark4$correctly_specified_groups/n_rep)
  
  cat("Analysis of the results for the benchmark method with bandwidth h5\n")
  results_benchmark5 <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                                         results_matrix_ = clustering_results_benchmark5,
                                         correct_specification_ = correct_specification)
  
  error_count_benchmark5[[j]]  <- table(factor(results_benchmark5$num_of_errors, levels = 0:8))
  correct_structure_benchmark5 <- c(correct_structure_benchmark5, results_benchmark5$correctly_specified_groups/n_rep)
  
  cat("Analysis of the results for the benchmark method with bandwidth h6\n")
  results_benchmark6 <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                                         results_matrix_ = clustering_results_benchmark6,
                                         correct_specification_ = correct_specification)
  
  error_count_benchmark6[[j]]  <- table(factor(results_benchmark6$num_of_errors, levels = 0:8))
  correct_structure_benchmark6 <- c(correct_structure_benchmark6, results_benchmark6$correctly_specified_groups/n_rep)
  j <- j + 1
}


#######################
#Output of the results#
#######################

tmp <- matrix(NA, nrow = length(different_T), ncol = 7)
tmp[, 1] <- as.vector(correct_structure)
tmp[, 2] <- as.vector(correct_structure_benchmark)
tmp[, 3] <- as.vector(correct_structure_benchmark2)
tmp[, 4] <- as.vector(correct_structure_benchmark3)
tmp[, 5] <- as.vector(correct_structure_benchmark4)
tmp[, 6] <- as.vector(correct_structure_benchmark5)
tmp[, 7] <- as.vector(correct_structure_benchmark6)

row.names(tmp) <- paste0("$T = ", different_T, "$")

filename = paste0("output/revision/clustering_comparison.tex")
output_matrix(matrix_ = tmp, filename = filename, numcols = 8)
line <- paste0("%This simulation was done for the following values of the parameters: n_ts = ",
               n_ts, ", with ", n_rep,
               " simulations for calculating size and power. Furthermore, for the error process we have a = ",
               a, " and sigma = ", sigma,
               ". There are no fixed effects. The grid is standard.")
write(line, file = filename, append = TRUE)   
