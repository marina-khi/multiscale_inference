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

different_T <- c(100, 250, 500) #Different lengths of time series
alpha       <- 0.05 #Different confidence levels

a     <- 0.25 
sigma <- 0.25

numCores = round(parallel::detectCores() * .80)

seed <- 135792468

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
  grid   <- construct_grid(t = t_len)
  

  m_matrix <- matrix(0, nrow = t_len, ncol = n_ts)
  for (i in (floor(n_ts / 3) + 1):(floor(2 * n_ts / 3))){
    m_matrix[, i] <- 0.35 * b_function((1:t_len)/t_len, 0.25, 0.25) - 0.35 * b_function((1:t_len)/t_len, 0.75, 0.25)
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
                    gaussian_sim = FALSE, comparison = TRUE)#Loop one-by-one using foreach
  } -> simulated_statistic
  stopCluster(cl)
  toc()
  
  cat("Performing HAC\n")
  
  tic()
  groups_mat <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_mat) <- paste0("rep_", 1:n_rep)
  rownames(groups_mat) <- paste0("ts_", 1:n_ts)
  
  groups_mat2 <- matrix(NA, ncol = n_rep, nrow = n_ts)
  colnames(groups_mat2) <- paste0("rep_", 1:n_rep)
  rownames(groups_mat2) <- paste0("ts_", 1:n_ts)

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
    statistic_vector <- simulated_statistic[1:(nrow(simulated_statistic)/8), i]
    statistic_matrix <- matrix(statistic_vector, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix <- forceSymmetric(statistic_matrix, uplo = "U")
    statistic_matrix <- as.dist(statistic_matrix)
    clustering       <- hclust(statistic_matrix, method = "complete")
    groups           <- cutree(clustering, k = 3)
    groups_mat[, i]  <- groups

    #Multiscale method with true lrv
    statistic_vector2 <- simulated_statistic[(nrow(simulated_statistic)/8 + 1):(2 * nrow(simulated_statistic)/8), i]
    statistic_matrix2 <- matrix(statistic_vector2, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix2 <- forceSymmetric(statistic_matrix2, uplo = "U")
    statistic_matrix2 <- as.dist(statistic_matrix2)
    clustering2       <- hclust(statistic_matrix2, method = "complete")
    groups2           <- cutree(clustering2, k = 3)
    groups_mat2[, i]  <- groups2
    
    #Benchmark method (h1)
    statistic_vector_benchmark <- simulated_statistic[(2 * nrow(simulated_statistic)/8 + 1):(3 * nrow(simulated_statistic) / 8), i]
    statistic_matrix_benchmark <- matrix(statistic_vector_benchmark, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark <- forceSymmetric(statistic_matrix_benchmark, uplo = "U")
    statistic_matrix_benchmark <- as.dist(statistic_matrix_benchmark)
    clustering_benchmark       <- hclust(statistic_matrix_benchmark, method = "complete")
    groups_benchmark           <- cutree(clustering_benchmark, k = 3)
    groups_benchmark_mat[, i]  <- groups_benchmark
    
    #Benchmark method (h2)
    statistic_vector_benchmark2 <- simulated_statistic[(3 * nrow(simulated_statistic)/8 + 1):(4 * nrow(simulated_statistic) / 8), i]
    statistic_matrix_benchmark2 <- matrix(statistic_vector_benchmark2, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark2 <- forceSymmetric(statistic_matrix_benchmark2, uplo = "U")
    statistic_matrix_benchmark2 <- as.dist(statistic_matrix_benchmark2)
    clustering_benchmark2       <- hclust(statistic_matrix_benchmark2, method = "complete")
    groups_benchmark2           <- cutree(clustering_benchmark2, k = 3)
    groups_benchmark_mat2[, i]  <- groups_benchmark2

    #Benchmark method (h3)
    statistic_vector_benchmark3 <- simulated_statistic[(4 * nrow(simulated_statistic)/8 + 1):(5 * nrow(simulated_statistic) / 8), i]
    statistic_matrix_benchmark3 <- matrix(statistic_vector_benchmark3, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark3 <- forceSymmetric(statistic_matrix_benchmark3, uplo = "U")
    statistic_matrix_benchmark3 <- as.dist(statistic_matrix_benchmark3)
    clustering_benchmark3       <- hclust(statistic_matrix_benchmark3, method = "complete")
    groups_benchmark3           <- cutree(clustering_benchmark3, k = 3)
    groups_benchmark_mat3[, i]  <- groups_benchmark3

    #Benchmark method (h4)
    statistic_vector_benchmark4 <- simulated_statistic[(5 * nrow(simulated_statistic)/8 + 1):(6 * nrow(simulated_statistic) / 8), i]
    statistic_matrix_benchmark4 <- matrix(statistic_vector_benchmark4, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark4 <- forceSymmetric(statistic_matrix_benchmark4, uplo = "U")
    statistic_matrix_benchmark4 <- as.dist(statistic_matrix_benchmark4)
    clustering_benchmark4       <- hclust(statistic_matrix_benchmark4, method = "complete")
    groups_benchmark4           <- cutree(clustering_benchmark4, k = 3)
    groups_benchmark_mat4[, i]  <- groups_benchmark4

    #Benchmark method (h5)
    statistic_vector_benchmark5 <- simulated_statistic[(6 * nrow(simulated_statistic)/8 + 1):(7 * nrow(simulated_statistic) / 8), i]
    statistic_matrix_benchmark5 <- matrix(statistic_vector_benchmark5, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark5 <- forceSymmetric(statistic_matrix_benchmark5, uplo = "U")
    statistic_matrix_benchmark5 <- as.dist(statistic_matrix_benchmark5)
    clustering_benchmark5       <- hclust(statistic_matrix_benchmark5, method = "complete")
    groups_benchmark5           <- cutree(clustering_benchmark5, k = 3)
    groups_benchmark_mat5[, i]  <- groups_benchmark5
    
    #Benchmark method (h6)
    statistic_vector_benchmark6 <- simulated_statistic[(7 * nrow(simulated_statistic)/8 + 1):nrow(simulated_statistic), i]
    statistic_matrix_benchmark6 <- matrix(statistic_vector_benchmark6, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
    statistic_matrix_benchmark6 <- forceSymmetric(statistic_matrix_benchmark6, uplo = "U")
    statistic_matrix_benchmark6 <- as.dist(statistic_matrix_benchmark6)
    clustering_benchmark6       <- hclust(statistic_matrix_benchmark6, method = "complete")
    groups_benchmark6           <- cutree(clustering_benchmark6, k = 3)
    groups_benchmark_mat6[, i]  <- groups_benchmark6
  }
  
  clustering_results            <- rbind(rep(3, n_rep), groups_mat)
  clustering_results2            <- rbind(rep(3, n_rep), groups_mat2)
  clustering_results_benchmark  <- rbind(rep(3, n_rep), groups_benchmark_mat)
  clustering_results_benchmark2 <- rbind(rep(3, n_rep), groups_benchmark_mat2)
  clustering_results_benchmark3 <- rbind(rep(3, n_rep), groups_benchmark_mat3)
  clustering_results_benchmark4 <- rbind(rep(3, n_rep), groups_benchmark_mat4)
  clustering_results_benchmark5 <- rbind(rep(3, n_rep), groups_benchmark_mat5)
  clustering_results_benchmark6 <- rbind(rep(3, n_rep), groups_benchmark_mat6)
  
  filename = paste0("output/revision/misc/results_for_T_", t_len, "_comparison.RData")
  save(clustering_results, clustering_results2, clustering_results_benchmark,
       clustering_results_benchmark2, clustering_results_benchmark3,
       clustering_results_benchmark4, clustering_results_benchmark5,
       clustering_results_benchmark6, file = filename)
}


######################################################
#Analysis of the clustering results for our procedure#
######################################################

cat("Analysis of the results for the multiscale method with unknown number of clusters\n")
correct_groups   <- c()
correct_structure <- c()

group_count <- list()
error_count <- list()

correct_groups2    <- c()
correct_structure2 <- c()

group_count2 <- list()
error_count2 <- list()

correct_groups_benchmark   <- c()
correct_structure_benchmark <- c()

group_count_benchmark <- list()
error_count_benchmark <- list()

correct_groups_benchmark2   <- c()
correct_structure_benchmark2 <- c()

group_count_benchmark2 <- list()
error_count_benchmark2 <- list()

j <- 1

for (t_len in different_T){
  filename = paste0("output/revision/misc/results_for_T_", t_len, "_comparison.RData")
  load(file = filename)
  cat("Analysis of the results for the multiscale method with known number of clusters\n")
  results <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                              results_matrix_ = clustering_results,
                              correct_specification_ = correct_specification)
  
  group_count[[j]] <- table(factor(clustering_results[1, ], levels = 1:5))
  error_count[[j]] <- table(factor(results$num_of_errors, levels = 0:8))
  
  correct_groups    <- c(correct_groups, results$correct_number_of_groups/n_rep)
  correct_structure <- c(correct_structure, results$correctly_specified_groups/n_rep)
  
  cat("Analysis of the results for the benchmark method with L2 distance and small bandwidth\n")
  results_benchmark <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                                        results_matrix_ = clustering_results_benchmark,
                                        correct_specification_ = correct_specification)
  
  group_count_benchmark[[j]] <- table(factor(clustering_results_benchmark[1, ], levels = 1:5))
  error_count_benchmark[[j]] <- table(factor(results_benchmark$num_of_errors, levels = 0:8))
  
  correct_groups_benchmark    <- c(correct_groups_benchmark, results_benchmark$correct_number_of_groups/n_rep)
  correct_structure_benchmark <- c(correct_structure_benchmark, results_benchmark$correctly_specified_groups/n_rep)

  cat("Analysis of the results for the benchmark method with L2 distance and large bandwidth\n")
  results_benchmark2 <- cluster_analysis(t_len_ = t_len, n_rep_ = n_rep, alpha_ = alpha,
                                         results_matrix_ = clustering_results_benchmark2,
                                         correct_specification_ = correct_specification)
  
  group_count_benchmark2[[j]] <- table(factor(clustering_results_benchmark2[1, ], levels = 1:5))
  error_count_benchmark2[[j]] <- table(factor(results_benchmark2$num_of_errors, levels = 0:8))
  
  correct_groups_benchmark2    <- c(correct_groups_benchmark2, results_benchmark2$correct_number_of_groups/n_rep)
  correct_structure_benchmark2 <- c(correct_structure_benchmark2, results_benchmark2$correctly_specified_groups/n_rep)
  j <- j + 1
}

#######################
#Output of the results#
#######################



tmp <- matrix(NA, nrow = length(different_T), ncol = 3)
tmp[, 1] <- as.vector(correct_structure)
tmp[, 2] <- as.vector(correct_structure_benchmark)
tmp[, 3] <- as.vector(correct_structure_benchmark2)
row.names(tmp) <- paste0("$T = ", different_T, "$")

filename = paste0("output/revision/clustering_comparison.tex")
output_matrix(tmp, filename)

addtorow     <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c("& \\multicolumn{6}{c}{nominal size $\\alpha$} \\\\\n",
                        "$T$ & 0.01 & 0.05 & 0.1 & 0.01 & 0.05 & 0.1 \\\\\n") 
print.xtable(xtable(matrix_, digits = c(3), align = "ccccccc"), type = "latex",
             file = filename, add.to.row = addtorow, include.colnames = FALSE)


line <- paste0("%This simulation was done for the following values of the parameters: n_ts = ", n_ts,
               ", with ", n_rep, " simulations for calculating size and power and ", sim_runs,
               " simulations to calculate the Gaussian quantiles. Furthermore, for the error process we have a = ",
               a, " and sigma = ", sigma,
               ". There are no fixed effects. The grid is standard.")
write(line, file = filename, append = TRUE)

# filename = paste0("output/tables/", n_ts, "_ts_correct_group_number.tex")
# rownames(correct_groups) <- different_T
# output_matrix(correct_groups, filename)
# 
# filename2 = paste0("output/tables/", n_ts, "_ts_correct_group_structure.tex")
# rownames(correct_structure) <- different_T
# output_matrix(correct_structure, filename2)



# filename = paste0("output/tables/", n_ts, "_ts_correct_group_number_benchmark.tex")
# rownames(correct_groups_benchmark) <- different_T
# output_matrix(correct_groups_benchmark, filename)
# 
# filename2 = paste0("output/tables/", n_ts, "_ts_correct_group_structure_benchmark.tex")
# rownames(correct_structure_benchmark) <- different_T
# output_matrix(correct_structure_benchmark, filename)