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

a_hat <- 0.25 
sigma <- 0.25
q     <- 25 #Parameters for the estimation of long-run-variance
r     <- 10

numCores  = round(parallel::detectCores() * .70)

###################################################
#Simulating the data and performing the clustering#
###################################################

#Constructing the set of pairwise comparisons
ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
ijset <- ijset[ijset$i < ijset$j, ]

for (t_len in different_T){
  #Constructing the grid
  u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
  h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len)
  
  m1 <- numeric(t_len)
  m2 <- numeric(t_len)
  m1 <- 0.35 * b_function((1:t_len)/t_len, 0.25, 0.25)
  m2 <- 2 * b_function((1:t_len)/t_len, 0.25, 0.025)

  cat("Calculating the clustering results for t = ", t_len,"\n")
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:n_rep, .combine = "cbind") %dopar% { 
    repl_clustering_revision(val, t_len, n_ts, sigma, a_hat, q, r, grid, m1, m2, h_grid[1]) # Loop one-by-one using foreach
  } -> simulated_statistic
  stopCluster(cl)
  toc()

  cat("Calculating the Gaussian quantiles\n")
  quantiles <- compute_quantiles(t_len = t_len, grid = grid, n_ts = n_ts,
                                 ijset = ijset, sigma = 1,
                                 sim_runs = sim_runs,
                                 deriv_order = 0,
                                 correction = TRUE, epidem = FALSE)
  probs  <- as.vector(quantiles$quant[1, ])
  quants <- as.vector(quantiles$quant[2, ])

  cat("Assesing the results\n")
  
  for (alpha in different_alpha){
    if (sum(probs == (1 - alpha)) == 0)
      pos <- which.min(abs(probs - (1 - alpha)))
    if (sum(probs == (1 - alpha)) != 0)
      pos <- which.max(probs == (1 - alpha))    
    quant <- quants[pos]
    
    groups_mat <- matrix(NA, ncol = n_rep, nrow = n_ts)
    colnames(groups_mat) <- paste0("rep_", 1:n_rep)
    rownames(groups_mat) <- paste0("ts_", 1:n_ts)
    
    groups_benchmark_mat <- matrix(NA, ncol = n_rep, nrow = n_ts)
    colnames(groups_benchmark_mat) <- paste0("rep_", 1:n_rep)
    rownames(groups_benchmark_mat) <- paste0("ts_", 1:n_ts)
    
    number_of_groups_vec <- c()
    for (i in 1:n_rep){
      #Multiscale method
      statistic_vector <- simulated_statistic[1:(nrow(simulated_statistic)/2), i]
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
      
      #Benchmark method
      statistic_vector_benchmark <- simulated_statistic[(nrow(simulated_statistic)/2 + 1):nrow(simulated_statistic), i]
      statistic_matrix_benchmark <- matrix(statistic_vector_benchmark, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
      statistic_matrix_benchmark <- forceSymmetric(statistic_matrix_benchmark, uplo = "U")
      statistic_matrix_benchmark <- as.dist(statistic_matrix_benchmark)
      clustering_benchmark       <- hclust(statistic_matrix_benchmark, method = "complete")
      groups_benchmark           <- cutree(clustering_benchmark, k = 3)
      groups_benchmark_mat[, i]  <- groups_benchmark
    }
    clustering_results <- rbind(number_of_groups_vec, groups_mat)
    clustering_results_benchmark <- rbind(rep(3, n_rep), groups_benchmark_mat)
    
    filename = paste0("output/revision/misc/results_for_T_", t_len, "_and_alpha_", alpha * 100, ".RData")
    save(clustering_results, file = filename)
    filename_benchmark = paste0("output/revision/misc/results_for_T_", t_len, "_and_alpha_", alpha * 100, "_benchmark.RData")
    save(clustering_results_benchmark, file = filename_benchmark)  
  }
}


####################################
#Analysis of the clustering results#
####################################

correct_groups   <- c()
correct_structure <- c()

alpha <- 0.05

group_count <- list()
error_count <- list()

j <- 0

for (t_len in different_T){
  filename = paste0("output/misc/results_for_T_", t_len, "_and_alpha_",
                    alpha*100, "_revision.RData")
  load(file = filename)
  correct_specification      <- c(rep(1, (floor(n_ts / 3))),
                                  rep(2, (floor(2 * n_ts / 3) - floor(n_ts / 3))),
                                  rep(3, n_ts - floor(2 * n_ts / 3)))
  correct_number_of_groups   <- 0 #Starting the counter from zero
  correctly_specified_groups <- 0
  
  num_of_errors     <- c()
  
  for (i in 1:n_rep){
    if (clustering_results[1, i] == 3) {
      correct_number_of_groups = correct_number_of_groups + 1
    }
    if ((clustering_results[1, i] == 2) | (clustering_results[1, i] == 3))
      groups123  <- clustering_results[2:(n_ts + 1), i]
    groups132  <- recode(groups123, "2=3;3=2")
    groups213  <- recode(groups123, "1=2;2=1")
    groups231  <- recode(groups123, "1=2;2=3;3=1")
    groups312  <- recode(groups123, "1=3;2=1;3=2")
    groups321  <- recode(groups123, "1=3;3=1")
    difference <- min(sum(correct_specification != groups132),
                      sum(correct_specification != groups213),
                      sum(correct_specification != groups231),
                      sum(correct_specification != groups312),
                      sum(correct_specification != groups321),
                      sum(correct_specification != groups123))
    if (clustering_results[1, i] == 4) {
      groups1234  <- clustering_results[2:(n_ts + 1), i]
      groups1243  <- recode(groups1234, "4=3;3=4")
      groups1342  <- recode(groups1234, "2=3;3=4;4=2")
      groups1324  <- recode(groups1234, "2=3;3=2")
      groups1423  <- recode(groups1234, "2=4;3=2;4=3")
      groups1432  <- recode(groups1234, "2=4;4=2")
      
      groups2134  <- recode(groups1234, "1=2;2=1")
      groups2143  <- recode(groups1234, "1=2;2=1;3=4;4=3")
      groups2314  <- recode(groups1234, "1=2;2=3;3=1")
      groups2341  <- recode(groups1234, "1=2;2=3;3=4;4=1")
      groups2413  <- recode(groups1234, "1=2;2=4;3=1;4=3")
      groups2431  <- recode(groups1234, "1=2;2=4;4=1")
      
      groups3124  <- recode(groups1234, "1=3;2=1;3=2")
      groups3142  <- recode(groups1234, "1=3;2=1;3=4;4=2")
      groups3214  <- recode(groups1234, "1=3;3=1")
      groups3241  <- recode(groups1234, "1=3;3=4;4=1")
      groups3412  <- recode(groups1234, "1=3;2=4;3=1;4=2")
      groups3421  <- recode(groups1234, "1=3;2=4;3=2;4=1")
      
      groups4123  <- recode(groups1234, "1=4;2=1;3=2;4=3")
      groups4132  <- recode(groups1234, "1=4;2=1;4=2")
      groups4213  <- recode(groups1234, "1=4;3=1;4=3")
      groups4231  <- recode(groups1234, "1=4;4=1")
      groups4312  <- recode(groups1234, "1=4;2=3;3=1;4=2")
      groups4321  <- recode(groups1234, "1=4;2=3;3=2;4=1")
      
      difference <- min(sum(correct_specification != groups1234),
                        sum(correct_specification != groups1243),
                        sum(correct_specification != groups1342),
                        sum(correct_specification != groups1324),
                        sum(correct_specification != groups1423),
                        sum(correct_specification != groups1432),
                        
                        sum(correct_specification != groups2134),
                        sum(correct_specification != groups2143),
                        sum(correct_specification != groups2314),
                        sum(correct_specification != groups2341),
                        sum(correct_specification != groups2413),
                        sum(correct_specification != groups2431),
                        
                        sum(correct_specification != groups3124),
                        sum(correct_specification != groups3142),
                        sum(correct_specification != groups3214),
                        sum(correct_specification != groups3241),
                        sum(correct_specification != groups3412),
                        sum(correct_specification != groups3421),
                        
                        sum(correct_specification != groups4123),
                        sum(correct_specification != groups4132),
                        sum(correct_specification != groups4213),
                        sum(correct_specification != groups4231),
                        sum(correct_specification != groups4312),
                        sum(correct_specification != groups4321))
    }
    if (difference == 0){
      correctly_specified_groups = correctly_specified_groups + 1
    }
    num_of_errors <- c(num_of_errors, difference)
  }
  
  if (alpha == 0.05){
    j <- j + 1
    group_count[[j]] <- table(factor(clustering_results[1, ], levels = 1:5))
    error_count[[j]] <- table(factor(num_of_errors, levels = 0:8))
  }
  
  correct_groups    <- c(correct_groups, correct_number_of_groups/n_rep)
  correct_structure <- c(correct_structure, correctly_specified_groups/n_rep)
  cat("Percentage of detecting true number of clusters",
      correct_number_of_groups/n_rep, "with alpha = ", alpha,
      ", T = ", t_len, "\n")
  cat("Percentage of detecting true clustering",
      correctly_specified_groups/n_rep, "with alpha = ", alpha,
      ", T = ", t_len, "\n")
  cat("Maximum number of errors is ", max(num_of_errors), "\n")
}

#######################
#Output of the results#
#######################

labels_ <- c("(a)", "(b)", "(c)")

pdf(paste0("output/plots/sim/histograms_groups.pdf"), width=8, height=2.9, paper="special")
par(mfrow = c(1,3))
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

pdf(paste0("output/plots/sim/histograms_errors.pdf"), width=8, height=2.9, paper="special")
par(mfrow = c(1,3))
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


filename = paste0("output/tables/", n_ts, "_ts_correct_group_number.tex")
creating_matrix_and_texing(correct_groups, different_T, different_alpha, filename)
filename2 = paste0("output/tables/", n_ts, "_ts_correct_group_structure.tex")
creating_matrix_and_texing(correct_structure, different_T, different_alpha, filename2)
