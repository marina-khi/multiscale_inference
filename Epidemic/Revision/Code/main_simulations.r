#This is the main file for producing the simulation results for size and power, that are reported in Section 4.1.
rm(list=ls())

library(multiscale)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

#The following file contains all the necessary auxiliary functions that produce all the results.
#All the necessary arguments for these functions are described in detail in the file.
source("functions/functions.r")

#Random generation of the seed. The seed is necessary for computing all the different
#specifications on comparable datasets.
seed <- sample(1:100000, 1)


##################
#Calculating size#
##################

n_sim     <- 5000               # number of simulation runs for power and size
sim_runs  <- 5000               # number of simulation runs to produce critical values
alpha_vec <- c(0.01, 0.05, 0.1) # different significance levels
n_ts_vec  <- c(5, 10, 50)       # different number of time series
t_len_vec <- c(100, 250, 500)   # different time series lengths
sigma_vec <- c(15, 10, 20)      # different overdispersion parameter

number_of_cols <- length(n_ts_vec) * length(alpha_vec) #Needed for the output

#As the mean function in the size simulations, we take the following function:
#lambda(u) = 5000 * exp(-(10 * u - 3) ^ 2 / 2) + 1000 for 0 <= u <= 1.
#Here is the plot of this function:

lambda_vec <- lambda_fct((1:100) / 100)

pdf(paste0("plots/lambda_fct.pdf"), width=5, height=5, paper="special")
par(mar = c(3, 1.5, 1.5, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
plot((1:100) / 100, lambda_vec,  ylim = c(0, max(lambda_vec) + 100), xlab="",
     ylab = "", mgp=c(2,0.5,0), type = "l", cex.lab = 0.85, cex.axis = 0.85)
title(xlab="u", line = 0.9, cex.lab = 0.85)
title(main = expression(Plot ~ of ~ the ~ "function" ~ lambda), line = 0.7, cex.main = 0.85)
title(sub = "(a)", cex.sub = 0.85, line = 2.0)
dev.off()

#However, you can change this function to whatever you like.
#In order to do this, change the definition of lambda_vec further in the code.
#Specifically, you will need to redefine lambda_vec to be a vector of length t_len
#with the values that your preferred function takes at time points 1/t_len, 2/t_len, ..., 1.

for (sigma in sigma_vec){
  size_matrix            <- matrix(NA, nrow = length(t_len_vec), ncol = number_of_cols)
  rownames(size_matrix)  <- paste0("$T = ", t_len_vec, "$")

  k <- 1
  for (n_ts in n_ts_vec){
    i <- 1
    for (t_len in t_len_vec){
      #Here you can change the functions for the size calculations
      lambda_vec <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10)
      
      set.seed(321)
      size <- calculate_size(t_len = t_len, n_ts = n_ts, alpha_vec = alpha_vec,
                             lambda_vec = lambda_vec, sigma = sigma,
                             n_sim = n_sim, sim_runs = sim_runs)
      size_matrix[i, (k * length(alpha_vec) - (length(alpha_vec) - 1)):(k * length(alpha_vec))] <- size
      i <- i + 1
    }
    k <- k + 1
  }
  print.xtable(xtable(size_matrix, digits = c(3),
                      align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
               type="latex", file=paste0("plots/size_overdispersion_", sigma, ".tex"),
               include.colnames = FALSE, sanitize.text.function = function(x) {x})
}

#Now the results of size simulations are stored as the tex tables in the folder ./plots/


#############################################################
#Calculating power for mean functions with different heights#
#############################################################

n_sim     <- 5000               # number of simulation runs for power and size
sim_runs  <- 5000               # number of simulation runs to produce critical values
alpha_vec <- c(0.01, 0.05, 0.1) # different significance levels
n_ts_vec  <- c(5, 10, 50)       # different number of time series
t_len_vec <- c(100, 250, 500)   # different time series lengths
sigma_vec <- c(15, 10, 20)      # different overdispersion parameter

number_of_cols <- length(n_ts_vec) * length(alpha_vec) #Needed for the output

#As the mean functions in these power simulations, we take the following functions:
#lambda_1(u) = 6000 * exp(-(10 * u - 3) ^ 2 / 2) + 1000 for 0 <= u <= 1.
#lambda(u) = 5000 * exp(-(10 * u - 3) ^ 2 / 2) + 1000 for 0 <= u <= 1.

#Here are the plots of these functions:

lambda_vec_1 <- lambda_fct((1:100) / 100, c = 1000, height = 6000, position = 10)
lambda_vec   <- lambda_fct((1:100) / 100, c = 1000, height = 5000, position = 10)


pdf(paste0("plots/lambda_fcts_height.pdf"), width=5, height=4, paper="special")
par(mar = c(3.3, 1.5, 1.5, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
par(mgp = c(3, 0.5, 0))
plot((1:100) / 100, lambda_vec_1,  ylim = c(0, max(lambda_vec_1, lambda_vec) + 100),
     xlab="", ylab = "", type = "l", cex.axis = 0.85)
lines((1:100) / 100, lambda_vec, type = "l", col = "red")
title(main = expression(Plot ~ of ~ the ~ "functions" ~ lambda[1] ~ and ~ lambda), line = 0.7, cex.main = 0.85)
title(xlab="u", line=0.9, cex.lab = 0.85)
title(sub = "(a) Scenario A", cex.sub = 0.85, line = 2.0)
legend("topright", inset = 0.02, legend=c(expression(lambda[1](u) ~" "), expression(lambda(u) ~" ")),
       col = c("black", "red"), lty = 1, cex = 0.75, ncol = 1)
dev.off()

#However, you can change the function to whatever you like.
#In order to do this, change the definition of lambda_vec_1 and lambda_vec further in the code.
#Specifically, you will need to redefine lambda_vec_1 and lambda_vec to be vectors
#of length t_len with the values that your preferred functions take at time points 1/t_len, 2/t_len, ..., 1.


for (sigma in sigma_vec){
  power_matrix           <- matrix(NA, nrow = length(t_len_vec), ncol = number_of_cols)
  rownames(power_matrix) <- paste0("$T = ", t_len_vec, "$")
  k <- 1
  for (n_ts in n_ts_vec){
    i <- 1
    for (t_len in t_len_vec){
      #Here you can change the functions for the power calculations
      lambda_vec_1 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 6000, position = 10)
      lambda_vec   <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10)

      set.seed(321)
      power <- calculate_power(t_len = t_len, n_ts = n_ts, alpha_vec = alpha_vec,
                               lambda_vec_1 = lambda_vec_1, lambda_vec_2 = lambda_vec,
                               sigma = sigma, n_sim = n_sim, sim_runs = sim_runs)
      power_matrix[i, (k * length(alpha_vec) - (length(alpha_vec) - 1)):(k * length(alpha_vec))] <- power
      i <- i + 1
    }
    k <- k + 1
  }
    
  print.xtable(xtable(power_matrix, digits = c(3),
                      align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
               type="latex", file=paste0("plots/power_sigma_", sigma, "_higher_peak.tex"),
               include.colnames = FALSE, sanitize.text.function = function(x) {x})
}

#Now the results of power simulations for this scenario are stored as the tex tables in the folder ./plots/


####################################################################
#Calculating power for mean functions with different peak locations#
####################################################################

n_sim     <- 5000               # number of simulation runs for power and size
sim_runs  <- 5000               # number of simulation runs to produce critical values
alpha_vec <- c(0.01, 0.05, 0.1) # different significance levels
n_ts_vec  <- c(5, 10, 50)       # different number of time series
t_len_vec <- c(100, 250, 500)   # different time series lengths
sigma_vec <- c(15, 10, 20)      # different overdispersion parameter

number_of_cols <- length(n_ts_vec) * length(alpha_vec) #Needed for the output

#As the mean function in these power simulations, we take the following functions:
#lambda_1(u) = 5000 * exp(-(9 * u - 3) ^ 2 / 2) + 1000 for 0 <= u <= 1.
#lambda(u) = 5000 * exp(-(10 * u - 3) ^ 2 / 2) + 1000 for 0 <= u <= 1.

#Here are the plots of these functions:

lambda_vec_1 <- lambda_fct((1:100) / 100, c = 1000, height = 5000, position = 9)
lambda_vec   <- lambda_fct((1:100) / 100, c = 1000, height = 5000, position = 10)

pdf(paste0("plots/lambda_fcts_shift.pdf"), width=5, height=4, paper="special")
par(mar = c(3.3, 1.5, 1.5, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
par(mgp = c(3, 0.5, 0))
plot((1:100) / 100, lambda_vec_1,  ylim = c(0, max(lambda_vec_1, lambda_vec) + 100),
     xlab="", ylab = "", type = "l", cex.axis = 0.85)
lines((1:100) / 100, lambda_vec, type = "l", col = "red")
title(main = expression(Plot ~ of ~ the ~ "functions" ~ lambda[1] ~ and ~ lambda), line = 0.7, cex.main = 0.85)
title(xlab="u", line=0.9, cex.lab = 0.85)
title(sub = "(b) Scenario B", cex.sub = 0.85, line = 2.0)
legend("topright", inset = 0.02, legend=c(expression(lambda[1](u) ~" "), expression(lambda(u) ~" ")),
       col = c("black", "red"), lty = 1, cex = 0.75, ncol = 1)
dev.off()


#However, you can change the function to whatever you like.
#In order to do this, change the definition of lambda_vec_1 and lambda_vec further in the code.
#Specifically, you will need to redefine lambda_vec_1 and lambda_vec to be vectors
#of length t_len with the values that your preferred functions take at time points 1/t_len, 2/t_len, ..., 1.


for (sigma in sigma_vec){
  power_matrix2           <- matrix(NA, nrow = length(t_len_vec), ncol = number_of_cols)
  rownames(power_matrix2) <- paste0("$T = ", t_len_vec, "$")

  k <- 1
  for (n_ts in n_ts_vec){
    i <- 1
    for (t_len in t_len_vec){
      #Here you can change the functions for the power calculations
      lambda_vec_1 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 9)
      lambda_vec   <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10)
      
      set.seed(321) # This is for calculating power for different specifications on comparable datasets
      power <- calculate_power(t_len = t_len, n_ts = n_ts, alpha_vec = alpha_vec,
                               lambda_vec_1 = lambda_vec_1, lambda_vec_2 = lambda_vec,
                               sigma = sigma, n_sim = n_sim, sim_runs = sim_runs)
      power_matrix2[i, (k * length(alpha_vec) - (length(alpha_vec) - 1)):(k * length(alpha_vec))] <- power
      i <- i + 1
    }
    k <- k + 1
  }
  print.xtable(xtable(power_matrix2, digits = c(3),
                      align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
               type="latex", file=paste0("plots/power_sigma_", sigma, "_shifted_peak.tex"),
               include.colnames = FALSE, sanitize.text.function = function(x) {x})
}

#Now the results of power simulations for this scenario are stored as the tex tables in the folder ./plots/
pdf(paste0("plots/figure_2.pdf"), width=10.5, height=4, paper="special")
lambda_vec_1 <- lambda_fct((1:100) / 100, c = 1000, height = 6000, position = 10)
lambda_vec   <- lambda_fct((1:100) / 100, c = 1000, height = 5000, position = 10)

layout(matrix(c(1, 2),ncol=2), widths=c(5, 5.5),
       heights=c(4, 4), TRUE)
par(mar = c(3.3, 1.5, 1.5, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
par(mgp = c(3, 0.5, 0))
plot((1:100) / 100, lambda_vec_1,  ylim = c(0, max(lambda_vec_1, lambda_vec) + 100),
     xlab="", ylab = "", type = "l", cex.axis = 0.85)
lines((1:100) / 100, lambda_vec, type = "l", col = "red")
title(main = expression(Plot ~ of ~ the ~ "functions" ~ lambda[1] ~ and ~ lambda), line = 0.7, cex.main = 0.85)
title(xlab="u", line=0.9, cex.lab = 0.85)
title(sub = "(a) Scenario A", cex.sub = 0.85, line = 2.0)
legend("topright", inset = 0.02, legend=c(expression(lambda[1](u) ~" "), expression(lambda(u) ~" ")),
       col = c("black", "red"), lty = 1, cex = 0.75, ncol = 1)

lambda_vec_1 <- lambda_fct((1:100) / 100, c = 1000, height = 5000, position = 9)
lambda_vec   <- lambda_fct((1:100) / 100, c = 1000, height = 5000, position = 10)

par(mar = c(3.3, 3, 1.5, 0)) #Margins for each plot
plot((1:100) / 100, lambda_vec_1,  ylim = c(0, max(lambda_vec_1, lambda_vec) + 100),
     xlab="", ylab = "", yaxt = "n", type = "l", cex.axis = 0.85)
lines((1:100) / 100, lambda_vec, type = "l", col = "red")
axis(2, at = seq(from = 0, to = 6000, by = 1000),labels=c("0", "1000", "", "3000", "", "5000", ""), las=0, cex.axis=0.85)
title(main = expression(Plot ~ of ~ the ~ "functions" ~ lambda[1] ~ and ~ lambda), line = 0.7, cex.main = 0.85)
title(xlab="u", line=0.9, cex.lab = 0.85)
title(sub = "(b) Scenario B", cex.sub = 0.85, line = 2.0)
legend("topright", inset = 0.02, legend=c(expression(lambda[1](u) ~" "), expression(lambda(u) ~" ")),
       col = c("black", "red"), lty = 1, cex = 0.75, ncol = 1)
dev.off()
