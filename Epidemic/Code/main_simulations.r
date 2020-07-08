rm(list=ls())
library(multiscale)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

n_sim     <- 5000               # number of simulation runs for power and size
sim_runs  <- 5000               # number of simulation runs for calculating quantiles
alpha_vec <- c(0.01, 0.05, 0.1) # significance level
n_ts_vec  <- c(5, 10, 50)           # number of time series
t_len_vec <- c(100, 250, 500)   # time series length
sigma_vec <- c(10, 15, 20)      # overdispersion parameter

number_of_cols <- length(n_ts_vec) * length(alpha_vec) #Needed for the output

#Auxiliary functions
source("functions/functions.r")

for (sigma in sigma_vec){
  number_of_cols         <- length(n_ts_vec) * length(alpha_vec)
  size_matrix            <- matrix(NA, nrow = length(t_len_vec), ncol = number_of_cols)
  rownames(size_matrix)  <- paste0("$T = ", t_len_vec, "$")

  k <- 1
  for (n_ts in n_ts_vec){
    i <- 1
    for (t_len in t_len_vec){
      lambda_vec <- lambda_fct((1:t_len) / t_len, c = 1000)
      pdf(paste0("plots/lambda_fct.pdf"), width=5, height=5, paper="special")
      par(mar = c(3, 2, 2, 0)) #Margins for each plot
      par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
      plot((1:t_len) / t_len, lambda_vec,  ylim = c(0, max(lambda_vec) + 100), xlab="u", ylab = "", mgp=c(2,0.5,0), type = "l")
      title(main = expression(Plot ~ of ~ the ~ "function" ~ lambda), line = 1)
      dev.off()

      set.seed(10) # This is for calculating size for different specifications on comparable datasets
      size <- calculate_size(t_len = t_len, n_ts = n_ts, alpha_vec = alpha_vec,
                             lambda_vec = lambda_vec, sigma = sigma,
                             n_sim = n_sim, sim_runs = sim_runs, iid = FALSE)
      size_matrix[i, (k * length(alpha_vec) - (length(alpha_vec) - 1)):(k * length(alpha_vec))] <- size
      i <- i + 1
    }
    k <- k + 1
  }
  print.xtable(xtable(size_matrix, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
               type="latex", file=paste0("plots/size_overdispersion_", sigma, ".tex"),
               include.colnames = FALSE, sanitize.text.function = function(x) {x})
}


for (sigma in sigma_vec){
  power_matrix           <- matrix(NA, nrow = length(t_len_vec), ncol = number_of_cols)
  rownames(power_matrix) <- paste0("$T = ", t_len_vec, "$")
  k <- 1
  for (n_ts in n_ts_vec){
    i <- 1
    for (t_len in t_len_vec){
      #Here you can change the functions for the power calculations
      lambda_vec_1 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 6000, position = 10)
      lambda_vec_2 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10)

      pdf(paste0("plots/lambda_fcts_height.pdf"), width=5, height=5, paper="special")
      par(mar = c(3, 2, 0, 0)) #Margins for each plot
      par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
      plot((1:t_len) / t_len, lambda_vec_1,  ylim = c(0, max(lambda_vec_1, lambda_vec_2) + 100), xlab="u", ylab = "", type = "l")
      lines((1:t_len) / t_len, lambda_vec_2, type = "l", col = "red")
      dev.off()
      
      set.seed(12) # This is for calculating power for different specifications on comparable datasets
      power <- calculate_power(t_len = t_len, n_ts = n_ts, alpha_vec = alpha_vec,
                               lambda_vec_1 = lambda_vec_1, lambda_vec_2 = lambda_vec_2,
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



for (sigma in sigma_vec){
  power_matrix           <- matrix(NA, nrow = length(t_len_vec), ncol = number_of_cols)
  rownames(power_matrix) <- paste0("$T = ", t_len_vec, "$")

  k <- 1
  for (n_ts in n_ts_vec){
    i <- 1
    for (t_len in t_len_vec){
      lambda_vec_1 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 9)
      lambda_vec_2 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10)
      
      pdf(paste0("plots/lambda_fcts_shift.pdf"), width=5, height=5, paper="special")
      par(mar = c(3, 2, 0, 0)) #Margins for each plot
      par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
      plot((1:t_len) / t_len, lambda_vec_1,  ylim = c(0, max(lambda_vec_1, lambda_vec_2) + 100), xlab="u", ylab = "", type = "l")
      lines((1:t_len) / t_len, lambda_vec_2, type = "l", col = "red")
      dev.off()
      
      set.seed(123) # This is for calculating power for different specifications on comparable datasets
      power <- calculate_power(t_len = t_len, n_ts = n_ts, alpha_vec = alpha_vec,
                               lambda_vec_1 = lambda_vec_1, lambda_vec_2 = lambda_vec_2,
                               sigma = sigma, n_sim = n_sim, sim_runs = sim_runs)
      power_matrix[i, (k * length(alpha_vec) - (length(alpha_vec) - 1)):(k * length(alpha_vec))] <- power
      i <- i + 1
    }
    k <- k + 1
  }
  print.xtable(xtable(power_matrix, digits = c(3),
                      align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
               type="latex", file=paste0("plots/power_sigma_", sigma, "_shifted_peak.tex"),
               include.colnames = FALSE, sanitize.text.function = function(x) {x})
}
