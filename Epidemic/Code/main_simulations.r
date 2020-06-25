rm(list=ls())
library(multiscale)
library(tictoc)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


n_sim     <- 5000               # number of simulation runs for power and size
sim_runs  <- 5000               # number of simulation runs for calculating quantiles
alpha_vec <- c(0.01, 0.05, 0.1) # significance level
n_ts_vec  <- c(5, 10, 50)       # number of time series
t_len_vec <- c(100, 250, 500)   # time series length
sigma_vec <- c(10, 15, 20)      # overdispersion parameter
const_c   <- c(500, 1000)       # constant in the lambda function

source("functions/functions.r")

for (c in const_c){
  for (sigma in sigma_vec){
    number_of_cols         <- length(n_ts_vec) * length(alpha_vec)
    size_matrix            <- matrix(NA, nrow = length(t_len_vec), ncol = number_of_cols)
    rownames(size_matrix)  <- paste0("T = ", t_len_vec)
    
    k <- 1
    for (n_ts in n_ts_vec){
      i <- 1
      for (t_len in t_len_vec){
        lambda_vec <- lambda_fct((1:t_len) / t_len, c)
        # if(t_len == 100) {
        #   pdf(paste0("plots/lambda_fct_", c, ".pdf"), width=5, height=5, paper="special")
        #   par(mar = c(3, 2, 2, 0)) #Margins for each plot
        #   par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
        #   plot(lambda_vec,  ylim = c(0, max(lambda_vec) + 100), xlab="u", ylab = "", mgp=c(2,0.5,0), type = "l")
        #   title(main = expression(Plot ~ of ~ the ~ "function" ~ lambda), line = 1)
        #   dev.off()
        # }
        
        #set.seed(1234) # This is for calculating size for different specifications on comparable datasets
        size <- calculate_size(t_len = t_len, n_ts = n_ts, alpha_vec = alpha_vec,
                               lambda_vec = lambda_vec, sigma = sigma,
                               n_sim = n_sim, sim_runs = sim_runs, iid = FALSE)
        size_matrix[i, (k * length(alpha_vec) - (length(alpha_vec) - 1)):(k * length(alpha_vec))] <- size
        i <- i + 1
      }
      k <- k + 1
    }
    
    print.xtable(xtable(size_matrix, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
                 type="latex", file=paste0("plots/size_overdispersion_", sigma, "_const_c_", c, "neg_binom.tex"), include.colnames = FALSE)
  }
}

for (c in const_c){
  for (sigma in sigma_vec){
    number_of_cols         <- length(n_ts_vec) * length(alpha_vec)
    size_matrix            <- matrix(NA, nrow = length(t_len_vec), ncol = number_of_cols)
    rownames(size_matrix)  <- paste0("T = ", t_len_vec)
    
    k <- 1
    for (n_ts in n_ts_vec){
      i <- 1
      for (t_len in t_len_vec){
        lambda_vec <- lambda_fct((1:t_len) / t_len, c)

        #set.seed(1234) # This is for calculating size for different specifications on comparable datasets
        size <- calculate_size(t_len = t_len, n_ts = n_ts, alpha_vec = alpha_vec,
                               lambda_vec = lambda_vec, sigma = sigma,
                               n_sim = n_sim, sim_runs = sim_runs, iid = TRUE)
        size_matrix[i, (k * length(alpha_vec) - (length(alpha_vec) - 1)):(k * length(alpha_vec))] <- size
        i <- i + 1
      }
      k <- k + 1
    }
    
    print.xtable(xtable(size_matrix, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
                 type="latex", file=paste0("plots/size_overdispersion_", sigma, "_const_c_", c, "iid_normal.tex"), include.colnames = FALSE)
  }
}
