# This is the main file for the analysis of both applications which is reported in Section 6.
rm(list=ls())

library(multiscale)
#library(xtable)
#options(xtable.floating = FALSE)
#options(xtable.timestamp = "")

##############################
#Defining necessary constants#
##############################

alpha    <- 0.05 #confidence level for application
sim_runs <- 1000


#################################
#Loading the exchange rates data#
#################################

exchange_rates <- read.csv("data/exchange_rates.csv", sep = ",", dec = ".", quote = '"', stringsAsFactors = FALSE)
exchange_rates[['exvzus']] <- NULL

#exchange_rates <- as.matrix(exchange_rates)
colSums(is.na(exchange_rates))

#exchange_rates <- exchange_rates[,colSums(is.na(exchange_rates)) <= 100] #Omitting the time series with too sparse data
exchange_rates <- na.omit(exchange_rates)#Deleting the rows with omitted variables
exchange_rates[["date"]] <- as.Date(as.character(exchange_rates[["date"]]), format = "%d/%m/%Y")

#exchange_rates <- exchange_rates[1:300, ]


t_len         <- nrow(exchange_rates)
n_ts          <- ncol(exchange_rates) - 1 #Updating the number of time series because of dropped stations

column_names  <- names(exchange_rates[, 2:(n_ts + 1)])

#exchange_rates <-matrix(unlist(exchange_rates), nrow = t_len)

#colnames(exchange_rates) <- column_names
exchange_rates[, 2:(n_ts + 1)] <- scale(exchange_rates[, 2:(n_ts + 1)], scale = FALSE)

pdf("plots/exchange_rates.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(x = exchange_rates[[1]], y = exchange_rates[[2]],
     ylim=c(min(exchange_rates[[2]], exchange_rates[[4]]),
            max(exchange_rates[[2]], exchange_rates[[4]] + 0.5)), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
lines(x = exchange_rates[[1]], y = exchange_rates[[4]], col="#604c38")
legend("topleft", inset = 0.02, legend=c("Australian Dollar", "New Zealand Dollar"),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()

pdf("plots/exchange_rates_1.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins


plot(x = exchange_rates[[1]], y = exchange_rates[[2]],
     ylim=c(min(exchange_rates[[2]], exchange_rates[[4]]),
            max(exchange_rates[[2]], exchange_rates[[4]] + 0.5)), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
rect(xleft=as.Date("2018-01-01", format = "%Y-%m-%d"),
     xright = as.Date("2018-12-31", format = "%Y-%m-%d"),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
lines(x = exchange_rates[[1]], y = exchange_rates[[2]], col="#EB811B")
lines(x = exchange_rates[[1]], y = exchange_rates[[4]], col="#604c38")
segments(as.Date("2018-01-01", format = "%Y-%m-%d"), par("usr")[3],
         as.Date("2018-12-31", format = "%Y-%m-%d"), par("usr")[3],
         col="#604c38", lwd = 3)
legend("topleft", inset = 0.02, legend=c("Australian Dollar", "New Zealand Dollar"),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()


pdf("plots/exchange_rates_2.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(x = exchange_rates[[1]], y = exchange_rates[[2]],
     ylim=c(min(exchange_rates[[2]], exchange_rates[[4]]),
            max(exchange_rates[[2]], exchange_rates[[4]] + 0.5)), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
rect(xleft=as.Date("2018-01-01", format = "%Y-%m-%d"),
     xright = as.Date("2018-12-31", format = "%Y-%m-%d"),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
rect(xleft=as.Date("2011-08-01", format = "%Y-%m-%d"),
     xright = as.Date("2013-06-30", format = "%Y-%m-%d"),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
lines(x = exchange_rates[[1]], y = exchange_rates[[2]], col="#EB811B")
lines(x = exchange_rates[[1]], y = exchange_rates[[4]], col="#604c38")
segments(as.Date("2018-01-01", format = "%Y-%m-%d"), par("usr")[3],
         as.Date("2018-12-31", format = "%Y-%m-%d"), par("usr")[3],
         col="#604c38", lwd = 3)
segments(as.Date("2011-08-01", format = "%Y-%m-%d"), par("usr")[3],
         as.Date("2013-06-30", format = "%Y-%m-%d"), par("usr")[3], col="#604c38", lwd = 3)
legend("topleft", inset = 0.02, legend=c("Australian Dollar", "New Zealand Dollar"),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()

#####################
#Estimating variance#
#####################

#Order selection
q <- 30:60
r <- 10:15
order_results <- c()

for (j in 2:(n_ts + 1)){
  criterion_matrix <- expand.grid(q = q, r = r)
  
  criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$SIC  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))
  
  for (i in 1:nrow(criterion_matrix)){
    FPE <- c()
    AIC <- c()
    AICC <- c()
    SIC <- c()
    HQ <- c()
    
    different_orders <- (1:9)
    
    for (order in different_orders){
      AR.struc      <- estimate_lrv(data=exchange_rates[[j]], q=criterion_matrix$q[[i]], r_bar=criterion_matrix$r[[i]], p=order)
      sigma_eta_hat <- sqrt(AR.struc$vareta)
      FPE <- c(FPE, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
      AIC <- c(AIC, t_len * log(sigma_eta_hat^2) + 2 * order)
      AICC <- c(AICC, t_len * log(sigma_eta_hat^2) + t_len * (1 + order / t_len)/(1 - (order +2)/t_len))
      SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
      HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
    }
    criterion_matrix$FPE[[i]]  <- which.min(FPE)
    criterion_matrix$AIC[[i]]  <- which.min(AIC)
    criterion_matrix$AICC[[i]] <- which.min(AICC)
    criterion_matrix$SIC[[i]]  <- which.min(SIC)
    criterion_matrix$HQ[[i]]   <- which.min(HQ)
  }
  maxim <- max(criterion_matrix[, 3:7])
  order_results <- c(order_results, maxim)
  cat("For exchange rate ", names(exchange_rates)[j], " the results are as follows: ", max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")
}


#Setting tuning parameters for testing
q     <- 55
r_bar <- 10


#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in 2:(n_ts+1)){
  AR.struc        <- estimate_lrv(data = exchange_rates[[i]], q = q, r_bar = r_bar, p=order_results[i-1])
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}

#Constructing the grid
u_grid <- seq(from = 100 / t_len, to = 1, by = 100 / t_len)
h_grid <- seq(from = 100 / t_len, to = 1 / 4, by = 100 / t_len)
h_grid <- h_grid[h_grid > log(t_len) / t_len]
grid <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len) #for plotting

#Calculating the statistic for real data
result <- multiscale_test(data = matrix(unlist(exchange_rates[, 2:(n_ts + 1)]), ncol = n_ts, byrow = FALSE),
                          sigma_vec = sigmahat_vector,
                          alpha = alpha,
                          n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)

source("functions.R")

for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]
  filename = paste0("plots/", column_names[i], "_vs_", column_names[j], "_presentation.pdf")
  produce_plots_talk(results = result, l = l, data_i = exchange_rates[[i + 1]],
                     data_j = exchange_rates[[j + 1]], dates = exchange_rates[[1]],
                     name_i = column_names[i], name_j = column_names[j],
                     filename = filename)
}
# 
# 
# #for the distance matrix we need a symmetrical one
# Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
# for (i in 1:(n_ts - 1)){
#   for (j in (i + 1):n_ts){
#     Delta_hat[i, j] <- result$stat_pairwise[i, j]
#     Delta_hat[j, i] <- result$stat_pairwise[i, j]
#   }
# }
# 
# currencies <- colnames(exchange_rates[, 2:(n_ts + 1)])
# 
# colnames(Delta_hat) <- currencies
# rownames(Delta_hat) <- currencies
# 
# delta_dist  <- as.dist(Delta_hat)
# res         <- hclust(delta_dist)
# n_cl        <- 5
# grid_points <- seq(1/t_len, 1, by = 1/t_len)
# 
# #Plotting dendrogram
# plot(res, cex = 0.8, xlab = "", ylab = "")
# rect.hclust(res, k = n_cl, border = 2:(n_cl + 1))
# 
# subgroups <- cutree(res, n_cl)
# 
# 
# 
# 
# 
# 
# for (cl in 1:n_cl){
#   countries_cluster <- colnames(Delta_hat)[subgroups == cl]
# 
#   if (length(countries_cluster) == 1){
#     m_hat_vec <- m_hat(grid_points, b = 1, covid_mat[, countries_cluster],
#                        grid_points, bw = bw_abs/t_len)
#     norm      <- integrate1_cpp(b = 1, data_points = covid_mat[, countries_cluster],
#                                   grid_points = grid_points,
#                                   bw = bw_abs/t_len, subdiv = 2000)$res
#       plot((1:t_len) / t_len, m_hat_vec/norm, yaxt = "n",
#            ylim = c(0, max(m_hat_vec/norm) + 1), xlab="u",
#            ylab = "", mgp = c(2, 0.5, 0), type = "l", col = "red")
#       title(main = paste("Cluster", cl), line = 1)
#       legend("topleft", inset = 0.02, legend=countries_cluster,
#              lty = 1, cex = 0.7, ncol = 1)
#     } else {
#       b_res_cl     <- b_res[subgroups == cl, subgroups == cl]
#       inds         <- which.max(apply(b_res_cl, 1, function(x) sum(x == 1, na.rm = TRUE)))
#       repr_country <- rownames(b_res_cl)[inds]
#       m_hat_vec    <- m_hat(grid_points, b = 1, covid_mat[, repr_country],
#                             grid_points, bw = bw_abs/t_len)
#       norm         <- integrate1_cpp(b = 1, data_points = covid_mat[, repr_country],
#                                      grid_points = grid_points,
#                                      bw = bw_abs/t_len, subdiv = 2000)$res
#       #cat("Country", repr_country, ", cluster", cl, " - success \n")
#       if (cl == 2) {height <- 8} else {height <- 3} #This should be manually adjusted for nice plots
#       plot(grid_points, m_hat_vec/norm,
#            ylim = c(0, max(m_hat_vec/norm) + height), xlab="u", yaxt = "n",
#            ylab = "m_hat(b * u)", mgp = c(2, 0.5, 0), type = "l", col = "red")
#       countries_cluster_1 <- countries_cluster[countries_cluster != repr_country]
#       for (country in countries_cluster_1){
#         b           <- max(1, b_res_cl[country, repr_country] / b_res_cl[repr_country, country])
#         m_hat_vec_1 <- m_hat(grid_points, b = b, covid_mat[, country],
#                              grid_points, bw = bw_abs/t_len)
#         m_hat_vec_1[(m_hat_vec_1 == 0 | is.nan(m_hat_vec_1))] <- NA
#         norm_1      <- integrate1_cpp(b = b, data_points = covid_mat[, country],
#                                       grid_points = grid_points,
#                                       bw = bw_abs/t_len, subdiv = 2000)$res
#         #cat("Country", country, " - success \n")
#         lines((1:length(m_hat_vec_1)) / t_len, m_hat_vec_1/(norm_1/(1/b)))
#       }
#       title(main = paste("Cluster", cl), line = 1)
#       legend("topleft", inset = 0.02, legend = countries_cluster,
#              lty = 1, cex = 0.7, ncol = 4)
#     }
#     dev.off()
#   }
# }
# 
# 
# 
# 
# n_cl       <- 12
# #results_output(res, covid_mat, Delta_hat, b_res, n_cl, countries,
# #               path = "plots/", bw_abs, grid_points, t_len)
#