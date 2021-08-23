# This is the main file for the analysis of both applications which is reported in Section 6.
rm(list=ls())

library(multiscale)
library(tictoc)
#library(xtable)
#options(xtable.floating = FALSE)
#options(xtable.timestamp = "")

alpha         <- 0.05 #confidence level for application
sim_runs      <- 1000

###########################################
#Loading the real station data for England#
###########################################

dow_jones <- read.csv("data/dow-jones.csv", skip = 15, header = TRUE)
nasdaq    <- read.csv("data/nasdaq.csv", skip = 15, header = TRUE)

indices <- merge(dow_jones, nasdaq, by = c("date"))
colnames(indices) <- c('date', 'dow-jones', 'nasdaq')

rm(dow_jones, nasdaq)

n_ts   <- ncol(indices) - 1
ind_names <- c("dow-jones", "nasdaq")
colSums(is.na(indices))

indices[, 2:(n_ts + 1)] <- scale(indices[, 2:(n_ts + 1)], scale = TRUE)
indices <- na.omit(indices)#Deleting the rows with ommitted variables
indices[, "date"] <- as.Date(as.character(indices[, "date"]), format = "%Y-%m-%d")
t_len  <- nrow(indices)

pdf("plots/indices.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(x = indices[, "date"], y = indices[, "dow-jones"], ylim=c(min(indices[, "dow-jones"], indices[, "nasdaq"]),
                          max(indices[, "dow-jones"], indices[, "nasdaq"]) + 2), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
lines(x = indices[, "date"], y = indices[, "nasdaq"], col="#604c38")
legend("topleft", inset = 0.02, legend=c("Dow Jones index", "NASDAQ index"),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()

pdf("plots/indices_1.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins


plot(x = indices[, "date"], y = indices[, "dow-jones"],
     ylim=c(min(indices[, "dow-jones"], indices[, "nasdaq"]),
            max(indices[, "dow-jones"], indices[, "nasdaq"]) + 2),
     type="l", col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
rect(xleft=as.Date("2012-01-01", format = "%Y-%m-%d"),
     xright = as.Date("2016-12-01", format = "%Y-%m-%d"),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
lines(x = indices[, "date"], y = indices[, "nasdaq"], col="#604c38")
lines(x = indices[, "date"], y = indices[, "dow-jones"], col="#EB811B")
segments(as.Date("2012-01-01", format = "%Y-%m-%d"), par("usr")[3],
         as.Date("2016-12-01", format = "%Y-%m-%d"), par("usr")[3], col="#604c38", lwd = 3)
legend("topleft", inset = 0.02, legend=c("Dow Jones index", "NASDAQ index"),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()


pdf("plots/indices_2.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(x = indices[, "date"], y = indices[, "dow-jones"],
     ylim=c(min(indices[, "dow-jones"], indices[, "nasdaq"]),
            max(indices[, "dow-jones"], indices[, "nasdaq"]) + 2),
     type="l", col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
rect(xleft=as.Date("2012-01-01", format = "%Y-%m-%d"),
     xright = as.Date("2016-12-01", format = "%Y-%m-%d"),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
rect(xleft=as.Date("2020-03-01", format = "%Y-%m-%d"),
     xright = as.Date("2021-08-12", format = "%Y-%m-%d"),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
lines(x = indices[, "date"], y = indices[, "nasdaq"], col="#604c38")
lines(x = indices[, "date"], y = indices[, "dow-jones"], col="#EB811B")
segments(as.Date("2012-01-01", format = "%Y-%m-%d"), par("usr")[3],
         as.Date("2016-12-01", format = "%Y-%m-%d"), par("usr")[3], col="#604c38", lwd = 3)
segments(as.Date("2020-03-01", format = "%Y-%m-%d"), par("usr")[3],
         as.Date("2021-08-12", format = "%Y-%m-%d"), par("usr")[3], col="#604c38", lwd = 3)
legend("topleft", inset = 0.02, legend=c("Dow Jones index", "NASDAQ index"),
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
                        AR.struc      <- estimate_lrv(data=indices[[j]], q=criterion_matrix$q[[i]], r_bar=criterion_matrix$r[[i]], p=order)
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
        cat("For stock ", names(indices)[j], " the results are as follows: ", max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")
}


#Setting tuning parameters for testing
q     <- 55
r_bar <- 15


#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in 2:(n_ts+1)){
        AR.struc        <- estimate_lrv(data = indices[[i]], q = q, r_bar = r_bar, p=order_results[i-1])
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
result <- multiscale_test(data = matrix(unlist(indices[, 2:(n_ts + 1)]), ncol = n_ts, byrow = FALSE),
                          sigma_vec = sigmahat_vector,
                          alpha = alpha,
                          n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)

source("functions.R")

for (l in seq_len(nrow(result$ijset))){
        i <- result$ijset[l, 1]
        j <- result$ijset[l, 2]
        if (result$stat_pairwise[i, j] > result$quant){
                filename = paste0("plots/", ind_names[i], "_vs_", ind_names[j], ".pdf")
                pdf(filename, width=5.5, height=10.5, paper="special")
                layout(matrix(c(1, 2, 3),ncol=1), widths=c(2.2, 2.2, 2.2),
                       heights=c(1.5, 1.5, 1.8), TRUE)
                
                #Setting the layout of the graphs
                
                par(cex = 1, tck = -0.025)
                par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
                par(oma = c(0.2, 1.5, 2, 0.2)) #Outer margins
                
                 
                plot(x = indices[, 1], y = indices[, i + 1],
                     ylim=c(min(indices[, i + 1], indices[, j + 1]),
                            max(indices[, i + 1], indices[, j + 1])), type="l",
                     col="blue", ylab="", xlab="", mgp=c(1, 0.5, 0))
                lines(x = indices[, 1], y = indices[, j + 1], col="red")
                
                title(main = "(a) observed normalised indices", font.main = 1, line = 0.5)
                legend("topright", inset = 0.02, legend=c(ind_names[i], ind_names[j]),
                       col = c("blue", "red"), lty = 1, cex = 0.95, ncol = 1)
                
                par(mar = c(0.5, 0.5, 3, 0)) #Margins for each plot
                
                #Plotting the smoothed version of the time series that we have
                smoothed_i  <- mapply(nadaraya_watson_smoothing, grid_points,
                                      MoreArgs = list(indices[, i + 1], grid_points, bw = 3.5 / t_len))
                smoothed_j  <- mapply(nadaraya_watson_smoothing, grid_points,
                                      MoreArgs = list(indices[, j + 1], grid_points, bw = 3.5 / t_len))
                
                plot(x = indices[, 1], y = smoothed_i, ylim=c(min(indices[, i + 1], indices[, j + 1]),
                                        max(indices[, i + 1], indices[, j + 1])), type="l",
                     col="black", ylab="", xlab = "", mgp=c(1,0.5,0))
                title(main = "(b) smoothed curves from (a)", font.main = 1, line = 0.5)
                lines(x = indices[, 1], y = smoothed_j, col="red")
                
                par(mar = c(2.7, 0.5, 3, 0)) #Margins for each plot
                gset    <- result$gset_with_values[[l]]
                a_t_set <- subset(gset, test == TRUE, select = c(u, h))
                if (nrow(a_t_set) > 0){
                        p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h) * t_len + 0.5,
                                              'endpoint' = (a_t_set$u + a_t_set$h) * t_len - 0.5, 'values' = 0)
                        p_t_set$values <- (1:nrow(p_t_set))/nrow(p_t_set)
                        
                        #Produce minimal intervals
                        p_t_set2  <- compute_minimal_intervals(p_t_set)
                        
                        plot(NA, xlim=c(0, t_len),  ylim = c(0, 1 + 1 / nrow(p_t_set)), xlab="", mgp=c(2, 0.5, 0), yaxt = "n")
                        title(main = "(c) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
                        #title(xlab = "day", line = 1.7, cex.lab = 0.9)
                        segments(p_t_set2$startpoint, p_t_set2$values, p_t_set2$endpoint, p_t_set2$values, lwd = 2)
                        segments(p_t_set$startpoint, p_t_set$values, p_t_set$endpoint, p_t_set$values, col = "gray")
                } else {
                        #If there are no intervals where the test rejects, we produce empty plots
                        
                        plot(NA, xlim=c(0, tlen),  ylim = c(0, 1), xlab="", ylab = "", mgp=c(2,0.5,0), yaxt = "n")
                        title(main = "(c) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
                        #title(xlab = "days since the hundredth case", line = 1.7, cex.lab = 0.9)
                }
                mtext(paste0("Comparison of ", ind_names[i], " and ", ind_names[j]), side = 3,
                      line = 0, outer = TRUE, font = 1, cex = 1.2)
                dev.off()
                
        }
}
