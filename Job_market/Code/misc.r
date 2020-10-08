pdf("plots/Germany_and_Italy.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(covid_mat[, 1], ylim=c(min(covid_mat[, 1], covid_mat[, 5]), max(covid_mat[, 1], covid_mat[, 5])), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
lines(covid_mat[, 5], col="#604c38")
title(ylab = "New daily cases", line = 1.7, cex.lab = 0.9)
title(xlab = "days since first Monday after reaching 100 cases", line = 1.7, cex.lab = 0.9)
legend("topright", inset = 0.02, legend=c("Germany", "Italy"),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()

pdf("plots/Germany_and_Italy_1.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins


plot(covid_mat[, 1], ylim=c(min(covid_mat[, 1], covid_mat[, 5]), max(covid_mat[, 1], covid_mat[, 5])), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
rect(xleft=18, xright = 27, ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
lines(covid_mat[, 5], col="#604c38")
lines(covid_mat[, 1], col="#EB811B")
segments(18, par("usr")[3], 27, par("usr")[3], col="#604c38", lwd = 3)
title(ylab = "New daily cases", line = 1.7, cex.lab = 0.9)
title(xlab = "days since first Monday after reaching 100 cases", line = 1.7, cex.lab = 0.9)
legend("topright", inset = 0.02, legend=c("Germany", "Italy"),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()


pdf("plots/Germany_and_Italy_2.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(covid_mat[, 1], ylim=c(min(covid_mat[, 1], covid_mat[, 5]), max(covid_mat[, 1], covid_mat[, 5])), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
rect(xleft=18, xright = 27, ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
rect(xleft=41, xright = 88, ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
lines(covid_mat[, 5], col="#604c38")
lines(covid_mat[, 1], col="#EB811B")
segments(18, par("usr")[3], 27, par("usr")[3], col="#604c38", lwd = 3)
segments(41, par("usr")[3], 88, par("usr")[3], col="#604c38", lwd = 3)
title(ylab = "New daily cases", line = 1.7, cex.lab = 0.9)
title(xlab = "days since first Monday after reaching 100 cases", line = 1.7, cex.lab = 0.9)
legend("topright", inset = 0.02, legend=c("Germany", "Italy"),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()


# #Plots of lambda functions
# n_sim     <- 5000               # number of simulation runs for power and size
# sim_runs  <- 5000               # number of simulation runs to produce critical values
# alpha_vec <- c(0.01, 0.05, 0.1) # different significance levels
# n_ts_vec  <- c(5, 10, 50)       # different number of time series
# t_len_vec <- c(100, 250, 500)   # different time series lengths
# sigma_vec <- c(15, 10, 20)      # different overdispersion parameter
# 
# number_of_cols <- length(n_ts_vec) * length(alpha_vec) #Needed for the output
# 
# #As the mean function in the size simulations, we take the following function:
# #lambda(u) = 5000 * exp(-(10 * u - 3) ^ 2 / 2) + 1000 for 0 <= u <= 1.
# #Here is the plot of this function:
# 
# lambda_vec <- lambda_fct((1:100) / 100)
# 
# pdf(paste0("plots_new/lambda_fct.pdf"), width=5, height=3, paper="special")
# par(mar = c(3, 2, 2, 0)) #Margins for each plot
# par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
# plot((1:100) / 100, lambda_vec,  ylim = c(0, max(lambda_vec) + 100), xlab="u",
#      ylab = "", mgp=c(2,0.5,0), type = "l")
# title(main = expression(Plot ~ of ~ the ~ "function" ~ lambda), line = 1)
# dev.off()
# 
# 
# lambda_vec_1 <- lambda_fct((1:100) / 100, c = 1000, height = 6000, position = 10)
# lambda_vec   <- lambda_fct((1:100) / 100, c = 1000, height = 5000, position = 10)
# 
# pdf(paste0("plots/lambda_fcts_height.pdf"), width=4.5, height=3, paper="special")
# par(mar = c(3, 2, 2, 0)) #Margins for each plot
# par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
# par(mgp = c(3, 0.5, 0))
# plot((1:100) / 100, lambda_vec_1,  ylim = c(0, max(lambda_vec_1, lambda_vec) + 100),
#      xlab="", ylab = "", type = "l", col = "#604c38")
# lines((1:100) / 100, lambda_vec, type = "l", col = "#EB811B")
# title(main = expression(Plot ~ of ~ the ~ "functions" ~ lambda[1] ~ and ~ lambda), line = 1)
# title(xlab="u", line=1.5)
# legend("topright", inset = 0.02, legend=c(expression(lambda[1](u) ~" "), expression(lambda(u) ~" ")),
#        col = c("#604c38", "#EB811B"), lty = 1, cex = 0.95, ncol = 1)
# dev.off()
# 
# lambda_vec_1 <- lambda_fct((1:100) / 100, c = 1000, height = 5000, position = 9)
# lambda_vec   <- lambda_fct((1:100) / 100, c = 1000, height = 5000, position = 10)
# 
# pdf(paste0("plots/lambda_fcts_shift.pdf"), width=4.5, height=3, paper="special")
# par(mar = c(3, 2, 2, 0)) #Margins for each plot
# par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
# par(mgp = c(3, 0.5, 0))
# plot((1:100) / 100, lambda_vec_1,  ylim = c(0, max(lambda_vec_1, lambda_vec) + 100),
#      xlab="", ylab = "", type = "l", col = "#604c38")
# lines((1:100) / 100, lambda_vec, type = "l", col = "#EB811B")
# title(main = expression(Plot ~ of ~ the ~ "functions" ~ lambda[1] ~ and ~ lambda), line = 1)
# title(xlab="u", line=1.5)
# legend("topright", inset = 0.02, legend=c(expression(lambda[1](u) ~" "), expression(lambda(u) ~" ")),
#        col = c("#604c38", "#EB811B"), lty = 1, cex = 0.95, ncol = 1)
# dev.off()

