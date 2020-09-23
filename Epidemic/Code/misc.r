pdf("plots_new/Germany_and_Italy.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 0.5, 0)) #Margins for each plot
par(oma = c(2, 1.5, 0.2, 0.2)) #Outer margins

plot(covid_mat[, 1], ylim=c(min(covid_mat[, 1], covid_mat[, 5]), max(covid_mat[, 1], covid_mat[, 5])), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
lines(covid_mat[, 5], col="#604c38")
legend("topright", inset = 0.02, legend=c("Germany", "Italy"),
       col = c("#604c38", "#EB811B"), lty = 1, cex = 0.95, ncol = 1)
dev.off()