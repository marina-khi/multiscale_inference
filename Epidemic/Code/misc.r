pdf("plots_new/Germany_and_Italy.pdf", width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(2.5, 2.5, 0, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(covid_mat[, 1], ylim=c(min(covid_mat[, 1], covid_mat[, 5]), max(covid_mat[, 1], covid_mat[, 5])), type="l",
     col="#604c38", ylab="Number of cases", xlab="Number of days since 100th case", mgp=c(1.3, 0.3, 0), cex.lab = 0.8, cex.axis = 0.5)
#title(ylab="Number of cases", line=3, cex.lab=1.2)
#title(xlab="Number of days since 100th case", line=3, cex.lab=1.2)
lines(covid_mat[, 5], col="#EB811B")
legend("topright", inset = 0.02, legend=c("Germany", "Italy"),
       col = c("#604c38", "#EB811B"), lty = 1, cex = 0.8, ncol = 1)

dev.off()
