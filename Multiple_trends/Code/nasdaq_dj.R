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

t_len  <- nrow(indices)
n_ts   <- ncol(indices) - 1
ind_names <- c("dow-jones", "nasdaq")
colSums(is.na(indices))

indices[, 2:(n_ts + 1)] <- scale(indices[, 2:(n_ts + 1)], scale = TRUE)
indices <- na.omit(indices)#Deleting the rows with ommitted variables
indices[, "date"] <- as.Date(as.character(indices[, "date"]), format = "%Y-%m-%d")

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
