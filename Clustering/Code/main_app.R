########################
#Analysis of covid data#
########################
rm(list=ls())

library(tidyr)
library(multiscale)
library(xtable)
library(aweek)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

library(dendextend)
library(tictoc)
library(Rcpp)

Rcpp::sourceCpp("example.cpp")

#Defining necessary constants
b_bar <- 1.05
bw_abs <- 7

#Loading the world coronavirus data
covid         <- read.csv("data/covid.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "")
covid         <- covid[!is.na(covid$countryterritoryCode), ]
covid$dateRep <- as.Date(covid$dateRep, format = "%d/%m/%Y")
covid         <- complete(covid, dateRep = seq.Date(min(dateRep), max(dateRep), by='day'),
                          countryterritoryCode, fill = list(cases = 0, deaths = 0))

#Now we "normalize" the data counting only the countries with more than 1000 deaths overall
#and taking the day of 100th case as the starting point
covid$weekday         <- weekdays(covid$dateRep)
covid$cumcases        <- 0
covid$cumdeaths       <- 0

covid_list <- list()
for (country in unique(covid$countryterritoryCode)){
  covid[covid$countryterritoryCode == country, "cumcases"]  <- cumsum(covid[covid$countryterritoryCode == country, "cases"])
  covid[covid$countryterritoryCode == country, "cumdeaths"] <- cumsum(covid[covid$countryterritoryCode == country, "deaths"])
  tmp <- max(covid[covid$countryterritoryCode == country, "cumdeaths"])
  if (tmp >= 1000){
    tmp_df <- covid[(covid$countryterritoryCode == country & covid$cumcases >= 100),
                        c("dateRep", "cases", "deaths", "cumcases", "cumdeaths", "weekday")]
    tmp_index <- match("Monday", tmp_df$weekday)
    #tmp_index = 1 #If we do not want to normalize by Mondays
    covid_list[[country]] <- tmp_df[tmp_index:nrow(tmp_df), ]
  }
}


#Calculate the number of days that we have data for all fivecountries.
#We are not considering CHN = China as it has too long dataset.
t_len     <- min(sapply(covid_list[names(covid_list) != "CHN"], NROW))
#t_len     <- 150 #We consider the first five months of the pandemic
countries <- names(covid_list)
dates     <- unique(covid$dateRep)
n_ts      <- length(covid_list) #Number of time series

#In order not to work with lists, we construct a matrix
#with number of cases for all countries and.
#It is not necessary, but it is more convenient to work with.
covid_mat           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(covid_mat) <- countries

i = 1
for (country in countries) {
  covid_mat[, i]        <- covid_list[[country]]$cases[1:t_len]
  i = i + 1
}

#Cleaning the data: there are weird cases in the dataset when the number of new cases is negative! 
sum(covid_mat < 0)
covid_mat[covid_mat < 0] <- 0


m_hat <- function(vect_u, b, data_p, grid_p, bw){
  m_hat_vec <- c()
  for (u in vect_u){
    result = sum((abs((grid_p - u * b) / bw) <= 1) * data_p)
    #norm = sum((abs((grid_p - u * b) / bw) <= 1))
    norm = min(floor((u * b + bw) * t_len), t_len) - max(ceiling((u * b - bw) * t_len), 1) + 1
    m_hat_vec <- c(m_hat_vec, result/norm)
  }
  return(m_hat_vec)
}

#Grid for b and for smoothing
b_grid      <- seq(1, b_bar, by = 0.05)
grid_points <- seq(1/t_len, 1, by = 1/t_len)


Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
b_res     <- matrix(data = rep(NA, n_ts * n_ts), nrow = n_ts, ncol = n_ts)

tic("Starting")
for (b in b_grid){
  norm_b <- c()
  norm   <- c()
  for (k in 1:n_ts){
    norm_b <- c(norm_b, integrate1_cpp(b = b, data_points = covid_mat[, k],
                                       grid_points = grid_points,
                                       bw = bw_abs/t_len, subdiv = 2000)$res)
    norm <- c(norm_b, integrate1_cpp(b = 1.0, data_points = covid_mat[, k],
                                     grid_points = grid_points,
                                     bw = bw_abs/t_len, subdiv = 2000)$res)
  }
  for (i in 1:(n_ts - 1)){
    for (j in (i + 1):n_ts){
      delta_ij <- 1/b * integrate2_cpp(b = b, data_points_1 = covid_mat[, i],
                                  data_points_2 = covid_mat[, j],
                                  norm_1 = norm_b[i], norm_2 = norm[j],
                                  grid_points = grid_points, bw = bw_abs/t_len,
                                  subdiv=2000)$res
      delta_ji <- 1/b * integrate2_cpp(b = b, data_points_1 = covid_mat[, j],
                                  data_points_2 = covid_mat[, i],
                                  norm_1 = norm_b[j], norm_2 = norm[i],
                                  grid_points = grid_points, bw = bw_abs/t_len,
                                  subdiv=2000)$res
      if (b == 1) {
        Delta_hat[i, j] <- min(delta_ij, delta_ji)
        b_res[i, j]     <- 1
        b_res[j, i]     <- 1
      } else {
        if (min(delta_ij, delta_ji) < Delta_hat[i, j]) {
          Delta_hat[i, j] <- min(delta_ij, delta_ji)
          if (delta_ij <= delta_ji) {
            b_res[i, j] <- b
            b_res[j, i] <- 1
          } else {
            b_res[j, i] <- b
            b_res[i, j] <- 1
          }
        }
      }
      Delta_hat[j, i] <- Delta_hat[i, j]
      #cat("b = ", b, ", Delta_hat = ", Delta_hat[i, j], "\n")
      }
  }  
}
toc()

colnames(Delta_hat) <- countries
rownames(Delta_hat) <- countries

colnames(b_res) <- countries
rownames(b_res) <- countries

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)
pdf("plots/dendrogram.pdf", width=15, height=6, paper="special")
plot(res, cex = 0.8)
rect.hclust(res, k = 6, border = 2:7)
dev.off()

subgroups <- cutree(res, 6)

for (cl in 1:6){
  countries_cluster <- colnames(Delta_hat)[subgroups == cl]
  pdf(paste0("plots/results_cluster_", cl, ".pdf"), width=7, height=6, paper="special")

  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
  
  if (length(countries_cluster) == 1){
    m_hat_vec <- m_hat(grid_points, b = 1, covid_mat[, countries_cluster],
                       grid_points, bw = bw_abs/t_len)
    norm <- integrate(m_hat, lower = 0, upper = 1, b = 1,
                      data_p = covid_mat[, countries_cluster], grid_p = grid_points,
                      bw = bw_abs/t_len, subdivisions=2000)$value
    plot((1:t_len) / t_len, m_hat_vec,
         ylim = c(0, max(m_hat_vec) + 10), xlab="u",
         ylab = "", mgp = c(2, 0.5, 0), type = "l")
    title(main = paste("Representative of cluster", cl), line = 1)
    
  } else {
    b_res_cl <- b_res[subgroups == cl, subgroups == cl]
    #colnames(b_res_cl) <- countries_cluster
    #rownames(b_res_cl) <- countries_cluster
    inds               <- which.min(rowSums(b_res_cl, na.rm = TRUE)) #dim(b_res_cl))
    repr_country       <- rownames(b_res_cl)[inds]
    m_hat_vec <- m_hat(grid_points, b = 1, covid_mat[, repr_country],
                       grid_points, bw = bw_abs/t_len)
    norm <- integrate(m_hat, lower = 0, upper = 1, b = 1,
              data_p = covid_mat[, repr_country], grid_p = grid_points,
              bw = bw_abs/t_len, subdivisions=2000)$value
    #cat("Country", repr_country, " - success \n")
    plot((1:t_len) / t_len, m_hat_vec/norm,
         ylim = c(0, max(m_hat_vec/norm) + 10), xlab="u",
         ylab = "m_hat(b * u)", mgp = c(2, 0.5, 0), type = "l")
    countries_cluster_1 <- countries_cluster[countries_cluster != repr_country]
    for (country in countries_cluster_1){
      b <- b_res_cl[country, repr_country] / b_res_cl[repr_country, country]
      m_hat_vec_1 <- m_hat(grid_points, b = b, covid_mat[, country],
                           grid_points, bw = bw_abs/t_len)
      m_hat_vec_1[(m_hat_vec_1 == 0 | is.nan(m_hat_vec_1))] <- NA
      norm_1 <- integrate(m_hat, lower = 0, upper = 1 / b, b = b,
                          data_p = covid_mat[, country], grid_p = grid_points,
                          bw = bw_abs/t_len, subdivisions=2000)$value
      #cat("Country", country, " - success \n")
      lines((1:length(m_hat_vec_1)) / t_len, m_hat_vec_1/norm_1)
    }
    title(main = paste("Representatives of cluster", cl), line = 1)
  }
  legend("topright", inset = 0.02, legend=countries_cluster,
         lty = 1, cex = 0.7, ncol = 1)
  dev.off()
}