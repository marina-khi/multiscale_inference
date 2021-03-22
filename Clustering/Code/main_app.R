########################
#Analysis of covid data#
########################
rm(list=ls())

library(tidyr)
library(aweek)
library(dendextend)
library(Rcpp)

require(rworldmap)

Rcpp::sourceCpp("example.cpp")

#Defining necessary constants
b_bar  <- 2
bw_abs <- 3.5

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
#covid$cumdeaths       <- 0

covid_list <- list()
for (country in unique(covid$countryterritoryCode)){
  covid[covid$countryterritoryCode == country, "cumcases"]  <- cumsum(covid[covid$countryterritoryCode == country, "cases"])
#  covid[covid$countryterritoryCode == country, "cumdeaths"] <- cumsum(covid[covid$countryterritoryCode == country, "deaths"])
  tmp <- max(covid[covid$countryterritoryCode == country, "cumcases"])
  if (tmp >= 20000){
    tmp_df <- covid[(covid$countryterritoryCode == country & covid$cumcases >= 100),
                        c("dateRep", "cases", "cumcases", "weekday")]
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
  covid_mat[, i] <- covid_list[[country]]$cases[1:t_len]
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
b_grid      <- seq(1, b_bar, by = 0.01)
grid_points <- seq(1/t_len, 1, by = 1/t_len)


Delta_hat_tmp <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
b_res         <- matrix(data = rep(NA, n_ts * n_ts), nrow = n_ts, ncol = n_ts)

for (b in b_grid){
  norm_b <- c()
  norm   <- c()
  for (k in 1:n_ts){
    norm_b <- c(norm_b, integrate1_cpp(b = b, data_points = covid_mat[, k],
                                       grid_points = grid_points,
                                       bw = bw_abs/t_len, subdiv = 2000)$res)
    norm <- c(norm, integrate1_cpp(b = 1.0, data_points = covid_mat[, k],
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
        Delta_hat_tmp[i, j] <- delta_ij
        Delta_hat_tmp[j, i] <- delta_ji
        b_res[i, j] <- 1
        b_res[j, i] <- 1
      } else {
        if (delta_ij < Delta_hat_tmp[i, j]) {
          Delta_hat_tmp[i, j] <- delta_ij
          b_res[i, j] <- b
          b_res[j, i] <- 1
        } 
        if (delta_ji < Delta_hat_tmp[j, i]) {
          Delta_hat_tmp[j, i] <- delta_ji
          b_res[j, i] <- b
          b_res[i, j] <- 1          
        }
      }
      #cat("b = ", b, ", Delta_hat = ", Delta_hat[i, j], "\n")
      }
  }  
}

#Delta_hat_tmp was a temporary non-symmetrical matrix,
#for the distance matrix we need a symmetrical one
Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
for (i in 1:(n_ts - 1)){
  for (j in (i + 1):n_ts){
    Delta_hat[i, j] <- min(Delta_hat_tmp[i, j], Delta_hat_tmp[j, i])
    Delta_hat[j, i] <- Delta_hat[i, j]
  }
}

colnames(Delta_hat) <- countries
rownames(Delta_hat) <- countries

colnames(b_res) <- countries
rownames(b_res) <- countries

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)

#Plotting world map
covid_map         <- data.frame(countries)
covid_map$cluster <- cutree(res, 6)
covid_map[covid_map$countries == 'XKX', "countries"] <- "KOS"

covidMap <- joinCountryData2Map(covid_map, 
                                nameJoinColumn="countries", 
                                joinCode="ISO3",
                                verbose = TRUE)

mapDevice('x11') #create a world shaped window

#plot the map
mapCountryData(covidMap, 
               nameColumnToPlot='cluster', 
               catMethod='categorical', 
               colourPalette = 2:7,
               numCats=6)

pdf("plots/dendrogram.pdf", width = 15, height = 6, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot(res, cex = 0.8, xlab = "", ylab = "")
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
    plot((1:t_len) / t_len, m_hat_vec,
         ylim = c(0, max(m_hat_vec) + 10), xlab="u",
         ylab = "", mgp = c(2, 0.5, 0), type = "l")
    title(main = paste("Representative of cluster", cl), line = 1)
  } else {
    b_res_cl     <- b_res[subgroups == cl, subgroups == cl]
    inds         <- which.max(apply(b_res_cl, 1, function(x) sum(x == 1, na.rm = TRUE)))
    #inds        <- which.min(rowSums(b_res_cl, na.rm = TRUE))
    repr_country <- rownames(b_res_cl)[inds]
    m_hat_vec    <- m_hat(grid_points, b = 1, covid_mat[, repr_country],
                          grid_points, bw = bw_abs/t_len)
    norm         <- integrate1_cpp(b = 1, data_points = covid_mat[, repr_country],
                                   grid_points = grid_points,
                                   bw = bw_abs/t_len, subdiv = 2000)$res
    #cat("Country", repr_country, " - success \n")
    plot(grid_points, m_hat_vec/norm,
         ylim = c(0, max(m_hat_vec/norm) + 10), xlab="u",
         ylab = "m_hat(b * u)", mgp = c(2, 0.5, 0), type = "l")
    countries_cluster_1 <- countries_cluster[countries_cluster != repr_country]
    for (country in countries_cluster_1){
      b <- b_res_cl[country, repr_country] / b_res_cl[repr_country, country]
      m_hat_vec_1 <- m_hat(grid_points, b = b, covid_mat[, country],
                           grid_points, bw = bw_abs/t_len)
      m_hat_vec_1[(m_hat_vec_1 == 0 | is.nan(m_hat_vec_1))] <- NA
      norm_1         <- integrate1_cpp(b = b, data_points = covid_mat[, country],
                                     grid_points = grid_points,
                                     bw = bw_abs/t_len, subdiv = 2000)$res
      #cat("Country", country, " - success \n")
      lines((1:length(m_hat_vec_1)) / t_len, m_hat_vec_1/norm_1)
    }
    title(main = paste("Representatives of cluster", cl), line = 1)
  }
  legend("topright", inset = 0.02, legend=countries_cluster,
         lty = 1, cex = 0.7, ncol = 1)
  dev.off()
}