#####################################################
#Analysis of covid data (provided by John Hopkins U)#
#####################################################
rm(list=ls())

library(tidyr)
library(aweek)
library(dendextend)
library(Rcpp)

require(rworldmap)

Rcpp::sourceCpp("integration.cpp")

#Defining necessary constants
b_bar  <- 2
bw_abs <- 7

#Loading the world coronavirus data
covid_deaths_tmp <- read.csv("data/time_series_covid19_deaths_global.csv", sep = ",", 
                             stringsAsFactors = FALSE, na.strings = "", check.names = FALSE)
covid_cases_tmp  <- read.csv("data/time_series_covid19_confirmed_global.csv", sep = ",", 
                             stringsAsFactors = FALSE, na.strings = "", check.names = FALSE)

names(covid_deaths_tmp)[names(covid_deaths_tmp) == "Country/Region"] <- 'CountryName'
covid_deaths_tmp <- covid_deaths_tmp[, -c(1, 3, 4)]
names(covid_cases_tmp)[names(covid_cases_tmp) == "Country/Region"] <- 'CountryName'
covid_cases_tmp  <- covid_cases_tmp[, -c(1, 3, 4)]


new_covid_cases  <- aggregate(. ~ CountryName, covid_cases_tmp, sum)
new_covid_deaths <- aggregate(. ~ CountryName, covid_deaths_tmp, sum)

covid_cases  <- gather(new_covid_cases, key = "dateRep",
                       value = "cumcases", 2:480)
covid_deaths <- gather(new_covid_deaths, key = "dateRep",
                       value = "cumdeaths", 2:480)

covid <- merge(covid_cases, covid_deaths, by = c("CountryName", "dateRep"))
rm(covid_cases_tmp, covid_deaths_tmp, new_covid_deaths, new_covid_cases,
   covid_cases, covid_deaths)

covid$dateRep <- as.Date(covid$dateRep, format = "%m/%d/%y")
covid$cases   <- 0
covid$deaths  <- 0
covid$weekday <- weekdays(covid$dateRep)

covid <- covid[order(covid$CountryName, covid$dateRep), ]

covid_list <- list()
for (country in unique(covid$CountryName)){
  cumdeaths_column <- covid[covid$CountryName == country, "cumdeaths"]
  cumcases_column  <- covid[covid$CountryName == country, "cumcases"]
  time_range       <- length(cumdeaths_column)
  covid[covid$CountryName == country, "deaths"] <- c(0, cumdeaths_column[2:time_range] - cumdeaths_column[1:(time_range - 1)])
  covid[covid$CountryName == country, "cases"]  <- c(0, cumcases_column[2:time_range] - cumcases_column[1:(time_range - 1)])
  tmp <- max(covid[covid$CountryName == country, "cumcases"])
  tmp_deaths <- max(cumdeaths_column)
  if (tmp >= 1000 & tmp_deaths >= 100 & country != "Cambodia"){
    #We restrict our attention only to the contries with more than 1000 cases and only starting from 100th case
    tmp_df <- covid[(covid$CountryName == country & covid$cumcases >= 100),
                    c("dateRep", "deaths", "cumdeaths", "cases", "cumcases", "weekday")]
    tmp_index <- match("Monday", tmp_df$weekday)
    if (nrow(tmp_df) > 300) {
      covid_list[[country]] <- tmp_df[tmp_index:nrow(tmp_df), ]
    }
  }
}

#Calculate the number of days that we have data for all countries.
#We are not considering CHN = China as it has too long dataset.
t_len     <- min(sapply(covid_list, NROW))
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
  covid_mat[, i] <- covid_list[[country]]$deaths[1:t_len]
  i = i + 1
}

#Cleaning the data: there are weird cases in the dataset when the number of new cases is negative! 
sum(covid_mat < 0)
covid_mat[covid_mat < 0] <- 0

#Nadaraya-Watson estimator
m_hat <- function(vect_u, b, data_p, grid_p, bw){
  m_hat_vec <- c()
  for (u in vect_u){
    result = sum((((grid_p - u * b) / bw < 1) & ((grid_p - u * b) / bw >= -1)) * data_p)
    norm = sum((((grid_p - u * b) / bw < 1) & ((grid_p - u * b) / bw >= -1)))
    m_hat_vec <- c(m_hat_vec, result/norm)
  }
  return(m_hat_vec)
}


#Grid for b and for smoothing
b_grid      <- seq(1, b_bar, by = 0.05)
grid_points <- seq(1/t_len, 1, by = 1/t_len)

#And finally calculating the distance matrix
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
    }
  }
  cat("b = ", b, " - success\n")
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

save(Delta_hat, b_res, file = "results_14days_deaths.RData")
#load("results_14days_deaths.RData")

n_cl       <- 15
delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)

#Plotting world map
covid_map         <- data.frame(countries)
covid_map$cluster <- cutree(res, n_cl)
covid_map[covid_map$countries == 'Czechia', "countries"] <- "Czech Republic"
covid_map[covid_map$countries == 'Taiwan*', "countries"] <- "Taiwan"

covidMap <- joinCountryData2Map(covid_map, 
                                nameJoinColumn="countries", 
                                joinCode="NAME",
                                verbose = TRUE)

mapDevice('x11') #create a world shaped window

#plot the map
mapCountryData(covidMap, 
               nameColumnToPlot='cluster', 
               catMethod='categorical', 
               colourPalette = rainbow(n_cl),
                 #c("#999999", "#E69F00", "#56B4E9", "#009E73",
                #               "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "red", "white"),
               numCats = n_cl,
               mapTitle = "")

pdf(paste0("plots/14days/dendrogram.pdf"), width = 15, height = 6, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot(res, cex = 0.8, xlab = "", ylab = "")
rect.hclust(res, k = n_cl, border = 2:(n_cl + 1))
dev.off()

subgroups <- cutree(res, n_cl)

for (cl in 1:n_cl){
  countries_cluster <- colnames(Delta_hat)[subgroups == cl]
  pdf(paste0("plots/14days/results_cluster_", cl, ".pdf"), width=7, height=6, paper="special")
  
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(1.5, 0.2, 0.2, 0.2)) #Outer margins
  
  if (length(countries_cluster) == 1){
    m_hat_vec <- m_hat(grid_points, b = 1, covid_mat[, countries_cluster],
                       grid_points, bw = bw_abs/t_len)
    norm      <- integrate1_cpp(b = 1, data_points = covid_mat[, countries_cluster],
                                grid_points = grid_points,
                                bw = bw_abs/t_len, subdiv = 2000)$res
    plot((1:t_len) / t_len, m_hat_vec/norm, yaxt = "n",
         ylim = c(0, max(m_hat_vec/norm) + 1), xlab="u",
         ylab = "", mgp = c(2, 0.5, 0), type = "l", col = "red")
    title(main = paste("Cluster", cl), line = 1)
    legend("topleft", inset = 0.02, legend=countries_cluster,
           lty = 1, cex = 0.7, ncol = 1)
  } else {
    b_res_cl     <- b_res[subgroups == cl, subgroups == cl]
    inds         <- which.max(apply(b_res_cl, 1, function(x) sum(x == 1, na.rm = TRUE)))
    repr_country <- rownames(b_res_cl)[inds]
    m_hat_vec    <- m_hat(grid_points, b = 1, covid_mat[, repr_country],
                          grid_points, bw = bw_abs/t_len)
    norm         <- integrate1_cpp(b = 1, data_points = covid_mat[, repr_country],
                                   grid_points = grid_points,
                                   bw = bw_abs/t_len, subdiv = 2000)$res
    #cat("Country", repr_country, ", cluster", cl, " - success \n")
    if (cl == 2) {height <- 8} else {height <- 3} #This should be manually adjusted for nice plots
    plot(grid_points, m_hat_vec/norm,
         ylim = c(0, max(m_hat_vec/norm) + height), xlab="u", yaxt = "n",
         ylab = "m_hat(b * u)", mgp = c(2, 0.5, 0), type = "l", col = "red")
    countries_cluster_1 <- countries_cluster[countries_cluster != repr_country]
    for (country in countries_cluster_1){
      b           <- max(1, b_res_cl[country, repr_country] / b_res_cl[repr_country, country])
      m_hat_vec_1 <- m_hat(grid_points, b = b, covid_mat[, country],
                           grid_points, bw = bw_abs/t_len)
      m_hat_vec_1[(m_hat_vec_1 == 0 | is.nan(m_hat_vec_1))] <- NA
      norm_1      <- integrate1_cpp(b = b, data_points = covid_mat[, country],
                                    grid_points = grid_points,
                                    bw = bw_abs/t_len, subdiv = 2000)$res
      #cat("Country", country, " - success \n")
      lines((1:length(m_hat_vec_1)) / t_len, m_hat_vec_1/(norm_1/(1/b)))
    }
    title(main = paste("Cluster", cl), line = 1)
    legend("topleft", inset = 0.02, legend = countries_cluster,
           lty = 1, cex = 0.7, ncol = 4)
  }
  dev.off()
}

