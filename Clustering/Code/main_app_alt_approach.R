###############################################
#Analysis of covid data - alternative approach#
###############################################
rm(list=ls())

library(tidyr)
library(aweek)
library(dendextend)
require(rworldmap)

#Defining necessary constants
b_bar  <- 2
bw_abs <- 14

#Loading the world coronavirus data
covid_tmp <- read.csv("data/time_series_covid19_confirmed_global.csv",
                      sep = ",", stringsAsFactors = FALSE, na.strings = "",
                      check.names = FALSE)
names(covid_tmp)[names(covid_tmp) == "Country/Region"] <- 'CountryName'
covid_tmp <- covid_tmp[, -c(1, 3, 4)]
new_covid <- aggregate(. ~ CountryName, covid_tmp, sum)
covid     <- gather(new_covid, key = "dateRep", value = "cumcases", 2:442)

rm(covid_tmp, new_covid)

covid$dateRep <- as.Date(covid$dateRep, format = "%m/%d/%y")
covid$cases   <- 0
covid$weekday <- weekdays(covid$dateRep)

covid_list <- list()
for (country in unique(covid$CountryName)){
  cumcases_column <- covid[covid$CountryName == country, "cumcases"]
  time_range      <- length(cumcases_column)
  covid[covid$CountryName == country, "cases"] <- c(0, cumcases_column[2:time_range] - cumcases_column[1:(time_range - 1)])
  tmp <- max(covid[covid$CountryName == country, "cumcases"])
  if (tmp >= 1000){
    #We restrict our attention only to the contries with more than 1000 cases\
    #and only starting from 100th case
    tmp_df <- covid[(covid$CountryName == country & covid$cumcases >= 100),
                    c("dateRep", "cases", "cumcases", "weekday")]
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
  covid_mat[, i] <- covid_list[[country]]$cases[1:t_len]
  i = i + 1
}

#Cleaning the data: there are weird cases when the number of new cases is negative! 
sum(covid_mat < 0)
covid_mat[covid_mat < 0] <- 0

m_hat <- function(vect_u, data_p, grid_p, bw){
  m_hat_vec <- c()
  for (u in vect_u){
    result <- sum((((u - grid_p) / bw < 1) & ((u - grid_p) / bw >= -1)) * data_p)
    norm   <- sum((((u - grid_p) / bw < 1) & ((u - grid_p) / bw >= -1)))
    if (norm == 0){
      m_hat_vec <- c(m_hat_vec, 0)
    } else {
      m_hat_vec <- c(m_hat_vec, result/norm)
    }
  }
  return(m_hat_vec)
}

grid_points <- (1:t_len)/sqrt(t_len)

#Step 2
norm   <- c()
a_vec  <- c()
b_vec  <- c()
c_vec  <- c()
norm_p <- c()
for (k in 1:n_ts) {
  norm <- c(norm, integrate(m_hat, lower = - Inf, upper = Inf,
                            data_p = covid_mat[, k], grid_p = grid_points,
                            bw = bw_abs/sqrt(t_len), subdivisions = 2000)$value)
  
  integrand1 <- function(x) {x * m_hat(x, data_p = covid_mat[, k],
                                       grid_p = grid_points,
                                       bw = bw_abs/sqrt(t_len)) / norm[k]}
  a_vec      <- c(a_vec, integrate(integrand1, lower = - Inf, upper = Inf,
                                   subdivisions = 2000)$value)
  
  integrand2 <- function(x) {x * x * m_hat(x, data_p = covid_mat[, k],
                                           grid_p = grid_points,
                                           bw = bw_abs/sqrt(t_len)) / norm[k]}
  tmp        <- integrate(integrand2, lower = - Inf, upper = Inf,
                          subdivisions = 2000)$value
  b_vec      <- c(b_vec, sqrt(tmp - a_vec[k]^2))
  c_vec      <- c(c_vec, norm[k] / b_vec[k])
  integrand3 <- function(x) {m_hat(a_vec[k] + b_vec[k] * x,
                                   data_p = covid_mat[, k], grid_p = grid_points,
                                   bw = bw_abs/sqrt(t_len)) / c_vec[k]}
  norm_p <- c(norm_p, integrate(integrand3, lower = - Inf, upper = Inf,
                                subdivisions = 2000)$value)
}

# 
# #Matrix with the distances: Step 3
# Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
# 
# for (i in 1:(n_ts - 1)){
#   p_i_star <- function(x) {(m_hat(a_vec[i] + b_vec[i] * x, data_p = covid_mat[, i],
#                                   grid_p = grid_points,
#                                   bw = bw_abs/sqrt(t_len)) / c_vec[i]) / norm_p[i]}
#   for (j in (i + 1):n_ts){
#     p_j_star <- function(x) {(m_hat(a_vec[j] + b_vec[j] * x, data_p = covid_mat[, j],
#                                     grid_p = grid_points,
#                                     bw = bw_abs/sqrt(t_len)) / c_vec[j]) / norm_p[j]}
#     integrand <- function(x) {(sqrt(p_i_star(x)) - sqrt(p_j_star(x)))^2}
#     # if (((i == 1) & ((j == 68) | (j == 117))) | ((i == 2) & (j == 24)) |
#     #     ((i == 9) & (j == 71))| ((i == 11) & (j == 52)) |
#     #     ((i == 14) & ((j == 45) | (j == 73))) | ((i == 15) & (j == 24)) |
#     #     ((i == 19) & ((j == 72) | (j == 89))) | ((i == 20) & (j == 62))){
#     #   tmp <- integrate(integrand, lower = -Inf, upper = Inf,
#     #                    subdivisions=3000)$value
#     # } else {
#       tmp <- integrate(integrand, lower = -Inf, upper = Inf,
#                        subdivisions=3000)$value
#     # }
#     # 
#     Delta_hat[i, j] <- tmp
#     Delta_hat[j, i] <- tmp
#     cat("i = ", i, ", j = ", j, " - success\n")
#   }
# }
# 
# #And now the clustering itself
# colnames(Delta_hat) <- countries
# rownames(Delta_hat) <- countries
# 
#save(Delta_hat, file = "results_alt_approach_28days.RData")
load("results_alt_approach_28days.RData")

n_cl       <- 12

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)
plot(res, cex = 0.8, xlab = "", ylab = "")
rect.hclust(res, k = n_cl, border = 2:(n_cl + 1))


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
               numCats = n_cl,
               mapTitle = "")

pdf(paste0("plots/28days/dendrogram_alt.pdf"), width = 15, height = 6, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot(res, cex = 0.8, xlab = "", ylab = "")
rect.hclust(res, k = n_cl, border = 2:(n_cl + 1))
dev.off()

plotting_grid <- seq(-5, 5, by = 1 / t_len)
subgroups <- cutree(res, n_cl)

for (cl in 1:n_cl){
  countries_cluster <- colnames(Delta_hat)[subgroups == cl]
  pdf(paste0("plots/28days/results_cluster_", cl, "_alt.pdf"), width=7, height=6, paper="special")
  
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(1.5, 0.2, 0.2, 0.2)) #Outer margins
  
  if (length(countries_cluster) == 1){
    ind <- match(countries_cluster, countries)
    m_hat_vec <- m_hat(a_vec[ind] + b_vec[ind] * plotting_grid,
                       covid_mat[, countries_cluster],
                       grid_points,
                       bw = bw_abs/sqrt(t_len)) / (c_vec[ind] * norm_p[ind])
    plot(plotting_grid, m_hat_vec, yaxt = "n",
         ylim = c(0, max(m_hat_vec) + 1), xlab="u",
         ylab = "", mgp = c(2, 0.5, 0), type = "l", col = "red")
    title(main = paste("Cluster", cl), line = 1)
    legend("topleft", inset = 0.02, legend=countries_cluster,
           lty = 1, cex = 0.7, ncol = 1)
  } else {
    inds <- match(countries_cluster, countries)
    i    <- inds[1]
    m_hat_vec <- m_hat(a_vec[i] + b_vec[i] * plotting_grid, covid_mat[, i],
                       grid_points,
                       bw = bw_abs/sqrt(t_len)) / (c_vec[i] * norm_p[i])
    plot(plotting_grid, m_hat_vec,
         ylim = c(0, max(m_hat_vec) + 1), xlab="u", yaxt = "n",
        mgp = c(2, 0.5, 0), type = "l", col = "red")
    for (j in inds[2:length(inds)]){
      m_hat_vec <- m_hat(a_vec[j] + b_vec[j] * plotting_grid, covid_mat[, j],
                         grid_points,
                         bw = bw_abs/sqrt(t_len)) / (c_vec[j] * norm_p[j])
      #cat("Country", country, " - success \n")
      lines(plotting_grid, m_hat_vec)
    }
    title(main = paste("Cluster", cl), line = 1)
    legend("topleft", inset = 0.02, legend = countries_cluster,
           lty = 1, cex = 0.7, ncol = 4)
  }
  dev.off()
}