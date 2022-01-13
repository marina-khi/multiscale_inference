data_load <- function(aligning = TRUE){
  #Loading the world coronavirus data
  covid_tmp <- read.csv("data/time_series_covid19_confirmed_global.csv",
                        sep = ",", stringsAsFactors = FALSE, na.strings = "",
                        check.names = FALSE)
  names(covid_tmp)[names(covid_tmp) == "Country/Region"] <- 'CountryName'
  covid_tmp <- covid_tmp[, -c(1, 3, 4)]
  
  new_covid <- aggregate(. ~ CountryName, covid_tmp, sum)
  covid     <- gather(new_covid, key = "dateRep", value = "cumcases",
                      #2:ncol(new_covid))
                      2:442)

  rm(covid_tmp, new_covid)

  covid$dateRep <- as.Date(covid$dateRep, format = "%m/%d/%y")
  covid$cases   <- 0
  covid$weekday <- weekdays(covid$dateRep)

  covid_list <- list()
  for (country in unique(covid$CountryName)){
    cumcases_column <- covid[covid$CountryName == country, "cumcases"]
    time_range      <- length(cumcases_column)
    covid[covid$CountryName == country,
          "cases"] <- c(0, cumcases_column[2:time_range] - cumcases_column[1:(time_range - 1)])
    tmp <- max(covid[covid$CountryName == country, "cumcases"])
    if (tmp >= 1000){
      #We restrict our attention only to the contries with more than 1000 cases
      #and only starting from 100th case
      tmp_df <- covid[(covid$CountryName == country & covid$cumcases >= 100),
                      c("dateRep", "cases", "cumcases", "weekday")]
      if (aligning){
        tmp_index <- match("Monday", tmp_df$weekday)
      } else {
        tmp_index <- 1
      }
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
  
  #Cleaning the data: there are weird cases in the dataset when the number of new cases is negative! 
  cat("Smaller than zero: ", sum(covid_mat < 0), " objects\n")
  covid_mat[covid_mat < 0] <- 0
  
  return(list(covid_mat = covid_mat, t_len = t_len, n_ts = n_ts, covid_list = covid_list))
}

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

m_hat_standard <- function(vect_u, data_p, grid_p, bw){
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


results_output <- function(res, covid_mat, Delta_hat, b_res, n_cl, countries,
                           path, bw_abs, grid_points, t_len){
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
  
  pdf(paste0(path, "dendrogram.pdf"), width = 15, height = 6, paper = "special")
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
  plot(res, cex = 0.8, xlab = "", ylab = "")
  rect.hclust(res, k = n_cl, border = 2:(n_cl + 1))
  dev.off()
  
  subgroups <- cutree(res, n_cl)
  
  for (cl in 1:n_cl){
    countries_cluster <- colnames(Delta_hat)[subgroups == cl]
    pdf(paste0(path, "results_cluster_", cl, ".pdf"), width=7, height=6, paper="special")
    
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
}


results_output_alt <- function(res, covid_mat, Delta_hat, n_cl, countries,
                               path, bw_abs, grid_points, t_len, a_vec, b_vec,
                               c_vec, norm_p){
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
  
  pdf(paste0(path, "dendrogram_alt.pdf"), width = 15, height = 6, paper = "special")
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
  plot(res, cex = 0.8, xlab = "", ylab = "")
  rect.hclust(res, k = n_cl, border = 2:(n_cl + 1))
  dev.off()
  
  plotting_grid <- seq(-5, 5, by = 1 / t_len)
  subgroups     <- cutree(res, n_cl)
  
  for (cl in 1:n_cl){
    countries_cluster <- colnames(Delta_hat)[subgroups == cl]
    pdf(paste0(path, "results_cluster_", cl, "_alt.pdf"), width=7, height=6, paper="special")
    
    #Setting the layout of the graphs
    par(cex = 1, tck = -0.025)
    par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
    par(oma = c(1.5, 0.2, 0.2, 0.2)) #Outer margins
    
    if (length(countries_cluster) == 1){
      ind <- match(countries_cluster, countries)
      m_hat_vec <- m_hat_standard(a_vec[ind] + b_vec[ind] * plotting_grid,
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
      m_hat_vec <- m_hat_standard(a_vec[i] + b_vec[i] * plotting_grid, covid_mat[, i],
                                  grid_points,
                                  bw = bw_abs/sqrt(t_len)) / (c_vec[i] * norm_p[i])
      plot(plotting_grid, m_hat_vec,
           ylim = c(0, max(m_hat_vec) + 1), xlab="u", yaxt = "n",
           mgp = c(2, 0.5, 0), type = "l", col = "red")
      for (j in inds[2:length(inds)]){
        m_hat_vec <- m_hat_standard(a_vec[j] + b_vec[j] * plotting_grid, covid_mat[, j],
                                    grid_points,
                                    bw = bw_abs/sqrt(t_len)) / (c_vec[j] * norm_p[j])
        lines(plotting_grid, m_hat_vec)
      }
      title(main = paste("Cluster", cl), line = 1)
      legend("topleft", inset = 0.02, legend = countries_cluster,
             lty = 1, cex = 0.7, ncol = 4)
    }
    dev.off()
  }
}