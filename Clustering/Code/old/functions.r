nadaraya_watson_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    for (i in 1:T_size){
      if (abs((grid_p[i] - u) / bw) <= 1){
        result = result + data_p[i]
        norm = norm + 1
      }
    }
    return(result/norm)
  }
}

integral_points <- seq(1/t_len, 1, by = 0.01/t_len)
tmp_grid  <- seq(1 / (2 * t_len), 1, by = 2 / (2 * t_len))


m_hat_vec <- m_hat(grid_points, b = 1, covid_mat[, 5], grid_points, bw = bw_abs/t_len)
plot(grid_points, covid_mat[, 5], type = "l")
lines(grid_points, m_hat_vec, col = "red")

m_hat_vec2 <- m_hat_2(grid_points, b = 1, covid_mat[, 5], grid_points, bw = bw_abs/t_len)
plot(grid_points, covid_mat[, 5], type = "l")
lines(grid_points, m_hat_vec2, col = "red")

m_hat_vec <- m_hat(integral_points, b = 1, covid_mat[, 2], grid_points, bw = bw_abs/t_len)
plot(grid_points, covid_mat[, 2], type = "l")
lines(integral_points, m_hat_vec, col = "blue")

m_hat_vec <- m_hat(integral_points, b = 1, covid_mat[, 80], grid_points, bw = bw_abs/t_len)
plot(grid_points, covid_mat[, 80], type = "l")
lines(integral_points, m_hat_vec, col = "green")

integrand <- function(vect_u, b, data_points_i, data_points_j, 
                      norm_b, norm, grid_points, bw) {
  tmp <- m_hat(vect_u, b, data_points_i, grid_points, bw)/(norm_b * 1/b) - m_hat(vect_u, b = 1, data_points_j, grid_points, bw) / (norm * 1/b)
  return(tmp^2)
}

integ_r <- c()
for (k in 1:n_ts){
  true_val <- integrate(m_hat, lower = 0, upper = 1/b, b = b,
                        data_p = covid_mat[, k], grid_p = grid_points,
                        bw = bw_abs/t_len, subdivisions=2000)$value
  #res <- integrate_test(b = b, data_points = covid_mat[, k], grid_points = grid_points, bw = bw_abs/t_len,
  #                      true_val = true_val)
  #cat(res$approximate - true_val, "\n")
  integ_r <- c(integ_r, true_val)
}


tmp_b <- m_hat(tmp_grid, b = b, data_p = covid_mat[, k], grid_points, bw = bw_abs/t_len)
norm_b <- c(norm_b, sum(tmp_b * (1 / t_len) * (1 / b), na.rm = TRUE))
norm_b <- c(norm_b, integrate(m_hat, lower = 0, upper = 1/b, b = b,
                              data_p = covid_mat[, k], grid_p = grid_points,
                              bw = bw_abs/t_len, subdivisions=1000)$value)
len <- floor(1 / b * t_len)
tmp <- m_hat(tmp_grid, b = 1, data_p = covid_mat[, k], grid_points, bw = bw_abs/t_len)
norm <- c(norm, sum(tmp[1:len] * (1 / t_len) * (1 / b), na.rm = TRUE))
norm <- c(norm, integrate(m_hat, lower = 0, upper = 1/b, b = 1,
                          data_p = covid_mat[, k], grid_p = grid_points, 
                          bw = bw_abs/t_len, subdivisions=1000)$value)

library("rjson")
covid_russia_tmp <- fromJSON(file = "data/data.json")
covid_russia_matrix <- matrix(data = rep(0, 286 * 85), ncol = 85, nrow = 286)
names <- c()
for (i in 1:85){
  covid_russia_matrix[, i] <- covid_russia_tmp[['data']][[i]]$confirmed
  names <- c(names, covid_russia_tmp[['data']][[i]]$name)
}
colnames(covid_russia_matrix) <- names

dates <- seq(from = as.Date(covid_russia_tmp$startDate, format = "%m/%d/%Y"), by = 1, length.out = 286)
covid_russia <- data.frame(dateRep = dates, cases_ru = rowSums(covid_russia_matrix))
covid_russia <- merge(covid_russia, covid_list$Russia, by = "dateRep")
