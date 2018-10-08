# Epanechnikov kernel function, which is defined as f(x) = 3/4(1-x^2)
# for |x|<1 and 0 elsewhere
epanechnikov_kernel <- function(x)
{
  if (abs(x)<1)
  {
    result = 3/4 * (1 - x*x)
  } else {
    result = 0
  }
  return(result)
}

# Biweight kernel function, which is defined as f(x) = 15/16(1-u^2)^2
# for |u|<1 and 0 elsewhere
biweight_kernel <- function(x)
{
  if (abs(x)<1)
  {
    result = (15/16) * (1 - x*x)^2
  } else {
    result = 0
  }
  return(result)
}

#Functions needed for local linear smoothing
s_t_1 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

s_t_2 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

s_t_0 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}


#Local Linear estimator using the epanechnikov kernel. 
local_linear_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    s_t_2_value = s_t_2(u, bw, T_size)
    s_t_1_value = s_t_1(u, bw, T_size) 
    for (i in 1:T_size){
      k = (s_t_2_value - s_t_1_value * ((grid_p[i] - u) / bw)) * epanechnikov_kernel((grid_p[i] - u) / bw)
      result = result + k * data_p[i]
      norm = norm + k 
    }
    return(result/norm)
  }
}

#Nadaraya-Watson estimator using the epanechnikov kernel. 
epanechnikov_smoothing <-function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result = 0
    norm = 0
    for (i in 1:length(data_p)){
      result = result + epanechnikov_kernel((u - grid_p[i])/bw) * data_p[i]
      norm = norm + epanechnikov_kernel((u - grid_p[i])/bw) 
    }
    return(result/norm)
  }
}


#Additive correction tern \lambda(h) that depends only on the bandwidth h
lambda <- function(h)
{
  result = tryCatch(sqrt(2*log(1/(2*h))), warning = function(w) print("h is exceeding h_max"))
  return(result)
}


#This functions finds minimal intervals as described in Duembgen(2002)
choosing_minimal_intervals <- function(dataset){
  set_cardinality <- nrow(dataset) 
  if (set_cardinality > 1) {
    dataset <- dataset[order(dataset$startpoint, -dataset$endpoint),] #Ordering such that we don't need to check previous intervals
    rownames(dataset) <- 1:nrow(dataset) #restoring the indices after ordering
    dataset[['contains']] <- numeric(set_cardinality)
    for (i in 1:(set_cardinality-1)){
      for (j in (i+1):set_cardinality){
        if ((dataset$startpoint[i] <= dataset$startpoint[j]) & (dataset$endpoint[i] >= dataset$endpoint[j])) {
          dataset[['contains']][i] <- 1 #We are marking all the intervals that contain at least one another interval
          break
        }
      }
    }
    p_t_set <- subset(dataset, contains == 0, select = c(startpoint, endpoint, values))#Subsetting everything not marked
  }
}


#Creating g_t_set over which we are taking the maximum (from Section 2.1)
creating_g_set <- function(T, kernel_method){
  u <- seq(from = 5/T, to = 1, by = 5/T)
  h <- seq(from = 3/T, to = 1/4+3/T, by = 5/T)

  #This is for full grid
  #u <- seq(from = 1/T, to = 1, by = 1/T)
  #h <- seq(from = 3/T, to = 1/4+1/T, by = 1/T)
  
  g_t_set_temp                  <- expand.grid(u = u, h = h) #Creating a dataframe with all possible combination of u and h
  g_t_set_temp$values           <-numeric(nrow(g_t_set_temp)) # Setting the values of the statistic to be zero
  g_t_set_temp$values_with_sign <-numeric(nrow(g_t_set_temp)) # Setting the values of the statistic to be zero

  if (kernel_method == "nw"){
    g_t_set <- subset(g_t_set_temp, u - h >= 0 & u + h <= 1, select = c(u, h, values, values_with_sign)) #Subsetting u and h such that [u-h, u+h] lies in [0,1]
  } else if (kernel_method == "ll"){
    g_t_set <- subset(g_t_set_temp, u >= 0 & u <= 1, select = c(u, h, values, values_with_sign)) #Subsetting u and h such that [u-h, u+h] lies in [0,1]
  } else {
    print('Given method is currently not supported')
  }
  
  g_t_set$lambda <- lambda(g_t_set[['h']]) #Calculating the lambda(h) in order to speed up the function psistar_statistic
  return(g_t_set)
}


#If we have already calculated quantiles and stored them in a file 'distribution.RData'
#then no need to calculate them once more, we just load them from this file.
#Ohterwise simulate the \Psi^star statistic 1000 times in order to calculate the quantiles
calculating_gaussian_quantile <- function(T, g_t_set, test_problem, kernel_ind, alpha = 0.05){
  filename = paste0("Shape/distribution/distr_T_", T,"_testing_", test_problem, ".RData")
  if(!file.exists(filename)) {
    gaussian_statistic_distribution <- replicate(1000, {
      z = rnorm(T, 0, 1)
      psistar_statistic(z, g_t_set, kernel_ind, 1)
    })
    save(gaussian_statistic_distribution, file = filename)
  } else {
    load(filename)
  }
  #Calculate the quantiles for gaussian statistic defined in the previous step
  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)
  return(gaussian_quantile)
}


#Everything as in previous function but using local linear estimate
calculating_gaussian_quantile_ll <- function(T, g_t_set, test_problem, kernel_ind, alpha = 0.05){
  filename = paste0("Shape/distribution/distr_T_", T,"_testing_", test_problem, "_type_ll.RData")
  if(!file.exists(filename)) {
    gaussian_statistic_distribution <- replicate(1000, {
      z = rnorm(T, 0, 1)
      psistar_statistic_ll(z, g_t_set, kernel_ind, 1)
    })
    save(gaussian_statistic_distribution, file = filename)
  } else {
    load(filename)
  }
  #Calculate the quantiles for gaussian statistic defined in the previous step
  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)
  return(gaussian_quantile)
}


#Create a matrix (for size and power table for example) and write them in the tex file
creating_matrix_and_texing <- function(vect, vect_T, vect_alpha, filename){
  matrix_ <- matrix(vect, nrow = length(vect_T), ncol = length(vect_alpha), byrow = TRUE)
  rownames(matrix_) <- vect_T
  colnames(matrix_) <- vect_alpha
  
  addtorow <- list()
  addtorow$pos <- list(0, 0)
  addtorow$command <- c("& \\multicolumn{3}{c}{nominal size $\\alpha$} \\\\\n",
                            "$T$ & 0.01 & 0.05 & 0.1 \\\\\n") 
  
  print.xtable(xtable(matrix_, digits = c(3), align = "cccc"), type="latex",  file=filename, add.to.row = addtorow, include.colnames = FALSE)
}

#Create a matrix (for size and power table for example) and write them in the tex file for comparison with SiZer
creating_matrix_and_texing_for_SiZer <- function(vect, vect_T, vect_alpha, filename){
  matrix_ <- matrix(vect, nrow = length(vect_T), ncol = 2 * length(vect_alpha), byrow = TRUE)
  rownames(matrix_) <- vect_T

  addtorow <- list()
  addtorow$pos <- list(0, 0, 0)
  addtorow$command <- c("& \\multicolumn{6}{c}{nominal size $\\alpha$} \\\\\n",
                        "$T$ & \\multicolumn{2}{c}{$0.01$} & \\multicolumn{2}{c}{$0.05$} & \\multicolumn{2}{c}{$0.10$} \\\\\n",
                        " & MT & SiZer & MT & SiZer & MT & SiZer \\\\\n") 
  
  print.xtable(xtable(matrix_, digits = c(3), align = "ccccccc"), type="latex",  file=filename, add.to.row = addtorow, include.colnames = FALSE)
}

#Estimate autocovariance function \gamma_q(l) for a given time series y_data by a sample autocovariance
sample_autocovariance <- function(l, y_data, q){
  # computes autocovariances at lags 0 to p 
  # for the ell-th differences of y_data

  if (l%%1==0)
  {
    T_size = length(y_data)
    y_data_diff = c(rep(0, q), diff(y_data, q))
    if (abs(l) >= T_size - q - 1) {
      print("Cannot estimate autocovariance from the data: sample size is too small")
    } else {
      result = sum(y_data_diff[(q + 1):(T_size - abs(l))] * y_data_diff[(q + abs(l) + 1):T_size])/(T_size - abs(l) - q)
    }
  } else {
    print('Check the input: l is not integer')
  }
  return(result)
}

AR_coefficients <- function(y_data, L1, L2, correction, p){
  pos <- 0
  a_mat <- matrix(0,ncol=L2-L1+1,nrow=p)
  for(q in L1:L2) {
    gamma_q_vector <- sapply(0:p, FUN = sample_autocovariance, y_data, q)
    Gamma_q_matrix <- matrix(0, ncol=p, nrow=p)
    for(i in 1:p){
      for(j in 1:p)
        Gamma_q_matrix[i,j] <- gamma_q_vector[abs(i - j)+1]
      }
    cov_vec <- gamma_q_vector[2:(p + 1)] + correction[(q - 1 + 1):(q - p + 1)]
    #cat("q = ", q,", p = ", p, ", L1 = ", L1, ", L2 = ", L2, "\n")
    a_mat[, q - L1 + 1] <- solve(Gamma_q_matrix) %*% cov_vec 
  }    
  a_hat <- rowMeans(a_mat)
  a_hat <- as.vector(a_hat)
  return(a_hat)
}

calculating_sigma_eta <- function(y_data, coefficients, p){
  y_diff <- c(0, diff(y_data))
  T_size <- length(y_diff)
  
  res    <- rep(0, T_size)

  for(i in (p+1):T_size){
    res[i] <- y_diff[i] - sum(coefficients * y_diff[(i-1):(i-p)])
  }
  
  var_eta <- mean(res^2)/2
  return(sqrt(var_eta))
}

corrections <- function(coefficients, sigma_eta, len){
  p        <- length(coefficients)  
  c_vec    <- rep(0,len)
  c_vec[1] <- 1
  for(j in 2:len){
    lags <- (j-1):max(j-p,1)
    c_vec[j] <- sum(coefficients[1:length(lags)] * c_vec[lags])  
  }
  c_vec <- c_vec * sigma_eta * sigma_eta
  return(c_vec)
}

#Calculate autocovariance function for AR(1) model: \varepsilon_t = a_1 \varepsilon_{t-1} + \eta_t based on true coefficients of the model
autocovariance_function_AR1 <- function(k, a_1, sigma_eta){
  if (k%%1==0)
  {
    result = sigma_eta * sigma_eta * a_1^(abs(k)) / (1 - a_1 * a_1)
  } else {
    print('Check the input: k is not integer')
  }
  return(result)
}

#Function that compares the minimal intervals for both SiZer and our method.
#For each run of this function the result is a pdf file with three plots:
#1. Possible time series together with underlying time trend
#2. Union of minimal intervals for our method for N_rep number of simulations
#3. Union of minimal intervals for SiZer for N_rep number of simulations
plotting_many_minimal_intervals <- function(trend_height, trend_width, T_size, SiZer_matrix, N_rep, kernel_ind, sigmahat, gaussian_quantile, a_1, sigma_eta){
  matrix_our_results   <- data.frame('startpoint' = SiZer_matrix$u - SiZer_matrix$h, 'endpoint' = SiZer_matrix$u + SiZer_matrix$h)
  matrix_their_results <- data.frame('startpoint' = SiZer_matrix$u - SiZer_matrix$h, 'endpoint' = SiZer_matrix$u + SiZer_matrix$h)
  
  biweight_trend  <- numeric(T_size)
  for (i in 1:T_size) {biweight_trend[i] = trend_height * biweight_kernel(trend_width *(i/T_size - 0.5))}
  
  for (col in 1:N_rep){
    y_data_ar_1_with_trend <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta)) + biweight_trend
    
    g_t_set     <- psihat_statistic_ll(y_data_ar_1_with_trend, SiZer_matrix, kernel_ind, sigmahat)[[1]]
    
    results_our   <- c()
    results_their <- c()
    
    for (row in 1:nrow(g_t_set)){
      i                = g_t_set[row, 'u']
      h                = g_t_set[row, 'h']
      q_h              = g_t_set[row, 'q_h']
      sd_m_hat_prime   = g_t_set[row, 'sd']
      XtWX_inverse_XtW = g_t_set$XtWX_inv_XtW[[row]]
      
      m_hat_prime <- (XtWX_inverse_XtW %*% y_data_ar_1_with_trend)[2]
      
      if (m_hat_prime - q_h * sd_m_hat_prime > 0){
        results_their = c(results_their, 1)
      } else if (m_hat_prime + q_h * sd_m_hat_prime < 0) {
        results_their = c(results_their, -1)
      } else {
        results_their = c(results_their, 0)
      }
      
      if (g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
        results_our = c(results_our, 1)
      } else if (-g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
        results_our = c(results_our, -1)
      } else {
        results_our = c(results_our, 0)
      }
    }
    matrix_our_results <- cbind(matrix_our_results, results_our)
    matrix_their_results <- cbind(matrix_their_results, results_their)
  }
  
  grid_points <- seq(from = 1/T_size, to = 1, length.out = T_size) #grid points for estimating
  
  pdffilename = paste0("Paper/Plots/min_int_with_T_", T_size, "_a1_", a_1*100, ".pdf")
  pdf(pdffilename, width=8, height=10, paper="special")
  
  par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  par(mar = c(1.5, 0.5, 0, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  
  plot(grid_points, y_data_ar_1_with_trend, ylim = c(min(y_data_ar_1_with_trend) - 0.2, max(y_data_ar_1_with_trend)+0.2), type = "l")
  lines(grid_points, biweight_trend)

  par(mar = c(1.5, 0.5, 2, 0)) #Margins for each plot
  
  plot(NA, xlim=c(0,1), ylim = c(-1, N_rep +1), main = "Our test")
  for (col in 3:(N_rep+2)){
    a_t_set <- subset(matrix_our_results, matrix_our_results[,col] != 0, select = c(startpoint, endpoint, col))
    colnames(a_t_set) <- c('startpoint', 'endpoint', 'values')
    p_t_set <- choosing_minimal_intervals(a_t_set)
    if (!is.null(p_t_set)) {segments(p_t_set$startpoint, col, p_t_set$endpoint, col)}
  }
  
  plot(NA, xlim=c(0,1), ylim = c(-1, N_rep +1), main = "SiZer")
  for (col in 3:(N_rep+2)){
    a_t_set <- subset(matrix_their_results, matrix_their_results[,col] != 0, select = c(startpoint, endpoint, col))
    colnames(a_t_set) <- c('startpoint', 'endpoint', 'values')
    p_t_set <- choosing_minimal_intervals(a_t_set)
    if (!is.null(p_t_set)){segments(p_t_set$startpoint, col, p_t_set$endpoint, col)}
  }
  dev.off()
}

calculating_SiZer_matrix <- function(different_i, different_h, T_size, T_star, alpha, gamma){
  
  SiZer_matrix              <- expand.grid(u = different_i, h = different_h) #Creating a dataframe with all possible combination of i and h
  SiZer_matrix$values       <- numeric(nrow(SiZer_matrix)) # Setting the values of the statistic to be zero
  SiZer_matrix$sd           <- numeric(nrow(SiZer_matrix)) # Setting the values of standard deviation to be zero
  SiZer_matrix$ESS_star     <- numeric(nrow(SiZer_matrix)) # Setting the values of ESS* to be zero
  SiZer_matrix$q_h          <- numeric(nrow(SiZer_matrix)) # Setting the values of the gaussian quantile to be zero
  SiZer_matrix$small_ESS    <- numeric(nrow(SiZer_matrix)) # Later we will delete all the row such that ESS* is too small
  SiZer_matrix$lambda       <- lambda(SiZer_matrix[['h']]) # Calculating the lambda(h) in order to speed up the function psistar_statistic
  SiZer_matrix$XtWX_inv_XtW <- I(vector(mode="list", length=nrow(SiZer_matrix)))
  
  for (row in 1:nrow(SiZer_matrix)){
    
    i = SiZer_matrix[row, 'u']
    h = SiZer_matrix[row, 'h']
    
    ESS      <- sum(sapply((i - seq(1/T_size, 1, by = 1/T_size))/h, epanechnikov_kernel))/epanechnikov_kernel(0)
    ESS_star <- (T_star/T_size) * ESS 
    l        <- T_size / ESS_star
    q_h      <- qnorm((1 + (1 - alpha)^(1 / l))/2)
    
    if (ESS_star <= 5){
      SiZer_matrix[row, 'small_ESS'] <- 1
    } else {
      x_matrix                     <- matrix(c(rep(1, T_size), seq(1/T_size - i, 1-i, by = 1/T_size)), nrow=T_size, byrow=FALSE)
      w_vector                     <- sapply(seq(1/T_size - i, 1-i, by = 1/T_size)/h, FUN = epanechnikov_kernel)/h
      XtW                          <- t(apply(x_matrix,2,function(x){x*w_vector}))   #t(X) %*% diag(w)   computed faster.  :)
      XtWX_inverse                 <- tryCatch({solve(XtW %*% x_matrix)},  
                                               error = function(e) {print("Something is wrong, the matrix can not be inverted")})
      SiZer_matrix$XtWX_inv_XtW[[row]] <- XtWX_inverse %*% XtW
      
      XtSigma <- matrix(data = NA, nrow =2, ncol =T_size)
      for (j in 1:T_size){      #t(X) %*% Sigma computed faster in order not to save the whole Sigma matrix (T_size * T_size) in the memory
        XtSigma[1, j] = 0
        XtSigma[2, j] = 0
        for (k in 1:T_size){
          XtSigma[1, j] = XtSigma[1, j] + gamma[abs(k - j) + 1] * w_vector[k] * w_vector[j]
          XtSigma[2, j] = XtSigma[2, j] + gamma[abs(k - j) + 1] * w_vector[k] * w_vector[j] * (k/T_size - i)   
        }
      }
      sd = sqrt((XtWX_inverse %*% XtSigma %*% x_matrix %*% XtWX_inverse)[2,2])
      
      SiZer_matrix[row, 'sd']       = sd
      SiZer_matrix[row, 'q_h']      = q_h
      SiZer_matrix[row, 'ESS_star'] = ESS_star
    }
  }
  
  SiZer_matrix           <- SiZer_matrix[SiZer_matrix$small_ESS != 1,]
  SiZer_matrix$small_ESS <- NULL
  SiZer_matrix$ESS_star  <- NULL
  
  return(SiZer_matrix)
}
