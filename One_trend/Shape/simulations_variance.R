histograms_for_variance_estimators <- function(a_1, sigma_eta, T_size, p, slope, N_rep, K1, K2, L1, L2){
  line_trend  <- numeric(T_size)
  line_trend = (1:T_size - 0.5*T_size) * slope/T_size
    
  a_hat_hall_vect    <- c() 
  a_hat_method1_vect <- c()
  a_hat_method2_vect <- c()
  a_hat_oracle_vect  <- c()
  
  for (i in 1:N_rep){
    set.seed(i)
    
    eps               <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))
    y_data_with_trend <- eps + line_trend
    
    a_hat_method1         <- AR_coefficients(y_data_with_trend, L1, L2, rep(0,L2), p)
    sigma_eta_hat_method1 <- calculating_sigma_eta(y_data_with_trend, a_hat_method1, p)
    sigma_hat_method1     <- sqrt(sigma_eta_hat_method1^2 / (1 - sum(a_hat_method1))^2) 
    
    corrections_value     <- corrections(a_hat_method1, sigma_eta_hat_method1, K2+1)
    a_hat_method2         <- AR_coefficients(y_data_with_trend, K1, K2, corrections_value, p)
    sigma_eta_hat_method2 <- calculating_sigma_eta(y_data_with_trend, a_hat_method2, p)
    sigma_hat_method2     <- sqrt(sigma_eta_hat_method2^2 / (1 - sum(a_hat_method2))^2) 
    
    a_hat_method1_vect <- c(a_hat_method1_vect, a_hat_method1)
    a_hat_method2_vect <- c(a_hat_method2_vect, a_hat_method2)
    
    result          <- estimating_sigma_for_AR1(y_data_with_trend, L1, L2)
    a_hat_hall_vect <- c(a_hat_hall_vect, result[[2]])

    gamma_1 <- sum(eps[1:(T_size - 1)] * eps[2:T_size])/T_size
    gamma_0 <- sum(eps * eps)/T_size
    a_hat_oracle_vect <- c(a_hat_oracle_vect, gamma_1 / gamma_0)
  }

  smallest_value <- min(a_hat_hall_vect, a_hat_method1_vect, a_hat_method2_vect, a_hat_oracle_vect)
  biggest_value <- max(a_hat_hall_vect, a_hat_method1_vect, a_hat_method2_vect, a_hat_oracle_vect)

  hist1 <- hist(a_hat_oracle_vect, breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), plot = FALSE)
  hist2 <- hist(a_hat_method1_vect, breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), plot = FALSE)
  hist3 <- hist(a_hat_method2_vect, breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), plot = FALSE)
  hist4 <- hist(a_hat_hall_vect, breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), plot = FALSE)

  highestCount <- max(hist1$counts, hist2$counts, hist3$counts, hist4$counts)
  
  pdfname = paste0("Paper/Plots/variance_histograms_a1_", a_1, "_slope_", slope*10, ".pdf")
  pdf(pdfname, width=10, height=5, paper="special")
  par(mfrow=c(1,4))
  par(mar = c(1.0, 1.0, 1.0, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.5, 0.2)) #Outer margins
  
  hist(a_hat_oracle_vect, main="Oracle", breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0))
  hist(a_hat_method1_vect, main = paste0("Method 1"), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0))
  hist(a_hat_method2_vect, main = paste0("Method 2"), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0))
  hist(a_hat_hall_vect, main = paste0("Hall and van Keilegom"), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0))
  dev.off()
}