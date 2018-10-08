histograms_for_variance_estimators <- function(a_1_star, sigma_eta_star, T_size, p, slope, N_rep, pdfname_a_hat, pdfname_lrv, L1, L2, K1, K2, M1, M2){
  lrv_star <- sigma_eta_star^2/(1-a_1_star)^2 
  
  line_trend  <- numeric(T_size)
  line_trend = (1:T_size - 0.5*T_size) * slope/T_size
    
  a_hat_hall_vect    <- c() 
  a_hat_method2_vect <- c()
  a_hat_oracle_vect  <- c()
  
  sigma_eta_hat_hall_vect    <- c() 
  sigma_eta_hat_method2_vect <- c()
  sigma_eta_hat_oracle_vect  <- c()
  
  lrv_hat_hall_vect    <- c() 
  lrv_hat_method2_vect <- c()
  lrv_hat_oracle_vect  <- c()
  
  for (i in 1:N_rep){

    eps               <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))
    y_data_with_trend <- eps + line_trend
    
    results_our_method         <- estimating_variance_new(y_data_with_trend, L1, L2, p, K1, K2)
    a_hat_method2_vect         <- c(a_hat_method2_vect, results_our_method[[2]])
    sigma_eta_hat_method2_vect <- c(sigma_eta_hat_method2_vect, results_our_method[[3]])
    lrv_hat_method2_vect       <- c(lrv_hat_method2_vect, results_our_method[[1]]^2)
    
    #cat("Our method:", results_our_method[[1]]^2)
        
    result_HvK         <- estimating_sigma_for_AR1(y_data_with_trend, M1, M2)
    a_hat_hall_vect    <- c(a_hat_hall_vect, result_HvK[[2]])
    sigma_eta_hat_hall <- calculating_sigma_eta(y_data_with_trend, result_HvK[[2]], p)
    lrv_hat_hall       <- sigma_eta_hat_hall^2/(1 - sum(result_HvK[[2]]))^2
    
    sigma_eta_hat_hall_vect <- c(sigma_eta_hat_hall_vect, sigma_eta_hat_hall)
    lrv_hat_hall_vect       <- c(lrv_hat_hall_vect, lrv_hat_hall)

    #cat("a_hat:", result_HvK[[2]], ", sigma_eta:", sigma_eta_hat_hall, ", LRV:", lrv_hat_hall, "\n")
    
    res_oracle               <- lm(eps[2:T_size] ~ eps[1:(T_size-1)] - 1)
    a_hat_oracle_vect         <- c(a_hat_oracle_vect, as.vector(res_oracle$coef))
    sigma_eta_hat_oracle_vect <- c(sigma_eta_hat_oracle_vect, sqrt(mean(res_oracle$residuals^2)))
    lrv_hat_oracle_vect       <- c(lrv_hat_oracle_vect, mean(res_oracle$residuals^2)/(1-sum(as.vector(res_oracle$coef)))^2)
    #cat("ORacle:",  mean(res_oracle$residuals^2)/(1-sum(as.vector(res_oracle$coef)))^2)
  }

  smallest_value <- min(a_hat_hall_vect, a_hat_method2_vect, a_hat_oracle_vect)
  biggest_value <- max(a_hat_hall_vect, a_hat_method2_vect, a_hat_oracle_vect)

  step_len       <- (biggest_value-smallest_value)/50
  if (smallest_value < -1.4){
    steps <- ceiling(50*(biggest_value-smallest_value)/(biggest_value+1.4))
    step_len <- (biggest_value-smallest_value)/steps  
  }
  
  breaks_grid <- seq(smallest_value, biggest_value, by = step_len)
  breaks_grid[length(breaks_grid)] <- biggest_value 

  hist1 <- hist(a_hat_oracle_vect, breaks = breaks_grid, plot = FALSE)
  hist3 <- hist(a_hat_method2_vect, breaks = breaks_grid, plot = FALSE)
  hist4 <- hist(a_hat_hall_vect, breaks = breaks_grid, plot = FALSE)

  highestCount <- max(hist1$counts, hist3$counts, hist4$counts)

  pdf(pdfname_a_hat, width=8, height=2.9, paper="special")
  par(mfrow=c(1,3))
  par(mar = c(3, 2, 0.5, 1)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.5, 0.2)) #Outer margins 
  hist(a_hat_method2_vect, main=NULL, breaks = breaks_grid, freq=TRUE, xlim=c(max(c(-1.4,smallest_value)),biggest_value), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0), cex.lab=1.1)
  mtext(side=1,text=expression(widehat(a)),line=2.5)
  segments(x0=a_1_star,y0=0,x1=a_1_star,y1=highestCount,col="red",lwd=1.5)
  hist(a_hat_hall_vect, main=NULL, breaks = breaks_grid, freq=TRUE, xlim=c(max(c(-1.4,smallest_value)),biggest_value), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0), cex.lab=1.1)
  mtext(side=1,text=expression(widehat(a)[HvK]),line=2.5)
  segments(x0=a_1_star,y0=0,x1=a_1_star,y1=highestCount,col="red",lwd=1.5)
  hist(a_hat_oracle_vect, main=NULL, breaks = breaks_grid, freq=TRUE, xlim=c(max(c(-1.4,smallest_value)),biggest_value), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0), cex.lab=1.1)
  mtext(side=1,text=expression(widehat(a)[oracle]),line=2.5)
  segments(x0=a_1_star,y0=0,x1=a_1_star,y1=highestCount,col="red",lwd=1.5)
  dev.off()
  

  
  smallest_value <- min(c(lrv_hat_method2_vect,lrv_hat_hall_vect,lrv_hat_oracle_vect))
  biggest_value  <- max(c(lrv_hat_method2_vect,lrv_hat_hall_vect,lrv_hat_oracle_vect))
  step_len       <- (biggest_value-smallest_value)/50

  breaks_grid <- seq(smallest_value, biggest_value, by = step_len)
  breaks_grid[length(breaks_grid)] <- biggest_value 
  
  hist1 <- hist(lrv_hat_method2_vect, breaks = breaks_grid, plot = FALSE)
  hist3 <- hist(lrv_hat_hall_vect, breaks = breaks_grid, plot = FALSE)
  hist4 <- hist(lrv_hat_oracle_vect, breaks = breaks_grid, plot = FALSE)
  
  highestCount <- max(hist1$counts, hist3$counts, hist4$counts)
  
  pdf(pdfname_lrv, width=8, height=2.9, paper="special")
  par(mfrow=c(1,3))
  par(mar = c(3, 2, 0.5, 1)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.5, 0.2)) #Outer margins 
  hist(lrv_hat_method2_vect, main = NULL, breaks = breaks_grid, freq=TRUE, xlim=c(smallest_value,biggest_value), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0), cex.lab = 1.1)
  mtext(side=1,text=expression(widehat(sigma)[]^2),line=2.75)
  segments(x0=lrv_star,y0=0,x1=lrv_star,y1=highestCount,col="red",lwd=1.5)
  hist(lrv_hat_hall_vect, main = NULL, breaks = breaks_grid, freq=TRUE, xlim=c(smallest_value,biggest_value), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0), cex.lab = 1.1)
  mtext(side=1,text=expression(widehat(sigma)[HvK]^2),line=2.75)
  segments(x0=lrv_star,y0=0,x1=lrv_star,y1=highestCount,col="red",lwd=1.5)
  hist(lrv_hat_oracle_vect, main = NULL, breaks = breaks_grid, freq=TRUE, xlim=c(smallest_value,biggest_value), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0), cex.lab = 1.1)
  mtext(side=1,text=expression(widehat(sigma)[oracle]^2),line=2.75)   
  segments(x0=lrv_star,y0=0,x1=lrv_star,y1=highestCount,col="red",lwd=1.5)
  dev.off()
  
  # compute MSE values
  a_mse2       <- mean((a_hat_method2_vect-a_1_star)^2)
  a_mse_HvK    <- mean((a_hat_hall_vect-a_1_star)^2)
  a_mse_oracle <- mean((a_hat_oracle_vect-a_1_star)^2)
  
  lrv_mse2       <- mean((lrv_hat_method2_vect-lrv_star)^2)
  lrv_mse_HvK    <- mean((lrv_hat_hall_vect-lrv_star)^2)
  lrv_mse_oracle <- mean((lrv_hat_oracle_vect-lrv_star)^2)
  return(list(a_mse2, a_mse_HvK, a_mse_oracle, lrv_mse2, lrv_mse_HvK, lrv_mse_oracle))
}