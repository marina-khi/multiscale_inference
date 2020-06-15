#' Epanechnikov kernel function.
#' @param x A number.
#' @return 3/4(1-x^2) for |x|<=1 and 0 elsewhere.
#' @example 
#' epanechnikov_kernel(1)

epanechnikov_kernel <- function(x)
{
  if (abs(x)<=1)
  {
    result = 3/4 * (1 - x*x)
  } else {
    result = 0
  }
  return(result)
}

#' Function needed for local linear smoothing
#' @param u      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
s_t_1 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

#' Function needed for local linear smoothing
#' @param u      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
s_t_2 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

#' Function needed for local linear smoothing
#' @param u      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
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

#Local Linear estimator using the epanechnikov kernel. 
nadaraya_watson_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    result = sum((abs((grid_points - u) / bw) <= 1) * data_p)
    norm = sum((abs((grid_points - u) / bw) <= 1))
    #for (i in 1:T_size){
    #  result = result + epanechnikov_kernel((grid_p[i] - u) / bw) * data_p[i]
    #  norm = norm + epanechnikov_kernel((grid_p[i] - u) / bw)
    #}
    return(result/norm)
  }
}


# correction factor for error variance

correct <- function(Y, bw=0.025)
{  Y <- as.matrix(Y)
nn <- dim(Y)[2]
TT <- dim(Y)[1]
X <- (1:TT)/TT
const <- rep(0,nn)
for(i in 1:nn)
{  lambda.fct <- nw(X,Y[,i],bw,TT)
resid <- Y[,i] - lambda.fct
pos <- (lambda.fct > 0)
resid <- resid[pos]/sqrt(lambda.fct[pos])
const[i] <- var(resid)
}   
const <- mean(const)
return(const)
}   

produce_plots <- function (results, l, data_i, data_j, smoothed_i, smoothed_j,
                           gov_resp_i, gov_resp_j, lagged_gov_resp_i, lagged_gov_resp_j,
                           country_i, country_j, text_){
  Tlen <- length(data_i)
  gset <- results$gset_with_vals[[l]]

  layout(matrix(c(1, 2, 3, 4),ncol=1), widths=c(3, 3, 3, 3),
         heights=c(1, 0.8, 1, 1), TRUE)
  #Setting the layout of the graphs

  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 1.5, 0.2)) #Outer margins

  plot(data_i, ylim=c(min(data_i, data_j), max(data_i, data_j)), type="l",
      col="blue", ylab="", xlab="", mgp=c(1, 0.5, 0))
  lines(data_j, col="red")
  title(main = expression((a) ~ observed ~ cases ~ per ~ day), line = 1)
  legend(x = 0, y = max(data_i, data_j) - 1, legend=c(country_i, country_j),
         col = c("blue", "red"), lty = 1, cex = 0.95, ncol = 1)

  par(mar = c(0.5, 0.5, 3.5, 0)) #Margins for each plot

  plot(smoothed_i, ylim=c(min(data_i, data_j), max(data_i, data_j)), type="l",
     col="blue", ylab="", xlab = "", mgp=c(1,0.5,0))
  lines(smoothed_j, col="red")
  title(main = expression((b) ~ smoothed ~ curve ~ from ~ (a)), line = 1)

  plot(gov_resp_i, ylim=c(0, 100), type="l",
       col="blue", ylab="", xlab = "", mgp=c(1,0.5,0))
  lines(gov_resp_j, col="red")
  lines(lagged_gov_resp_i, col="blue", lty = "dashed", lwd = 3)
  lines(lagged_gov_resp_j, col="red", lty = "dashed", lwd = 3)
  
  title(main = expression((c) ~ index ~ of ~ government ~ response), line = 1)
  
  par(mar = c(0.5, 0.5, 3, 0)) #Margins for each plot

  a_t_set <- subset(gset, test == 1, select = c(u, h))
  if (nrow(a_t_set) > 0){
    p_t_set <- data.frame('startpoint' = a_t_set$u - a_t_set$h, 'endpoint' = a_t_set$u + a_t_set$h, 'values' = 0)
    p_t_set$values <- (1:nrow(p_t_set))/nrow(p_t_set)
    
    #plot(NA, xlim=c(0,Tlen),  ylim = c(0, 1 + 1/nrow(p_t_set)), xlab="days", mgp=c(2,0.5,0))
    #title(main = expression((c) ~ all ~ intervals ~ where ~ the ~ test ~ rejects), line = 1)
    #segments(p_t_set[['startpoint']], p_t_set$plottingindex, p_t_set$endpoint, p_t_set$plottingindex, lwd = 2)
  
    #Produce minimal intervals
    p_t_set2               <- compute_minimal_intervals(p_t_set)

    plot(NA, xlim=c(0,Tlen),  ylim = c(0, 1 + 1/nrow(p_t_set)), xlab="days", mgp=c(2,0.5,0))
    title(main = expression((d) ~ minimal ~ intervals ~ produced ~ by ~ our ~ test), line = 1)
    segments(p_t_set2$startpoint * Tlen, p_t_set2$values, p_t_set2$endpoint *Tlen, p_t_set2$values, lwd = 2)
    segments(p_t_set$startpoint * Tlen, p_t_set$values, p_t_set$endpoint *Tlen, p_t_set$values, col = "gray")
  } else {
    #If there are no intervals where the test rejects, we produce empty plots
    #plot(NA, xlim=c(0,Tlen),  ylim = c(0, 1), xlab="days", mgp=c(2,0.5,0))
    #title(main = expression((c) ~ all ~ intervals ~ where ~ the ~ test ~ rejects), line = 1)
    
    plot(NA, xlim=c(0,Tlen),  ylim = c(0, 1), xlab="days", mgp=c(2,0.5,0))
    title(main = expression((d) ~ minimal ~ intervals ~ produced ~ by ~ our ~ test), line = 1)
  }
  mtext(text_, side = 3, line = 0, outer = TRUE)
  
  #SiZer 
  #par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  #plot_SiZer_map(sort(unique(gset[,1])), sort(unique(gset[,2])), test.results = result$test_matrices[[l]], plot.title = expression((c) ~ SiZer ~ map ~ 'for' ~ italic(T)[MS]))
  
  #axis_at = seq(4/Tlen, 724/Tlen, by = 30/Tlen)
  #axis_labels = seq(as.Date("2018/4/4"), as.Date("2020/4/4"), by = 30)
  #axis(1, at=axis_at, labels = axis_labels, mgp=c(1,0.5,0))
  return(a_t_set)
}