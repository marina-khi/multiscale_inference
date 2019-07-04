simulating_data <- function(T, a1, sigma_eta, sim.design, slope.fac = 0){
  # compute AR error terms
  eps <- rep(0,T+50)
  eta <- rnorm(T+50, mean=0, sd=sigma_eta)
  for(i in 2:(T+50)) {eps[i] <- a1*eps[i-1] + eta[i]}
  eps <- eps[-(1:50)]

  #true long-run error variance (sqrt thereof)
  sigma  <- sqrt(sigma_eta^2/((1 - a1)^2))
  
  int.plus <- NULL
  int.minus <- NULL

  if(sim.design == "constant"){
    #constant trend function
    trend <- rep(0,T)
  } else if (sim.design == "linear"){
    #linear trend function
    slope <- slope.fac * sqrt(sigma_eta^2/(1-a1^2))
    trend <- slope * (1:T/T)
    int.plus <- c(0, 1)
  } else if (sim.design == "line"){
    trend <- slope * (1:T/T)
    int.plus <- c(0, 1)
  } else if (sim.design == "blocks"){
    #trend function of blocks example
    h_for_blocks <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
    t_for_blocks <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
    trend <- numeric(T)
    for(i in 1:T) {trend[i] <- (0.6 / 9.2) * (sum(h_for_blocks*(1 + sign(i/T - t_for_blocks))/2) + 2) + 0.2}
  } else if (sim.design == "brokenline") {
    #trend function which is a broken line
    u.lower <- 0.45
    u.upper <- 0.55 
  
    int.plus  <- c(u.lower,u.upper)
    int.minus <- "empty" 
  
    trend.fct <- function(u,u.lower,u.upper){ as.double(u >= u.lower & u <= u.upper) * 0.5 * (u-u.lower) / (u.upper-u.lower) + 0.5 * as.double(u > u.upper) }
    trend <- trend.fct(seq(1/T,1,by=1/T),u.lower,u.upper) 
  } else if (sim.design == "spike"){
    # trend function which is a spike
    u.lower <- 0.45
    u.upper <- 0.55 
  
    int.plus <- c(u.lower,0.5)
    int.minus <- c(0.5,u.upper)
  
    trend.fct <- function(u,u.lower,u.upper){ as.double(u >= u.lower & u <= 0.5) * (u-u.lower) / (0.5-u.lower) + as.double(u > 0.5 & u <= u.upper) * (u.upper-u) / (u.upper-0.5)}
    trend <- 0.75 * trend.fct(seq(1/T,1,by=1/T),u.lower,u.upper) 
  } else if (sim.design == "bump"){
    #trend function which is a bump
    u.lower <- 0.45
    u.upper <- 0.55 
  
    int.plus <- c(u.lower,0.5)
    int.minus <- c(0.5,u.upper)
  
    bump  <- function(u){
      arg <- (u-0.5)/(u.upper-0.5)
      return(0.5 * as.double(u >= u.lower & u <= u.upper) * (1-arg^2)^2)
    }
    trend <- bump((1:T)/T) * slope.fac
  } else if (sim.design == 'sine'){
    trend <- sin(6*pi*(1:T)/T)
  } else {cat("Unrecognized sim.design \n")}
  
  #simulated time series
  data <- trend + eps
  return(list(trend = trend, data = data, sigma = sigma, int.plus = int.plus, int.minus = int.minus))
}