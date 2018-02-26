psihat_statistic <- function(y_data, g_t_set, kernel_function = epanechnikov, sigmahat){
  #Wrapper of a C function from psihat_statistic.C
  #Arguments:
  # y_data          vector of length T 
  # g_t_set         dataframe with columns "u", "h", "values" (equal to 0), "lambda"
  # kernel_function type of kernel function used. Currently only epanechnikov and biweight kernels are supported
  # sigmahat        appropriate estimate of sqrt(long-run variance)
  
  T <- as.integer(length(y_data))

  #Passing an indicator to C functions that determines a necessary kernel function
  if (kernel_function = epanechnikov) {kernel_ind = 1}
  else {kernel_ind = 2}

  
  
    
  d <- as.integer(dim(x)[2])

if (!n == length(y))
  stop("Predictor and Response of different length")

if (!length(bw) == d)
  stop("Wrong dimension of bandwidth vector")

storage.mode(x)       <- "double"
storage.mode(y)       <- "double"
storage.mode(bw)      <- "double"
storage.mode(N)       <- "integer"
storage.mode(iterate) <- "integer"

if(grids==TRUE)
  Grid <- Grid1
else
  Grid <- rep(seq(0,1,1/(N-1)),d)

L <- ceiling(0.5*d*(d-1)*N)
L <- as.integer(L)


xvec     <- as.vector(x)
gridvec  <- as.vector(Grid)
m0       <- vector(mode="double",length=1)
kde      <- vector(mode="double",length=d*N)
regbf    <- vector(mode="double",length=d*N)
regnw    <- regbf       
conv     <- vector(mode="double",length=d)
fh_int   <- vector(mode="double",length=d)
fh_int2a <- vector(mode="double",length=L)
fh_int2b <- vector(mode="double",length=L)  
reg_int  <- vector(mode="double",length=d) 


result <- .C("sbfnw", xvec, y, n, d, gridvec, N, bw, kde, m0, regbf, regnw,
             iterate, conv, fh_int, fh_int2a, fh_int2b, reg_int)

list(grid=matrix(result[[5]], ncol=d),bw=result[[7]], kde=matrix(result[[8]], ncol=d),
     m0=result[[9]], sbf=matrix(result[[10]], ncol=d), nwr=matrix(result[[11]], ncol=d),
     iteration=result[[12]], conv=result[[13]], fh_int=matrix(result[[14]],ncol=d),
     fh_int2a=matrix(result[[15]],ncol=ceiling(0.5*d*(d-1))),
     fh_int2b=matrix(result[[16]],ncol=ceiling(0.5*d*(d-1))), 
     nwr_int=matrix(result[[17]],ncol=d)) 
}




# Function that calculates the multiscale statistic \hat{\Psi}_T.
# It takes the following entities as arguments.
#   y-data: the data
#   g_t_set: range of different locations u and bandwidths h
#   k_function: type of kernel function
#   sigmahat: the estimator of the square root of the long-run error variance \sigma^2
# It produces the value of the test statistic as an output
psihat_statistic <- function(y_data, g_t_set, kernel_function = epanechnikov_kernel, sigmahat) {
  g_t_set$values <- abs(mapply(psi_average, u = g_t_set$u, h = g_t_set$h,
                               MoreArgs = list(data = y_data, k_function = kernel_function))/sigmahat) - .subset2(g_t_set, 'lambda')
  result = max(g_t_set$values)
  return(list(g_t_set, result))
}

#old version
#psihat_statistic <- function(y_data, g_t_set, k_function = epanechnikov_kernel, sigmahat) {
#  g_t_set_card = nrow(g_t_set)
#  for (i in 1:g_t_set_card) {
#    g_t_set[['values']][i] <- abs(psi_average(y_data, .subset2(g_t_set,'u')[i], .subset2(g_t_set,'h')[i], k_function)/sigmahat) - .subset2(g_t_set, 'lambda')[i]
#  }
#  result = max(g_t_set$values)
#  return(result)
#}


# Function that calculates the auxiliary statistic \Psi^star_T.
# The only difference with the previus function is in the return values.
psistar_statistic <- function(y_data, g_t_set, kernel_function = epanechnikov_kernel, sigmahat) {
  g_t_set$values <- abs(mapply(psi_average, u = g_t_set$u, h = g_t_set$h,
                               MoreArgs = list(data = y_data, k_function = kernel_function))/sigmahat) - .subset2(g_t_set, 'lambda')
  result = max(g_t_set$values)
  return(result)
}


#This functions chooses minimal intervals as described in Dumbgen()
choosing_minimal_intervals <- function(dataset){
  set_cardinality <- nrow(dataset) 
  if (set_cardinality > 1) {
    dataset <- dataset[order(dataset$startpoint, -dataset$endpoint),]
    rownames(dataset) <- 1:nrow(dataset)
    dataset[['contains']] <- numeric(set_cardinality)
    for (i in 1:(set_cardinality-1)){
      for (j in (i+1):set_cardinality){
        if ((dataset$startpoint[i] <= dataset$startpoint[j])&(dataset$endpoint[i] >= dataset$endpoint[j])) {
          #          cat("[",dataset$startpoint[i], ", ", dataset$endpoint[i],"], [", dataset$startpoint[j], ", ", dataset$endpoint[j], "]", sep="")
          dataset[['contains']][i] <- 1
          break
        }
      }
    }
    p_t_set <- subset(dataset, contains == 0, select = c(startpoint, endpoint, values))
  }
}


#Creating g_t_set over which we are taking the maximum (from Section 2.1)
creating_g_set <- function(T){
  u <- seq(4/T, 1, length.out = T/4)
  h <- seq(3/T, 1/4+3/T, length.out = T/20)
  
  g_t_set_temp <- expand.grid(u = u, h = h) #Creating a dataframe with all possible combination of u and h
  g_t_set_temp$values <-numeric(nrow(g_t_set_temp)) # Setting the values of the statistic to be zero
  
  g_t_set <- subset(g_t_set_temp, u - h >= 0 & u + h <= 1, select = c(u, h, values)) #Subsetting u and h such that [u-h, u+h] lies in [0,1]
  g_t_set$lambda <- lambda(g_t_set[['h']]) #Calculating the lambda(h) in order to speed up the function psistar_statistic
  return(g_t_set)
}