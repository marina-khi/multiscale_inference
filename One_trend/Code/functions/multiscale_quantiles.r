multiscale_quantiles <- function(T, grid, weights, kappa=0.1, SimRuns=1000, probs=seq(0.5,0.995,by=0.005)){


  # Function that computes the quantiles of the multiscale test statistics under the null.
  #
  # Inputs:
  # T            sample size 
  # grid         grid of location-bandwidth points as produced by the function 'grid_construction',
  #              list with the elements 'gset', 'bws', 'gtype' 
  # weights      kernel weights matrix produced by the function 'kernel_weights'
  # kappa        specifies which order statistic to use for the version of the multiscale test
  #              defined in Section ??
  # SimRuns      number of simulations used to compute the critical value 
  #              Default is SimRuns=1000.
  # probs        vector of probability levels (1-alpha) for which the quantiles are computed 
  #              Default is probs=seq(0.5,0.995,by=0.005).
  #           
  # Outputs: list with the elements  
  # quant_ms     vector of quantiles for all probability levels in the vector 'probs' 
  #              for our version of the multiscale statistic defined in Section 3
  # quant_uncor  vector of quantiles for all probability levels in the vector 'probs' 
  #              for a version of the multiscale statistic without the additive 
  #              correction term 'lambda'
  # quant_order  vector of quantiles for all probability levels in the vector 'probs' 
  #              for the version of the multiscale statistic defined in Section ??
  # quant_rows   matrix of quantiles for the row-wise version of the test,
  #              each row of the matrix contains the quantiles for the probability
  #              levels in 'probs' and a specific bandwidth h        


  gset    <- grid$gset
  h.vec   <- grid$bws   
  N       <- dim(gset)[1]
  N.h     <- length(h.vec)
 
  Psi.ms    <- rep(NA,SimRuns)
  Psi.uncor <- rep(NA,SimRuns)
  Psi.order <- rep(NA,SimRuns)
  Psi.rows  <- matrix(NA,ncol=SimRuns,nrow=N.h)

  cat("","\n")
  cat("Computing critical value of the multiscale test","\n")
  progbar <- txtProgressBar(min = 1, max = SimRuns, style = 3, char = ".")
 
  for(pos in 1:SimRuns)
  {  zeta <- rnorm(T, mean=0, sd=1)
     res  <- multiscale_statistics(data=zeta, weights=weights, sigmahat=1, grid=grid) 
    
     Psi.ms[pos]    <- res$stat_ms
     Psi.uncor[pos] <- res$stat_uncor
     Psi.rows[,pos] <- res$stat_rows
    
     temp.order     <- res$stat_order
     temp.order     <- sort(temp.order)
     Psi.order[pos] <- temp.order[max(1,floor((1-kappa)*N.h))]       

     setTxtProgressBar(progbar, pos)
  } 
  close(progbar)

  quant.ms    <- as.vector(quantile(Psi.ms, probs=probs))
  quant.uncor <- as.vector(quantile(Psi.uncor, probs=probs))
  quant.order <- as.vector(quantile(Psi.order, probs=probs))
  quant.rows  <- matrix(NA,nrow=N.h,ncol=length(probs))
  for(pos in 1:N.h)
     quant.rows[pos,] <- as.vector(quantile(Psi.rows[pos,], probs=probs))

  quant.ms    <- rbind(probs, quant.ms)
  quant.uncor <- rbind(probs, quant.uncor)
  quant.order <- rbind(probs, quant.order)
  quant.rows  <- rbind(probs, quant.rows)

  colnames(quant.ms)    <- NULL
  rownames(quant.ms)    <- NULL
  colnames(quant.uncor) <- NULL
  rownames(quant.uncor) <- NULL
  colnames(quant.order) <- NULL
  rownames(quant.order) <- NULL
  colnames(quant.rows)  <- NULL
  rownames(quant.rows)  <- NULL 

  return(list(quant_ms=quant.ms, quant_uncor=quant.uncor, quant_order=quant.order, quant_rows=quant.rows))

}
