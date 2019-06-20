critical_values <- function(T, grid, weights, kappa=0.1, SimRuns=1000, probs=seq(0.5,0.995,by=0.005)){


  # Function that computes the critical value of the multiscale test.
  #
  # Inputs:
  # T           sample size 
  # grid        grid of location-bandwidth points as produced by the function 'grid_construction',
  #             list with the elements 'gset', 'bws', 'gtype' 
  # weights     kernel weights matrix produced by the function 'kernel_weights'
  # kappa       specifies which order statistic to use for the version of the multiscale test
  #             defined in Section ??
  # SimRuns     number of simulations used to compute the critical value 
  #             Default is SimRuns=1000.
  # probs       vector of significance levels for which the critical value is computed
  #             Default is probs=seq(0.5,0.995,by=0.005).
  #           
  # Outputs: list with the elements  
  # quant_ms     vector of critical values for all significance levels in the vector 'probs' 
  #              for our version of the multiscale test defined in Section 3
  # quant_uncor  vector of critical values for all significance levels in the vector 'probs' 
  #              for a version of the multiscale test without the additive correction term 
  #              lambda
  # quant_order  vector of critical values for all significance levels in the vector 'probs' 
  #              for the version of the multiscale test defined in Section ??
  # quant_rows   matrix of critical values for the rowwise version of the test,
  #              each row of the matrix contains the critical values for the
  #              significance levels in 'probs' and a specific bandwidth h        


  gset    <- grid$gset
  gtype   <- grid$gtype
  h.vec   <- grid$bws   
  N       <- dim(gset)[1]
  N.h     <- length(h.vec)
  

  filename.ms    <- paste0("distributions/distr_T_", T, "_multiscale.txt")
  filename.uncor <- paste0("distributions/distr_T_", T, "_uncorrected.txt")
  filename.order <- paste0("distributions/distr_T_", T, "_orderstat.txt")
  filename.rows  <- paste0("distributions/distr_T_", T, "_rowwise.txt")


  if(gtype == "default" & file.exists(filename.ms)){ 
    quant.ms    <- as.matrix(read.table(filename.ms, header=FALSE))
    quant.ms    <- as.vector(quant.ms)
    quant.uncor <- as.matrix(read.table(filename.uncor, header=FALSE))
    quant.uncor <- as.vector(quant.uncor)
    quant.order <- as.matrix(read.table(filename.order, header=FALSE))
    quant.order <- as.vector(quant.order)
    quant.rows  <- as.matrix(read.table(filename.rows, header=FALSE))
  }


  if(gtype != "default" | !file.exists(filename.ms)){
 
    Psi.ms    <- rep(NA,SimRuns)
    Psi.uncor <- rep(NA,SimRuns)
    Psi.order <- rep(NA,SimRuns)
    Psi.rows  <- matrix(NA,ncol=SimRuns,nrow=N.h)

    cat("","\n")
    cat("Computing critical value of the test","\n")
    progbar <- txtProgressBar(min = 1, max = SimRuns, style = 3, char = ".")
  
    wghts <- kernel_weights(T, grid)
    
    for(pos in 1:SimRuns)
    {  zeta <- rnorm(T, mean=0, sd=1)
       res  <- multiscale_statistics(data=zeta, weights=wghts, sigmahat=1, grid=grid) 
    
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
       quant.rows[pos,] <- quantile(Psi.rows[pos,], probs=probs)
  }


  if(gtype == "default" & !file.exists(filename.ms)) 
  { write.table(quant.ms, file=filename.ms, col.names=FALSE, row.names=FALSE)
    write.table(quant.uncor, file=filename.uncor, col.names=FALSE, row.names=FALSE)
    write.table(quant.order, file=filename.order, col.names=FALSE, row.names=FALSE)
    write.table(quant.rows, file=filename.rows, col.names=FALSE, row.names=FALSE)
  }


  return(list(probs=probs, quant_ms=quant.ms, quant_uncor=quant.uncor, quant_order=quant.order, quant_rows=quant.rows))

}
