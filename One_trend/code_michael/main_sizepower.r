rm(list=ls())

source("functions/grid_construction.r")
dyn.load("functions/kernel_weights.dll")
source("functions/kernel_weights.r")
source("functions/multiscale_statistics.r")
source("functions/critical_value.r")
source("functions/multiscale_test.r")
source("functions/inputs_for_plots.r")


# Parameters

T          <- 1000            # sample size
a1         <- -0.25          # AR parameter 
sigma_eta  <- 1              # standard deviation of the innovation term in the AR model
sim.design <- "spike"        # trend specification: "constant", "blocks", ...

alpha      <- 0.05                  # significance level
SimRuns    <- 1000                  # number of simulation runs to produce critical values
kappa      <- 0.1                   # parameter to determine order statistic for the version
                                    # of the multiscale statistic from Section ?? 
grid       <- grid_construction(T)  # grid of location-bandwidth points 
                                    # the grid should be constructed via the function 
                                    # 'grid construction'. See the corresponding R-file 
                                    # the function definition for details 

Nsim       <- 1000   # number of simulation runs for size/power calculations


# Size/power calculations

sizepower.ms    <- 0
sizepower.uncor <- 0
sizepower.order <- 0
sizepower.rows  <- 0

power.plus.ms    <- 0
power.plus.uncor <- 0
power.plus.order <- 0
power.plus.rows  <- 0

power.minus.ms    <- 0
power.minus.uncor <- 0
power.minus.order <- 0
power.minus.rows  <- 0

plus.left   <- 0
plus.right  <- 1
minus.left  <- 0
minus.right <- 1

source("simulations/sim.r")  # implicitly specifies plus.left, plus.right, ...
                             # for the simulation design under consideration 
gset <- grid$gset
N    <- dim(gset)[1]
pos.power.plus  <- rep(1,N)
pos.power.minus <- rep(1,N)  

for(i in 1:N)
{  if(gset[i,1] - gset[i,2] > plus.right | gset[i,1] + gset[i,2] < plus.left)
     pos.power.plus[i] <- 0
   if(gset[i,1] - gset[i,2] > minus.right | gset[i,1] + gset[i,2] < minus.left)
     pos.power.minus[i] <- 0
}   

wghts <- kernel_weights(T=T,grid=grid)

for(loops in 1:Nsim)

{  # Simulate data

   source("simulations/sim.r") 

   # Compute test results 

   res <- multiscale_test(alpha=alpha, data=data, weights=wghts, sigmahat=sigma, grid=grid, kappa=kappa, SimRuns=SimRuns)
  
   res.ms    <- res$test_ms
   res.uncor <- res$test_uncor
   res.order <- res$test_order 
   res.rows  <- res$test_rows 

   res.plus.ms    <- res$test_ms * pos.power.plus
   res.plus.uncor <- res$test_uncor * pos.power.plus
   res.plus.order <- res$test_order * pos.power.plus
   res.plus.rows  <- res$test_rows * pos.power.plus

   res.minus.ms    <- res$test_ms * pos.power.minus
   res.minus.uncor <- res$test_uncor * pos.power.minus
   res.minus.order <- res$test_order * pos.power.minus
   res.minus.rows  <- res$test_rows * pos.power.minus

   if(sum(res.ms == 1) + sum(res.ms == -1) > 0)
     sizepower.ms <- sizepower.ms + 1 
   if(sum(res.uncor == 1) + sum(res.uncor == -1) > 0)
     sizepower.uncor <- sizepower.uncor + 1 
   if(sum(res.order == 1) + sum(res.order == -1) > 0)
     sizepower.order <- sizepower.order + 1 
   if(sum(res.rows == 1) + sum(res.rows == -1) > 0)
     sizepower.rows <- sizepower.rows + 1 

   if(sum(res.plus.ms == 1) > 0)
     power.plus.ms <- power.plus.ms + 1 
   if(sum(res.plus.uncor == 1) > 0)
     power.plus.uncor <- power.plus.uncor + 1 
   if(sum(res.plus.order == 1) > 0)
     power.plus.order <- power.plus.order + 1 
   if(sum(res.plus.rows == 1) > 0)
     power.plus.rows <- power.plus.rows + 1 

   if(sum(res.minus.ms == -1) > 0)
     power.minus.ms <- power.minus.ms + 1 
   if(sum(res.minus.uncor == -1) > 0)
     power.minus.uncor <- power.minus.uncor + 1 
   if(sum(res.minus.order == -1) > 0)
     power.minus.order <- power.minus.order + 1 
   if(sum(res.minus.rows == -1) > 0)
     power.minus.rows <- power.minus.rows + 1 

   print(loops)
}

sizepower.ms    <- sizepower.ms/Nsim
sizepower.uncor <- sizepower.uncor/Nsim
sizepower.order <- sizepower.order/Nsim
sizepower.rows  <- sizepower.rows/Nsim

power.plus.ms    <- power.plus.ms/Nsim
power.plus.uncor <- power.plus.uncor/Nsim
power.plus.order <- power.plus.order/Nsim
power.plus.rows  <- power.plus.rows/Nsim

power.minus.ms    <- power.minus.ms/Nsim
power.minus.uncor <- power.minus.uncor/Nsim
power.minus.order <- power.minus.order/Nsim
power.minus.rows  <- power.minus.rows/Nsim

cat("","\n")
cat(paste0("Size/power Multiscale = ", sizepower.ms),"\n")
cat(paste0("Size/power Uncorrected = ", sizepower.uncor),"\n")
cat(paste0("Size/power Orderstat = ", sizepower.order),"\n")
cat(paste0("Size/power Rowwise = ", sizepower.rows),"\n")

cat("","\n")
cat(paste0("Power plus Multiscale = ", power.plus.ms),"\n")
cat(paste0("Power plus Uncorrected = ", power.plus.uncor),"\n")
cat(paste0("Power plus Orderstat = ", power.plus.order),"\n")
cat(paste0("Power plus Rowwise = ", power.plus.rows),"\n")

cat("","\n")
cat(paste0("Power minus Multiscale = ", power.minus.ms),"\n")
cat(paste0("Power minus Uncorrected = ", power.minus.uncor),"\n")
cat(paste0("Power minus Orderstat = ", power.minus.order),"\n")
cat(paste0("Power minus Rowwise = ", power.minus.rows),"\n")