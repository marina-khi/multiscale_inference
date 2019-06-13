rm(list=ls())

source("functions/grid_construction.r")
dyn.load("functions/kernel_weights.so")
source("functions/kernel_weights.r")
source("functions/multiscale_statistics.r")
source("functions/critical_value.r")
source("functions/multiscale_test.r")
source("functions/inputs_for_plots.r")


# Parameters

T          <- 500             # sample size
a1         <- -0.25           # AR parameter 
sigma_eta  <- 1               # standard deviation of the innovation term in the AR model
sim.design <- "spike"         # trend specification: "constant", "blocks", ...

## spike: 0.49-0.51, a1=-0.25,sigma_eta=0.3
## spike: 0-1, a1=-0.25, sigma_eta=1.2

alpha      <- 0.05                  # significance level
SimRuns    <- 1000                  # number of simulation runs to produce critical values
kappa      <- 0.1                   # parameter to determine order statistic for the version
                                    # of the multiscale statistic from Section ?? 
grid       <- grid_construction(T)  # grid of location-bandwidth points 
                                    # the grid should be constructed via the function 
                                    # 'grid construction'. See the corresponding R-file 
                                    # the function definition for details 


# Simulate data

source("simulations/sim.r") 


# Compute test results 

wghts     <- kernel_weights(T=T,grid=grid)
res       <- multiscale_test(alpha=alpha, data=data, weights=wghts, sigmahat=sigma, grid=grid, kappa=kappa, SimRuns=SimRuns)
res.ms    <- res$test_ms
res.uncor <- res$test_uncor
res.order <- res$test_order
res.rows  <- res$test_rows


# Construct SiZer plots

res.ms    <- inputs_for_plots(res.ms, grid)
res.uncor <- inputs_for_plots(res.uncor, grid)
res.order <- inputs_for_plots(res.order, grid)
res.rows  <- inputs_for_plots(res.rows, grid)

u.grid <- res.ms$ugrid
h.grid <- res.ms$hgrid
t.pts  <- (1:T)/T

dev.new()
par(mfrow=c(3,2))
plot(t.pts,data,xlab="",ylab="",type="l")
lines(t.pts,trend,col="red",lwd=2)
plot(t.pts,data,xlab="",ylab="",type="l")
lines(t.pts,trend,col="red",lwd=2)
image(x=u.grid, y=log(h.grid,10), z=t(res.ms$mat), col=res.ms$cols, xlab="u", ylab=expression(log[10](h)), main = "Multiscale")
image(x=u.grid, y=log(h.grid,10), z=t(res.uncor$mat), col=res.uncor$cols, xlab="u", ylab=expression(log[10](h)), main = "Uncorrected")
image(x=u.grid, y=log(h.grid,10), z=t(res.order$mat), col=res.order$cols, xlab="u", ylab=expression(log[10](h)), main = "Order")
image(x=u.grid, y=log(h.grid,10), z=t(res.rows$mat), col=res.rows$cols, xlab="u", ylab=expression(log[10](h)), main = "Rowwise")
