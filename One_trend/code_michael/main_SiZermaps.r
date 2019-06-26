rm(list=ls())

source("functions/grid_construction.r")
source("functions/multiscale_test.r")

dyn.load("functions/SiZer_functions.so")
source("functions/SiZer_functions.r")

source("functions/long_run_variance.r")

source("functions/SiZermap.r")


# Parameters

T          <- 500                   # sample size
a1         <- 0.5                   # AR parameter 
sigma_eta  <- 1                     # standard deviation of the innovation term in the AR model
sim.design <- "bump"                # trend specification: "constant", "blocks", "spike", ...
lrv.type   <- "true"                # long-run error variance: "true" or "estimated"

## spike: 0.49-0.51, a1=-0.25,sigma_eta=0.3
## spike: 0-1, a1=-0.25, sigma_eta=1.2

alpha      <- 0.05                  # significance level
SimRuns    <- 1000                  # number of simulation runs to produce critical values
kappa      <- 0.1                   # parameter to determine order statistic for the version
                                    # of the multiscale statistic from Section ?? 


# Simulate data

source("simulations/sim.r") 


# Construct grid with effective sample size ESS.star(u,h) >= 5 for any (u,h)

grid      <- grid_construction(T)
gset      <- grid$gset
u.grid    <- sort(unique(gset[,1]))
h.grid    <- sort(unique(gset[,2]))
ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
deletions <- ess$del
grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)


# Compute long-run error variance

if(lrv.type == "estimated")
{ AR.struc <- AR_lrv(data=data, q=25, r.bar=10, p=1)    
  a.hat    <- AR.struc$ahat 
  vareta   <- AR.struc$vareta   
  sigma    <- sqrt(AR.struc$lrv)
  autocov  <- AR_acf(coefs=a.hat, var.eta=vareta, len=T)   
}


# Compute multiscale test results 

res       <- multiscale_test(alpha=alpha, data=data, sigmahat=sigma, grid=grid, kappa=kappa, SimRuns=SimRuns)
res.ms    <- res$test_ms
res.uncor <- res$test_uncor
res.order <- res$test_order
res.rows  <- res$test_rows


# Compute SiZer test results

sizer.wghts  <- SiZer_weights(T=T, grid=grid)
sizer.stats  <- SiZer_statistics(data=data, weights=sizer.wghts, autocov=autocov)
sizer.vals   <- sizer.stats$vals
sizer.std    <- sizer.stats$std 
sizer.quants <- SiZer_quantiles(alpha=alpha, T=T, grid=grid, autocov=autocov)
res.sizer    <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)
res.sizer    <- res.sizer$test


# Construct SiZer plots

u.grid <- res$ugrid
h.grid <- res$hgrid
t.pts  <- (1:T)/T

dev.new()
par(mfrow=c(3,2))
plot(t.pts,data,xlab="",ylab="",type="l")
lines(t.pts,trend,col="red",lwd=2)
SiZermap(u.grid=u.grid, h.grid=h.grid, test.results=res.ms, plot.title="Multiscale")
SiZermap(u.grid=u.grid, h.grid=h.grid, test.results=res.uncor, plot.title="Uncorrected")
SiZermap(u.grid=u.grid, h.grid=h.grid, test.results=res.order, plot.title="Order")
SiZermap(u.grid=u.grid, h.grid=h.grid, test.results=res.rows, plot.title="Row-wise")
SiZermap(u.grid=u.grid, h.grid=h.grid, test.results=res.sizer, plot.title="Row-wise SiZer")
