rm(list=ls())

source("functions/grid_construction.r")
dyn.load("functions/kernel_weights.so")
source("functions/kernel_weights.r")
source("functions/multiscale_statistics.r")
source("functions/multiscale_quantiles.r")
source("functions/multiscale_testing.r")

dyn.load("functions/SiZer_functions.so")
source("functions/SiZer_functions.r")

source("functions/long_run_variance.r")


# Parameters

T          <- 500                   # sample size
a1         <- -0.5                  # AR parameter 
sigma_eta  <- 1                     # standard deviation of the innovation term in the AR model
sim.design <- "bump"                # trend specification: "constant", "line", "blocks", 
                                    # "spike", "bump", ...
lrv.type   <- "estimated"           # long-run error variance: "true" or "estimated"

## line, slope=1, a1=-0.5, sigma_eta=1 (T=500)
## spike: 0.45-0.55 (height=), a1=-0.5, sigma_eta=1 (T=500)

alpha      <- 0.05                  # significance level
SimRuns    <- 1000                  # number of simulation runs to produce critical values
kappa      <- 0.1                   # parameter to determine order statistic for the version
                                    # of the multiscale statistic from Section ?? 

Nsim       <- 1000                  # number of simulation runs for size/power calculations


# Construct grid with effective sample size ESS.star(u,h) >= 5 for any (u,h)

# run simulation file in order to compute autocovariance function "autocov" 
# of error process which is needed to calculate ESS.star and some of the parameters
# of the simulation setup
int.plus  <- NULL
int.minus <- NULL
source("simulations/sim.r")  

grid      <- grid_construction(T)
gset      <- grid$gset
u.grid    <- sort(unique(gset[,1]))
h.grid    <- sort(unique(gset[,2]))
ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
deletions <- ess$del
grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)


# Compute kernel weights and critical value for multiscale test  

wghts  <- kernel_weights(T=T, grid=grid)
quants <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=kappa, SimRuns=SimRuns)


# Compute kernel weights for row-wise SiZer

sizer.wghts  <- SiZer_weights(T=T, grid=grid)


# Construct subgrids for size/power calculations

gset        <- grid$gset
N           <- dim(gset)[1]
gset.full   <- grid$gset_full
h.grid.full <- sort(unique(gset.full[,2]))
h.len       <- length(h.grid.full)
u.len       <- length(unique(gset.full[,1]))
pos.full    <- grid$pos_full

u.left  <- (gset[,1] - gset[,2])
u.right <- (gset[,1] + gset[,2])

if(is.null(int.plus)){
  pos.plus <- matrix(1, ncol=u.len, nrow=h.len)
} else { 
  temp <- rep(1,N)
  temp[u.left > int.plus[2] | u.right < int.plus[1]] <- 100 
  pos.plus <- rep(NA,length(pos.full)) 
  pos.plus[!is.na(pos.full)] <- temp
  pos.plus <- matrix(pos.plus, ncol=u.len, byrow=TRUE)
}
if(is.null(int.minus)){
  pos.minus <- matrix(1, ncol=u.len, nrow=h.len) 
} else {
  temp <- rep(1,N)
  temp[u.left > int.minus[2] | u.right < int.minus[1]] <- 100 
  pos.minus <- rep(NA,length(pos.full)) 
  pos.minus[!is.na(pos.full)] <- temp
  pos.minus <- matrix(pos.minus, ncol=u.len, byrow=TRUE)
}


# Size/power simulations

power.rows.pm    <- matrix(0,ncol=h.len,nrow=5)
power.rows.plus  <- matrix(0,ncol=h.len,nrow=5)
power.rows.minus <- matrix(0,ncol=h.len,nrow=5)

power.global.pm    <- rep(0,5)
power.global.plus  <- rep(0,5)
power.global.minus <- rep(0,5)

spurious.rows.pm    <- matrix(0,ncol=h.len,nrow=5)
spurious.rows.plus  <- matrix(0,ncol=h.len,nrow=5)
spurious.rows.minus <- matrix(0,ncol=h.len,nrow=5)

spurious.global.pm    <- rep(0,5)
spurious.global.plus  <- rep(0,5)
spurious.global.minus <- rep(0,5)

for(loops in 1:Nsim)

{  source("simulations/sim.r")

   if(lrv.type == "estimated")
   { AR.struc <- AR_lrv(data=data, q=25, r.bar=10, p=1)    
     a.hat    <- AR.struc$ahat 
     vareta   <- AR.struc$vareta   
     sigma    <- sqrt(AR.struc$lrv)
     autocov  <- AR_acf(coefs=a.hat, var.eta=vareta, len=T)
   }

   stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma, grid=grid) 
   vals     <- stats$values
   test.res <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
 
   sizer.stats  <- SiZer_statistics(data=data, weights=sizer.wghts, autocov=autocov)
   sizer.vals   <- sizer.stats$vals
   sizer.std    <- sizer.stats$std
   sizer.quants <- SiZer_quantiles(alpha=alpha, T=T, grid=grid, autocov=autocov)
   sizer.res    <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)
   
   test.res[[7]] <- sizer.res$test

   for(k in 1:5)
   {  temp.plus  <- test.res[[k+2]]
      temp.plus  <- temp.plus * pos.plus
      temp.plus  <- (temp.plus == 1)
      temp.minus <- test.res[[k+2]]
      temp.minus <- temp.minus * pos.minus
      temp.minus <- (temp.minus == -1)
    
      temp.rows.plus    <- rowSums(temp.plus,na.rm=TRUE)
      temp.global.plus  <- sum(as.vector(temp.plus),na.rm=TRUE)
      temp.rows.minus   <- rowSums(temp.minus,na.rm=TRUE)
      temp.global.minus <- sum(as.vector(temp.minus),na.rm=TRUE)

      temp <- (temp.rows.minus > 0 | temp.rows.plus > 0)
      power.rows.pm[k,] <- power.rows.pm[k,] + temp 
      temp <- (temp.global.minus > 0 | temp.global.plus > 0)
      power.global.pm[k] <- power.global.pm[k] + temp 

      temp <- (temp.rows.plus > 0)
      power.rows.plus[k,] <- power.rows.plus[k,] + temp 
      temp <- (temp.global.plus > 0)
      power.global.plus[k] <- power.global.plus[k] + temp 

      temp <- (temp.rows.minus > 0)
      power.rows.minus[k,] <- power.rows.minus[k,] + temp 
      temp <- (temp.global.minus > 0)
      power.global.minus[k] <- power.global.minus[k] + temp 
   }

   for(k in 1:5)
   {  temp.plus  <- test.res[[k+2]]
      temp.plus  <- temp.plus * pos.plus
      temp.plus  <- (temp.plus == 100)
      temp.minus <- test.res[[k+2]]
      temp.minus <- temp.minus * pos.minus
      temp.minus <- (temp.minus == -100)
    
      temp.rows.plus    <- rowSums(temp.plus,na.rm=TRUE)
      temp.global.plus  <- sum(as.vector(temp.plus),na.rm=TRUE)
      temp.rows.minus   <- rowSums(temp.minus,na.rm=TRUE)
      temp.global.minus <- sum(as.vector(temp.minus),na.rm=TRUE)

      temp <- (temp.rows.minus > 0 | temp.rows.plus > 0)
      spurious.rows.pm[k,] <- spurious.rows.pm[k,] + temp 
      temp <- (temp.global.minus > 0 | temp.global.plus > 0)
      spurious.global.pm[k] <- spurious.global.pm[k] + temp 

      temp <- (temp.rows.plus > 0)
      spurious.rows.plus[k,] <- spurious.rows.plus[k,] + temp 
      temp <- (temp.global.plus > 0)
      spurious.global.plus[k] <- spurious.global.plus[k] + temp 

      temp <- (temp.rows.minus > 0)
      spurious.rows.minus[k,] <- spurious.rows.minus[k,] + temp 
      temp <- (temp.global.minus > 0)
      spurious.global.minus[k] <- spurious.global.minus[k] + temp 
   }

   print(loops)  
}

power.rows.pm    <- power.rows.pm/Nsim
power.rows.plus  <- power.rows.plus/Nsim
power.rows.minus <- power.rows.minus/Nsim

power.global.pm    <- power.global.pm/Nsim
power.global.plus  <- power.global.plus/Nsim
power.global.minus <- power.global.minus/Nsim

spurious.rows.pm    <- spurious.rows.pm/Nsim
spurious.rows.plus  <- spurious.rows.plus/Nsim
spurious.rows.minus <- spurious.rows.minus/Nsim

spurious.global.pm    <- spurious.global.pm/Nsim
spurious.global.plus  <- spurious.global.plus/Nsim
spurious.global.minus <- spurious.global.minus/Nsim

# delete bandwidth rows for which EES.star < 5 at any location
keep.rows <- rowSums(is.na(matrix(pos.full,ncol=u.len,byrow=TRUE)))
keep.rows <- (keep.rows < u.len)

power.rows.pm    <- power.rows.pm[,keep.rows]
power.rows.plus  <- power.rows.plus[,keep.rows]
power.rows.minus <- power.rows.minus[,keep.rows]

spurious.rows.pm    <- spurious.rows.pm[,keep.rows]
spurious.rows.plus  <- spurious.rows.plus[,keep.rows]
spurious.rows.minus <- spurious.rows.minus[,keep.rows]

h.grid <- h.grid.full[keep.rows]


# Plot results

cols <- c("red","blue","yellow","darkgreen","purple")
plot.min <- min(c(as.vector(power.rows.pm),as.vector(power.rows.plus),as.vector(power.rows.minus)))
plot.max <- max(c(as.vector(power.rows.pm),as.vector(power.rows.plus),as.vector(power.rows.minus)))

pdf(file="plots/power_comparison.pdf", width=8, height=10, paper="special") 
par(mfcol=c(3,2))
#
#pdf(file="plots/power_comparison.pdf", width=6, height=6, paper="special") 
plot(h.grid, power.rows.pm[1,], type="l", lty=1, col=cols[1], lwd=1.5, ylim=c(plot.min,plot.max), main="power", xlab="bandwidth h", ylab="")
points(h.grid, power.rows.pm[1,] ,col=cols[1], pch=19)
for(k in c(2,4,5))
{  lines(h.grid, power.rows.pm[k,], lty=1, col=cols[k], lwd=1.5)
   points(h.grid, power.rows.pm[k,], col=cols[k], pch=19)  
}
legend(0.01,plot.max-0.01,c(expression(T[MS]), expression(T[UC]), expression(T[RW]), expression(T[SiZer])), col = cols[c(1,2,4,5)], lty = rep(1,4), lwd=2, bg="white")
#dev.off()
#
plot(h.grid, power.rows.plus[1,], type="l", lty=1, col=cols[1], lwd=1.5, ylim=c(plot.min,plot.max), main="power (increases)", xlab="bandwidth h", ylab="")
points(h.grid, power.rows.plus[1,] ,col=cols[1], pch=19)
for(k in c(2,4,5))
{  lines(h.grid, power.rows.plus[k,], lty=1, col=cols[k], lwd=1.5)
   points(h.grid, power.rows.plus[k,], col=cols[k], pch=19)  
}
legend(0.01,plot.max-0.01,c(expression(T[MS]), expression(T[UC]), expression(T[RW]), expression(T[SiZer])), col = cols[c(1,2,4,5)], lty = rep(1,4), lwd=2, bg="white")
#
plot(h.grid, power.rows.minus[1,], type="l", lty=1, col=cols[1], lwd=1.5, ylim=c(plot.min,plot.max), main="power (decreases)", xlab="bandwidth h", ylab="")
points(h.grid, power.rows.minus[1,] ,col=cols[1], pch=19)
for(k in c(2,4,5))
{  lines(h.grid, power.rows.minus[k,], lty=1, col=cols[k], lwd=1.5)
   points(h.grid, power.rows.minus[k,], col=cols[k], pch=19)  
}
legend(0.01,plot.max-0.01,c(expression(T[MS]), expression(T[UC]), expression(T[RW]), expression(T[SiZer])), col = cols[c(1,2,4,5)], lty = rep(1,4), lwd=2, bg="white")
#
plot.min <- min(c(as.vector(spurious.rows.pm),as.vector(spurious.rows.plus),as.vector(spurious.rows.minus)))
plot.max <- max(c(as.vector(spurious.rows.pm),as.vector(spurious.rows.plus),as.vector(spurious.rows.minus)))
#
plot(h.grid, spurious.rows.pm[1,], type="l", lty=1, col=cols[1], lwd=1.5, ylim=c(plot.min,plot.max), main="spurious power", xlab="bandwidth h", ylab="")
points(h.grid, spurious.rows.pm[1,] ,col=cols[1], pch=19)
for(k in c(2,4,5))
{  lines(h.grid, spurious.rows.pm[k,], lty=1, col=cols[k], lwd=1.5)
   points(h.grid, spurious.rows.pm[k,], col=cols[k], pch=19)  
}
legend(0.01,plot.max-0.01,c(expression(T[MS]), expression(T[UC]), expression(T[RW]), expression(T[SiZer])), col = cols[c(1,2,4,5)], lty = rep(1,4), lwd=2, bg="white")
#
plot(h.grid, spurious.rows.plus[1,], type="l", lty=1, col=cols[1], lwd=1.5, ylim=c(plot.min,plot.max), main="spurious power (increases)", xlab="bandwidth h", ylab="")
points(h.grid, spurious.rows.plus[1,] ,col=cols[1], pch=19)
for(k in c(2,4,5))
{  lines(h.grid, spurious.rows.plus[k,], lty=1, col=cols[k], lwd=1.5)
   points(h.grid, spurious.rows.plus[k,], col=cols[k], pch=19)  
}
legend(0.01,plot.max-0.01,c(expression(T[MS]), expression(T[UC]), expression(T[RW]), expression(T[SiZer])), col = cols[c(1,2,4,5)], lty = rep(1,4), lwd=2, bg="white")
#
plot(h.grid, spurious.rows.minus[1,], type="l", lty=1, col=cols[1], lwd=1.5, ylim=c(plot.min,plot.max), main="spurious power (decreases)", xlab="bandwidth h", ylab="")
points(h.grid, spurious.rows.minus[1,] ,col=cols[1], pch=19)
for(k in c(2,4,5))
{  lines(h.grid, spurious.rows.minus[k,], lty=1, col=cols[k], lwd=1.5)
   points(h.grid, spurious.rows.minus[k,], col=cols[k], pch=19)  
}
legend(0.01,plot.max-0.01,c(expression(T[MS]), expression(T[UC]), expression(T[RW]), expression(T[SiZer])), col = cols[c(1,2,4,5)], lty = rep(1,4), lwd=2, bg="white")
#
dev.off()


































