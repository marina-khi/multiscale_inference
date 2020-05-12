rm(list=ls())
source("functions.r")
source("data_read.r")


alpha <- 0.05  # significance level

Y <- covid_mat
nn <- dim(Y)[2]
TT <- dim(Y)[1]

correction <- correct(Y)
res <- statistics(Y, correction)
test.stat <- max(res$ms)
crit.val <- critical.value(nn, TT, alpha)
test.res <- as.numeric(test.stat > crit.val)
   
print(test.res)


if(nn=2)
{ dev.new()
  par(mfcol=c(3,1))
  
  plot(Y[,1],ylim=c(min(Y),max(Y)),type="l",col="blue", ylab="", xlab="time")
  lines(Y[,2],col="red")
  
  X <- (1:TT)/TT
  lambda.hat <- nw(X,Y[,1],bw=0.025,TT)
  plot(lambda.hat,ylim=c(min(Y),max(Y)),type="l",col="blue", ylab="", xlab="time")
  lambda.hat <- nw(X,Y[,2],bw=0.025,TT)
  lines(lambda.hat,col="red")
  
  pos <- (res$stats[[1]] > crit.val)
  if(sum(pos) > 0)
  { ints <- as.matrix(intervals(TT)[pos, ]) 
    nb.ints <- dim(ints)[1]
    ints <- as.vector(ints)
    ints[ints==0] <- NA
    ints <- matrix(ints,nrow=nb.ints)
    plot(ints[1,], type="l", ylim=c(0,dim(ints)[1]+1), ylab="", xlab="time")
    for(jj in 2:dim(ints)[1])
       lines(jj*ints[jj,])
  }
}
  
  
  
  
  
  