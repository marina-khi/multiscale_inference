rm(list=ls())
source("functions.r")
source("data_read.r")


alpha <- 0.05  # significance level

Y <- covid_mat
nn <- dim(Y)[2]
TT <- dim(Y)[1]

scaling <- sigma.square(Y)$homo
res <- statistics(Y, scaling)
crit.val <- critical.value(nn, TT, alpha)

pos <- 1
for(i in 1:(nn-1))
{  for(j in (i+1):nn)
   {  pdftitle <- paste0("Results/Results_",countries[i],"_",countries[j],".pdf")
      pdf(pdftitle,width=4,height=9)
      layout(matrix(c(1,2,3),nrow=3),heights=c(0.35,0.325,0.325))
      par(mar=c(5,3,5,1),mgp=c(2.5,1,0))
      # plot raw data
      min.plot <- min(c(Y[,i],Y[,j]))
      max.plot <- max(c(Y[,i],Y[,j]))
      plot(Y[,i], ylim=c(min.plot,max.plot), type="l", col="blue", ylab="", xlab="time (days)")
      title(main=paste0("Comparison of ",countries[i], " and ", countries[j]), line=3.5)
      title(main="(a) time series of new cases per day", font.main=1,line=1.25)
      lines(Y[,j], col="red") 
      legend(x=0, y=0.975*max(c(Y[,i],Y[,j])), legend=c(countries[i],countries[j]), col=c("blue","red"), lty=1)
      # plot smoothed data
      par(mar=c(5,3,2.5,1))
      X <- (1:TT)/TT
      lambda.hat <- nw(X, Y[,i], bw=7/TT, TT)
      plot(lambda.hat, ylim=c(min.plot,max.plot), type="l", col="blue", ylab="", xlab="time (days)")
      title(main="(b) smoothed time series", font.main=1, line=1.25)
      lambda.hat <- nw(X, Y[,j], bw=7/TT, TT)
      lines(lambda.hat, col="red")
      # plot intervals
      pos.ints <- (res$stats[[pos]] > crit.val)
      if(sum(pos.ints) > 0)
      { ints <- intervals(TT)[pos.ints, ]
        if(sum(pos.ints)==1) 
          ints <- matrix(as.vector(intervals(TT)[pos.ints, ]),nrow=1)
        nb.ints <- dim(ints)[1]
        ints <- as.vector(ints)
        ints[ints==0] <- NA
        ints <- matrix(ints,nrow=nb.ints)
        pos.minimal <- rep(0,dim(ints)[1])
        for(ii in 1:(dim(ints)[1]-1))
        {  for(jj in (ii+1):dim(ints)[1])
           {  if((sum(ints[ii,] * ints[jj,], na.rm=TRUE) == sum(ints[ii,], na.rm=TRUE)) & (pos.minimal[ii] != -1))
              { pos.minimal[ii] <- 1
                pos.minimal[jj] <- -1
              }  
           }   
        }  
        linewidth <- ifelse(pos.minimal[1] == 1, 2, 1)
        linecol <- ifelse(pos.minimal[1] == 1, "black", "gray")
        plot(ints[1,], type="l", ylim=c(0,dim(ints)[1]+1), ylab="", xlab="time (days)", yaxt='n', lwd=linewidth, col=linecol)
        title(main="(c) (minimal) intervals produced by the test", font.main=1, line=1.25)
        for(jj in 2:dim(ints)[1])
        {  linewidth <- ifelse(pos.minimal[jj] == 1, 2, 1) 
           linecol <- ifelse(pos.minimal[jj] == 1, "black", "gray") 
           lines(jj*ints[jj,], lwd=linewidth, col=linecol)
        }   
      }
      if(sum(pos.ints) == 0)
      { plot(1:TT, 1:TT, col="white", ylab="", xlab="time", yaxt='n') 
        title(main="(c) (minimal) intervals produced by the test", font.main=1, line=1.25)
      }
      pos <- pos+1
      dev.off()
   }
}  

pos = 1
for(i in 1:(nn-1)) {
  for(j in (i+1):nn) {
    cat(countries[i], " vs ", countries[j], ", stat = ", max(res$stats[[pos]]), "\n")
    pos <- pos+1
  }
}