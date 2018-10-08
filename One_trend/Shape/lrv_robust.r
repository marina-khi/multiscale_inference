rm(list=ls())

p          <- 1    # AR order
a_star_vec <- c(-0.95,-0.75,-0.5,-0.25,0.25,0.5,0.75,0.95)  # AR parameters
slope_fac  <- 10   # slope of linear trend = slope_fac * sqrt(1/(1-a_star^2))

T_vec <- c(500)   # sample sizes 
Nsim  <- 1000     # number of simulation runs

tpars <- list()
tpars[[1]] <- c(20,20,5,10,15,25)
tpars[[2]] <- c(20,20,5,15,15,25)
tpars[[3]] <- c(25,25,5,10,20,30)
tpars[[4]] <- c(25,25,5,15,20,30)
tpars[[5]] <- c(30,30,5,10,25,35)
tpars[[6]] <- c(30,30,5,15,25,35)
tpars[[7]] <- c(35,35,5,10,30,40)
tpars[[8]] <- c(35,35,5,15,30,40)
#tpars[[9]] <- c(40,40,5,10,35,45)
#tpars[[10]] <- c(40,40,5,15,35,45)


for(ii in 1:8)
{  temp <- tpars[[ii]]
   L1   <- temp[1]
   L2   <- temp[2]
   K1   <- temp[3]
   K2   <- temp[4]
   M1   <- temp[5]
   M2   <- temp[6]


# Functions

source("lrv_functions.r")


# Simulations

for(T in T_vec)

{  a_mat1 <- a_mat2 <- a_mat_HvK <- a_mat_oracle <- numeric()
   sig_mat1 <- sig_mat2 <- sig_mat_HvK <- sig_mat_oracle <- numeric()
   lrv_mat1 <- lrv_mat2 <- lrv_mat_HvK <- lrv_mat_oracle <- numeric()

   for(a_star in a_star_vec)

   {  slope <- slope_fac * sqrt(1/(1-a_star^2))

      a_vec1 <- a_vec2 <- a_vec_HvK <- a_vec_oracle <- rep(NA,Nsim)
      sig_vec1 <- sig_vec2 <- sig_vec_HvK <- sig_vec_oracle <- rep(NA,Nsim)
      lrv_vec1 <- lrv_vec2 <- lrv_vec_HvK <- lrv_vec_oracle <- rep(NA,Nsim)
      lrv_star <- 1/(1-a_star)^2 
 
      for(nb in 1:Nsim)

      {  set.seed(nb)

         # simulate data
         source("lrv_sim.r")  

         # compute estimators
         a_hat1       <- AR_coef(y_data,L1,L2,rep(0,L2+1),p)
         sig_eta_hat1 <- sigma_eta(y_data,a_hat1,p)
         lrv_hat1     <- sig_eta_hat1/(1-sum(a_hat1))^2 
         correct      <- corrections(a_hat1,sig_eta_hat1,K2+1)
         a_hat2       <- AR_coef(y_data,K1,K2,correct,p)
         sig_eta_hat2 <- sigma_eta(y_data,a_hat2,p)
         lrv_hat2     <- sig_eta_hat2/(1-sum(a_hat2))^2 
         a_vec1[nb]   <- a_hat1
         a_vec2[nb]   <- a_hat2
         sig_vec1[nb] <- sig_eta_hat1
         sig_vec2[nb] <- sig_eta_hat2
         lrv_vec1[nb] <- lrv_hat1
         lrv_vec2[nb] <- lrv_hat2
 
         # compute HvK estimators
         a_vec_HvK[nb]   <- AR_coef_HvK(y_data,M1,M2,p)
         sig_vec_HvK[nb] <- sigma_eta(y_data,a_vec_HvK[nb],p)
         lrv_vec_HvK[nb] <- sig_vec_HvK[nb]/(1-sum(a_vec_HvK[nb]))^2

         # compute oracle estimators
         res_oracle         <- lm(eps[2:T] ~ eps[1:(T-1)] - 1)
         a_vec_oracle[nb]   <- as.vector(res_oracle$coef)
         sig_vec_oracle[nb] <- mean(res_oracle$residuals^2)
         lrv_vec_oracle[nb] <- sig_vec_oracle[nb]/(1-sum(a_vec_oracle[nb]))^2 

         print(nb)
      } 

      # compute MSE values
      a_mat1 <- c(a_mat1,mean((a_vec1-a_star)^2))
      a_mat2 <- c(a_mat2,mean((a_vec2-a_star)^2))
      a_mat_HvK <- c(a_mat_HvK,mean((a_vec_HvK-a_star)^2))
      a_mat_oracle <- c(a_mat_oracle,mean((a_vec_oracle-a_star)^2))

      lrv_mat1 <- c(lrv_mat1,mean((lrv_vec1-lrv_star)^2))
      lrv_mat2 <- c(lrv_mat2,mean((lrv_vec2-lrv_star)^2))
      lrv_mat_HvK <- c(lrv_mat_HvK,mean((lrv_vec_HvK-lrv_star)^2))
      lrv_mat_oracle <- c(lrv_mat_oracle,mean((lrv_vec_oracle-lrv_star)^2))

      # Plot results

      # (a) summary plots

      bounds <- "yes"

      a_min <- min(c(a_vec1,a_vec2,a_vec_HvK,a_vec_oracle))
      a_max <- max(c(a_vec1,a_vec2,a_vec_HvK,a_vec_oracle))
      sig_min <- min(c(sig_vec1,sig_vec2,sig_vec_HvK,sig_vec_oracle))
      sig_max <- max(c(sig_vec1,sig_vec2,sig_vec_HvK,sig_vec_oracle))
      lrv_min <- min(c(lrv_vec1,lrv_vec2,lrv_vec_HvK,lrv_vec_oracle))
      lrv_max <- max(c(lrv_vec1,lrv_vec2,lrv_vec_HvK,lrv_vec_oracle))

      if(bounds == "yes")
      { a_min <- max(c(a_min,-1))
        a_max <- min(c(a_max,1))
        sig_min <- min(c(sig_vec1,sig_vec2,sig_vec_oracle))
        sig_max <- max(c(sig_vec1,sig_vec2,sig_vec_oracle))
        lrv_min <- min(c(lrv_vec1,lrv_vec2,lrv_vec_oracle))
        lrv_max <- max(c(lrv_vec1,lrv_vec2,lrv_vec_oracle))
      }

      sig_star <- 1

      plots <- "no"
      if(plots=="yes")
      { dev.new()
        par(mfrow=c(3,4))
        hist(a_vec_oracle,breaks=100,xlim=c(a_min,a_max),main="oracle")
        abline(v=a_star,col="red",lwd=2)
        hist(a_vec1,breaks=100,xlim=c(a_min,a_max),main="our method 1")
        abline(v=a_star,col="red",lwd=2)
        hist(a_vec2,breaks=100,xlim=c(a_min,a_max),main="our method 2")
        abline(v=a_star,col="red",lwd=2)
	hist(a_vec_HvK,breaks=100,xlim=c(a_min,a_max),main="HvK method")
	abline(v=a_star,col="red",lwd=2)
	hist(sig_vec_oracle,breaks=100,xlim=c(sig_min,sig_max),main="oracle")
	abline(v=sig_star,col="red",lwd=2)
	hist(sig_vec1,breaks=100,xlim=c(sig_min,sig_max),main="our method 1")
	abline(v=sig_star,col="red",lwd=2)
	hist(sig_vec2,breaks=100,xlim=c(sig_min,sig_max),main="our method 2")
	abline(v=sig_star,col="red",lwd=2)
	hist(sig_vec_HvK,breaks=100,xlim=c(sig_min,sig_max),main="HvK method")
	abline(v=sig_star,col="red",lwd=2)
	hist(lrv_vec_oracle,breaks=100,xlim=c(lrv_min,lrv_max),main="oracle")
	abline(v=lrv_star,col="red",lwd=2)
	hist(lrv_vec1,breaks=100,xlim=c(lrv_min,lrv_max),main="our method 1")
	abline(v=lrv_star,col="red",lwd=2)
	hist(lrv_vec2,breaks=100,xlim=c(lrv_min,lrv_max),main="our method 2")
	abline(v=lrv_star,col="red",lwd=2)
	hist(lrv_vec_HvK,breaks=100,xlim=c(lrv_min,lrv_max),main="HvK method")
	abline(v=lrv_star,col="red",lwd=2)
      } 
   }

   # (b) MSE plots
   
   a_mat_min <- min(c(a_mat2,a_mat_HvK,a_mat_oracle))
   a_mat_max <- max(c(a_mat2,a_mat_HvK,a_mat_oracle))

   name_spec <- paste("T=",T,"_slope=",slope_fac,"_(L1,L2,K1,K2,M1,M2)=(",L1,",",L2,",",K1,",",K2,",",M1,",",M2,")",sep="")
   pdfname = paste("Results/MSE_a_",name_spec,".pdf",sep="")
   pdf(pdfname, width=5.25, height=4.16, paper="special")
   par(mar = c(4, 4, 3, 0)) #Margins for each plot
   par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins 
   plot(a_mat2,type="l",lty=1,xaxt='n',ylab="MSE",xlab=expression(a[1]),ylim=c(a_mat_min,0.04))
   points(a_mat2,pch=19)
   lines(a_mat_HvK,lty="dashed",lwd=1.5)
   points(a_mat_HvK,pch=19)
   lines(a_mat_oracle,lty="dotted",lwd=1.5)
   points(a_mat_oracle,pch=19)
   axis(1, at=1:length(a_star_vec), labels=a_star_vec)
   legend( "topright", cex = 1, bty = "n", legend = c(expression(widehat(a)), expression(widehat(a)[HvK]), expression(widehat(a)[oracle])), lty=c("solid","dashed","dotted"),lwd=1.5,y.intersp=1.25)
   #title(main = bquote("T = " ~ .(T) ~ "and" ~ s[beta] == .(slope_fac)),line = 1.5)
   title(main = bquote("q = " *.(L1)* ", (" *underline(r)* ", " *bar(r)* ") = (" *.(K1)* "," *.(K2)* "), (" *m[1]* ", " *m[2]* ") = (" *.(M1)* "," *.(M2)* ")"),line = 1.5)
   dev.off()

   lrv_mat1 <- log(lrv_mat1)
   lrv_mat2 <- log(lrv_mat2)
   lrv_mat_HvK <- log(lrv_mat_HvK)
   lrv_mat_oracle <- log(lrv_mat_oracle)

   lrv_mat_min <- min(c(lrv_mat2,lrv_mat_HvK,lrv_mat_oracle))
   lrv_mat_max <- max(c(lrv_mat2,lrv_mat_HvK,lrv_mat_oracle))

   pdfname = paste("Results/MSE_lrv_",name_spec,".pdf",sep="")
   pdf(pdfname, width=5.25, height=4.16, paper="special")
   par(mar = c(4, 4, 3, 0)) #Margins for each plot
   par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins 
   plot(lrv_mat2,type="l",lty=1,xaxt='n',ylab="log(MSE)",xlab=expression(a[1]),ylim=c(lrv_mat_min,lrv_mat_max))
   points(lrv_mat2,pch=19)
   lines(lrv_mat_HvK,lty="dashed",lwd=1.5)
   points(lrv_mat_HvK,pch=19)
   lines(lrv_mat_oracle,lty="dotted",lwd=1.5)
   points(lrv_mat_oracle,pch=19)
   axis(1, at=1:length(a_star_vec), labels=a_star_vec)
   legend( "topleft", cex = 1, bty = "n", legend = c(expression(widehat(sigma)^2), expression(widehat(sigma)[HvK]^2), expression(widehat(sigma)[oracle]^2)), lty=c("solid","dashed","dotted"),lwd=1.5,y.intersp=1.25)
   #title(main = bquote("T = " ~ .(T) ~ "and" ~ s[beta] == .(slope_fac)),line = 1.5)
   title(main = bquote("q = " *.(L1)* ", (" *underline(r)* ", " *bar(r)* ") = (" *.(K1)* "," *.(K2)* "), (" *m[1]* ", " *m[2]* ") = (" *.(M1)* "," *.(M2)* ")"),line = 1.5)
   dev.off()

}
}
