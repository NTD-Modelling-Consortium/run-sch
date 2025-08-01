## Function that takes r0 and k, and fixed parameter z (=e^-gamma)
## Returns unnormalised density for prior
## Hyperparameters
## Mplus - worm burden for which q-th quantile should not exceed
## delta controls smoothness of logistic cut-off
# d5PYtKsLEzm7RJy!
set.seed(NULL)
set.seed(9)
library(mvtnorm)
prior_mu<-log(c(1.5,0.1))#c(1.1,-1.6)
rho<--0.5
s1<-0.5 #0.5
s2<-0.7 #0.7
#prior_acceptances<-0
#prior_samples<-0
z.val <- exp(-0.0007) # Mansoni

#### for breakpoint
f.function <- function(M, z, k){ #eqn 3
  return((1+(1-z)*M/k)^(-(k+1)))
}

theta.function <- function(M, z, k){ # eqn 4
  ans <- 1- ( (1+(1-z)*M/k) / (1+(2-z)*M/2/k) )^(k+1)
  return(ans)
}

critical.M <- function(z,k){ # eqn A.3
  numerator <- (k*( (2-z) / (2*(1-z)) )^(1/(k+2)) - k)
  denominator <- (z-1) * ( (2-z) / (2*(1-z)) )^(1/(k+2)) + (1-z/2)
  ans <- numerator/denominator
  return(ans)
}

critical.R0 <- function(z, k){
  M.crit <- critical.M(z,k)
  theta.val <- theta.function(M.crit, z, k)
  f.val <- f.function(M.crit, z, k)
  ans <- (theta.val*f.val)^(-1)
  return(ans)
}

prior_Sigma<-matrix(c(s1^2,rho*s1*s2,rho*s1*s2,s2^2),2,2)
rprior<-function(n) {
  ret<-exp(mvtnorm::rmvnorm(n,prior_mu,prior_Sigma))
  colnames(ret)<-c("R_0","k")
  wh<-which(ret[,1]<critical.R0(z.val,ret[,2]))
  if (length(wh)>0) {
    #cat("MVN prior resamples:",length(wh),"\n")
    #prior_acceptances<<-prior_acceptances-length(wh)
    ret[wh,]<-rprior(length(wh))
  }
  #prior_acceptances<<-prior_acceptances+n
  #prior_samples<<-prior_samples+n
  return(ret)
}
dprior<-function(x,log=FALSE) {
  #cat("Test conformability:",x,",",prior_mu,";",dim(x),",",dim(prior_mu),"\n")
  if (length(x)!=2) {
    stop("Error: Prior density expects length 2 input.\n")
  } else if (prior_samples==0) {
    stop("Prior samples must be generated to estimate density normalising constant.\n")
  } else if (x[2]<=0 || x[1]<critical.R0(z.val,x[2])) {
    if (log) {
      return(-Inf)
    } else {
      return(0)
    }
  } else {
    if (log) {
      return(mvtnorm::dmvnorm(log(c(x)),prior_mu,prior_Sigma,log=T)-sum(log(x))+log(prior_samples)-log(prior_acceptances))
    } else {
      return(mvtnorm::dmvnorm(log(c(x)),prior_mu,prior_Sigma)/prod(x)*prior_samples/prior_acceptances)
    }
  }
}
prior_mvn<-list(rprior=rprior,dprior=dprior)

prior_acceptances<<-10^7
prior_samples<<-13280021

# samp<-rprior(10000)
# plot_samp = cbind(samp, group=1)
# 
# # Simon Spencer's code to provide proper prior densities.
# r0_max<-20
# k_max<-1.5
# # If in doubt, apply trapezium rule!
# L<-301
# x<-seq(0.01,r0_max,length.out=L)
# y<-seq(0.0001,k_max,length.out=L)
# den<-matrix(1,L,L)
# for (i in 1:L) {
#  for (j in 1:L) {
#    den[i,j]<-dprior(c(x[i],y[j]))
#  }
# }
# nc<-(sum(den)-0.5*sum(den[1,]+den[L,]+den[,1]+den[,L])+0.25*sum(den[1,1]+den[1,L]+den[L,1]+den[L,L]))*(x[2]-x[1])*(y[2]-y[1])
# print(nc)
# 
# 
# library(scales)
# png(paste0("~/Documents/SCH-endgame/r0_k_prior_sch.png"),width=10,height=6,units="in",res=600)
# par(mfrow=c(1,2))
# plot(plot_samp, col=alpha(plot_samp[,3], 0.05), xlim=c(0,max(plot_samp[,1])), ylim=c(0,max(plot_samp[,2])))
# abline(v=1,col="grey")
# contour(x,y,den,xlab="R0",ylab="k", xlim=c(0,max(plot_samp[,1])), ylim=c(0,max(plot_samp[,2])))
# abline(v=1,col="grey")
# dev.off()
