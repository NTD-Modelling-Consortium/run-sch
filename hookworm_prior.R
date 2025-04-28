### prior for STH hookworm
# from Bayesian results of Truscott et al. (2019) Heterogeneity in transmission parameters of hookworm infection within the baseline data from the TUMIKIA study in Kenya
# assume prevalence is uniform bw 0, 1
# linear relationship bw k and prevalence 
# power law relationship bw R0 and k

inflation.factor <- sqrt(116)
hessianmat <- matrix(c(723.3789, -1436.099, -1436.099,3328.839), 2,2 ) 

#### for breakpoint
f.function <- function(M, z, k){ #eqn 3
  return((1+(1-z)*M/k)^(-(k+1)))
}

theta.function <- function(M, z, k){ # eqn 4
  ans <- 1- ( (1+(1-z)*M/k) / (1+(2-z)*M/k) )^(k+1) 
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

#z.val <- exp(-0.07)  #ascaris
#z.val <- exp(-0.0035) #trichuris
z.val <- exp(-0.02) #hookworm
k_lb <- 0.01
k_ub <- 1
#R_lb <- critical.R0(z.val,k_ub)   # close to breakpoint - this is not actually used
R_lb = 1
R_ub <- 30
#prob_prior <- 0.75

########### functions for sampling from prior
rprop0 <- function(m){
  prev_samples <- rep(NA,m)
  k_samples <- rep(NA,m)
  r_samples <- rep(NA,m)

  i<-0
  while(i<m){
    newprev <- runif(1,0.01,1)
    # error: newk <- rnorm(1, mean=0.77163*newprev, sd=0.4038265*newprev)
    newk <- rnorm(1, mean=0.77163*newprev, sd=0.2019133*newprev)
    if(newk>=k_lb & newk<=k_ub){
      explogr <- log(0.5711742) - 0.3818534*log(newk)
      r.sd <- inflation.factor*sqrt(c(1,log(newk))%*%solve(hessianmat,c(1, log(newk))))
      logr <- rnorm(1, mean=explogr, sd=r.sd)
      newr <- exp(logr) 
      critR <- critical.R0(z.val,newk)
      if (newr>=critR & newr<=R_ub){
        i<-i+1
        prev_samples[i] <- newprev
        k_samples[i] <- newk
        r_samples[i] <- newr
      } 
    }
  }
  ret<-cbind(r_samples,k_samples)
  colnames(ret)<-c("R0","k")
  #hist(prev_samples)
  #Sys.sleep(10)
  return(ret)
}

dprop0 <- function(x,log=FALSE){
  if (length(x)!=2) {stop("dprop0 is expecting a vector of length 2.\n")}
  a<-x[1]
  b<-x[2]
  # a is R0, b is k
  Rcrit <- critical.R0(z.val,b)
  #p(a, b) = p(a|b)p(b) = p(a|b) sum p(b|prev)p(prev)
  #p(a,b|A) = P(a|b) sum p(b|p) p(p) / p(A|b) sum p(B|p)p(p) # A = a accepted, B= b accepted
  if(a<Rcrit | a>R_ub | b>k_ub | b<k_lb) {
    integrate_value=0
  }else{
    #Integrate - Riemann sum over small intervals in P
    L<-1001
    p.samples<-seq(0.01,1,length.out=L)
    k.fitted <- 0.77163*p.samples
    k.sd <- 0.2019133*p.samples #0.4038265*p.samples
    probk <- (p.samples[2]-p.samples[1])*(sum(dnorm(b,mean=k.fitted,sd=k.sd))-0.5*sum(dnorm(b,k.fitted[c(1,L)],k.sd[c(1,L)])))/0.99
    nc_k <- (p.samples[2]-p.samples[1])*(sum(pnorm(k_ub,k.fitted,k.sd)-pnorm(k_lb,k.fitted,k.sd))-
                     0.5*sum(pnorm(k_ub,k.fitted[c(1,L)],k.sd[c(1,L)])-pnorm(k_lb,k.fitted[c(1,L)],k.sd[c(1,L)])))/0.99
    # Draw r|k
    logr <- log(a)
    explogr <- log(0.5711742) - 0.3818534*log(b)
    r.sd <- inflation.factor*sqrt(c(1,log(b))%*%solve(hessianmat,c(1, log(b))))
    problogr <- dnorm(logr, mean=explogr, sd=r.sd)/a
    nc_r <- pnorm(log(R_ub), explogr, r.sd) - pnorm(log(Rcrit),explogr, r.sd)
    integrate_value <- probk*problogr/nc_k/nc_r
  } 
  if (log) {
    return(log(integrate_value))
  } else {
    return(integrate_value)
  }
}
Prior<-list(rprior=rprop0,dprior=dprop0)
#
# Test code - not usually run
#
#sam<-rprop0(500)
#plot(sam)
## If in doubt, apply trapezium rule!
#L<-201
#x<-seq(R_lb,R_ub,length.out=L)
#y<-seq(k_lb,k_ub,length.out=L)
#den<-matrix(NA,L,L)
#for (j in 1:L) {
#  for (i in 1:L) {
#    den[i,j]<-dprop0(c(x[i],y[j]))
#  }
#  #nc_r <- sum(den[,j]-0.5*(den[1,j]+den[L,j]))*(x[2]-x[1])
#  #print(nc_r)
#}
#nc<-(sum(den)-0.5*sum(den[1,]+den[L,]+den[,1]+den[,L])+0.25*sum(den[1,1]+den[1,L]+den[L,1]+den[L,L]))*(x[2]-x[1])*(y[2]-y[1])
#print(nc)

####### Plot samples for SI 
# get data from DW3_Forecasting_Paper
#   plot(df$Prev, df$k, xlim=c(0,1), ylim=c(0,1.1), xlab="Prevalence", ylab="k")
#   m <- 200
#   # run through code of rprop0
#    points(prev_samples, k_samples, col="blue", pch=2)
#    points(df$Prev, df$k)
# # #
# # #
#    plot(df$k, df$R0, xlim=c(0,0.5), xlab="k", ylab="R0")
#    points(k_samples, r_samples, col="blue", pch=2)
#    points(df$k, df$R0)
#
#
# r_samples <- c()
# r.sd.store <- c()
# i <- 1
# while(length(r_samples) < samplesize){
#  explogr <- log(0.5711742) - 0.3818534*log(k_samples[i])
#  r.sd <- inflation.factor*sqrt(c(1,log(k_samples[i]))%*%solve(hessianmat,c(1, log(k_samples[i]))))
#  logr <- rnorm(1, mean=explogr, sd=r.sd)
#  newr <- exp(logr)
#  if(newr > 1){
#  r.sd.store <- c(r.sd.store, r.sd)
#   r_samples <- c(r_samples, newr)
#   i <- i+1
#  }
# }
#

# # points(k_samples, r.sd.store, col="red")
#
# #plot(log(df$k), log(df$R0))
# #points(log(k_samples), log(r_samples), col="blue")
