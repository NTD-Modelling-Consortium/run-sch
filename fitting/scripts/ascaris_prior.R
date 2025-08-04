### Based on equations from Truscott et al. (2019) Heterogeneity in transmission parameters of hookworm infection within the baseline data from the TUMIKIA study in Kenya

#d5PYtKsLEzm7RJy! 
 # Get path from environment variable
 kPathToFittingInputs <- Sys.getenv("PATH_TO_FITTING_INPUTS")
 raw_data_file <- file.path(kPathToFittingInputs, "RawDataForPrior.csv")
 
 if (file.exists(raw_data_file)) {
   d<-read.csv(raw_data_file)
 } else {
   stop("RawDataForPrior.csv not found in PATH_TO_FITTING_INPUTS directory")
 }
 log.likelihood <- function(pars,pos,neg,mu) {
    a<-pars[1]
    b<-pars[2]
    prob.zero <- (1+mu/(a+b*mu))^(-a-b*mu)
    return(-sum(pos*log(1-prob.zero)+neg*log(prob.zero)))
  }
  R0<-function(W,k,z) {
   # return((1+(1-z)*W/k)^(k+1))
    f <- 1/(1+W*(1-z)/k)^(k+1) * (1- ( (1+W*(1-z)/k) / (1+W*(2-z)/k) )^(k+1))
    return(1/f)
  }
  R0.Simon<-function(W,k,z) {
    return((1+W*(1-z)/k)^(k+1)/(1- ( (1+W*(1-z)/k) / (1+W*(2-z)/k) )^(k+1)))
  } 
  sigma.R0<-function(W,k,z) { # derivative of R wrt k
    R.val <- R0(W,k,z)
    h <- (1 + W*(1-z)/k)
    h.prime <- -W*(1-z)/k^2
    g <- (1 + W*(2-z)/k)
    g.prime <- -W*(2-z)/k^2
    ans <- abs(log(h) + (k+1)*h.prime/h - log(1-h/g) - (k+1)/(1-h/g)*(g*h.prime - h*g.prime)/g^2)
    return(R.val*ans)
  }

  #### for breakpoint
  f.function <- function(M, z, k){ #eqn 3
    return((1+(1-z)*M/k)^(-(k+1)))
  }
  
theta.function <- function(M, z, k){ # eqn 4
  #return(1- ( (1+(1-z)*M/k) / (1+(2-z)*M/2/k) )^(k+1)) 
  return(1- ( (1+(1-z)*M/k) / (1+(2-z)*M/k) )^(k+1)) 
}
  
critical.M <- function(z,k){ # eqn A.3
  numerator <- (k*( (2-z) / (2*(1-z)) )^(1/(k+2)) - k)
  denominator <- (z-1) * ( (2-z) / (2*(1-z)) )^(1/(k+2)) + (1-z/2)
  return(numerator/denominator)
}  
critical.R0 <- function(z, k){
  M.crit <- critical.M(z,k)
  theta.val <- theta.function(M.crit, z, k)
  f.val <- f.function(M.crit, z, k)
  return((theta.val*f.val)^(-1))
}    
z.val <- exp(-0.07)  #ascaris
#z.val <- exp(-0.0035) #trichuris
#z.val <- exp(-0.02) #hookworm
#k_lb <- 0.001 original value used to fit all IUs
k_lb <- 0.01 
k_ub <- 0.5
#R_lb <- critical.R0(z.val,k_ub)   # close to breakpoint - this is not actually used
R_lb = 1
R_ub <- 12
prob_prior <- 0.75

# pre-calculate
prior_fit<-optim(c(0.333,0.017),log.likelihood,pos=d$Positive,neg=d$Negative,mu=d$MeanIntensity,hessian=T)  

####### to draw n samples from prior- using raw data d
rprop0 <-function(m) {
  # number of samples from fitted prior 
  n <- rbinom(1, m, prob_prior)
  inflation.factor<-sqrt(28) # sqrt(n-p)
  # maximum likelihood for a and b
  v<-prior_fit$`par`
  hess<-prior_fit$hessian
  k.samples <- rep(NA,m)
  R.samples <- rep(NA,m)
  i<-0
  rejections<-0
  while(i<n) {
    log.W.sample<-runif(1,log(0.01),log(40))
    W.sample<-exp(log.W.sample)
    k.fitted<-v[1]+v[2]*W.sample
    R.fitted<-R0(W.sample,k.fitted,z.val)
    k.sd.sample <- inflation.factor*sqrt(c(1,W.sample)%*%solve(hess,c(1,W.sample)))
    R0.sd.sample<- sigma.R0(W.sample,k.fitted,z.val)*k.sd.sample
    # Draw k|W ~ N(k.mean(W),k.sd^2), where k.sd is inflated by a factor of sqrt(28).    
    k.sample<-rnorm(1,k.fitted,k.sd.sample)
    # Draw R|W ~ N(R.mean(W),R.sd^2)
    R.sample<-rnorm(1,R.fitted,R0.sd.sample)
    if (is.na(k.sample) | is.na(R.sample)) {cat("NA error",W.sample[j],k.fitted,R.fitted,k.sd.sample,R0.sd.sample,k.sample,R.sample,"\n")} 
    #if (k.sample>=k_lb && k.sample<=k_ub && R.sample>=critical.R0(z.val,k.sample) && R.sample<=R_ub) {
    if (k.sample>=k_lb && k.sample<=k_ub && R.sample>=R_lb && R.sample<=R_ub) {
      i<-i+1
      k.samples[i] <- k.sample
      R.samples[i] <- R.sample
    }
  }
  ##### Uniform draws ###############
  if (m-n>0) {
    k.samples[(n+1):m] <- runif(m-n, k_lb, k_ub)
    #critR <- critical.R0(z.val,k.samples[(n+1):m])
    #R.samples[(n+1):m] <- runif(m-n, critR, R_ub)
    R.samples[(n+1):m] <- runif(m-n, R_lb, R_ub)
  }
  ret<-cbind(R.samples,k.samples)
  colnames(ret)<-c("R0","k")
  return(ret)
}


####### gives the joint density of R0 and k
dprop0<-function(x,log=TRUE){
  if (length(x)!=2) {stop("dprop0 is expecting a vector of length 2.\n")}
  a<-x[1]
  b<-x[2]
	# a is R0, b is k
  #if (b<k_lb || b>k_ub || a < critical.R0(z.val,b) || a>R_ub) { # enforce bounds on both parts of prior
  if (b<k_lb || b>k_ub || a < R_lb || a>R_ub) { # enforce bounds on both parts of prior
	  integrate_value<-0
    unif_den<-0
	} else {
 	  # this is an approximation
    inflation.factor<-sqrt(28) # sqrt(n-p)
    v<-prior_fit$`par`
    hess<-prior_fit$hessian
  	#Integrate - Reimann sum over small intervals in M
    L<-1001
    log.W.samples<-seq(log(0.01),log(40),length.out=L)
    W.samples<-exp(log.W.samples)
  	#W.samples<-seq(0.01,40,length.out=L)
 	  k.fitted<-v[1]+v[2]*W.samples
    R.fitted<-R0(W.samples,k.fitted,z.val)
    k.sd.samples<-rep(NA,L)
    R.sd.samples<-rep(NA,L)
    for (i in 1:L) {
      k.sd.samples[i] <- inflation.factor*sqrt(c(1,W.samples[i])%*%solve(hess,c(1,W.samples[i])))
      R.sd.samples[i]<- sigma.R0(W.samples[i],k.fitted[i],z.val)*k.sd.samples[i]
    }
    integrand<-dnorm(b, mean=k.fitted, sd=k.sd.samples)*dnorm(a,R.fitted,R.sd.samples)/W.samples/(log(40)-log(0.01))
    #nc_parts<-(pnorm(R_ub,R.fitted,R.sd.samples)-pnorm(critical.R0(z.val,b),R.fitted,R.sd.samples))*
    nc_parts<-(pnorm(R_ub,R.fitted,R.sd.samples)-pnorm(R_lb,R.fitted,R.sd.samples))*
              (pnorm(k_ub,k.fitted,k.sd.samples)-pnorm(k_lb,k.fitted,k.sd.samples))/W.samples/(log(40)-log(0.01))
    nc <- sum((nc_parts[1:(L-1)]+nc_parts[2:L])/2*(W.samples[2:L]-W.samples[1:(L-1)]))
    integrate_value <- sum((integrand[1:(L-1)]+integrand[2:L])/2*(W.samples[2:L]-W.samples[1:(L-1)]))/nc
    #(W.samples[2]-W.samples[1])*(sum(integrand)-(integrand[1]+integrand[L])/2)/nc
    #unif_den<-dunif(a,critical.R0(z.val,b),R_ub)*dunif(b,k_lb,k_ub)
    unif_den<-dunif(a,R_lb,R_ub)*dunif(b,k_lb,k_ub)
	}	
  if (is.nan(integrate_value) | is.na(integrate_value)) {integrate_value<-0; cat("Prior is nan",a,b,"\n")}
  if (is.nan(unif_den) | is.na(unif_den)) {unif_den<-0; cat("Prior unif is nan",a,b,"\n")}  
	if (log) {
    return(log(prob_prior*integrate_value + (1-prob_prior)*unif_den))
  } else {
    return(prob_prior*integrate_value + (1-prob_prior)*unif_den)
  }
}

dprior<-function(x,log=TRUE) { # Must take a single vector (currently)
  if (log) {
    return(dunif(x[1],R_lb,R_ub,log=TRUE)+dunif(x[2],k_lb,k_ub,log=TRUE))
  } else {
    return(dunif(x[1],R_lb,R_ub)*dunif(x[2],k_lb,k_ub))
  }
}
rprior<-function(n) { # Must produce a matrix (currently)
  ret<-cbind(runif(n,R_lb,R_ub),runif(n,k_lb,k_ub))
  colnames(ret) <- c("R_0","k")
  return(ret)
}
#UniformPrior<-list(rprior=rprior,dprior=dprior)
Prior<-list(rprior=rprop0,dprior=dprop0)
#
# Test code - not usually run
#
#sam<-rprop0(500)
#plot(sam)
## If in doubt, apply trapezium rule!
#L<-101
#x<-seq(R_lb,R_ub,length.out=L)
#y<-seq(k_lb,k_ub,length.out=L)
#den<-matrix(1,L,L)
#for (i in 1:L) {
#  for (j in 1:L) {
#    den[i,j]<-dprop0(c(x[i],y[j]))
#  }
#}
#nc_asc<-(sum(den)-0.5*sum(den[1,]+den[L,]+den[,1]+den[,L])+0.25*sum(den[1,1]+den[1,L]+den[L,1]+den[L,L]))*(x[2]-x[1])*(y[2]-y[1])
#nc_asc2<-sum((den[1:(L-1),1:(L-1)]+den[2:L,1:(L-1)]+den[1:(L-1),2:L]+den[2:L,2:L])/4*(x[2:L]-x[1:(L-1)])*(y[2:L]-y[1:(L-1)]))
#print(nc_asc)
#print(nc_asc2)

#
# Even testier code to check I can integrate over a log uniform distribution
#

#L<-1001
#rW<-function(n) {exp(runif(n,log(0.01),log(40)))}
#dW<-function(W) {
#  return(1/W/(log(40)-log(0.01)))
#}
#W<-rW(100)
#hist(W)
#Wseq<-seq(0.01,40,length.out=L)
#den<-dW(Wseq)
#nc<-(sum(den)-0.5*(den[1]+den[L]))*(Wseq[2]-Wseq[1])
#print(nc)
#logWseq<-seq(log(0.01),log(40),length.out=L)
#Wseq<-exp(logWseq)
#den<-rep(1/(log(40)-log(0.01)),L)/Wseq
#nc<-sum((den[1:(L-1)]+den[2:L])/2*(Wseq[2:L]-Wseq[1:(L-1)]))
#print(nc)
