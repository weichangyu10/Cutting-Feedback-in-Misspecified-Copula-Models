#Additional supporting functions, e.g., log pdfs
expit <- function(z){
  
  1/(1+exp(-z))
  
}

dLogLogHalfNormal<- function(x, s.sq){
  
  return(-0.5*log(pi*s.sq/2) - 0.5*exp(2*x)/s.sq + x)
  
}

dLogExpHalfCauchy<- function(x, s.sq){
  
  return(-log(pi*sqrt(s.sq)/2) - log(exp(-x) + exp(x)/s.sq))
  
}

dLogOmegaPDF <- function(x){
  
  return( log(2) + x + exp(x) - 2*log(exp(exp(x)) + 1) )
  
}

LogPDFGumbelCop <- function(F1, F2, tau){
  
  psi <- 1/(1-tau)
  F1tilde <- -log(F1)
  F2tilde <- -log(F2)
  ans <- F1tilde + F2tilde + (psi - 1)* log(F1tilde) + (psi - 1)* log(F2tilde) - (F1tilde^psi + F2tilde^psi)^(1/psi) + log( (psi - 1)*( F1tilde^psi + F2tilde^psi )^(1/psi - 2) + ( F1tilde^psi + F2tilde^psi )^(2/psi - 2) )
  return(ans)

}

LogPDFTCop <- function(F1, F2, tau){
  
  rho <- sin(0.5*pi*tau)
  R <- matrix(c(1,rho,rho,1),byrow=TRUE, ncol=2)
  StartConst <-  -0.5*log(1 - rho^2) + 2*( lgamma(0.5) - lgamma(2/2)) + lgamma(3/2) - lgamma(1/2) 
  q1 <- qt(F1, df=1); q2 <- qt(F2, df=1)
  dataPart <- -(3/2)*log(1 + (q1^2 + q2^2 - 2*rho*q1*q2)/(1-rho^2)) + log(1 + q1^2) + log(1 + q2^2)
  return( dataPart + StartConst )
  
}

LogNormal.Objective <- function(x, y){
  
  mu <- x[1]; v <- x[2]
  return( sum(dlnorm(y, meanlog = mu, sdlog = sqrt(exp(v)), log = TRUE )) + dnorm(x = mu, mean=0, sd = 100, log = TRUE) + dLogLogHalfNormal(x = v, s.sq=10000 ) )
  
}

Gamma.Objective <- function(x, y){
  
  a <- x[1]; b <- x[2]
  return( sum(dgamma(y, shape = exp(a), rate = exp(b), log = TRUE )) + sum(dLogExpHalfCauchy(x = c(a,b), s.sq=5)))
  
}

LogNormal.Gamma.Gumbel.Objective <- function(x, y){
  
  mu <- x[1]; v <- x[2]; a <- x[3]; b <- x[4]; omega <- x[5]
  tau <- exp(log(exp(exp(omega)) - 1) - log(exp(exp(omega)) + 1))
  LogLik <- sum(LogPDFGumbelCop( F1 = plnorm(y[,1], meanlog = mu, sdlog = sqrt(exp(v)) ), F2 = pgamma(y[,2], shape = exp(a), rate = exp(b) ), tau = tau ) + dlnorm(y[,1], meanlog = mu, sdlog = sqrt(exp(v)), log = TRUE ) + dgamma(y[,2], shape = exp(a), rate = exp(b), log = TRUE ))
  #browser()
  ans <- LogLik + dnorm(x = mu, mean=0, sd = 100, log = TRUE) + dLogLogHalfNormal(v, s.sq=10000) + sum(dLogExpHalfCauchy(x = c(a,b), s.sq=5)) + dLogOmegaPDF(omega)
  
  return( ans )
  
}


#VB with Correctly Specified Marginals and misspecified gumbel copula. Note that parameterisation is (lognormal mean, lognormal log-variance, log-shape, log-rate, logit(0.5*(tau+1))).
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#m.init denote initial values of VB posterior mean
#C.init denote initial values of VB posterior covariance Cholesky
#maxRuns denotes the number of MCMC draws
VB.LogNormal.Gamma.Gumbel.uncut <- function(y, m.init, C.init, maxRuns=10000){
  
  n <- nrow(y); p <- length(m.init)
  zMAT <- matrix(rnorm( (maxRuns*p) ), nrow=maxRuns, ncol=p)
  ELBOstore <- rep(0, maxRuns)
  gradstore <- matrix(0,nrow=maxRuns,ncol=p)
  mstore <- matrix(0,nrow=maxRuns,ncol=p)
  Cstore <- array(0,dim=c(p,p,maxRuns))
  
  Accugrad.m <- rep(0,p)
  Change.m <- rep(0,p)
  Accuchange.m <- rep(0,p)
  Accugrad.C <- matrix(0,nrow = p, ncol=p)
  Accuchange.C <- matrix(0,nrow = p, ncol=p)
  Change.C <- matrix(0,nrow = p, ncol=p)
  
  LowerID <- matrix(1,nrow=p,ncol=p)
  LowerID[row(LowerID) < col(LowerID)] <- 0
  
  WindowSize <- 5000
  
  redraw <- TRUE
  it <- 0
  maxLhat <- -9999999999999
  meanLhatcriterion <- 0
  rho <- 0.85
  epsilon <- 10^(-6)
  m.curr <- m.init
  C.curr <- C.init
  
  
  
  for(it in 1:maxRuns){
    
    zdraw <- zMAT[it,]
    theta.curr <- c(C.curr %*% zdraw) + m.curr
    grad.curr <- grad(func=LogNormal.Gamma.Gumbel.Objective, x=theta.curr, y=y) + c(t(solve(C.curr)) %*% zdraw)
    Accugrad.m <- rho*Accugrad.m + (1-rho)*grad.curr^2
    Change.m <- sqrt(Accuchange.m +  epsilon)/sqrt(Accugrad.m + epsilon) *grad.curr
    m.new <- m.curr + Change.m
    Accuchange.m <- rho*Accuchange.m + (1-rho)*Change.m^2
    
    grad.C <- matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr))
    grad.C <- LowerID*grad.C
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*( grad.C )^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * ( grad.C )
    C.new <- C.curr + Change.C
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    m.curr <- m.new
    C.curr <- C.new
    gradstore[it,] <- grad.curr
    Cstore[,,it] <- C.curr
    mstore[it,] <- m.curr
    
  }
  
  return(list(mstore = mstore, Cstore = Cstore, gradstore=gradstore))
  
}

#VB for correctly specified lognormal marginal. Used for Cut I posterior inference. Parameterisation is (lognormal mean, lognormal log-variance)
#y is a user-input vector for observed Y_1 values.
#m.init denote initial values of VB posterior mean
#C.init denote initial values of VB posterior covariance Cholesky
#maxRuns denotes the number of MCMC draws
VB.Marginal.LogNormal.cut1 <- function(y, m.init, C.init, maxRuns=10000){
  
  n <- length(y); p <- length(m.init)
  zMAT <- matrix(rnorm( (maxRuns*p) ), nrow=maxRuns, ncol=p)
  ELBOstore <- rep(0, maxRuns)
  gradstore <- matrix(0,nrow=maxRuns,ncol=p)
  mstore <- matrix(0,nrow=maxRuns,ncol=p)
  Cstore <- array(0,dim=c(p,p,maxRuns))
  
  Accugrad.m <- rep(0,p)
  Change.m <- rep(0,p)
  Accuchange.m <- rep(0,p)
  Accugrad.C <- matrix(0,nrow = p, ncol=p)
  Accuchange.C <- matrix(0,nrow = p, ncol=p)
  
  it <- 0
  maxLhat <- -9999999999999
  meanLhatcriterion <- 0
  rho <- 0.85
  epsilon <- 10^(-6)
  m.curr <- m.init
  C.curr <- C.init
  
  J.C <- C.curr
  J.C[row(J.C) < col(J.C)] <- 0
  J.C[row(J.C) >= col(J.C)] <- 1
  
  for(it in 1:maxRuns){
    
    zdraw <- zMAT[it,]
    theta.curr <- c(C.curr %*% zdraw) + m.curr
    grad.curr <- grad(func=LogNormal.Objective, x=theta.curr, y=y)
    
    Accugrad.m <- rho*Accugrad.m + (1-rho)*grad.curr^2
    Change.m <- sqrt(Accuchange.m +  epsilon)/sqrt(Accugrad.m + epsilon) *grad.curr
    m.new <- m.curr + Change.m
    Accuchange.m <- rho*Accuchange.m + (1-rho)*Change.m^2
    
    grad.C <- matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr))
    grad.C <- grad.C*J.C
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*( grad.C )^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * ( grad.C )
    C.new <- C.curr + Change.C
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    
    m.curr <- m.new
    C.curr <- C.new
    gradstore[it,] <- grad.curr
    Cstore[,,it] <- C.curr
    mstore[it,] <- m.curr
  }
  
  return(list(mstore = mstore, Cstore = Cstore, gradstore=gradstore))
  
}

#VB for correctly specified gamma marginal. Used for Cut I posterior inference. Parameterisation is (log-shape, log-rate)
#y is a user-input vector for observed Y_2
#m.init denote initial values of VB posterior mean
#C.init denote initial values of VB posterior covariance Cholesky
#maxRuns denotes the number of MCMC draws
VB.Marginal.Gamma.cut1 <- function(y, m.init, C.init, maxRuns=10000){
  
  n <- length(y); p <- length(m.init)
  zMAT <- matrix(rnorm( (maxRuns*p) ), nrow=maxRuns, ncol=p)
  ELBOstore <- rep(0, maxRuns)
  gradstore <- matrix(0,nrow=maxRuns,ncol=p)
  mstore <- matrix(0,nrow=maxRuns,ncol=p)
  Cstore <- array(0,dim=c(p,p,maxRuns))
  
  Accugrad.m <- rep(0,p)
  Change.m <- rep(0,p)
  Accuchange.m <- rep(0,p)
  Accugrad.C <- matrix(0,nrow = p, ncol=p)
  Accuchange.C <- matrix(0,nrow = p, ncol=p)
  
  it <- 0
  maxLhat <- -9999999999999
  meanLhatcriterion <- 0
  rho <- 0.85
  epsilon <- 10^(-6)
  m.curr <- m.init
  C.curr <- C.init
  
  J.C <- C.curr
  J.C[row(J.C) < col(J.C)] <- 0
  J.C[row(J.C) >= col(J.C)] <- 1
  
  for(it in 1:maxRuns){
    
    zdraw <- zMAT[it,]
    theta.curr <- c(C.curr %*% zdraw) + m.curr
    grad.curr <- grad(func=Gamma.Objective, x=theta.curr, y=y)
    
    Accugrad.m <- rho*Accugrad.m + (1-rho)*grad.curr^2
    Change.m <- sqrt(Accuchange.m +  epsilon)/sqrt(Accugrad.m + epsilon) *grad.curr
    m.new <- m.curr + Change.m
    Accuchange.m <- rho*Accuchange.m + (1-rho)*Change.m^2
    
    grad.C <- matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr))
    grad.C <- grad.C*J.C
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*( grad.C )^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * ( grad.C )
    C.new <- C.curr + Change.C
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2

    m.curr <- m.new
    C.curr <- C.new
    gradstore[it,] <- grad.curr
    Cstore[,,it] <- C.curr
    mstore[it,] <- m.curr
  }
  
  return(list(mstore = mstore, Cstore = Cstore, gradstore=gradstore))
  
}

#VB for misspecified gumbel copula. Used for Cut I posterior inference. Parameterisation is (lognormal mean, lognormal log-variance, log-shape, log-rate, logit(0.5*(tau+1)))
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#m.init denote initial values of VB posterior mean
#C.init denote initial values of VB posterior covariance Cholesky
#maxRuns denotes the number of MCMC draws
VB.LogNormal.Gamma.Gumbel.cut1 <- function(y, m.init, C.init, maxRuns=10000){
  
  n <- nrow(y); p <- length(m.init)
  zMAT <- matrix(rnorm( (maxRuns*p) ), nrow=maxRuns, ncol=p)
  ELBOstore <- rep(0, maxRuns)
  gradstore <- matrix(0,nrow=maxRuns,ncol=p)
  mstore <- matrix(0,nrow=maxRuns,ncol=p)
  Cstore <- array(0,dim=c(p,p,maxRuns))
  
  Accugrad.m <- rep(0,p)
  Change.m <- rep(0,p)
  Accuchange.m <- rep(0,p)
  Accugrad.C <- matrix(0,nrow = p, ncol=p)
  Accuchange.C <- matrix(0,nrow = p, ncol=p)
  
  WindowSize <- 5000
  
  redraw <- TRUE
  it <- 0
  maxLhat <- -9999999999999
  meanLhatcriterion <- 0
  rho <- 0.85
  epsilon <- 10^(-6)
  m.curr <- m.init
  C.curr <- C.init
  
  J.C <- matrix(0,nrow=p,ncol=p)
  J.C[p,] <- 1
  
  for(it in 1:maxRuns){
    
    zdraw <- zMAT[it,]
    theta.curr <- c(C.curr %*% zdraw) + m.curr
    grad.curr <- grad(func=LogNormal.Gamma.Gumbel.Objective, x=theta.curr, y=y)
    grad.curr.m <- c(rep(0,4),grad.curr[5])
    Accugrad.m <- rho*Accugrad.m + (1-rho)*grad.curr.m^2
    Change.m <- sqrt(Accuchange.m +  epsilon)/sqrt(Accugrad.m + epsilon) *grad.curr.m
    m.new <- m.curr + Change.m
    Accuchange.m <- rho*Accuchange.m + (1-rho)*Change.m^2
    
    grad.C <- matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr))
    grad.C <- grad.C*J.C
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*( grad.C )^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * ( grad.C )
    C.new <- C.curr + Change.C
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2

    m.curr <- m.new
    C.curr <- C.new
    gradstore[it,] <- grad.curr
    Cstore[,,it] <- C.curr
    mstore[it,] <- m.curr
    
  }
  
  return(list(mstore = mstore, Cstore = Cstore, gradstore=gradstore))
  
}

#Computes KL-divergence between Gamma(shape1, rate1) and Gamma(shape2, rate2)
KL.GAMMA <- function(shape1, rate1, shape2, rate2){
  
  Integrand.Gamma <- function(x){
    
    a1 <- dgamma(x, shape=shape1, rate=rate1, log = TRUE)
    a2 <- dgamma(x, shape=shape2, rate=rate2, log = TRUE)
    return(exp(a1)*(a1-a2))
    
  }
  return(integrate(f=Integrand.Gamma, lower = 0, upper=30)$value)
  
}