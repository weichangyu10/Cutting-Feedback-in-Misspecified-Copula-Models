#Additional supporting functions, e.g., log pdfs
expit <- function(z){
  
  1/(1+exp(-z))
  
}

expit <- function(z){
  
  1/(1+exp(-z))
  
}

dLogLogHalfNormal<- function(x, s.sq){
  
  return(-0.5*log(pi*s.sq/2) - 0.5*exp(2*x)/s.sq + x)
  
}

dLogExpHalfCauchy<- function(x, s.sq){
  
  return(-log(pi*sqrt(s.sq)/2) - log(exp(-x) + exp(x)/s.sq))
  
}

dLogCorrectOmegaPDF <- function(x){
  
  return( log(2) + x + exp(x) - 2*log(exp(exp(x)) + 1) )
  
}



LogPDFGumbelCop <- function(F1, F2, tau){
  
  psi <- 1/(1-tau)
  F1tilde <- -log(F1)
  F2tilde <- -log(F2)
  ans <- F1tilde + F2tilde + (psi - 1)* log(F1tilde) + (psi - 1)* log(F2tilde) - (F1tilde^psi + F2tilde^psi)^(1/psi) + log( (psi - 1)*( F1tilde^psi + F2tilde^psi )^(1/psi - 2) + ( F1tilde^psi + F2tilde^psi )^(2/psi - 2) )
  return(ans)
  
}

KLPosNormalLogN <- function(Pos.m1, Pos.var1, LogN.m2, LogN.v2){
  
  KLIntegrand <- function(x){
    
    a1 <- dnorm(x=x, mean=Pos.m1, sd = sqrt(Pos.var1), log = TRUE) - log(1 - pnorm(0, mean = Pos.m1, sd = sqrt(Pos.var1)))
    a2 <- dlnorm(x=x, meanlog = LogN.m2, sdlog = sqrt(LogN.v2), log = TRUE)
    return(exp(a1)*(a1-a2))
    
  }
  grid.val <- seq(0.0001, 10, 0.01)
  Integrand.val <- KLIntegrand(grid.val)
  grid.length <- length(grid.val)
  ans <- sum(0.5*0.01*(Integrand.val[1:(grid.length-1)] + Integrand.val[2:(grid.length)]))
  return(ans)
}

KLPosNormalGamma <- function(Pos.m1, Pos.var1, Gamma.A, Gamma.B){
  
  KLIntegrand <- function(x){
    
    a1 <- dnorm(x=x, mean=Pos.m1, sd = sqrt(Pos.var1), log = TRUE) - log(1 - pnorm(0, mean = Pos.m1, sd = sqrt(Pos.var1)))
    a2 <- dgamma(x=x, shape = Gamma.A, rate = Gamma.B, log = TRUE)
    return(exp(a1)*(a1-a2))
    
  }
  grid.val <- seq(0.0001, 10, 0.01)
  Integrand.val <- KLIntegrand(grid.val)
  grid.length <- length(grid.val)
  ans <- sum(0.5*0.01*(Integrand.val[1:(grid.length-1)] + Integrand.val[2:(grid.length)]))
  return(ans)
}

KL.GAMMA <- function(shape1, rate1, shape2, rate2){
  
  Integrand.Gamma <- function(x){
    
    a1 <- dgamma(x, shape=shape1, rate=rate1, log = TRUE)
    a2 <- dgamma(x, shape=shape2, rate=rate2, log = TRUE)
    return(exp(a1)*(a1-a2))
    
  }
  return(integrate(f=Integrand.Gamma, lower = 0, upper=30)$value)
  
}


#MCMC with misspecified marginals and correctly-specified copula.
#y is a user-input matrix with two columns. First column for observed Y_1 values and second column for observed Y_2 values.
#mu1.init and v1.init denote initial values of mean and log-variance of first zero-truncated Gaussian marginal.
#mu2.init and v2.init denote initial values of mean and log-variance of second zero-truncated Gaussian marginal.
#omega init denotes initial value of log[log(1+tau) - log(1-tau)] , where tau is the Kendall's tau copula parameter.
#B denotes the number of MCMC draws
NewMCMCfull.Normal.Normal.Gumbel <- function(y, mu1.init, v1.init, mu2.init, v2.init, omega.init, B=10000){
  
  n <- nrow(y); d <-  ncol(y)
  toss.uni.vec <- runif(B)
  accept.count <- 0
  theta.guess <- c(mu1.init, v1.init, mu2.init, v2.init, omega.init)
  
  PosNormal.PosNormal.Gumbel.Objective <- function(x){
    
    mu1 <- x[1]; v1 <- x[2]; mu2 <- x[3]; v2 <- x[4]; omega <- x[5]
    tau <- exp(log(exp(exp(omega)) - 1) - log(exp(exp(omega)) + 1))
    mycop <- gumbelCopula(param = (1/(1-tau)) ,dim=d)
    normFac1 <- ( 1 - pnorm(0,mean=mu1,sd=sqrt(exp(v1))))
    normFac2 <- ( 1 - pnorm(0,mean=mu2,sd=sqrt(exp(v2))))
    pn1 <- (pnorm(y[,1],mean=mu1,sd=sqrt(exp(v1))) - pnorm(0,mean=mu1,sd=sqrt(exp(v1))))/normFac1
    pn2 <- (pnorm(y[,2],mean=mu2,sd=sqrt(exp(v2))) - pnorm(0,mean=mu2,sd=sqrt(exp(v2))))/normFac2
    pn1[which(pn1==1)] <- 0.99999; pn2[which(pn2==1)] <- 0.99999
    LogLik <- sum( dCopula( cbind( pn1,  pn2), mycop, log = TRUE) + dnorm(y[,1], mean = mu1, sd = sqrt(exp(v1)), log=TRUE) - log(normFac1) + dnorm(y[,2],mean = mu2, sd = sqrt(exp(v2)), log=TRUE ) - log(normFac2) )
    ans <- LogLik + sum(dnorm(x = c(mu1,mu2), mean=0, sd = 100, log = TRUE)) + sum(dLogLogHalfNormal(c(v1,v2), s.sq=10000)) + dLogCorrectOmegaPDF(omega)
    
  }
  theta.hat <- optim(par=theta.guess, fn = PosNormal.PosNormal.Gumbel.Objective, control=list(fnscale=-1))$par
  Proposal.Var <- -solve(pracma::hessian(PosNormal.PosNormal.Gumbel.Objective,theta.hat))
  if(sum(diag(Proposal.Var)<0)>0){
    
    #theta.hat <- optim(par=theta.guess, fn = PosNormal.PosNormal.Gumbel.Objective, control=list(fnscale=-1), method = "L-BFGS-B")$par
    #Proposal.Var <- -solve(pracma::hessian(PosNormal.PosNormal.Gumbel.Objective,theta.hat))
    Proposal.Var <- matrix(c(1.338165e-02,-2.217045e-03,3.703887e-03,-1.123939e-03,1.345817e-04,-2.217045e-03,8.850499e-04,-4.756793e-04 ,5.943735e-04,1.476397e-04,3.703887e-03,-4.756793e-04 ,1.304693e-03,-4.469364e-04,5.574192e-05,-1.123939e-03,5.943735e-04,-4.469364e-04,1.229792e-03,2.686280e-04,1.345817e-04,1.476397e-04 ,5.574192e-05,2.686280e-04,3.625042e-04),nrow=5,ncol=5,byrow=TRUE)
 
  }
  theta.proposed.MAT <- rmvnorm(B, mean= theta.hat, sigma=(Proposal.Var))
  theta.MAT <- matrix(0,nrow=B,ncol=5)
  theta.curr <- theta.hat
  logPostCurr <- PosNormal.PosNormal.Gumbel.Objective(theta.curr)
  for(t in 1:B){
    
    theta.proposed <- theta.proposed.MAT[t,]
    logPostNew <- PosNormal.PosNormal.Gumbel.Objective(theta.proposed)
    logQratio <- dmvnorm(theta.curr, mean= theta.hat, sigma=Proposal.Var, log = TRUE) - dmvnorm(theta.proposed, mean= theta.hat, sigma=Proposal.Var, log = TRUE)
    accept.prob <- exp(logPostNew - logPostCurr + logQratio)
    if(toss.uni.vec[t] <= accept.prob){
      
      theta.curr <- theta.proposed
      logPostCurr <- logPostNew
      accept.count <- accept.count + 1
      
    }
    theta.MAT[t,] <- theta.curr
    
  }
  
  return(list( mu1.vec = theta.MAT[,1], v1.vec = theta.MAT[,2], mu2.vec = theta.MAT[,3], v2.vec = theta.MAT[,4], tau.vec = exp(log(exp(exp(theta.MAT[,5])) - 1) - log(exp(exp(theta.MAT[,5])) + 1)) , accept.rate = (accept.count/B) ))

  
}

#MCMC for correctly-specified Gumbel.
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#omega init denotes initial value of log[log(1+tau) - log(1-tau)] , where tau is the Kendall's tau copula parameter.
#B denotes the number of MCMC draws
NewMCMCcut.Gumbelcop <-  function(y, omega.init, B=10000){
  
  n <- nrow(y); d <-  ncol(y)
  toss.uni.vec <- runif(B)
  accept.count <- 0
  theta.guess <- omega.init
  Rank.Gumbelcop.Objective <- function(x){
    
    omega <- x
    tau <- exp(log(exp(exp(omega)) - 1) - log(exp(exp(omega)) + 1))
    B1.rank <- rank(y[,1])/(n+1); B2.rank <- rank(y[,2])/(n+1); A1.rank <- B1.rank - 1/(n+1); A2.rank <- B2.rank - 1/(n+1)
    cc.local <- gumbelCopula(param = (1/(1-tau)) ,dim=d)
    CopComp <- pCopula(u=cbind(B1.rank, B2.rank), copula = cc.local) - pCopula(u=cbind(B1.rank, A2.rank), copula = cc.local) - pCopula(u=cbind(A1.rank, B2.rank), copula = cc.local) + pCopula(u=cbind(A1.rank, A2.rank), copula = cc.local)
    return(sum(log(CopComp[CopComp>0 ])) + (-99999)*(sum(CopComp==0)>0) + dLogCorrectOmegaPDF(omega))
    
  }
  theta.hat <- optim(par=theta.guess, fn = Rank.Gumbelcop.Objective, control=list(fnscale=-1), method = "Brent", lower = -20, upper=20)$par
  Proposal.Var <- -1/c(pracma::hessian(Rank.Gumbelcop.Objective,theta.hat))
  theta.proposed.vec <- rnorm(B, mean= theta.hat, sd=sqrt(Proposal.Var))
  theta.vec <- rep(0,B)
  theta.curr <- theta.hat
  logPostCurr <- Rank.Gumbelcop.Objective(theta.curr)
  
  for(t in 1:B){
    
    theta.proposed <- theta.proposed.vec[t]
    logPostNew <- Rank.Gumbelcop.Objective(theta.proposed)
    logQratio <- dnorm(theta.curr, mean= theta.hat, sd=sqrt(Proposal.Var), log = TRUE) - dnorm(theta.proposed, mean= theta.hat, sd=sqrt(Proposal.Var), log = TRUE)
    accept.prob <- exp(logPostNew - logPostCurr + logQratio)
    if(toss.uni.vec[t] <= accept.prob){
      
      theta.curr <- theta.proposed
      logPostCurr <- logPostNew
      accept.count <- accept.count + 1
      
    }
    theta.vec[t] <- theta.curr
    
    
  }
  
  return(list(  tau.vec = exp(log(exp(exp(theta.vec)) - 1) - log(exp(exp(theta.vec)) + 1)), Proposal.Var=Proposal.Var, accept.rate = (accept.count/B) ))
}

#MCMC cut II posterior inference with misspecified marginals and correctly-specified copula.
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#mu1.init and v1.init denote initial values of mean and log-variance of first zero-truncated Gaussian marginal.
#mu2.init and v2.init denote initial values of mean and log-variance of second zero-truncated Gaussian marginal.
#tau.input denotes conditional draw of tau, where tau is the Kendall's tau copula parameter.
#B denotes the number of MCMC draws
NewMCMCcut.Normal.Normal.Gumbel <- function(y, mu1.init, v1.init, mu2.init, v2.init, tau.input, B=10000){
  
  n <- nrow(y); d <-  ncol(y)
  toss.uni.vec <- runif(B)
  accept.count <- 0
  theta.guess <- c(mu1.init, v1.init, mu2.init, v2.init)
  
  PosNormal.PosNormal.Gumbel.Objective <- function(x){
    
    mu1 <- x[1]; v1 <- x[2]; mu2 <- x[3]; v2 <- x[4]
    tau <- tau.input
    mycop <- gumbelCopula(param = (1/(1-tau)) ,dim=d)
    normFac1 <- ( 1 - pnorm(0,mean=mu1,sd=sqrt(exp(v1))))
    normFac2 <- ( 1 - pnorm(0,mean=mu2,sd=sqrt(exp(v2))))
    pn1 <- (pnorm(y[,1],mean=mu1,sd=sqrt(exp(v1))) - pnorm(0,mean=mu1,sd=sqrt(exp(v1))))/normFac1
    pn2 <- (pnorm(y[,2],mean=mu2,sd=sqrt(exp(v2))) - pnorm(0,mean=mu2,sd=sqrt(exp(v2))))/normFac2
    pn1[which(pn1==1)] <- 0.99999; pn2[which(pn2==1)] <- 0.99999
    LogLik <- sum( dCopula( cbind( pn1,  pn2), mycop, log = TRUE) + dnorm(y[,1], mean = mu1, sd = sqrt(exp(v1)), log=TRUE) - log(normFac1) + dnorm(y[,2],mean = mu2, sd = sqrt(exp(v2)), log=TRUE ) - log(normFac2) )
    ans <- LogLik + sum(dnorm(x = c(mu1,mu2), mean=0, sd = 100, log = TRUE)) + sum(dLogLogHalfNormal(c(v1,v2), s.sq=10000))
    
  }
  theta.hat <- optim(par=theta.guess, fn = PosNormal.PosNormal.Gumbel.Objective, control=list(fnscale=-1))$par
  Proposal.Var <- -solve(pracma::hessian(PosNormal.PosNormal.Gumbel.Objective,theta.hat))
  if(sum(diag(Proposal.Var)<0)>0){

    theta.hat <- theta.guess
    Proposal.Var <- matrix(c(2.773107e-06,-4.112522e-09,8.683159e-07,8.092845e-08,-4.112522e-09,2.759551e-06,7.703463e-07,1.629308e-06,8.683159e-07,7.703463e-07,2.351099e-04,-1.992250e-04,8.092845e-08,1.629308e-06,-1.992250e-04 ,7.030746e-04),nrow=4,ncol=4,byrow=TRUE)

  }
  theta.proposed.MAT <- rmvnorm(B, mean= theta.hat, sigma=(Proposal.Var))
  theta.MAT <- matrix(0,nrow=B,ncol=4)
  theta.curr <- theta.hat
  logPostCurr <- PosNormal.PosNormal.Gumbel.Objective(theta.curr)
  for(t in 1:B){
    
    theta.proposed <- theta.proposed.MAT[t,]
    logPostNew <- PosNormal.PosNormal.Gumbel.Objective(theta.proposed)
    logQratio <- dmvnorm(theta.curr, mean= theta.hat, sigma=Proposal.Var, log = TRUE) - dmvnorm(theta.proposed, mean= theta.hat, sigma=Proposal.Var, log = TRUE)
    accept.prob <- exp(logPostNew - logPostCurr + logQratio)
    if(toss.uni.vec[t] <= accept.prob){
      
      theta.curr <- theta.proposed
      logPostCurr <- logPostNew
      accept.count <- accept.count + 1
      
    }
    theta.MAT[t,] <- theta.curr
    
  }
  
  return(list( mu1.vec = theta.MAT[,1], v1.vec = theta.MAT[,2], mu2.vec = theta.MAT[,3], v2.vec = theta.MAT[,4], accept.rate = (accept.count/B) ))
  
  
}

#MCMC cut II posterior inference with correctly-specified marginals and copula.
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#mu.init and v.init denote initial values of mean and log-variance of Lognormal marginal.
#a.init and b.init denote initial values of log-shape and log-rate of Gamma marginal.
#tau.input denotes conditional draw of tau, where tau is the Kendall's tau copula parameter.
#B denotes the number of MCMC draws
NewMCMCcut.LogNormal.Gamma.Gumbel <- function(y, mu.init, v.init, a.init, b.init, tau.input, B=10000){
  
  n <- nrow(y); d <-  ncol(y)
  toss.uni.vec <- runif(B)
  accept.count <- 0
  theta.guess <- c(mu.init, v.init, a.init, b.init)
  
  LogNormal.Gamma.Gumbel.Objective <- function(x){
    
    mu <- x[1]; v <- x[2]; a <- x[3]; b <- x[4]
    tau <- tau.input
    mycop <- gumbelCopula(param = (1/(1-tau)) ,dim=d)
    LogLik <- sum(dCopula( cbind( plnorm(y[,1], meanlog = mu, sdlog = sqrt(exp(v)) ),  pgamma(y[,2], shape = exp(a), rate = exp(b) ) ), mycop, log = TRUE) + dlnorm(y[,1], meanlog = mu, sdlog = sqrt(exp(v)), log = TRUE ) + dgamma(y[,2], shape = exp(a), rate = exp(b), log = TRUE ) )
    ans <- LogLik + dnorm(x = mu, mean=0, sd = 100, log = TRUE) + dLogLogHalfNormal(v, s.sq=10000) + sum(dLogExpHalfCauchy(x = c(a,b), s.sq=5))
    return(ans)
    
  }
  theta.hat <- optim(par=theta.guess, fn = LogNormal.Gamma.Gumbel.Objective, control=list(fnscale=-1))$par
  Proposal.Var <- -solve(pracma::hessian(LogNormal.Gamma.Gumbel.Objective,theta.hat))
  if(sum(diag(Proposal.Var)<0)>0){
    
    theta.hat <- theta.guess
    Proposal.Var <- matrix(c(2.773107e-06,-4.112522e-09,8.683159e-07,8.092845e-08,-4.112522e-09,2.759551e-06,7.703463e-07,1.629308e-06,8.683159e-07,7.703463e-07,2.351099e-04,-1.992250e-04,8.092845e-08,1.629308e-06,-1.992250e-04 ,7.030746e-04),nrow=4,ncol=4,byrow=TRUE)
    
  }
  theta.proposed.MAT <- rmvnorm(B, mean= theta.hat, sigma=(Proposal.Var))
  theta.MAT <- matrix(0,nrow=B,ncol=4)
  theta.curr <- theta.hat
  logPostCurr <- LogNormal.Gamma.Gumbel.Objective(theta.curr)
  for(t in 1:B){
    
    theta.proposed <- theta.proposed.MAT[t,]
    logPostNew <- LogNormal.Gamma.Gumbel.Objective(theta.proposed)
    logQratio <- dmvnorm(theta.curr, mean= theta.hat, sigma=Proposal.Var, log = TRUE) - dmvnorm(theta.proposed, mean= theta.hat, sigma=Proposal.Var, log = TRUE)
    accept.prob <- exp(logPostNew - logPostCurr + logQratio)
    if(toss.uni.vec[t] <= accept.prob){
      
      theta.curr <- theta.proposed
      logPostCurr <- logPostNew
      accept.count <- accept.count + 1
      
    }
    theta.MAT[t,] <- theta.curr
    
  }
  
  return(list( mu.vec = theta.MAT[,1], v.vec = theta.MAT[,2], a.vec = theta.MAT[,3], b.vec = theta.MAT[,4], accept.rate = (accept.count/B) ))
  
  
}

#MCMC posterior inference with correctly-specified marginals and copula.
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#mu.init and v.init denote initial values of mean and log-variance of Lognormal marginal.
#a.init and b.init denote initial values of log-shape and log-rate of Gamma marginal.
#omega init denotes initial value of log[log(1+tau) - log(1-tau)] , where tau is the Kendall's tau copula parameter.
#B denotes the number of MCMC draws
NewMCMCfull.LogNormal.Gamma.Gumbel <- function(y, mu.init, v.init, a.init, b.init, omega.init, B=10000){
  
  n <- nrow(y); d <-  ncol(y)
  toss.uni.vec <- runif(B)
  accept.count <- 0
  theta.guess <- c(mu.init, v.init, a.init, b.init, omega.init)
  
  LogNormal.Gamma.Gumbel.Objective <- function(x){
    
    mu <- x[1]; v <- x[2]; a <- x[3]; b <- x[4]; omega <- x[5]
    tau <- exp(log(exp(exp(omega)) - 1) - log(exp(exp(omega)) + 1))
    mycop <- gumbelCopula(param = (1/(1-tau)) ,dim=d)
    LogLik <- sum(dCopula( cbind( plnorm(y[,1], meanlog = mu, sdlog = sqrt(exp(v)) ),  pgamma(y[,2], shape = exp(a), rate = exp(b) ) ), mycop, log = TRUE) + dlnorm(y[,1], meanlog = mu, sdlog = sqrt(exp(v)), log = TRUE ) + dgamma(y[,2], shape = exp(a), rate = exp(b), log = TRUE ) )
    ans <- LogLik + dnorm(x = mu, mean=0, sd = 100, log = TRUE) + dLogLogHalfNormal(v, s.sq=10000) + sum(dLogExpHalfCauchy(x = c(a,b), s.sq=5)) + dLogCorrectOmegaPDF(omega)
    ans
    
  }
  theta.hat <- optim(par=theta.guess, fn = LogNormal.Gamma.Gumbel.Objective, control=list(fnscale=-1))$par
  Proposal.Var <- -solve(pracma::hessian(LogNormal.Gamma.Gumbel.Objective,theta.hat))
  theta.proposed.MAT <- rmvnorm(B, mean= theta.hat, sigma=(Proposal.Var))
  theta.MAT <- matrix(0,nrow=B,ncol=5)
  theta.curr <- theta.hat
  logPostCurr <- LogNormal.Gamma.Gumbel.Objective(theta.curr)
  
  for(t in 1:B){
    
    theta.proposed <- theta.proposed.MAT[t,]
    logPostNew <- LogNormal.Gamma.Gumbel.Objective(theta.proposed)
    logQratio <- dmvnorm(theta.curr, mean= theta.hat, sigma=Proposal.Var, log = TRUE) - dmvnorm(theta.proposed, mean= theta.hat, sigma=Proposal.Var, log = TRUE)
    accept.prob <- exp(logPostNew - logPostCurr + logQratio)
    if(toss.uni.vec[t] <= accept.prob){
      
      theta.curr <- theta.proposed
      logPostCurr <- logPostNew
      accept.count <- accept.count + 1
      
    }
    theta.MAT[t,] <- theta.curr
    
  }
  
  return(list( mu.vec = theta.MAT[,1], v.vec = theta.MAT[,2], a.vec = theta.MAT[,3], b.vec = theta.MAT[,4], tau.vec = exp(log(exp(exp(theta.MAT[,5])) - 1) - log(exp(exp(theta.MAT[,5])) + 1)) , accept.rate = (accept.count/B) ))
  
}