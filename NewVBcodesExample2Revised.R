#Additional supporting functions, e.g., log pdfs
expit <- function(z){
  
  1/(1+exp(-z))
  
}


PosNormal.PosNormal.Gumbel.Objective <- function(x, y1, y2){
  
  mu1 <- x[1]; v1 <- x[2]; mu2 <- x[3]; v2 <- x[4]; omega <- x[5]
  tau <- exp(log(exp(exp(omega)) - 1) - log(exp(exp(omega)) + 1))
  mycop <- gumbelCopula(param = (1/(1-tau)) ,dim=d)
  normFac1 <- ( 1 - pnorm(0,mean=mu1,sd=sqrt(exp(v1))))
  normFac2 <- ( 1 - pnorm(0,mean=mu2,sd=sqrt(exp(v2))))
  pn1 <- (pnorm(y1,mean=mu1,sd=sqrt(exp(v1))) - pnorm(0,mean=mu1,sd=sqrt(exp(v1))))/normFac1
  pn2 <- (pnorm(y2,mean=mu2,sd=sqrt(exp(v2))) - pnorm(0,mean=mu2,sd=sqrt(exp(v2))))/normFac2
  pn1[which(pn1==1)] <- 0.99999; pn2[which(pn2==1)] <- 0.99999
  LogLik <- sum( dCopula( cbind( pn1,  pn2), mycop, log = TRUE) + dnorm(y1, mean = mu1, sd = sqrt(exp(v1)), log=TRUE) - log(normFac1) + dnorm(y2,mean = mu2, sd = sqrt(exp(v2)), log=TRUE ) - log(normFac2) )
  ans <- LogLik + sum(dnorm(x = c(mu1,mu2), mean=0, sd = 100, log = TRUE)) + sum(dLogLogHalfNormal(c(v1,v2), s.sq=10000)) + dLogCorrectOmegaPDF(omega)
  
}

Rank.Gumbelcop.Objective <- function(x, y1, y2){

  n <- length(y1)
  omega <- x
  tau <- exp(log(exp(exp(omega)) - 1) - log(exp(exp(omega)) + 1))
  B1.rank <- rank(y1)/(n+1); B2.rank <- rank(y2)/(n+1); A1.rank <- B1.rank - 1/(n+1); A2.rank <- B2.rank - 1/(n+1)
  cc.local <- gumbelCopula(param = (1/(1-tau)) ,dim=d)
  CopComp <- pCopula(u=cbind(B1.rank, B2.rank), copula = cc.local) - pCopula(u=cbind(B1.rank, A2.rank), copula = cc.local) - pCopula(u=cbind(A1.rank, B2.rank), copula = cc.local) + pCopula(u=cbind(A1.rank, A2.rank), copula = cc.local)
  return(sum(log(CopComp[CopComp>0 ])) + (-99999)*(sum(CopComp==0)>0) + dLogCorrectOmegaPDF(omega))

}

#VB with misspecified gumbel and Correctly Specified Marginals copula. Note that parameterisation is (mean1, log-variance1, mean2, log-variance2, log[log(1+tau) - log(1-tau)]).
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#m.init denote initial values of VB posterior mean
#C.init denote initial values of VB posterior covariance Cholesky
#maxRuns denotes the number of VB draws
NewVBuncut.PosNormal.PosNormal.Gumbel <- function(y1, y2, mu.init, C.init, maxRuns=10000){
  
  
  n <- length(y1); p <- length(mu.init)
  zMAT <- matrix(rnorm( (maxRuns*p) ), nrow=maxRuns, ncol=p)
  ELBOstore <- rep(0, maxRuns)
  gradstore <- matrix(0, nrow=maxRuns, ncol=p)
  mustore <- matrix(0, nrow=maxRuns, ncol=p)
  Cstore <- array(0,dim=c(p,p,maxRuns))
  
  Accugrad.mu <- rep(0,p)
  Change.mu <- rep(0,p)
  Accuchange.mu <- rep(0,p)
  Accugrad.C <- matrix(0,nrow=p,ncol=p)
  Accuchange.C <- matrix(0,nrow=p,ncol=p)
  
  rho <- 0.9
  epsilon <- 10^(-4)
  mu.curr <- mu.init
  C.curr <- C.init
  for(it in 1:maxRuns){
    
    zdraw <- zMAT[it,]
    theta.curr <- c(C.curr %*% zdraw) + mu.curr
    grad.curr <- numDeriv::grad(func=PosNormal.PosNormal.Gumbel.Objective, x=theta.curr, y1=y1, y2=y2)
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon)*grad.curr
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    #browser()
    Key.C <- matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr)) 
    Key.C[row(Key.C) < col(Key.C)] <- 0
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*Key.C^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * Key.C
    C.new <- C.curr + Change.C
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    
    mu.curr <- mu.new
    C.curr <- C.new
    gradstore[it,] <- grad.curr
    Cstore[,,it] <- C.curr
    mustore[it,] <- mu.curr
    if((it%%100)==0){
      
      cat("done with it:", it, "\n")
      
    }
    
    
  }
  
  
  return(list(mustore = mustore, Cstore = Cstore, gradstore=gradstore))
  
}

#VB for correctly specified gumbel copula. Used for Cut II posterior inference. Parameterisation is log[log(1+tau) - log(1-tau)]
#y is a user-input vector for observed Y_1 values.
#m.init denote initial values of VB posterior mean
#c.init denote initial values of VB posterior SD
#maxRuns denotes the number of VB draws
NewVBcut.Gumbel <- function(y1, y2, mu.init, c.init, maxRuns=10000){
  
  
  n <- length(y1)
  zstore <- rnorm( maxRuns )
  ELBOstore <- rep(0, maxRuns)
  gradstore <- rep(0, maxRuns)
  mustore <- rep(0, maxRuns)
  Cstore <- rep(0, maxRuns)
  
  Accugrad.mu <- 0
  Change.mu <- 0
  Accuchange.mu <- 0
  Accugrad.C <- 0
  Accuchange.C <- 0
  
  rho <- 0.85
  epsilon <- 10^(-6)
  mu.curr <- mu.init
  c.curr <- c.init
  for(it in 1:maxRuns){
    
    zdraw <- zstore[it]
    theta.curr <- c.curr * zdraw + mu.curr
    grad.curr <- numDeriv::grad(func=Rank.Gumbelcop.Objective, x=theta.curr, y1=y1, y2=y2)
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon)*grad.curr
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    #browser()
    Key.C <- grad.curr * zdraw + 1/c.curr
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*Key.C^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * Key.C
    c.new <- c.curr + Change.C
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    
    mu.curr <- mu.new
    c.curr <- c.new
    gradstore[it] <- grad.curr
    Cstore[it] <- c.curr
    mustore[it] <- mu.curr
    if((it%%100)==0){
      
      cat("done with it:", it, "\n")
      
    }
    
    
  }
  
  
  return(list(mustore = mustore, Cstore = Cstore, gradstore=gradstore))
  
}

#VB cut II inference with misspecified gumbel and Correctly Specified Marginals copula. Note that parameterisation is (mean1, log-variance1, mean2, log-variance2, log[log(1+tau) - log(1-tau)]).
#y1 is a user-input vector for observed Y_1 values.
#y2 is a user-input vector for observed Y_2 values.
#omega.fix denotes the fixed value of 
#m.init denote initial values of VB posterior mean
#maxRuns denotes the number of VB draws
NewVBcut.PosNormal.PosNormal.Gumbel <- function(y1, y2, mu.init, maxRuns=10000){
  
  
  n <- length(y1); p <- length(mu.init) - 1
  zMAT <- matrix(rnorm( (maxRuns*p) ), nrow=maxRuns, ncol=p)
  ELBOstore <- rep(0, maxRuns)
  gradstore <- matrix(0, nrow=maxRuns, ncol=p)
  mustore <- matrix(0, nrow=maxRuns, ncol=p)
  Cstore <- array(0,dim=c(p,p,maxRuns))
  
  Accugrad.mu <- rep(0,p)
  Change.mu <- rep(0,p)
  Accuchange.mu <- rep(0,p)
  Accugrad.C <- matrix(0,nrow=p,ncol=p)
  Accuchange.C <- matrix(0,nrow=p,ncol=p)
  omega.f <-  mu.init[5]
  
  
  PosNormal.PosNormal.Gumbel.Objective.local <- function(x){
    
    mu1 <- x[1]; v1 <- x[2]; mu2 <- x[3]; v2 <- x[4]
    tau <- exp(log(exp(exp(omega.f)) - 1) - log(exp(exp(omega.f)) + 1))
    mycop <- gumbelCopula(param = (1/(1-tau)) ,dim=d)
    normFac1 <- ( 1 - pnorm(0,mean=mu1,sd=sqrt(exp(v1))))
    normFac2 <- ( 1 - pnorm(0,mean=mu2,sd=sqrt(exp(v2))))
    pn1 <- (pnorm(y1,mean=mu1,sd=sqrt(exp(v1))) - pnorm(0,mean=mu1,sd=sqrt(exp(v1))))/normFac1
    pn2 <- (pnorm(y2,mean=mu2,sd=sqrt(exp(v2))) - pnorm(0,mean=mu2,sd=sqrt(exp(v2))))/normFac2
    pn1[which(pn1==1)] <- 0.99999; pn2[which(pn2==1)] <- 0.99999
    LogLik <- sum( dCopula( cbind( pn1,  pn2), mycop, log = TRUE) + dnorm(y1, mean = mu1, sd = sqrt(exp(v1)), log=TRUE) - log(normFac1) + dnorm(y2,mean = mu2, sd = sqrt(exp(v2)), log=TRUE ) - log(normFac2) )
    ans <- LogLik + sum(dnorm(x = c(mu1,mu2), mean=0, sd = 100, log = TRUE)) + sum(dLogLogHalfNormal(c(v1,v2), s.sq=10000))
    return(ans)
    
  }
  theta.hat <- optim(par = mu.init[1:4], fn = PosNormal.PosNormal.Gumbel.Objective.local, control=list(fnscale=-1))$par
  Proposal.Var <- -solve(pracma::hessian(PosNormal.PosNormal.Gumbel.Objective.local,theta.hat))
  if(sum(diag(Proposal.Var)<0)>0){
    
    Proposal.Var <- matrix(c(0.006145967,-0.0019272762,0.0018868778,-0.0011190089,-0.001927276,0.0020823533,-0.0002670583,0.0013281358,0.001886878,-0.0002670583,0.0008595655,-0.0003071021,-0.001119009,0.0013281358,-0.0003071021 ,0.0014544822),nrow=4,ncol=4,byrow=TRUE)
    
  }
  rho <- 0.95
  epsilon <- 10^(-4)
  mu.curr <- theta.hat
  C.curr <- t(chol(Proposal.Var))
  for(it in 1:maxRuns){
    
    zdraw <- zMAT[it,]
    theta.curr <- c(C.curr %*% zdraw) + mu.curr
    grad.curr <- numDeriv::grad(func=PosNormal.PosNormal.Gumbel.Objective.local, x=theta.curr)
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon)*grad.curr
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    #browser()
    Key.C <- matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr)) 
    Key.C[row(Key.C) < col(Key.C)] <- 0
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*Key.C^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * Key.C
    C.new <- C.curr + Change.C
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    
    mu.curr <- mu.new
    C.curr <- C.new
    gradstore[it,] <- grad.curr
    Cstore[,,it] <- C.curr
    mustore[it,] <- mu.curr
    if((it%%100)==0){
      
      cat("done with it:", it, "\n")
      
    }
    
    
  }
  
  
  return(list(mustore = mustore, Cstore = Cstore, gradstore=gradstore))
  
}