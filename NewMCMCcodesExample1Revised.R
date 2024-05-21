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

dLogOmegaCorrectPDF <- function(x){
  
  return( -log(1+exp(-x)) - log(1+exp(x)) )
  
}

LogPDFTCop <- function(F1, F2, tau){
  
  rho <- sin(0.5*pi*tau)
  R <- matrix(c(1,rho,rho,1),byrow=TRUE, ncol=2)
  StartConst <-  -0.5*log(1 - rho^2) + 2*( lgamma(0.5) - lgamma(2/2)) + lgamma(3/2) - lgamma(1/2) 
  q1 <- qt(F1, df=1); q2 <- qt(F2, df=1)
  dataPart <- -(3/2)*log(1 + (q1^2 + q2^2 - 2*rho*q1*q2)/(1-rho^2)) + log(1 + q1^2) + log(1 + q2^2)
  return( dataPart + StartConst )
  
}


#MCMC with Correctly Specified Marginals and Copula
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#mu.init and v.init denote initial values of mean and log-variance of lognormal marginal.
#a.init and b.init denote initial values of log-shape and log-rate of gamma marginal.
#omega init denotes initial value of logit(0.5*(tau+1)), where tau is the Kendall's tau copula parameter.
#B denotes the number of MCMC draws
NewMCMCfull.LogNormal.Gamma.T <- function(y, mu.init, v.init, a.init, b.init, omega.init, B=10000){
  
  n <- nrow(y); d <-  ncol(y)
  toss.uni.vec <- runif(B)
  accept.count <- 0
  theta.guess <- c(mu.init, v.init, a.init, b.init, omega.init)
  
  LogNormal.Gamma.T.Objective <- function(x){
    
    mu <- x[1]; v <- x[2]; a <- x[3]; b <- x[4]; omega <- x[5]
    tau <- 2*expit(omega) - 1
    mycop <- tCopula(param = sin(0.5*pi*tau), dim=d, df=1 )
    LogLik <- sum(dCopula( cbind( plnorm(y[,1], meanlog = mu, sdlog = sqrt(exp(v)) ),  pgamma(y[,2], shape = exp(a), rate = exp(b) ) ), mycop, log = TRUE) + dlnorm(y[,1], meanlog = mu, sdlog = sqrt(exp(v)), log = TRUE ) + dgamma(y[,2], shape = exp(a), rate = exp(b), log = TRUE ) )
    ans <- LogLik + dnorm(x = mu, mean=0, sd = 100, log = TRUE) + dLogLogHalfNormal(v, s.sq=10000) + sum(dLogExpHalfCauchy(x = c(a,b), s.sq=5)) + dLogOmegaCorrectPDF(omega)
    ans
    
  }
  theta.hat <- optim(par=theta.guess, fn = LogNormal.Gamma.T.Objective, control=list(fnscale=-1))$par
  Proposal.Var <- -solve(pracma::hessian(LogNormal.Gamma.T.Objective,theta.hat))
  theta.proposed.MAT <- rmvnorm(B, mean= theta.hat, sigma=(Proposal.Var))
  theta.MAT <- matrix(0,nrow=B,ncol=5)
  theta.curr <- theta.hat
  logPostCurr <- LogNormal.Gamma.T.Objective(theta.curr)
  
  for(t in 1:B){
    
    theta.proposed <- theta.proposed.MAT[t,]
    logPostNew <- LogNormal.Gamma.T.Objective(theta.proposed)
    logQratio <- dmvnorm(theta.curr, mean= theta.hat, sigma=Proposal.Var, log = TRUE) - dmvnorm(theta.proposed, mean= theta.hat, sigma=Proposal.Var, log = TRUE)
    accept.prob <- exp(logPostNew - logPostCurr + logQratio)
    if(toss.uni.vec[t] <= accept.prob){
      
      theta.curr <- theta.proposed
      logPostCurr <- logPostNew
      accept.count <- accept.count + 1
      
    }
    theta.MAT[t,] <- theta.curr
    
  }
  
  return(list( mu.vec = theta.MAT[,1], v.vec = theta.MAT[,2], a.vec = theta.MAT[,3], b.vec = theta.MAT[,4], tau.vec = (2*expit(theta.MAT[,5])-1) , accept.rate = (accept.count/B) ))
  
}

#MCMC with Correctly Specified Marginals and misspecified gumbel copula.
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#mu.init and v.init denote initial values of mean and log-variance of lognormal marginal.
#a.init and b.init denote initial values of log-shape and log-rate of gamma marginal.
#omega init denotes initial value of logit(0.5*(tau+1)), where tau is the Kendall's tau copula parameter.
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
    ans <- LogLik + dnorm(x = mu, mean=0, sd = 100, log = TRUE) + dLogLogHalfNormal(v, s.sq=10000) + sum(dLogExpHalfCauchy(x = c(a,b), s.sq=5)) + dLogOmegaPDF(omega)
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

#MCMC for correctly specified lognormal marginal. Used for Cut I posterior inference.
#y is a user-input vector for observed Y_1 values.
#mu.init and v.init denote initial values of mean and log-variance of lognormal marginal.
#B denotes the number of MCMC draws.
NewMCMCcutMarginal.LogNormal <- function(y, mu.init, v.init, B=10000){
  
  n <- length(y)
  y <- log(y)
  toss.uni.vec <- runif(B)
  accept.count <- 0
  theta.guess <- c(mu.init, v.init)
  Normal.Objective <- function(x){
    
    mu <- x[1]; v <- x[2]
    LogLik <- sum(dnorm(y, mean = mu, sd = sqrt(exp(v)), log = TRUE ) )
    ans <- LogLik + dnorm(mu, mean=0, sd=100, log = TRUE) + dLogLogHalfNormal(x = v, s.sq=10000)
    ans
    
  }
  theta.hat <- optim(par=theta.guess, fn = Normal.Objective, control=list(fnscale=-1))$par
  Proposal.Var <- -solve(pracma::hessian(Normal.Objective,theta.hat))
  theta.proposed.MAT <- rmvnorm(B, mean= theta.hat, sigma=(Proposal.Var))
  theta.MAT <- matrix(0,nrow=B,ncol=2)
  theta.curr <- theta.hat
  logPostCurr <- Normal.Objective(theta.curr)
  for(t in 1:B){
    
    theta.proposed <- theta.proposed.MAT[t,]
    logPostNew <- Normal.Objective(theta.proposed)
    logQratio <- dmvnorm(theta.curr, mean= theta.hat, sigma=Proposal.Var, log = TRUE) - dmvnorm(theta.proposed, mean= theta.hat, sigma=Proposal.Var, log = TRUE)
    accept.prob <- exp(logPostNew - logPostCurr + logQratio)
    if(toss.uni.vec[t] <= accept.prob){
      
      theta.curr <- theta.proposed
      logPostCurr <- logPostNew
      accept.count <- accept.count + 1
      
    }
    theta.MAT[t,] <- theta.curr
    
  }
  return(list( mu.vec = theta.MAT[,1], v.vec = theta.MAT[,2], accept.rate = (accept.count/B)))
  
  
}

#MCMC for correctly specified gamma marginal. Used for Cut I posterior inference.
#y is a user-input vector for observed Y_2 values.
#a.init and b.init denote initial values of log-shape and log-rate of gamma marginal.
#B denotes the number of MCMC draws
NewMCMCcutMarginal.Gamma <- function(y, a.init, b.init, B=10000){
  
  n <- length(y)
  toss.uni.vec <- runif(B)
  accept.count <- 0
  theta.guess <- c(a.init, b.init)
  Gamma.Objective <- function(x){
    
    a <- x[1]; b <- x[2]
    LogLik <- sum(dgamma(y, shape = exp(a), rate = exp(b), log = TRUE ) )
    ans <- LogLik + sum(dLogExpHalfCauchy(x = c(a,b), s.sq=5))
    ans
    
  }
  theta.hat <- optim(par=theta.guess, fn = Gamma.Objective, control=list(fnscale=-1))$par
  Proposal.Var <- -solve(pracma::hessian(Gamma.Objective,theta.hat))
  theta.proposed.MAT <- rmvnorm(B, mean= theta.hat, sigma=(Proposal.Var))
  theta.MAT <- matrix(0,nrow=B,ncol=2)
  theta.curr <- theta.hat
  logPostCurr <- Gamma.Objective(theta.curr)
  for(t in 1:B){
    
    theta.proposed <- theta.proposed.MAT[t,]
    logPostNew <- Gamma.Objective(theta.proposed)
    logQratio <- dmvnorm(theta.curr, mean= theta.hat, sigma=Proposal.Var, log = TRUE) - dmvnorm(theta.proposed, mean= theta.hat, sigma=Proposal.Var, log = TRUE)
    accept.prob <- exp(logPostNew - logPostCurr + logQratio)
    if(toss.uni.vec[t] <= accept.prob){
      
      theta.curr <- theta.proposed
      logPostCurr <- logPostNew
      accept.count <- accept.count + 1
      
    }
    theta.MAT[t,] <- theta.curr
    
  }
  return(list( a.vec = theta.MAT[,1], b.vec = theta.MAT[,2], accept.rate = (accept.count/B)))
  
  
}

#MCMC for misspecified gumbel copula. Used for Cut I posterior inference.
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#mu.val and v.val denote conditional draw from their respective Cut I marginal posteriors.
#a.init and b.init denote conditional draw from their respective Cut I marginal posteriors.
#omega init denotes initial value of logit(0.5*(tau+1)), where tau is the Kendall's tau copula parameter.
#B denotes the number of MCMC draws
NewMCMCcut1.LogNormal.Gamma.Gumbel <- function(y, mu.val, v.val, a.val, b.val, omega.init, B=10000){
  
  n <- length(y)
  toss.uni.vec <- runif(B)
  accept.count <- 0
  theta.guess <- c(omega.init)
  
  LogNormal.Gamma.Gumbel.Objective <- function(x){
    
    omega <- x
    tau <- exp(log(exp(exp(omega)) - 1) - log(exp(exp(omega)) + 1))
    mycop <- gumbelCopula(param = (1/(1-tau)) ,dim=d)
    LogLik <- sum(dCopula( cbind( plnorm(y[,1], meanlog = mu.val, sdlog = sqrt(exp(v.val)) ),  pgamma(y[,2], shape = exp(a.val), rate = exp(b.val) ) ), mycop, log = TRUE)  )
    ans <- LogLik + dLogOmegaPDF(omega)
    ans
    
  }
  theta.hat <- optim(par=theta.guess, fn = LogNormal.Gamma.Gumbel.Objective, control=list(fnscale=-1), method = "Brent", lower = -5, upper = 5)$par
  Proposal.Var <- -solve(pracma::hessian(LogNormal.Gamma.Gumbel.Objective,theta.hat))
  theta.proposed.vec <- rnorm(B, mean= theta.hat, sd=sqrt(c(Proposal.Var)) )
  theta.vec <- rep(0, B)
  theta.curr <- theta.hat
  logPostCurr <- LogNormal.Gamma.Gumbel.Objective(theta.curr)
  
  for(t in 1:B){
    
    theta.proposed <- theta.proposed.vec[t]
    logPostNew <- LogNormal.Gamma.Gumbel.Objective(theta.proposed)
    logQratio <- dnorm(theta.curr, mean= theta.hat, sd=sqrt(c(Proposal.Var)), log = TRUE) - dnorm(theta.proposed, mean= theta.hat, sd=sqrt(c(Proposal.Var)), log = TRUE)
    accept.prob <- exp(logPostNew - logPostCurr + logQratio)
    if(toss.uni.vec[t] <= accept.prob){
      
      theta.curr <- theta.proposed
      logPostCurr <- logPostNew
      accept.count <- accept.count + 1
      
    }
    theta.vec[t] <- theta.curr
    
  }
  
  return(list( tau.vec = exp(log(exp(exp(theta.vec)) - 1) - log(exp(exp(theta.vec)) + 1)),  accept.rate = (accept.count/B) ))
  
}

#MCMC for correctly specified T copula. Used for Correctly-specified Cut I posterior inference.
#y is a user-input matrix with two columns. First column for obserevd Y_1 values and second column for observed Y_2 values.
#mu.val and v.val denote conditional draw from their respective Cut I marginal posteriors.
#a.init and b.init denote conditional draw from their respective Cut I marginal posteriors.
#omega init denotes initial value of logit(0.5*(tau+1)), where tau is the Kendall's tau copula parameter.
#B denotes the number of MCMC draws
NewMCMCcut1.LogNormal.Gamma.T <- function(y, mu.val, v.val, a.val, b.val, omega.init, B=10000){
  
  n <- length(y)
  toss.uni.vec <- runif(B)
  accept.count <- 0
  theta.guess <- c(omega.init)
  
  LogNormal.Gamma.T.Objective <- function(x){
    
    omega <- x
    tau <- 2*expit(omega) - 1
    mycop <- tCopula(param = sin(0.5*pi*tau), dim=d, df=1 )
    LogLik <- sum(dCopula( cbind( plnorm(y[,1], meanlog = mu.val, sdlog = sqrt(exp(v.val)) ),  pgamma(y[,2], shape = exp(a.val), rate = exp(b.val) ) ), mycop, log = TRUE))
    ans <- LogLik + dLogOmegaCorrectPDF(omega)
    ans
    
  }
  theta.hat <- optim(par=theta.guess, fn = LogNormal.Gamma.T.Objective, control=list(fnscale=-1), method = "Brent", lower = -5, upper = 5)$par
  Proposal.Var <- -solve(pracma::hessian(LogNormal.Gamma.T.Objective,theta.hat))
  theta.proposed.vec <- rnorm(B, mean= theta.hat, sd=sqrt(c(Proposal.Var)) )
  theta.vec <- rep(0, B)
  theta.curr <- theta.hat
  logPostCurr <- LogNormal.Gamma.T.Objective(theta.curr)
  
  for(t in 1:B){
    
    theta.proposed <- theta.proposed.vec[t]
    logPostNew <- LogNormal.Gamma.T.Objective(theta.proposed)
    logQratio <- dnorm(theta.curr, mean= theta.hat, sd=sqrt(c(Proposal.Var)), log = TRUE) - dnorm(theta.proposed, mean= theta.hat, sd=sqrt(c(Proposal.Var)), log = TRUE)
    accept.prob <- exp(logPostNew - logPostCurr + logQratio)
    if(toss.uni.vec[t] <= accept.prob){
      
      theta.curr <- theta.proposed
      logPostCurr <- logPostNew
      accept.count <- accept.count + 1
      
    }
    theta.vec[t] <- theta.curr
    
  }
  
  return(list( tau.vec = (2*expit(theta.vec)-1),  accept.rate = (accept.count/B) ))
  
}

#Evaluates the log objective function at a given omega = logit(0.5*(tau+1)) value. Used for IFM copula parameter estimation.
#BigU is a user-input matrix with two columns. First column denotes the cumulative probabilities of observed Y_1 based on IFM estimates of marginal parameters. Second column denotes the cumulative probabilities of observed Y_2 based on IFM estimates of marginal parameters.
LogN.Gamma.Gumbel.Objective.ForIFM <- function(BigU, omega){
  
  tau <- exp(log(exp(exp(omega)) - 1) - log(exp(exp(omega)) + 1))
  mycop <- gumbelCopula(param = (1/(1-tau)) ,dim=d)
  LogLik <- sum(dCopula( BigU, mycop, log = TRUE)  )
  ans <- LogLik + dLogOmegaPDF(omega)
  ans
  
}