MakePhiArray <- function(x){
  
  xprime <- 2*pnorm(x[1:70]) - 1
  PhiArray.local <- array(0,dim=c(4,4,5))
  PhiArray.local[,,1] <- matrix(c(1,xprime[1],xprime[2],xprime[3],xprime[1],1,xprime[4],xprime[5],xprime[2],xprime[4],1,xprime[6],xprime[3],xprime[5],xprime[6],1),nrow=4,ncol=4)
  PhiArray.local[,,2] <- matrix(xprime[7:22],nrow=4,ncol=4)
  PhiArray.local[,,3] <- matrix(xprime[23:38],nrow=4,ncol=4)
  PhiArray.local[,,4] <- matrix(xprime[39:54],nrow=4,ncol=4)
  PhiArray.local[,,5] <- matrix(xprime[55:70],nrow=4,ncol=4)
  
  
  return(PhiArray.local)
}

USMacro.Cut.Marginal.ST.Single.Objective <- function(x, y){
  
  xi <- x[1]; omega <- exp(x[2]); alpha <- x[3]; nu <- exp(x[4])
  ans <- sum(dst(y, xi=xi, omega = omega, alpha = alpha, nu = nu, log=TRUE)) + dnorm(xi,mean=0,sd=100,log=TRUE) + sum(dnorm(x[2:4],mean=0,sd=2,log=TRUE))
  return(ans)
  
}
USMacro.Cut.Marginal.TruncST.Single.Objective <- function(x, y){
  
  xi <- x[1]; omega <- exp(x[2]); alpha <- x[3]; nu <- exp(x[4])
  normConst <- 1 - pst(0,dp=c(xi,omega,alpha,nu))
  ans <- sum(dst(y, xi=xi, omega = omega, alpha = alpha, nu = nu, log=TRUE)) + dnorm(xi,mean=0,sd=100,log=TRUE) + sum(dnorm(x[2:4],mean=0,sd=2,log=TRUE)) - length(y)*log(normConst)
  return(ans)
  
}

NewCopulaObjectiveFunctionUncut <- function(x, y1, y2, y3, y4){
  
  NumSamp <- length(y1)
  normConst3 <- 1 - pst(x = 0, xi = x[9], omega = exp(x[10]), alpha = x[11], nu = exp(x[12]) ); normConst4 <-  1 - pst(x = 0, xi = x[13], omega = exp(x[14]), alpha = x[15], nu = exp(x[16]) )
  U1 <- pst(x = y1, xi = x[1], omega = exp(x[2]), alpha = x[3], nu = exp(x[4]) ); U2 <- pst(x = y2, xi = x[5], omega = exp(x[6]), alpha = x[7], nu = exp(x[8]) ); U3 <- exp(log(pst(x = y3, xi = x[9], omega = exp(x[10]), alpha = x[11], nu = exp(x[12]) ) - pst(x = 0, xi = x[9], omega = exp(x[10]), alpha = x[11], nu = exp(x[12]) ) ) - log(normConst3) ); U4 <- exp(log(pst(x = y4, xi = x[13], omega = exp(x[14]), alpha = x[15], nu = exp(x[16]) ) - pst(x = 0, xi = x[13], omega = exp(x[14]), alpha = x[15], nu = exp(x[16]) ) ) - log(normConst4) )
  Umat.Local <- cbind(U1, U2, U3, U4)
  PhiArray.local <- MakePhiArray(x[17:86])
  Lik.val <- VineCopLikelihood(Umat = Umat.Local, PhiArray = PhiArray.local, p = 4) + sum(dst(x = y1, xi = x[1], omega = exp(x[2]), alpha = x[3], nu = exp(x[4]), log = TRUE)) + sum(dst(x = y2, xi = x[5], omega = exp(x[6]), alpha = x[7], nu = exp(x[8]), log = TRUE)) + sum(dst(x = y3, xi = x[9], omega = exp(x[10]), alpha = x[11], nu = exp(x[12]), log = TRUE)) - NumSamp*log(normConst3) + sum(dst(x = y4, xi = x[13], omega = exp(x[14]), alpha = x[15], nu = exp(x[16]), log = TRUE)) - NumSamp*log(normConst4)
  Prior.val <- sum(dnorm(x[17:22], sd=sqrt( exp(x[87]) ), log=TRUE)) + sum(dnorm(x[23:38], sd=sqrt( exp(x[88]) ), log=TRUE)) + sum(dnorm(x[39:54], sd=sqrt( exp(x[89]) ), log=TRUE)) + sum(dnorm(x[55:70], sd=sqrt( exp(x[90]) ), log=TRUE)) + sum(dnorm(x[71:86], sd=sqrt( exp(x[91]) ), log=TRUE))
  #browser()
  Hyperprior.val <- 5*log(2/pi) - 5*0.5*log(1) - sum(x[87:91]) -  sum(log(1 + exp(2*x[87:91])))
  return((Lik.val + Prior.val + Hyperprior.val))
}

NewCopulaObjectiveFunctionCut <- function(x, U1, U2, U3, U4){
  
  Umat.Local <- cbind(U1, U2, U3, U4)
  PhiArray.local <- MakePhiArray(x[1:70])
  Lik.val <- VineCopLikelihood(Umat = Umat.Local, PhiArray = PhiArray.local, p = 4)
  Prior.val <- sum(dnorm(x[1:6], sd=sqrt( exp(x[71]) ), log=TRUE)) + sum(dnorm(x[7:22], sd=sqrt( exp(x[72]) ), log=TRUE)) + sum(dnorm(x[23:38], sd=sqrt( exp(x[73]) ), log=TRUE)) + sum(dnorm(x[39:54], sd=sqrt( exp(x[74]) ), log=TRUE)) + sum(dnorm(x[55:70], sd=sqrt( exp(x[75]) ), log=TRUE))
  #browser()
  Hyperprior.val <- 5*log(2/pi) - 5*0.5*log(1) - sum(x[71:75]) -  sum(log(1 + exp(2*x[71:75])))
  return((Lik.val + Prior.val + Hyperprior.val))
}

NewCopulaObjectiveFunctionReverseCut <- function(x, U1, U2, U3, U4){
  
  Umat.Local <- cbind(U1, U2, U3, U4)
  PhiArray.local <- MakePhiArray(x[1:70])
  #browser()
  #browser()
  Lik.val <- VineCopLikelihood(Umat = Umat.Local, PhiArray = PhiArray.local, p = 4)
  Prior.val <- sum(dnorm(x[1:6], sd=sqrt( exp(x[71]) ), log=TRUE)) + sum(dnorm(x[7:22], sd=sqrt( exp(x[72]) ), log=TRUE)) + sum(dnorm(x[23:38], sd=sqrt( exp(x[73]) ), log=TRUE)) + sum(dnorm(x[39:54], sd=sqrt( exp(x[74]) ), log=TRUE)) + sum(dnorm(x[55:70], sd=sqrt( exp(x[75]) ), log=TRUE))
  #browser()
  Hyperprior.val <- 5*log(2/pi) - 5*0.5*log(1) - sum(x[71:75]) -  sum(log(1 + exp(2*x[71:75])))
  return((Lik.val + Prior.val + Hyperprior.val))
}

NewVB.cut.ST.USMacro <- function(y, mu.init, C.init, maxRuns=10000){
  
  n <- length(y); p <- length(mu.init)
  zMAT <- matrix(rnorm( (maxRuns*p) ), nrow=maxRuns, ncol=p)
  ELBOstore <- rep(0, maxRuns)
  gradstore <- matrix(0,nrow=maxRuns,ncol=p)
  mustore <- matrix(0,nrow=maxRuns,ncol=p)
  Cstore <- array(0,dim=c(p,p,maxRuns))
  
  Accugrad.mu <- rep(0,p)
  Change.mu <- rep(0,p)
  Accuchange.mu <- rep(0,p)
  Accugrad.C <- matrix(0,nrow = p, ncol=p)
  Accuchange.C <- matrix(0,nrow = p, ncol=p)
  
  LowerID <- matrix(1,nrow=p,ncol=p)
  LowerID[row(LowerID) < col(LowerID)] <- 0
  LowerID <- c(LowerID)
  
  WindowSize <- 5000
  
  redraw <- TRUE
  it <- 0
  maxLhat <- -9999999999999
  meanLhatcriterion <- 0
  rho <- 0.95
  epsilon <- 10^(-4)
  mu.curr <- mu.init
  C.curr <- C.init
  
  for(it in 1:maxRuns){
    
    zdraw <- zMAT[it,]
    theta.curr <- c(C.curr %*% zdraw) + mu.curr
    grad.curr <- numDeriv::grad(func=USMacro.Cut.Marginal.ST.Single.Objective, x=theta.curr, y=y)
    
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon) *grad.curr
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    
    Key.C <- (matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr)) )
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*Key.C^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * Key.C
    C.new <- C.curr + Change.C
    C.new[row(C.new) < col(C.new)] <- 0
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    Accuchange.C[row(Accuchange.C)< col(Accuchange.C)] <- 0
    
    mu.curr <- mu.new
    C.curr <- C.new
    gradstore[it,] <- grad.curr
    Cstore[,,it] <- C.curr
    mustore[it,] <- mu.curr
    if(it%%100==0){
      
      cat("Completed iter: ", it, "\n")
      
    }
    #browser()
    
    
  }
  
  return(list(mustore = mustore, Cstore = Cstore, gradstore=gradstore))
  
  
  
}

NewVB.cut.TruncST.USMacro <- function(y, mu.init, C.init, maxRuns=10000){
  
  n <- length(y); p <- length(mu.init)
  zMAT <- matrix(rnorm( (maxRuns*p) ), nrow=maxRuns, ncol=p)
  ELBOstore <- rep(0, maxRuns)
  gradstore <- matrix(0,nrow=maxRuns,ncol=p)
  mustore <- matrix(0,nrow=maxRuns,ncol=p)
  Cstore <- array(0,dim=c(p,p,maxRuns))
  
  Accugrad.mu <- rep(0,p)
  Change.mu <- rep(0,p)
  Accuchange.mu <- rep(0,p)
  Accugrad.C <- matrix(0,nrow = p, ncol=p)
  Accuchange.C <- matrix(0,nrow = p, ncol=p)
  
  LowerID <- matrix(1,nrow=p,ncol=p)
  LowerID[row(LowerID) < col(LowerID)] <- 0
  LowerID <- c(LowerID)
  
  WindowSize <- 5000
  
  redraw <- TRUE
  it <- 0
  maxLhat <- -9999999999999
  meanLhatcriterion <- 0
  rho <- 0.95
  epsilon <- 10^(-4)
  mu.curr <- mu.init
  C.curr <- C.init
  
  for(it in 1:maxRuns){
    
    zdraw <- zMAT[it,]
    theta.curr <- c(C.curr %*% zdraw) + mu.curr
    grad.curr <- numDeriv::grad(func=USMacro.Cut.Marginal.TruncST.Single.Objective, x=theta.curr, y=y)
    
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon) *grad.curr
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    
    Key.C <- (matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr)) )
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*Key.C^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * Key.C
    C.new <- C.curr + Change.C
    C.new[row(C.new) < col(C.new)] <- 0
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    Accuchange.C[row(Accuchange.C)< col(Accuchange.C)] <- 0
    
    mu.curr <- mu.new
    C.curr <- C.new
    gradstore[it,] <- grad.curr
    Cstore[,,it] <- C.curr
    mustore[it,] <- mu.curr
    if(it%%100==0){
      
      cat("Completed iter: ", it, "\n")
      
    }
    #browser()
    
    
  }
  
  return(list(mustore = mustore, Cstore = Cstore, gradstore=gradstore))
  
  
  
}

NewCopulaObjectiveFunction <- function(x, Umat.Local){
  
  PhiArray.local <- MakePhiArray(x[1:70])
  Lik.val <- VineCopLikelihood(Umat = Umat.Local, PhiArray = PhiArray.local, p = 4)
  Prior.val <- sum(dnorm(x[1:6], sd=sqrt( exp(x[71]) ), log=TRUE)) + sum(dnorm(x[7:22], sd=sqrt( exp(x[72]) ), log=TRUE)) + sum(dnorm(x[23:38], sd=sqrt( exp(x[73]) ), log=TRUE)) + sum(dnorm(x[39:54], sd=sqrt( exp(x[74]) ), log=TRUE)) + sum(dnorm(x[55:70], sd=sqrt( exp(x[75]) ), log=TRUE))
  #browser()
  Hyperprior.val <- 5*log(2/pi) - 5*0.5*log(1) - sum(x[71:75]) -  sum(log(1 + exp(2*x[71:75])))
  return((Lik.val + Prior.val + Hyperprior.val))
}

NewCopulaVB <- function(UMAT.input, maxRuns=500000, mu.init, C.init){
  
  p <- length(mu.init)
  CArray <- array(0,dim=c(p,p,maxRuns))
  muMAT <- matrix(0,nrow=maxRuns,ncol=p)
  thetaMAT <- matrix(0,nrow=maxRuns,ncol=p)
  gradMAT <- matrix(0,nrow=maxRuns,ncol=p)
  zMAT <- rmvnorm(maxRuns,mean = rep(0,p), sigma = diag(1,p))
  
  mu.curr <- mu.init
  C.curr <- C.init
  rho <- 0.85
  epsilon <- 10^(-6)
  
  Accugrad.mu <- rep(0,p)
  Accuchange.mu <- rep(0,p)
  Accugrad.C <- matrix(0,nrow=p,ncol=p)
  Accuchange.C <- matrix(0,nrow=p,ncol=p)
  
  
  
  for(t in 1:maxRuns){
    
    zdraw <- zMAT[t,]
    theta.curr <- c(C.curr %*% matrix(zdraw,ncol=1)) + mu.curr
    
    grad.curr <- pracma::grad(f=NewCopulaObjectiveFunction, x0 = theta.curr, Umat.Local=UMAT.input)
    #browser()
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon) *grad.curr
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    
    KeyC <- (matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr)) )
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*KeyC^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * KeyC
    C.new <- C.curr + Change.C
    C.new[row(C.new) < col(C.new)] <- 0
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    Accuchange.C[row(Accuchange.C)< col(Accuchange.C)] <- 0
    
    CArray[,,t] <- C.new
    muMAT[t,] <- mu.new
    gradMAT[t,] <- c(grad.curr)
    thetaMAT[t,] <- c(theta.curr)
    zMAT[t,] <- zdraw
    mu.curr <- mu.new
    C.curr <- C.new
    
    
    
    if( (t%%10)==0 ){
      
      cat("Completed iter", t, " \n")
      
    }
    
    
  }
  return(list(muMAT=muMAT,CArray=CArray))
  
  
}

NewCopulaVBUncut <- function(y1, y2, y3, y4, maxRuns=500000, mu.init, C.init){
  
  p <- length(mu.init)
  CArray <- array(0,dim=c(p,p,maxRuns))
  muMAT <- matrix(0,nrow=maxRuns,ncol=p)
  thetaMAT <- matrix(0,nrow=maxRuns,ncol=p)
  gradMAT <- matrix(0,nrow=maxRuns,ncol=p)
  zMAT <- rmvnorm(maxRuns,mean = rep(0,p), sigma = diag(1,p))
  
  mu.curr <- mu.init
  C.curr <- C.init
  rho <- 0.85
  epsilon <- 10^(-6)
  
  Accugrad.mu <- rep(0,p)
  Accuchange.mu <- rep(0,p)
  Accugrad.C <- matrix(0,nrow=p,ncol=p)
  Accuchange.C <- matrix(0,nrow=p,ncol=p)
  
  
  
  for(t in 1:maxRuns){
    
    zdraw <- zMAT[t,]
    theta.curr <- c(C.curr %*% matrix(zdraw,ncol=1)) + mu.curr
    
    grad.curr <- pracma::grad(f=NewCopulaObjectiveFunctionUncut , x0 = theta.curr, y1=y1,y2=y2,y3=y3,y4=y4)
    #browser()
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon) *grad.curr
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    
    KeyC <- (matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr)) )
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*KeyC^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * KeyC
    C.new <- C.curr + Change.C
    C.new[row(C.new) < col(C.new)] <- 0
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    Accuchange.C[row(Accuchange.C)< col(Accuchange.C)] <- 0
    
    CArray[,,t] <- C.new
    muMAT[t,] <- mu.new
    gradMAT[t,] <- c(grad.curr)
    thetaMAT[t,] <- c(theta.curr)
    zMAT[t,] <- zdraw
    mu.curr <- mu.new
    C.curr <- C.new
    
    
    
    if( (t%%10)==0 ){
      
      cat("Completed iter", t, " \n")
      
    }
    
    
  }
  return(list(muMAT=muMAT,CArray=CArray))
  
  
}

NewCopulaVBCut <- function(y1, y2, y3, y4, maxRuns=500000, mu.init, C.init){
  
  p <- length(mu.init)
  CArray <- array(0,dim=c(p,p,maxRuns))
  muMAT <- matrix(0,nrow=maxRuns,ncol=p)
  thetaMAT <- matrix(0,nrow=maxRuns,ncol=p)
  gradMAT <- matrix(0,nrow=maxRuns,ncol=p)
  zMAT <- rmvnorm(maxRuns,mean = rep(0,p), sigma = diag(1,p))
  
  mu.curr <- mu.init
  C.curr <- C.init
  rho <- 0.85
  epsilon <- 10^(-6)
  
  Accugrad.mu <- rep(0,p)
  Accuchange.mu <- rep(0,p)
  Accugrad.C <- matrix(0,nrow=p,ncol=p)
  Accuchange.C <- matrix(0,nrow=p,ncol=p)
  
  
  
  for(t in 1:maxRuns){
    
    zdraw <- zMAT[t,]
    theta.curr <- c(C.curr %*% matrix(zdraw,ncol=1)) + mu.curr
    normConst3.temp <- 1 - pst(x = 0, xi = theta.curr[9], omega = exp(theta.curr[10]), alpha = theta.curr[11], nu = exp(theta.curr[12]) ); normConst4.temp <-  1 - pst(x = 0, xi = theta.curr[13], omega = exp(theta.curr[14]), alpha = theta.curr[15], nu = exp(theta.curr[16]) )
    U1.temp <- pst(x = y1, xi = theta.curr[1], omega = exp(theta.curr[2]), alpha = theta.curr[3], nu = exp(theta.curr[4]) ); U2.temp <- pst(x = y2, xi = theta.curr[5], omega = exp(theta.curr[6]), alpha = theta.curr[7], nu = exp(theta.curr[8]) ); U3.temp <- exp(log(pst(x = y3, xi = theta.curr[9], omega = exp(theta.curr[10]), alpha = theta.curr[11], nu = exp(theta.curr[12]) ) - pst(x = 0, xi = theta.curr[9], omega = exp(theta.curr[10]), alpha = theta.curr[11], nu = exp(theta.curr[12]) ) ) - log(normConst3.temp) ); U4.temp <- exp(log(pst(x = y4, xi = theta.curr[13], omega = exp(theta.curr[14]), alpha = theta.curr[15], nu = exp(theta.curr[16]) ) - pst(x = 0, xi = theta.curr[13], omega = exp(theta.curr[14]), alpha = theta.curr[15], nu = exp(theta.curr[16]) ) ) - log(normConst4.temp) )
    grad.curr <- pracma::grad(f=NewCopulaObjectiveFunctionCut , x0 = theta.curr[17:91], U1=U1.temp, U2=U2.temp, U3=U3.temp, U4=U4.temp)
    grad.curr <- c(rep(0,16),grad.curr)
    #browser()
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon) *grad.curr
    Change.mu[1:16] <- rep(0,16)
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    
    KeyC <- (matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr)) )
    KeyC[row(KeyC) < col(KeyC)] <- 0
    KeyC[1:16,1:16] <- 0
    Accugrad.C <- rho*Accugrad.C  + (1-rho)*KeyC^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * KeyC
    C.new <- C.curr + Change.C
    Change.C[row(C.new) < col(C.new)] <- 0
    Change.C[1:16,1:16] <- 0
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    Accuchange.C[row(Accuchange.C)< col(Accuchange.C)] <- 0
    
    CArray[,,t] <- C.new
    muMAT[t,] <- mu.new
    gradMAT[t,] <- c(grad.curr)
    thetaMAT[t,] <- c(theta.curr)
    zMAT[t,] <- zdraw
    mu.curr <- mu.new
    C.curr <- C.new
    
    
    
    if( (t%%10)==0 ){
      
      cat("Completed iter", t, " \n")
      
    }
    
    
  }
  return(list(muMAT=muMAT,CArray=CArray))
  
  
}


LogScoreMetric.FinalPart <- function(y1, y2, y3, y4, VBmean.sub, VBcov.sub, numDraws = 5000){
  
  Para.draws <- rmvnorm(numDraws, mean = VBmean.sub, sigma = VBcov.sub)
  LogScore.draws <- matrix(0,nrow=numDraws,ncol=(length(y1) - 4))
  MarginalVec <- rep(0,numDraws)
  NumSamp <- length(y1)
  for(draw.id in 1:numDraws){
    
    PhiArray.temp <- MakePhiArray(Para.draws[draw.id,17:86])
    normConst3.temp <- 1 - pst(x = 0, xi = Para.draws[draw.id,9], omega = exp(Para.draws[draw.id,10]), alpha = Para.draws[draw.id,11], nu = exp(Para.draws[draw.id,12]) ); normConst4.temp <-  1 - pst(x = 0, xi = Para.draws[draw.id,13], omega = exp(Para.draws[draw.id,14]), alpha = Para.draws[draw.id,15], nu = exp(Para.draws[draw.id,16]) )
    U1.temp <- pst(x = y1, xi = Para.draws[draw.id,1], omega = exp(Para.draws[draw.id,2]), alpha = Para.draws[draw.id,3], nu = exp(Para.draws[draw.id,4]) ); U2.temp <- pst(x = y2, xi = Para.draws[draw.id,5], omega = exp(Para.draws[draw.id,6]), alpha = Para.draws[draw.id,7], nu = exp(Para.draws[draw.id,8]) ); U3.temp <- exp(log(pst(x = y3, xi = Para.draws[draw.id,9], omega = exp(Para.draws[draw.id,10]), alpha = Para.draws[draw.id,11], nu = exp(Para.draws[draw.id,12]) ) - pst(x = 0, xi = Para.draws[draw.id,9], omega = exp(Para.draws[draw.id,10]), alpha = Para.draws[draw.id,11], nu = exp(Para.draws[draw.id,12]) ) ) - log(normConst3.temp) ); U4.temp <- exp(log(pst(x = y4, xi = Para.draws[draw.id,13], omega = exp(Para.draws[draw.id,14]), alpha = Para.draws[draw.id,15], nu = exp(Para.draws[draw.id,16]) ) - pst(x = 0, xi = Para.draws[draw.id,13], omega = exp(Para.draws[draw.id,14]), alpha = Para.draws[draw.id,15], nu = exp(Para.draws[draw.id,16]) ) ) - log(normConst4.temp) )
    Umat.temp <- cbind(U1.temp, U2.temp, U3.temp, U4.temp)
    LogScore.draws[draw.id,] <- LogScoreMetricUPart(Umat = Umat.temp, PhiArray = PhiArray.temp, p = 4)
    MarginalVec[draw.id] <- sum(dst(x = y1, xi = Para.draws[draw.id,1], omega = exp(Para.draws[draw.id,2]), alpha = Para.draws[draw.id,3], nu = exp(Para.draws[draw.id,4]), log = TRUE)) + sum(dst(x = y2, xi = Para.draws[draw.id,5], omega = exp(Para.draws[draw.id,6]), alpha = Para.draws[draw.id,7], nu = exp(Para.draws[draw.id,8]), log = TRUE)) + sum(dst(x = y3, xi = Para.draws[draw.id,9], omega = exp(Para.draws[draw.id,10]), alpha = Para.draws[draw.id,11], nu = exp(Para.draws[draw.id,12]), log = TRUE)) - NumSamp*log(normConst3.temp) + sum(dst(x = y4, xi = Para.draws[draw.id,13], omega = exp(Para.draws[draw.id,14]), alpha = Para.draws[draw.id,15], nu = exp(Para.draws[draw.id,16]), log = TRUE)) - NumSamp*log(normConst4.temp)
    
  }
  
  
  
  return(list(LogScore.draws=LogScore.draws, MarginalVec=MarginalVec))
  
  
}

pst.modified <- function (x, xi = 0, omega = 1, alpha = 0, nu = Inf){
  

  dp.std <- c(0, 1, alpha, nu)
  dp.std <- cbind(0,1,alpha, nu)
  delta <- alpha/sqrt(1 + alpha^2)
  browser()
  pr <- rep(NA, length(x))
  z <- ((x - xi)/omega)
  nu0 <- (8.2 + 3.55 * log(log(length(z) + 1)))
  browser()
  zero.alpha.ind <- which(alpha==0)
  p <- rep(0,length(xi))
  p[zero.alpha.ind] <- pt(z[zero.alpha.ind], df = nu[zero.alpha.ind])
  fp <- function(v, alpha, nu, t.value) psn(sqrt(v) * t.value, 0, 1, alpha) * dchisq(v * nu, nu) * nu
  browser()
  p <- pst_int(z, 0, 1, alpha, nu)


}

DrawQu <- function(Amat, Bmat, etamat, wmat){
  
  TSlength <- nrow(Amat)
  zDrawMat <- matrix(rnorm(TSlength*4),nrow=TSlength, ncol=4)*wmat + etamat
  pnormzDrawMat <- pmax(pmin(pnorm(zDrawMat),matrix(0.99,nrow=TSlength,ncol=4)),matrix(0.01,nrow=TSlength,ncol=4))
  uDrawMat <- pnormzDrawMat*(Bmat - Amat) + Amat
  #browser()
  return(uDrawMat)
  
}

DrawQPhiLogTau <- function(mu, C){
  
  draw.len <- length(mu)
  return(mu+c(C%*%rnorm(draw.len)))
  
}

EstimateXi <- function(AMAT, BMAT, mu, C, etamat, wmat, numDraws){
  
  d <- length(mu)
  numPars <- d + 0.5*d*(d+1) + 2*length(c(etamat))
  E1.est <- rep(0,numPars)
  E2.est <- rep(0,numPars)
  E3.est <- rep(0,numPars)
  E4.est <- rep(0,numPars)
  
  for(b in 1:numDraws){
    
    DrawU.curr <- DrawQu(Amat = AMAT, Bmat = BMAT, etamat = etamat, wmat = wmat)
    DrawTauPhi.curr <- DrawQPhiLogTau(mu=mu, C=C)
    Logh.curr <- NewCopulaObjectiveFunctionReverseCut(x = DrawTauPhi.curr, U1 = DrawU.curr[,1], U2 = DrawU.curr[,2], U3 = DrawU.curr[,3], U4 = DrawU.curr[,4])
    Logq.curr <- LogQEval(UMAT = DrawU.curr, Amat = AMAT, Bmat = BMAT, vcopulaPar = DrawTauPhi.curr, vmu = mu, mC = C, veta = etamat, vwSq = wmat)
    GradLogQ.curr <- c(GradLogQ(UMAT = DrawU.curr, Amat = AMAT, Bmat = BMAT, vcopulaPar = DrawTauPhi.curr, vmu = mu, mC = C, veta = etamat, vwSq = wmat))
    E1.est <- E1.est + (Logh.curr - Logq.curr)*(GradLogQ.curr^2)
    E2.est <- E2.est + (Logh.curr - Logq.curr)*GradLogQ.curr
    E3.est <- E3.est + GradLogQ.curr
    E4.est <- E4.est + (GradLogQ.curr^2)

  }
  denom <- ( (E4.est/numDraws) - (E3.est/numDraws)^2  )
  pmax.denom <- pmax(abs(denom),rep(0.0000001,length(denom)))
  sign.denom <- sign(denom)
  sign.denom[sign.denom==0] <- 1
  denom <- sign.denom*pmax.denom
  vxi <- ((E1.est/numDraws) - (E2.est/numDraws)*(E3.est/numDraws))/denom
  
  
  return(vxi)

}

NewCopulaVBReverseCutRversion <- function(AMAT, BMAT, muInit, CInit, EtaMatInit, WMatInit, maxRuns){
  
  d <- length(muInit)
  S <- 100
  LengthTS = nrow(EtaMatInit)
  numPar <- d + d*(d+1)/2 + 8*LengthTS
  
  muCurr = muInit
  CCurr = CInit
  rho = 0.7
  eps = 1/1000000
  
  AccugradMu <- rep(0,d)
  AccuchangeMu <- rep(0,d)
  AccugradC <- matrix(0,nrow=d,ncol=d)
  AccuchangeC <- matrix(0,nrow=d,ncol=d)
  AccugradEta <- matrix(0,nrow=LengthTS,ncol=4)
  AccuchangeEta <- matrix(0,nrow=LengthTS,ncol=4)
  AccugradW <- matrix(0,nrow=LengthTS,ncol=4)
  AccuchangeW <- matrix(0,nrow=LengthTS,ncol=4)
  
  MuStore <- matrix(0,nrow=maxRuns, ncol=d)
  CStore <- array(0,dim = c(d,d,maxRuns))
  EtaStore <- array(0,dim = c(LengthTS,4,maxRuns))
  WStore <- array(0,dim = c(LengthTS,4,maxRuns))
  
  EtaMatCurr = EtaMatInit
  WMatCurr = WMatInit
  GradLogQVec <- rep(0,numPar)
  gEstCurr <- rep(0,numPar)
  
  UdrawCurr <- matrix(0,nrow=nrow(AMAT),ncol=ncol(AMAT))
  DrawTauPhi <- rep(0,d)
  
  Cz <- matrix(0,nrow=d,ncol=S)
  zDrawMat <- matrix(0,nrow=LengthTS,4)
  anstemp <- NULL
  logHcurr <- NULL
  logQcurr <- NULL
  
  gradMu <- rep(0,d)
  ChangeC <- matrix(0,nrow=d,ncol=d)
  CNew <- matrix(0,nrow=d,ncol=d)
  
  gradEta <- matrix(0,nrow=LengthTS,ncol=4)
  ChangeEta <- matrix(0,nrow=LengthTS,ncol=4)
  WNew <- matrix(0,nrow=LengthTS,ncol=4)
  
  MAT99 <- matrix(0.99,nrow=LengthTS,ncol=4)
  MAT01 <- matrix(0.01,nrow=LengthTS,ncol=4)
  
  for(t in 1:maxRuns){
    
    #XiEstCurr = EstimateXiCpp(AMAT, BMAT, EtaMatCurr, WMatCurr, muCurr, CCurr, 100)
    Cz <- matrix(rnorm(d*S),nrow=d,ncol=S)
    zDrawMat <- matrix(rnorm(LengthTS*4),nrow=LengthTS,ncol=4)
    gEstCurr <- rep(0,numPar)
    for(s in 1:S){
      
      DrawTauPhiCurr = muCurr + Cz[,s]
      anstemp <- pnorm(zDrawMat*exp(0.5*WMatCurr) + EtaMatCurr)
      anstemp <- pmax(pmin(anstemp,MAT99),MAT01)
      UdrawCurr <- anstemp*(BMAT - AMAT) + AMAT
      logHcurr = NewCopulaObjectiveFunctionReverseCutCpp(DrawTauPhiCurr, UdrawCurr[,1], UdrawCurr[,2], UdrawCurr[,3], UdrawCurr[,4])
      logQcurr <- LogQEval(UdrawCurr, AMAT, BMAT, DrawTauPhiCurr, muCurr, CCurr, EtaMatCurr, WMatCurr)
      GradLogQVec <- c(GradLogQ(UdrawCurr, AMAT, BMAT, DrawTauPhiCurr, muCurr, CCurr, EtaMatCurr, WMatCurr))
      #gEstCurr <- gEstCurr + (logHcurr - logQcurr  - XiEstCurr)*GradLogQVec
      gEstCurr <- gEstCurr + (logHcurr - logQcurr)*GradLogQVec
      
    }
    gEstCurr <- gEstCurr/S
    #browser()
    gradMu <- gEstCurr[1:d]
    AccugradMu <- rho*AccugradMu + (1-rho)*(gradMu^2)
    ChangeMu <- sqrt( AccuchangeMu + eps )/sqrt( AccugradMu + eps )*gradMu
    muNew = muCurr + ChangeMu
    AccuchangeMu = rho*AccuchangeMu + (1-rho)*(ChangeMu^2)
    #browser()
    gradC = invvechCpp( gEstCurr[(d+1):(  d+d*(d+1)/2 )] )
    AccugradC = rho*AccugradC + (1-rho)*(gradC^2)
    ChangeC = ( sqrt( AccuchangeC + eps ) )/sqrt( AccugradC + eps )*gradC
    CNew = CCurr + ChangeC
    AccuchangeC <- rho*AccuchangeC + (1-rho)*(ChangeC^2)
    
    gradEta <- matrix( gEstCurr[(d+d*(d+1)/2 + 1):(d+d*(d+1)/2 + 4*LengthTS)], nrow=LengthTS,ncol=4,byrow=TRUE)
    AccugradEta <- rho*AccugradEta + (1-rho)*(gradEta^2)
    ChangeEta <- sqrt( AccuchangeEta + eps )/sqrt( AccugradEta + eps)*gradEta
    EtaNew <- EtaMatCurr + ChangeEta
    AccuchangeEta <- rho*AccuchangeEta + (1-rho)*(ChangeEta^2)
    
    gradW <- matrix( gEstCurr[(d+d*(d+1)/2 + 4*LengthTS+1):(d+d*(d+1)/2 + 8*LengthTS)] ,nrow=LengthTS,ncol=4,byrow=TRUE)
    AccugradW <- rho*AccugradW + (1-rho)*(gradW^2)
    ChangeW <- sqrt(AccuchangeW + eps)/sqrt(AccugradW + eps)*gradW
    WNew <- WMatCurr + ChangeW
    AccuchangeW <- rho*AccuchangeW + (1-rho)*(ChangeW^2)
    
    muCurr <- muNew
    CCurr <- CNew
    EtaMatCurr <- EtaNew
    WMatCurr <- WNew
    MuStore[t,] <- muCurr
    CStore[,,t] <- CCurr
    EtaStore[,,t] <- EtaMatCurr
    WStore[,,t] <- WMatCurr
    if(is.nan(sum(muCurr))){
      
      browser()
      
    }
    if((t%%10)==0){
      
      cat("Completed iter ", "\n")
      
    }

  }
  
  return( list(MuStore=MuStore,CStore=CStore,EtaStore=EtaStore,WStore=WStore) )
  
}

NewCopulaVBReverseCutMarginalPar <- function(y1, y2, y3, y4, maxRuns=500000, mu.init, C.init){
  
  p <- length(mu.init)
  CArray <- array(0,dim=c(p,p,maxRuns))
  muMAT <- matrix(0,nrow=maxRuns,ncol=p)
  thetaMAT <- matrix(0,nrow=maxRuns,ncol=p)
  gradMAT <- matrix(0,nrow=maxRuns,ncol=p)
  zMAT <- rmvnorm(maxRuns,mean = rep(0,p), sigma = diag(1,p))
  
  mu.curr <- mu.init
  C.curr <- C.init
  rho <- 0.85
  epsilon <- 10^(-6)
  
  Accugrad.mu <- rep(0,p)
  Accuchange.mu <- rep(0,p)
  Accugrad.C <- matrix(0,nrow=p,ncol=p)
  Accuchange.C <- matrix(0,nrow=p,ncol=p)
  
  
  
  for(t in 1:maxRuns){
    
    zdraw <- zMAT[t,]
    theta.curr <- c(C.curr %*% matrix(zdraw,ncol=1)) + mu.curr
    theta.curr <- c(theta.curr[76:91],theta.curr[1:75])
    normConst3.temp <- 1 - pst(x = 0, xi = theta.curr[9], omega = exp(theta.curr[10]), alpha = theta.curr[11], nu = exp(theta.curr[12]) ); normConst4.temp <-  1 - pst(x = 0, xi = theta.curr[13], omega = exp(theta.curr[14]), alpha = theta.curr[15], nu = exp(theta.curr[16]) )
    U1.temp <- pst(x = y1, xi = theta.curr[1], omega = exp(theta.curr[2]), alpha = theta.curr[3], nu = exp(theta.curr[4]) ); U2.temp <- pst(x = y2, xi = theta.curr[5], omega = exp(theta.curr[6]), alpha = theta.curr[7], nu = exp(theta.curr[8]) ); U3.temp <- exp(log(pst(x = y3, xi = theta.curr[9], omega = exp(theta.curr[10]), alpha = theta.curr[11], nu = exp(theta.curr[12]) ) - pst(x = 0, xi = theta.curr[9], omega = exp(theta.curr[10]), alpha = theta.curr[11], nu = exp(theta.curr[12]) ) ) - log(normConst3.temp) ); U4.temp <- exp(log(pst(x = y4, xi = theta.curr[13], omega = exp(theta.curr[14]), alpha = theta.curr[15], nu = exp(theta.curr[16]) ) - pst(x = 0, xi = theta.curr[13], omega = exp(theta.curr[14]), alpha = theta.curr[15], nu = exp(theta.curr[16]) ) ) - log(normConst4.temp) )
    #browser()
    grad.curr <- pracma::grad(f=NewCopulaObjectiveFunctionUncut , x0 = theta.curr, y1=y1, y2=y2, y3=y3, y4=y4)
    grad.curr[17:91] <- 0
    grad.curr <- c(grad.curr[17:91],grad.curr[1:16])
    #browser()
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon) *grad.curr
    #browser()
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    
    KeyC <- (matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr)) )
    KeyC[row(KeyC) < col(KeyC)] <- 0
    KeyC[1:75,1:75] <- 0
    Accugrad.C <- rho*Accugrad.C  + (1-rho)*KeyC^2
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * KeyC
    Change.C[row(Change.C ) < col(Change.C)] <- 0
    Change.C[1:75,1:75] <- 0
    C.new <- C.curr + Change.C
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    Accuchange.C[row(Accuchange.C)< col(Accuchange.C)] <- 0
    
    CArray[,,t] <- C.new
    muMAT[t,] <- mu.new
    gradMAT[t,] <- c(grad.curr)
    thetaMAT[t,] <- c(theta.curr)
    zMAT[t,] <- zdraw
    mu.curr <- mu.new
    C.curr <- C.new
    
    
    
    if( (t%%10)==0 ){
      
      cat("Completed iter", t, " \n")
      
    }
    
    
  }
  return(list(muMAT=muMAT,CArray=CArray))
  
  
}

# NewCopulaVBReverseCut <- function(AMAT, BMAT, mu.init, C.init, etamat.init, wmat.init, maxRuns = 50000){
#   
#   d <- length(mu.init)
#   S <- 100
#   LengthTS <- nrow(etamat.init)
#   CArray <- array(0,dim=c(d,d,maxRuns))
#   muMAT <- matrix(0,nrow=maxRuns,ncol=d)
#   etaArray <- array(0,dim=c(LengthTS,4,maxRuns))
#   WArray <- array(0,dim=c(LengthTS,4,maxRuns))
#   
#   mu.curr <- mu.init
#   C.curr <- C.init
#   rho <- 0.85
#   epsilon <- 10^(-6)
#   
#   Accugrad.mu <- rep(0,d)
#   Accuchange.mu <- rep(0,d)
#   Accugrad.C <- matrix(0,nrow=d,ncol=d)
#   Accuchange.C <- matrix(0,nrow=d,ncol=d)
#   Accugrad.eta <- matrix(0,nrow=LengthTS,ncol=4)
#   Accuchange.eta <- matrix(0,nrow=LengthTS,ncol=4)
#   Accugrad.W <- matrix(0,nrow=LengthTS,ncol=4)
#   Accuchange.W <- matrix(0,nrow=LengthTS,ncol=4)
#   
#   etamat.curr <- etamat.init
#   wmat.curr <- wmat.init
#   mu.curr <- mu.init
#   C.curr <- C.init
#   
#   for(t in 1:maxRuns){
#     
#     XiEst.curr <- EstimateXi(AMAT = AMAT, BMAT = BMAT, mu = mu.curr, C = C.curr, etamat = etamat.curr, wmat = wmat.curr, numDraws = 100)
#     g.est <- rep(0,(d+0.5*d*(d+1)+8*LengthTS))
#     #browser()
#     for(s in 1:S){
#       
#       DrawU.curr <- DrawQu(Amat = AMAT, Bmat = BMAT, etamat = etamat.curr, wmat = wmat.curr)
#       DrawTauPhi.curr <- DrawQPhiLogTau(mu=mu.curr, C=C.curr)
#       Logh.curr <- NewCopulaObjectiveFunctionReverseCut(x = DrawTauPhi.curr, U1 = DrawU.curr[,1], U2 = DrawU.curr[,2], U3 = DrawU.curr[,3], U4 = DrawU.curr[,4])
#       Logq.curr <- LogQEval(UMAT = DrawU.curr, Amat = AMAT, Bmat = BMAT, vcopulaPar = DrawTauPhi.curr, vmu = mu.curr, mC = C.curr, veta = etamat.curr, vwSq = wmat.curr)
#       GradLogQ.curr <- c(GradLogQ(UMAT = DrawU.curr, Amat = AMAT, Bmat = BMAT, vcopulaPar = DrawTauPhi.curr, vmu = mu.curr, mC = C.curr, veta = etamat.curr, vwSq = wmat.curr))
#       g.est <- g.est + (Logh.curr - Logq.curr - XiEst.curr)*GradLogQ.curr
#       
#     }
#     g.est <- g.est/S
#     #browser()
#     grad.mu <- g.est[1:d]
#     Accugrad.mu <- rho*Accugrad.mu + (1-rho)*(grad.mu^2)
#     Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon)*grad.mu
#     mu.new <- mu.curr + Change.mu
#     Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
#     #browser()
#     grad.C <- invvech(g.est[(d+1):(d+0.5*d*(d+1))])
#     grad.C[row(grad.C)<col(grad.C)] <- 0
#     Accugrad.C <- rho*Accugrad.C + (1-rho)*(grad.C^2)
#     Change.C <- sqrt(Accuchange.C +  epsilon)/sqrt(Accugrad.C + epsilon)*grad.C
#     C.new <- C.curr + Change.C
#     Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
#     #browser()
#     grad.eta <- g.est[(d+0.5*d*(d+1)+1):(d+0.5*d*(d+1)+4*LengthTS)]
#     grad.eta <- matrix(grad.eta,nrow=LengthTS,ncol=4,byrow=TRUE)
#     Accugrad.eta <- rho*Accugrad.eta + (1-rho)*(grad.eta^2)
#     Change.eta <- sqrt(Accuchange.eta +  epsilon)/sqrt(Accugrad.eta + epsilon)*grad.eta
#     etamat.new <- etamat.curr + Change.eta
#     Accuchange.eta <- rho*Accuchange.eta + (1-rho)*Change.eta^2
#     #browser()
#     grad.W <- g.est[(d+0.5*d*(d+1)+4*LengthTS+1):(d+0.5*d*(d+1)+8*LengthTS)]
#     grad.W <- matrix(grad.W,nrow=LengthTS,ncol=4,byrow=TRUE)
#     Accugrad.W <- rho*Accugrad.W + (1-rho)*(grad.W^2)
#     Change.W <- sqrt(Accuchange.W +  epsilon)/sqrt(Accugrad.W + epsilon)*grad.W
#     wmat.new <- wmat.curr + Change.W
#     Accuchange.W <- rho*Accuchange.W + (1-rho)*Change.W^2
#     #browser()
#     mu.curr <- mu.new
#     C.curr <- C.new
#     etamat.curr <- etamat.new
#     wmat.curr <- wmat.new
#     muMAT[t,] <- mu.curr
#     CArray[,,t] <- C.curr
#     etaArray[,,t] <- etamat.curr
#     WArray[,,t] <- wmat.curr
#     
#     if( (t%%10)==0 ){
#       
#       cat("Completed iter", t, " \n")
#       
#     }
#     
#   }
#   
#   return(list(muMAT=muMAT,CArray=CArray,etaArray=etaArray,WArray=WArray))
#   
#   
# }

# LogQEval.R <- function( UMAT, Amat, Bmat, vcopulaPar, vmu, mCvech, veta, vwSq ){
#   
#   mC.input <- invvech(mCvech)
#   mC.input[row(mC.input)<col(mC.input)] <- 0
#   LogQEval(UMAT = UMAT, Amat = Amat, Bmat = Bmat, vcopulaPar = vcopulaPar, vmu = vmu, mC = mC.input, veta = veta, vwSq = vwSq)
#   
# }