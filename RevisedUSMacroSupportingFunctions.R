dLogExpHalfCauchy<- function(x, s.sq){
  
  return(-log(pi*sqrt(s.sq)/2) - log(exp(-x) + exp(x)/s.sq))
  
}

dLogExpExp<- function(x, lambda.rate){
  
  return(log(lambda.rate) + x - lambda.rate*exp(x))
  
}

UpperTailEst.NormProb <- function(z){
  
  return(-0.5*log(2*pi) + log( sqrt(z^2 + 2*log(2)) - z ) - z^2/2)
  
}

ApproxQnorm <- function(p){
  
  return(-sqrt( -2*log(p) - log(-2*log(p) ) - log(2*pi) ))
  
}

ApproxPnorm <- function(q){
  
  return( -log(2) - q^2/2 - 0.5*log(2*pi) - 0.5*log(q) )
  
}

ComputeConditionalUs <- function( M = 4, p = 4, U, PhiList ){
  
  N <- length(U)
  QnormU <- qnorm(U)
  ExtremelyHighID <- which( U > 0.9999683 )
  ExtremelyLowID <- which( U < 0.00003167 )
  QnormU[ExtremelyHighID] <-  -ApproxQnorm(1-U[ExtremelyHighID])
  QnormU[ExtremelyLowID] <-  ApproxQnorm(U[ExtremelyLowID])
  ConditionalQnormU.MAT <- diag(QnormU)
  ConditionalU.MAT <- diag(U)
  cores=detectCores()
  
  for(r in 1:( (p+1)*M-1)){
    
    ivec <- seq(r+1,N,1)
    jvec <- ivec - r
    svec <- ceiling(jvec/M)
    tvec <- ceiling(ivec/M)
    kvec <- tvec - svec
    l1vec <- ivec - M*(tvec-1)
    l2vec <- jvec - M*(svec-1)
    
    QNORM.u.i.jPlusOne <- diag(ConditionalQnormU.MAT[ivec,(jvec+1)])
    QNORM.u.j.iMinusOne <- diag(ConditionalQnormU.MAT[jvec,(ivec-1)])
    
    
    
    phi.temp <- c(apply(cbind(l1vec,l2vec,pmin(kvec+1, rep(p+1,length(kvec)) )), 1, function(s){ PhiList[s[1],s[2],s[3]] }))
    #browser()
    h.uijPlus.preNorm <- ( QNORM.u.i.jPlusOne - phi.temp*QNORM.u.j.iMinusOne )/sqrt( 1 - phi.temp^2  )
    h.ujiMinus.preNorm <- ( QNORM.u.j.iMinusOne - phi.temp*QNORM.u.i.jPlusOne )/sqrt( 1 - phi.temp^2  )

    ConditionalQnormU.MAT[row(ConditionalQnormU.MAT)==(col(ConditionalQnormU.MAT)+r)] <- diag(ConditionalQnormU.MAT[ivec,(jvec+1)])*as.numeric(kvec>p) + h.uijPlus.preNorm*as.numeric(kvec<=p)
    ConditionalQnormU.MAT[col(ConditionalQnormU.MAT)==(row(ConditionalQnormU.MAT)+r)] <- diag( ConditionalQnormU.MAT[jvec,(ivec-1)] )*as.numeric(kvec>p) + h.ujiMinus.preNorm*as.numeric(kvec<=p)
    if(sum(is.nan(ConditionalQnormU.MAT[row(ConditionalQnormU.MAT)==(col(ConditionalQnormU.MAT)+r)]))>0){
      
      browser()
      
    }
    if(sum(is.nan(ConditionalQnormU.MAT[col(ConditionalQnormU.MAT)==(row(ConditionalQnormU.MAT)+r)]))>0){
      
      browser()
      
    }
    tempQnorm.lower <- ConditionalQnormU.MAT[row(ConditionalQnormU.MAT)==(col(ConditionalQnormU.MAT)+r)]
    tempQnorm.upper <- ConditionalQnormU.MAT[col(ConditionalQnormU.MAT)==(row(ConditionalQnormU.MAT)+r)]
    LowerTri.update <- pnorm(tempQnorm.lower)
    UpperTri.update <- pnorm( tempQnorm.upper )
    Extreme.LowerTri <- which(abs(tempQnorm.lower)> 4)
    Extreme.UpperTri <- which(abs(tempQnorm.upper)> 4)
    tempLower <- LowerTri.update
    tempUpper <- UpperTri.update
    if(length(Extreme.LowerTri)> 0 ){
      
      
      tempLower[Extreme.LowerTri] <- (1-exp(ApproxPnorm( abs(tempQnorm.lower[Extreme.LowerTri]) )))*as.numeric(sign(tempQnorm.lower[Extreme.LowerTri])==1) + exp(ApproxPnorm( abs(tempQnorm.lower[Extreme.LowerTri]) ))*as.numeric(sign(tempQnorm.lower[Extreme.LowerTri])==-1)
      
    }
    if(length(Extreme.UpperTri)>0){
      
      tempUpper[Extreme.UpperTri] <- (1-exp(ApproxPnorm( abs(tempQnorm.upper[Extreme.UpperTri]) )))*as.numeric(sign(tempQnorm.upper[Extreme.UpperTri])==1) + exp(ApproxPnorm( abs(tempQnorm.upper[Extreme.UpperTri]) ))*as.numeric(sign(tempQnorm.upper[Extreme.UpperTri])==-1)
      
    }
    
    ConditionalU.MAT[row(ConditionalU.MAT)==(col(ConditionalU.MAT)+r)] <- tempLower
    ConditionalU.MAT[col(ConditionalU.MAT)==(row(ConditionalU.MAT)+r)] <- tempUpper

    
  }
  
  return(list(ConditionalQnormU.MAT= QConditionalQnormU.MAT, ConditionalU.MAT=ConditionalU.MAT))
  
}


Truncated.SkewT.logPDF <- function(theta, dat.y){
  
  dat.length <- length(dat.y)
  xi <- theta[1]; omega = exp(theta[2]); alpha = theta[3]; nu = exp(theta[4])
  normConst <-  pst(x = 0, xi = xi, omega = omega, alpha = alpha, nu = nu, lower.tail = FALSE, log.p = TRUE)
  #The extra penalty dexp(nu, rate = 0.2, log = TRUE) + LaplacesDemon::dlaplace(x = alpha, scale = 5, log = TRUE)
  #encourages numerical stability. Otherwise, optimal alpha and nu may be very large
  Ans <- sum(dst(x=dat.y, xi = xi, omega = omega, alpha = alpha, nu = nu, log = TRUE)) - dat.length*normConst + dexp(nu, rate = 0.05, log = TRUE) + LaplacesDemon::dlaplace(x = alpha, scale = 10, log = TRUE)
  return(Ans)
  
}

Truncated.SkewT.logPDF.Plotter <- function(x, xi, omega, alpha, nu){
  
  normConst <- 1-pst(x = 0,xi = xi, omega = omega, alpha = alpha, nu = nu)
  Ans <- (dst(x,xi = xi, omega = omega, alpha = alpha, nu = nu)/normConst)*as.numeric(x>=0)
  return(Ans)
  
}

Truncated.SkewT.logCDF.Plotter <- function(x, xi, omega, alpha, nu){
  
  normConst <- 1-pst(x = 0,xi = xi, omega = omega, alpha = alpha, nu = nu)
  Ans <- ((pst(x,xi = xi, omega = omega, alpha = alpha, nu = nu) - pst(0,xi = xi, omega = omega, alpha = alpha, nu = nu))/normConst)*as.numeric(x>=0)
  return(Ans)
  
}

EstimateCorrelationsFromU <- function(BigU, p=4){
  
  BigT <- nrow(BigU)
  m <- ncol(BigU)
  N <- m * BigT
  CovArray <- array(0,dim = c(m,m,p+1))
  CovArray[,,1] <- var(BigU)
  for(k in 1:p){
    
    TransformedU <- cbind( BigU[1:(BigT-k),], BigU[(1+k):BigT,]  )
    CovArray[,,(k+1)] <- var(TransformedU)[1:m,(m+1):(2*m)]
    
    
  }
  CovMAT <- matrix(0,nrow=N,ncol=N)
  for(t in 1:BigT){
    for(tprime in t:(min(t+p,BigT)) ){
      
      k.temp <- tprime - t
      a.t <- (t-1)*m+1; b.t <- t*m
      a.tprime <- (tprime-1)*m+1; b.tprime <- tprime*m
      CovMAT[a.t:b.t, a.tprime:b.tprime] <- CovArray[,,k.temp+1]
      CovMAT[a.tprime:b.tprime, a.t:b.t] <- t(CovArray[,,k.temp+1])
      
    }
  }
  
  CorMAT <- cov2cor(CovMAT)
  return(CorMAT)
  
}

ConvertCorToPCor <- function( M, BigT, m, p=4 ){
  
  D <- ncol(M)
  PcorMAT <- matrix(0,nrow=D,ncol=D)
  diag(PcorMAT) <- rep(1,D)
  PcorMAT[ (col(PcorMAT))==(row(PcorMAT)+1) ] <- M[ (col(M))==(row(M)+1) ]
  PcorMAT[ (row(PcorMAT))==(col(PcorMAT)+1) ] <- M[ (row(M))==(col(M)+1) ]
  for(k in 2:(D-1)){
    for(j in 1:(D-k)){
      
      r1 <- matrix( M[j,(j+1):(j+k-1)], nrow=1)
      r3 <- matrix( M[j+k,(j+1):(j+k-1)], ncol=1)
      R2 <- M[(j+1):(j+k-1),(j+1):(j+k-1)]
      Denom <- sqrt(c(1 - r1%*%solve(R2)%*%t(r1) )*c(1 -  t(r3)%*%solve(R2)%*%r3 ) )
      PcorMAT[j,j+k] <- (M[j,j+k] - c(r1%*%solve(R2)%*%r3))/Denom
      PcorMAT[j+k,j] <- PcorMAT[j,j+k]
      
    }
  }
  
  PhiArray <- array(0,dim = c(m,m,p+1))
  for(k in 1:(p+1)){
    
    PhiArray[,,k] <- PcorMAT[1:m, ((k-1)*m+1):(k*m)]
    
  }
  
  N <- m*BigT
  PhiMAT <- matrix(0,nrow=N,ncol=N)
  for(t in 1:BigT){
    for(tprime in t:min(t+p,BigT)){
      
      k.temp <- tprime - t
      a.t <- (t-1)*m+1; b.t <- t*m
      a.tprime <- (tprime-1)*m+1; b.tprime <- tprime*m
      PhiMAT[a.t:b.t, a.tprime:b.tprime] <- PhiArray[,,k.temp+1]
      PhiMAT[a.tprime:b.tprime, a.t:b.t] <- t(PhiArray[,,k.temp+1])
      
    }
  }
  
  
  return(list(PhiArray=PhiArray,PhiMAT=PhiMAT))
  
}

ComputeLogKfuncSameT <- function(QNormUcondMat, tFrom, phiMat, m){
  
  a.t <- (tFrom-1)*m+1; b.t <- tFrom*m
  OuterRange <- seq(a.t+1, b.t, 1)
  #OuterRangeReps is synonymous with i
  OuterRangeReps <- rep(OuterRange, 1:(m-1))
  VV <- col(diag(1,m))
  #InnerRange is synonymous with i
  InnerRange <- (tFrom-1)*m + (t(VV))[upper.tri(VV)]
  QNorm.u.i.jPlusOne <- diag(QNormUcondMat[OuterRangeReps,InnerRange+1])
  QNorm.u.j.iMinusOne <- diag(QNormUcondMat[InnerRange, OuterRangeReps-1])
  l1 <- OuterRangeReps - m*(tFrom-1)
  l2 <- InnerRange - m*(tFrom-1)
  phi.l1.l2.k <- diag(phiMat[l1,l2])
  AnsK <- sum(-0.5*log(1-phi.l1.l2.k^2) + 0.5*( -phi.l1.l2.k^2*( QNorm.u.i.jPlusOne^2 + QNorm.u.j.iMinusOne^2 ) + 2*phi.l1.l2.k*QNorm.u.i.jPlusOne*QNorm.u.j.iMinusOne )/(1-phi.l1.l2.k^2))
  
  return(AnsK)
  
}

ComputeLogKfuncDifferentT <- function(QNormUcondMat, tFrom, tTo, phiMat, m){
  
  a.t <- (tFrom-1)*m+1; b.t <- tFrom*m
  a.s <- (tTo-1)*m+1; b.s <- tTo*m
  OuterRange <- seq(a.t, b.t, 1)
  #OuterRangeReps is synonymous with i
  #browser()
  OuterRangeReps <- rep(OuterRange, rep(m,m))
  #InnerRange is synonymous with i
  InnerRange <- rep(seq(a.s,b.s,1),m)
  QNorm.u.i.jPlusOne <- diag(QNormUcondMat[OuterRangeReps,InnerRange+1])
  QNorm.u.j.iMinusOne <- diag(QNormUcondMat[InnerRange, OuterRangeReps-1])
  l1 <- OuterRangeReps - m*(tFrom-1)
  l2 <- InnerRange - m*(tTo-1)
  phi.l1.l2.k <- diag(phiMat[l1,l2])
  AnsK <- sum(-0.5*log(1-phi.l1.l2.k^2) + 0.5*( -phi.l1.l2.k^2*( QNorm.u.i.jPlusOne^2 + QNorm.u.j.iMinusOne^2 ) + 2*phi.l1.l2.k*QNorm.u.i.jPlusOne*QNorm.u.j.iMinusOne )/(1-phi.l1.l2.k^2))

  return(AnsK)
  
}

ComputeLikelihoodCopulaPart <- function(QNormCondU.input, PhiArray.input, p =4, m){
  
  maxT <- ncol(QNormCondU.input)/m
  AnsLikCop <- ComputeLogKfuncSameT(QNormUcondMat = QNormCondU.input, tFrom = 1, phiMat = PhiArray.input[,,1], m = m)
  for(tInd in 2:maxT){
    
    AnsLikCop <- AnsLikCop + ComputeLogKfuncSameT(QNormUcondMat = QNormCondU.input, tFrom = tInd, phiMat = PhiArray.input[,,1], m = m)
    for(tPrime in (max(1,tInd-p)):(tInd-1) ){
      
      AnsLikCop <- AnsLikCop + ComputeLogKfuncDifferentT(QNormUcondMat = QNormCondU.input, tFrom= tInd, tTo =  tPrime, phiMat = PhiArray.input[,,tInd-tPrime+1], m = 4)

    }
    
  }
  return(AnsLikCop)
  
}


NewLogJointPosteriorUSMacro <- function( Y.mat, PAR.sub, p =4 ){
  
  BigT <- nrow(Y.mat)
  m <- ncol(Y.mat)
  N <- m * BigT
  
  xi.sub <- PAR.sub[1:m] 
  log.omega.sub <-  PAR.sub[(m+1):(2*m)]
  alpha.sub <-  PAR.sub[(2*m+1):(3*m)]
  log.nu.sub <-  PAR.sub[(3*m+1):(4*m)]
  
  #Need to isolate the first m(m-1)/2 entries as the first matrix in array
  PCorArray.sub <- array(0, dim = c(m,m,p+1))
  PCorArray.sub.first.slice <- matrix(0,nrow=m,ncol=m)
  PCorArray.sub.first.slice[lower.tri(PCorArray.sub.first.slice)] <- PAR.sub[(4*m+1):(4*m + m*(m-1)/2 )]
  PCorArray.sub.first.slice <- PCorArray.sub.first.slice + t(PCorArray.sub.first.slice)
  diag(PCorArray.sub.first.slice) <- Inf
  PCorArray.sub[,,1] <- PCorArray.sub.first.slice
  PCorArray.sub[,,2:(p+1)] <- array( PAR.sub[(4*m + m*(m-1)/2 + 1):( 4*m + m*(m-1)/2 + p*m^2 )], dim = c(m,m,p)  )
  Transformed.PCorArray.sub <- PCorArray.sub
  PCorArray.sub <- 2*pnorm(PCorArray.sub)-1
  log.Tau.sub <- PAR.sub[-(1:( 4*m + m*(m-1)/2 + p*m^2 ))]
  
  

  U1.sub <- pst(x = Y.mat[,1], xi = xi.sub[1], omega = exp(log.omega.sub[1]), alpha = alpha.sub[1], nu = exp(log.nu.sub[1]))
  U2.sub <- pst(x = Y.mat[,2], xi = xi.sub[2], omega = exp(log.omega.sub[2]), alpha = alpha.sub[2], nu = exp(log.nu.sub[2]))
  U3.sub <- Truncated.SkewT.logCDF.Plotter(x = Y.mat[,3], xi = xi.sub[3], omega = exp(log.omega.sub[3]), alpha = alpha.sub[3], nu = exp(log.nu.sub[3]))
  U4.sub <- Truncated.SkewT.logCDF.Plotter(x = Y.mat[,4], xi = xi.sub[4], omega = exp(log.omega.sub[4]), alpha = alpha.sub[4], nu = exp(log.nu.sub[4]))
 
  Ucombine.vec.sub <- c(rbind(U1.sub,U2.sub,U3.sub,U4.sub))
  
  ConditionalUObj.sub <- ComputeConditionalUs(U = Ucombine.vec.sub,PhiList = PCorArray.sub)
  LikCop <- ComputeLikelihoodCopulaPart(QNormCondU.input = ConditionalUObj.sub$ConditionalQnormU.MAT, PhiArray.input = PCorArray.sub, p = 4, m = 4)
  LikMarg1 <- sum(dst(x = Y.mat[,1], xi = xi.sub[1], omega = exp(log.omega.sub[1]), alpha = alpha.sub[1], nu = exp(log.nu.sub[1]), log = TRUE))
  LikMarg2 <- sum(dst(x = Y.mat[,2], xi = xi.sub[2], omega = exp(log.omega.sub[2]), alpha = alpha.sub[2], nu = exp(log.nu.sub[2]), log = TRUE))
  LikMarg3 <- Truncated.SkewT.logPDF(theta = c(xi.sub[3],log.omega.sub[3],alpha.sub[3],log.nu.sub[3]), dat.y = Y.mat[,3])
  LikMarg4 <- Truncated.SkewT.logPDF(theta = c(xi.sub[4],log.omega.sub[4],alpha.sub[4],log.nu.sub[4]), dat.y = Y.mat[,4])
  
  AnsLik <- LikCop + LikMarg1 + LikMarg2 + LikMarg3 + LikMarg4
  AnsPrior <- sum(dLogExpHalfCauchy(log.omega.sub, s.sq = 5)) + sum(dnorm(xi.sub, sd=100, log = TRUE)) + sum(LaplacesDemon::dlaplace(x = alpha.sub, scale = 10, log = TRUE)) + sum(dLogExpExp(x=log.nu.sub, lambda.rate=0.05)) 
  AnsPrior <- AnsPrior + sum(dnorm(PAR.sub[(4*m+1):(4*m + m*(m-1)/2 )],mean = 0, sd = exp(log.Tau.sub[1]), log = TRUE))  + sum(dnorm(c(Transformed.PCorArray.sub[,,2]),mean = 0, sd = exp(log.Tau.sub[2]), log = TRUE)) + sum(dnorm(c(Transformed.PCorArray.sub[,,3]),mean = 0, sd = exp(log.Tau.sub[3]), log = TRUE)) + sum(dnorm(c(Transformed.PCorArray.sub[,,4]),mean = 0, sd = exp(log.Tau.sub[4]), log = TRUE)) + sum(dnorm(c(Transformed.PCorArray.sub[,,5]),mean = 0, sd = exp(log.Tau.sub[5]), log = TRUE)) + sum(dLogExpHalfCauchy(log.Tau.sub, s.sq = 1))

  
  return(AnsLik+AnsPrior)
  
  
  
  
}

NegativeNewLogJointPosteriorUSMacro <- function( Y.mat, PAR.sub, p =4 ){
  
  BigT <- nrow(Y.mat)
  m <- ncol(Y.mat)
  N <- m * BigT
  
  xi.sub <- PAR.sub[1:m] 
  log.omega.sub <-  PAR.sub[(m+1):(2*m)]
  alpha.sub <-  PAR.sub[(2*m+1):(3*m)]
  log.nu.sub <-  PAR.sub[(3*m+1):(4*m)]
  
  #Need to isolate the first m(m-1)/2 entries as the first matrix in array
  PCorArray.sub <- array(0, dim = c(m,m,p+1))
  PCorArray.sub.first.slice <- matrix(0,nrow=m,ncol=m)
  PCorArray.sub.first.slice[lower.tri(PCorArray.sub.first.slice)] <- PAR.sub[(4*m+1):(4*m + m*(m-1)/2 )]
  PCorArray.sub.first.slice <- PCorArray.sub.first.slice + t(PCorArray.sub.first.slice)
  diag(PCorArray.sub.first.slice) <- Inf
  PCorArray.sub[,,1] <- PCorArray.sub.first.slice
  PCorArray.sub[,,2:(p+1)] <- array( PAR.sub[(4*m + m*(m-1)/2 + 1):( 4*m + m*(m-1)/2 + p*m^2 )], dim = c(m,m,p)  )
  Transformed.PCorArray.sub <- PCorArray.sub
  PCorArray.sub <- 2*pnorm(PCorArray.sub)-1
  log.Tau.sub <- PAR.sub[-(1:( 4*m + m*(m-1)/2 + p*m^2 ))]
  
  
  
  U1.sub <- pst(x = Y.mat[,1], xi = xi.sub[1], omega = exp(log.omega.sub[1]), alpha = alpha.sub[1], nu = exp(log.nu.sub[1]))
  U2.sub <- pst(x = Y.mat[,2], xi = xi.sub[2], omega = exp(log.omega.sub[2]), alpha = alpha.sub[2], nu = exp(log.nu.sub[2]))
  U3.sub <- Truncated.SkewT.logCDF.Plotter(x = Y.mat[,3], xi = xi.sub[3], omega = exp(log.omega.sub[3]), alpha = alpha.sub[3], nu = exp(log.nu.sub[3]))
  U4.sub <- Truncated.SkewT.logCDF.Plotter(x = Y.mat[,4], xi = xi.sub[4], omega = exp(log.omega.sub[4]), alpha = alpha.sub[4], nu = exp(log.nu.sub[4]))
  
  Ucombine.vec.sub <- c(rbind(U1.sub,U2.sub,U3.sub,U4.sub))
  
  ConditionalUObj.sub <- ComputeConditionalUs(U = Ucombine.vec.sub,PhiList = PCorArray.sub)
  LikCop <- ComputeLikelihoodCopulaPart(QNormCondU.input = ConditionalUObj.sub$ConditionalQnormU.MAT, PhiArray.input = PCorArray.sub, p = 4, m = 4)
  LikMarg1 <- sum(dst(x = Y.mat[,1], xi = xi.sub[1], omega = exp(log.omega.sub[1]), alpha = alpha.sub[1], nu = exp(log.nu.sub[1]), log = TRUE))
  LikMarg2 <- sum(dst(x = Y.mat[,2], xi = xi.sub[2], omega = exp(log.omega.sub[2]), alpha = alpha.sub[2], nu = exp(log.nu.sub[2]), log = TRUE))
  LikMarg3 <- sum(log(Truncated.SkewT.logPDF.Plotter(x = Y.mat[, 3], xi = xi.sub[3], omega = exp(log.omega.sub[3]), alpha = alpha.sub[3], nu = exp(log.nu.sub[3])) ))
  LikMarg4 <- sum(log(Truncated.SkewT.logPDF.Plotter(x = Y.mat[, 4], xi = xi.sub[4], omega = exp(log.omega.sub[4]), alpha = alpha.sub[4], nu = exp(log.nu.sub[4])) ))
  AnsLik <- LikCop + LikMarg1 + LikMarg2 + LikMarg3 + LikMarg4
  AnsPrior <- sum(dLogExpHalfCauchy(log.omega.sub, s.sq = 5)) + sum(dnorm(xi.sub, sd=100, log = TRUE)) + sum(LaplacesDemon::dlaplace(x = alpha.sub, scale = 10, log = TRUE)) + sum(dLogExpExp(x=log.nu.sub, lambda.rate=0.05)) 
  AnsPrior <- AnsPrior + sum(dnorm(PAR.sub[(4*m+1):(4*m + m*(m-1)/2 )],mean = 0, sd = exp(0.5*log.Tau.sub[1]), log = TRUE))  + sum(dnorm(c(Transformed.PCorArray.sub[,,2]),mean = 0, sd = exp(0.5*log.Tau.sub[2]), log = TRUE)) + sum(dnorm(c(Transformed.PCorArray.sub[,,3]),mean = 0, sd = exp(0.5*log.Tau.sub[3]), log = TRUE)) + sum(dnorm(c(Transformed.PCorArray.sub[,,4]),mean = 0, sd = exp(0.5*log.Tau.sub[4]), log = TRUE)) + sum(dnorm(c(Transformed.PCorArray.sub[,,5]),mean = 0, sd = exp(0.5*log.Tau.sub[5]), log = TRUE)) + sum(dLogExpHalfCauchy(log.Tau.sub, s.sq = 1))

  
  return(-AnsLik-AnsPrior)
  
  
  
  
}

NewUncutPosteriorUSMacro <- function(Y.DAT, xi.init, log.omega.init, alpha.init, log.nu.init, PCor.init, p=4, mu.init, C.init, maxITERS=2000){
  
  BigT <- nrow(Y.DAT)
  m <- ncol(Y.DAT)
  N <- m * BigT
  logTauInit <- rep(log(3)^2, p+1)
  
  numPars <- 4*m + m*(m-1)/2 + p*m^2 + p + 1
  zDrawDictionary <- matrix(rnorm(maxITERS*numPars),nrow=maxITERS, ncol=numPars)
  MuSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  CSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  
  AccuGrad.Mu <- rep(0, numPars)
  AccuChange.Mu <- rep(0, numPars)
  AccuGrad.C <- rep(0, numPars)
  AccuChange.C <- rep(0, numPars)
  FirstSlice.init <- PCor.init[,,1]
  
  #mu.curr <- c(xi.init, log.omega.init, alpha.init, log.nu.init, FirstSlice.init[lower.tri(FirstSlice.init)], c(PCor.init[,,2:(p+1)]), logTauInit )
  #C.curr <- rep(0.1,numPars)
  mu.curr <- mu.init
  C.curr <- C.init
  eps <- 10^(-5); delta <- 0.9
  
  for(it in 1:maxITERS){
    
    zDraw <- zDrawDictionary[it,]
    theta.draw <- mu.curr + C.curr * zDraw
    #browser()
    grad.mu <- numDeriv::grad(func = NewLogJointPosteriorUSMacrocpp, method = "simple", x = theta.draw, Ymat = Y.DAT, p=4)
    #browser()
    AccuGrad.Mu <- delta*AccuGrad.Mu + (1 - delta)*( grad.mu^2 )
    Change.Mu <- sqrt( AccuChange.Mu + eps )/sqrt( AccuGrad.Mu + eps )*grad.mu
    mu.new <- mu.curr + Change.Mu
    AccuChange.Mu <- delta*AccuChange.Mu + (1 - delta)*( Change.Mu^2 )
    #browser()
    grad.C <- grad.mu*zDraw + 1/C.curr
    AccuGrad.C <- delta*AccuGrad.C + (1 - delta)*( grad.C^2 )
    Change.C <- sqrt( AccuChange.C + eps )/sqrt( AccuGrad.C + eps )*grad.C
    C.new <- C.curr + Change.C
    AccuChange.C <- delta*AccuChange.C + (1 - delta)*( Change.C^2 )
    #browser()
    MuSTORE[it,] <- mu.new
    CSTORE[it,] <- C.new
    mu.curr <- mu.new
    C.curr <- C.new
    #browser()
    
    if((it%%10)==0){
      
      cat("Completed iter: ", it, "\n")
      
    }
    
    
  }
  
  return(list( MuSTORE=MuSTORE, CSTORE=CSTORE ))
  
}


NewReverseCutPosteriorUSMacro.StageI <- function(Y.DAT, p=4, mu.init, C.init, maxITERS=2000){
  
  BigT <- nrow(Y.DAT)
  m <- ncol(Y.DAT)
  N <- m * BigT
  logTauInit <- rep(log(3)^2, p+1)
  RankMAT.fix <- cbind(rank(Y.DAT[,1], ties.method = "first"), rank(Y.DAT[,2], ties.method = "first"), rank(Y.DAT[,3], ties.method = "first"), rank(Y.DAT[,4], ties.method = "first"))
  
  numPars <- length(mu.init)
  zDrawDictionary <- matrix(rnorm(maxITERS*numPars),nrow=maxITERS, ncol=numPars)
  MuSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  CSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  
  AccuGrad.Mu <- rep(0, numPars)
  AccuChange.Mu <- rep(0, numPars)
  AccuGrad.C <- rep(0, numPars)
  AccuChange.C <- rep(0, numPars)
  
  mu.curr <- mu.init
  C.curr <- C.init
  eps <- 10^(-5); delta <- 0.9
  
  for(it in 1:maxITERS){
    
    zDraw <- zDrawDictionary[it,]
    theta.draw <- mu.curr + C.curr * zDraw
    #browser()
    grad.mu <- numDeriv::grad(func = NewLogJointPosteriorUSMacroReverseStage1cpp, method = "simple", x = theta.draw, RankMAT = RankMAT.fix, p=4)
    #browser()
    AccuGrad.Mu <- delta*AccuGrad.Mu + (1 - delta)*( grad.mu^2 )
    Change.Mu <- sqrt( AccuChange.Mu + eps )/sqrt( AccuGrad.Mu + eps )*grad.mu
    mu.new <- mu.curr + Change.Mu
    AccuChange.Mu <- delta*AccuChange.Mu + (1 - delta)*( Change.Mu^2 )
    #browser()
    grad.C <- grad.mu*zDraw + 1/C.curr
    AccuGrad.C <- delta*AccuGrad.C + (1 - delta)*( grad.C^2 )
    Change.C <- sqrt( AccuChange.C + eps )/sqrt( AccuGrad.C + eps )*grad.C
    C.new <- C.curr + Change.C
    AccuChange.C <- delta*AccuChange.C + (1 - delta)*( Change.C^2 )
   #browser()
    MuSTORE[it,] <- mu.new
    CSTORE[it,] <- C.new
    mu.curr <- mu.new
    C.curr <- C.new
    #browser()
    
    if((it%%10)==0){
      
      cat("Completed iter: ", it, "\n")
      
    }
    
    
  }
  
  return(list( MuSTORE=MuSTORE, CSTORE=CSTORE ))
  
}

NewReverseCutPosteriorUSMacro.StageII <- function(Y.DAT, p=4, mu.init, mu.fix, C.init, C.fix, maxITERS=2000){
  
  BigT <- nrow(Y.DAT)
  m <- ncol(Y.DAT)
  N <- m * BigT
  logTauInit <- rep(log(3)^2, p+1)
  
  
  numPars <- 4*m + m*(m-1)/2 + p*m^2 + p + 1
  fix.inds <- seq(4*m+1, numPars,1)
  num.fix.inds <- length(fix.inds)
  zDrawDictionary <- matrix(rnorm(maxITERS*numPars),nrow=maxITERS, ncol=numPars)
  MuSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  CSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  
  AccuGrad.Mu <- rep(0, numPars)
  AccuChange.Mu <- rep(0, numPars)
  AccuGrad.C <- rep(0, numPars)
  AccuChange.C <- rep(0, numPars)
  
  mu.curr <- c(mu.init,mu.fix)
  C.curr <- c(C.init,C.fix)
  eps <- 10^(-5); delta <- 0.9
  
  for(it in 1:maxITERS){
    
    zDraw <- zDrawDictionary[it,]
    theta.draw <- mu.curr + C.curr * zDraw
    #browser()
    grad.mu <- numDeriv::grad(func = NewLogJointPosteriorUSMacrocpp, method = "simple", x = theta.draw, Ymat = Y.DAT, p=4)
    #browser()
    AccuGrad.Mu <- delta*AccuGrad.Mu + (1 - delta)*( grad.mu^2 )
    Change.Mu <- sqrt( AccuChange.Mu + eps )/sqrt( AccuGrad.Mu + eps )*grad.mu
    Change.Mu[fix.inds] <- 0
    mu.new <- mu.curr + Change.Mu
    AccuChange.Mu <- delta*AccuChange.Mu + (1 - delta)*( Change.Mu^2 )
    #browser()
    grad.C <- grad.mu*zDraw + 1/C.curr
    AccuGrad.C <- delta*AccuGrad.C + (1 - delta)*( grad.C^2 )
    Change.C <- sqrt( AccuChange.C + eps )/sqrt( AccuGrad.C + eps )*grad.C
    Change.C[fix.inds] <- 0
    C.new <- C.curr + Change.C
    AccuChange.C <- delta*AccuChange.C + (1 - delta)*( Change.C^2 )
    #browser()
    MuSTORE[it,] <- mu.new
    CSTORE[it,] <- C.new
    mu.curr <- mu.new
    C.curr <- C.new
    #browser()
    
    if((it%%10)==0){
      
      cat("Completed iter: ", it, "\n")
      
    }
    
    
  }
  
  return(list( MuSTORE=MuSTORE, CSTORE=CSTORE ))
  
}

NewCutIPosteriorSkewTMarginalUSMacro <- function(vy, mu.init, C.init, maxITERS=2000){
  
  BigT <- length(vy)
  numPars <- 4
  zDrawDictionary <- matrix(rnorm(maxITERS*numPars),nrow=maxITERS, ncol=numPars)
  MuSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  CSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  
  AccuGrad.Mu <- rep(0, numPars)
  AccuChange.Mu <- rep(0, numPars)
  AccuGrad.C <- rep(0, numPars)
  AccuChange.C <- rep(0, numPars)
  
  mu.curr <- mu.init
  C.curr <- C.init
  eps <- 10^(-5); delta <- 0.9
  
  for(it in 1:maxITERS){
    
    zDraw <- zDrawDictionary[it,]
    theta.draw <- mu.curr + C.curr * zDraw
    #browser()
    grad.mu <- numDeriv::grad(func = NewLogJointPosteriorSkewTUSMacroCutStage1cpp, method = "simple", x = theta.draw, Yvec = vy)
    #browser()
    AccuGrad.Mu <- delta*AccuGrad.Mu + (1 - delta)*( grad.mu^2 )
    Change.Mu <- sqrt( AccuChange.Mu + eps )/sqrt( AccuGrad.Mu + eps )*grad.mu
    mu.new <- mu.curr + Change.Mu
    AccuChange.Mu <- delta*AccuChange.Mu + (1 - delta)*( Change.Mu^2 )
    #browser()
    grad.C <- grad.mu*zDraw + 1/C.curr
    AccuGrad.C <- delta*AccuGrad.C + (1 - delta)*( grad.C^2 )
    Change.C <- sqrt( AccuChange.C + eps )/sqrt( AccuGrad.C + eps )*grad.C
    C.new <- C.curr + Change.C
    AccuChange.C <- delta*AccuChange.C + (1 - delta)*( Change.C^2 )
    #browser()
    MuSTORE[it,] <- mu.new
    CSTORE[it,] <- C.new
    mu.curr <- mu.new
    C.curr <- C.new
    #browser()
    
    if((it%%10)==0){
      
      cat("Completed iter: ", it, "\n")
      
    }
    
    
  }
  
  return(list( MuSTORE=MuSTORE, CSTORE=CSTORE ))
  
}

NewCutIPosteriorTruncatedSkewTMarginalUSMacro <- function(vy, mu.init, C.init, maxITERS=2000){
  
  BigT <- length(vy)
  numPars <- 4
  zDrawDictionary <- matrix(rnorm(maxITERS*numPars),nrow=maxITERS, ncol=numPars)
  MuSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  CSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  
  AccuGrad.Mu <- rep(0, numPars)
  AccuChange.Mu <- rep(0, numPars)
  AccuGrad.C <- rep(0, numPars)
  AccuChange.C <- rep(0, numPars)
  
  mu.curr <- mu.init
  C.curr <- C.init
  eps <- 10^(-5); delta <- 0.9
  
  for(it in 1:maxITERS){
    
    zDraw <- zDrawDictionary[it,]
    theta.draw <- mu.curr + C.curr * zDraw
    #browser()
    grad.mu <- numDeriv::grad(func = NewLogJointPosteriorTruncatedSkewTUSMacroCutStage1cpp, method = "simple", x = theta.draw, Yvec = vy)
    #browser()
    AccuGrad.Mu <- delta*AccuGrad.Mu + (1 - delta)*( grad.mu^2 )
    Change.Mu <- sqrt( AccuChange.Mu + eps )/sqrt( AccuGrad.Mu + eps )*grad.mu
    mu.new <- mu.curr + Change.Mu
    AccuChange.Mu <- delta*AccuChange.Mu + (1 - delta)*( Change.Mu^2 )
    #browser()
    grad.C <- grad.mu*zDraw + 1/C.curr
    AccuGrad.C <- delta*AccuGrad.C + (1 - delta)*( grad.C^2 )
    Change.C <- sqrt( AccuChange.C + eps )/sqrt( AccuGrad.C + eps )*grad.C
    C.new <- C.curr + Change.C
    AccuChange.C <- delta*AccuChange.C + (1 - delta)*( Change.C^2 )
    #browser()
    MuSTORE[it,] <- mu.new
    CSTORE[it,] <- C.new
    mu.curr <- mu.new
    C.curr <- C.new
    #browser()
    
    if((it%%10)==0){
      
      cat("Completed iter: ", it, "\n")
      
    }
    
    
  }
  
  return(list( MuSTORE=MuSTORE, CSTORE=CSTORE ))
  
}

NewCutIPosteriorUSMacro.StageII <- function(Y.DAT, p=4, mu.init, mu.fix, C.init, C.fix, maxITERS=2000){
  
  BigT <- nrow(Y.DAT)
  m <- ncol(Y.DAT)
  N <- m * BigT
  logTauInit <- rep(log(3)^2, p+1)
  
  
  numPars <- 4*m + m*(m-1)/2 + p*m^2 + p + 1
  fix.inds <- seq(1, 4*m,1)
  num.fix.inds <- length(fix.inds)
  zDrawDictionary <- matrix(rnorm(maxITERS*numPars),nrow=maxITERS, ncol=numPars)
  MuSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  CSTORE <- matrix(0,nrow=maxITERS, ncol=numPars)
  
  AccuGrad.Mu <- rep(0, numPars)
  AccuChange.Mu <- rep(0, numPars)
  AccuGrad.C <- rep(0, numPars)
  AccuChange.C <- rep(0, numPars)
  
  mu.curr <- c(mu.fix,mu.init)
  C.curr <- c(C.fix,C.init)
  eps <- 10^(-5); delta <- 0.9
  
  for(it in 1:maxITERS){
    
    zDraw <- zDrawDictionary[it,]
    
    theta.draw <- mu.curr + C.curr * zDraw
    #browser()
    grad.mu <- numDeriv::grad(func = NewLogJointPosteriorUSMacrocpp, method = "simple", x = theta.draw, Ymat = Y.DAT, p=4)
    #browser()
    AccuGrad.Mu <- delta*AccuGrad.Mu + (1 - delta)*( grad.mu^2 )
    Change.Mu <- sqrt( AccuChange.Mu + eps )/sqrt( AccuGrad.Mu + eps )*grad.mu
    Change.Mu[fix.inds] <- 0
    mu.new <- mu.curr + Change.Mu
    AccuChange.Mu <- delta*AccuChange.Mu + (1 - delta)*( Change.Mu^2 )
    #browser()
    grad.C <- grad.mu*zDraw + 1/C.curr
    AccuGrad.C <- delta*AccuGrad.C + (1 - delta)*( grad.C^2 )
    Change.C <- sqrt( AccuChange.C + eps )/sqrt( AccuGrad.C + eps )*grad.C
    Change.C[fix.inds] <- 0
    C.new <- C.curr + Change.C
    AccuChange.C <- delta*AccuChange.C + (1 - delta)*( Change.C^2 )
    MuSTORE[it,] <- mu.new
    CSTORE[it,] <- C.new
    mu.curr <- mu.new
    C.curr <- C.new
    #browser()
    
    if((it%%10)==0){
      
      cat("Completed iter: ", it, "\n")
      
    }
    
    
  }
  
  return(list( MuSTORE=MuSTORE, CSTORE=CSTORE ))
  
}
