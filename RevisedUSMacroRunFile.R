library(doParallel)
library(foreach)
library(Rcpp)
library(RcppArmadillo)
library(sn)
source("RevisedUSMacroSupportingFunctions.R")
sourceCpp("RevisedUSMacroSupportingFunctions.cpp")


USmacro.dat <- read.csv(file = "macro2022Q3.csv",header=TRUE)


SELMobj1 <- selm(GDP.Growth~1, family = "ST", data = USmacro.dat, method = "MLE")
Sobj1 <- summary(SELMobj1, "dp")
Marginal.MLE1 <- Sobj1@param.table[,1]

SELMobj2 <- selm(Inflation~1, family = "ST", data = USmacro.dat, method = "MLE")
Sobj2 <- summary(SELMobj2, "dp")
Marginal.MLE2 <- Sobj2@param.table[,1]

Marginal.MLE3 <- optim(par = c(-1,5,-3,7), fn = Truncated.SkewT.logPDF, method = "L-BFGS-B", control = list(fnscale=-1), dat.y = USmacro.dat$UR)$par
Marginal.MLE3 <- c(Marginal.MLE3[1], exp(Marginal.MLE3[2]), Marginal.MLE3[3], exp(Marginal.MLE3[4]))

Marginal.MLE4 <- optim(par = c(-1,5,-3,7), fn = Truncated.SkewT.logPDF, method = "L-BFGS-B", control = list(fnscale=-1), dat.y = USmacro.dat$IR)$par
Marginal.MLE4 <- c(Marginal.MLE4[1], exp(Marginal.MLE4[2]), Marginal.MLE4[3], exp(Marginal.MLE4[4]))

Y.input <- cbind(USmacro.dat$GDP.Growth, USmacro.dat$Inflation, USmacro.dat$UR, USmacro.dat$IR)
xi.input <- c(Marginal.MLE1[1],Marginal.MLE2[1], Marginal.MLE3[1], Marginal.MLE4[1])
log.omega.input <- log(c(Marginal.MLE1[2],Marginal.MLE2[2], Marginal.MLE3[2], Marginal.MLE4[2]))
alpha.input <- c(Marginal.MLE1[3],Marginal.MLE2[3], Marginal.MLE3[3], Marginal.MLE4[3])
log.nu.input <- log(c(Marginal.MLE1[4],Marginal.MLE2[4], Marginal.MLE3[4], Marginal.MLE4[4]))

U1.input <- pst(x = Y.input[,1], xi = xi.input[1], omega = exp(log.omega.input[1]), alpha = alpha.input[1], nu = exp(log.nu.input[1]))
U2.input <- pst(x = Y.input[,2], xi = xi.input[2], omega = exp(log.omega.input[2]), alpha = alpha.input[2], nu = exp(log.nu.input[2]))
U3.input <- Truncated.SkewT.logCDF.Plotter(x = Y.input[,3], xi = xi.input[3], omega = exp(log.omega.input[3]), alpha = alpha.input[3], nu = exp(log.nu.input[3]))
U4.input <- Truncated.SkewT.logCDF.Plotter(x = Y.input[,4], xi = xi.input[4], omega = exp(log.omega.input[4]), alpha = alpha.input[4], nu = exp(log.nu.input[4]))

UcombineMAT.input <- cbind(U1.input,U2.input,U3.input,U4.input)
Ucombine.vec.input <- c(t(UcombineMAT.input))
Cor.input <- EstimateCorrelationsFromU(BigU = cbind(U1.input,U2.input, U3.input, U4.input))
PCorArray.input <- ConvertCorToPCor(Cor.input[1:20,1:20], BigT = length(USmacro.dat$GDP.Growth), m = 4)$PhiArray
TransformedPCorArray.input <- qnorm((PCorArray.input+1)/2)
PAR.input <- unname(c(xi.input, log.omega.input, alpha.input, log.nu.input, (TransformedPCorArray.input[,,1])[lower.tri(TransformedPCorArray.input[,,1])],c(TransformedPCorArray.input[,,2:5])))
PAR.input <- c(PAR.input, rep(log(0.3^2), 5))

NewLogJointPosteriorUSMacro(Y.mat = Y.input, PAR.sub = PAR.input, p = 4)
NewLogJointPosteriorUSMacrocpp(Ymat = Y.input, PARsub = PAR.input, p = 4)

Gradout <- numDeriv::grad(func = NewLogJointPosteriorUSMacrocpp, method = "simple", x = PAR.input, Ymat = Y.input, p=4)

Optim.test <- optim(par = PAR.input, fn = NewLogJointPosteriorUSMacrocpp, method = "L-BFGS-B", hessian = TRUE, control = list(fnscale=-1),  Ymat = Y.input, p=4)
write.csv( x = data.frame( mu = Optim.test$par, C = sqrt(-1/diag(Optim.test$hessian))), file = "UncutPosteriorStartingValues.csv")
UNCUT.INIT = read.csv("UncutPosteriorStartingValues.csv",header=TRUE)
VB.USmacro.uncut <- NewUncutPosteriorUSMacro(Y.DAT = Y.input, xi.init = xi.input, log.omega.init = log.omega.input, alpha.init = alpha.input, log.nu.init = log.nu.input,  PCor.init = TransformedPCorArray.input, p = 4, maxITERS = 5000, mu.init =  UNCUT.INIT$mu, C.init = UNCUT.INIT$C )
min(which(apply(VB.USmacro.uncut$MuSTORE,1,function(s){ sum(is.nan(s))>0 })))


DrawSize <- 5000
lag.p <- 4
PredictiveDensityValue <- matrix(0,nrow=183,ncol=(8*4))
BANDWIDTH.fix <- scan(file = "BANDWIDTHs",sep=",")
Forward.t.fix <- floor((1:32 - 1)/4)+1

for(t.train.end in 84:266){
  
  Train.PART <- Y.input[1:t.train.end,]
  Test.PART <- Y.input[-(1:t.train.end),]
  VB.uncut.obj.temp <-NewUncutPosteriorUSMacro(Y.DAT = Train.PART, xi.init = xi.input, log.omega.init = log.omega.input, alpha.init = alpha.input, log.nu.init = log.nu.input, PCor.init = TransformedPCorArray.input, p = 4, maxITERS = 12500, mu.init =  UNCUT.INIT$mu, C.init = UNCUT.INIT$C)
  u1.draw.MAT <- matrix(0,nrow=lag.p,ncol=DrawSize)
  u2.draw.MAT <- matrix(0,nrow=lag.p,ncol=DrawSize)
  u3.draw.MAT <- matrix(0,nrow=lag.p,ncol=DrawSize)
  u4.draw.MAT <- matrix(0,nrow=lag.p,ncol=DrawSize)
  
  MuFINAL.temp <- colMeans(VB.uncut.obj.temp$MuSTORE[10001:12500,1:86])
  VarFINAL.temp <- colMeans(VB.uncut.obj.temp$CSTORE[10001:12500,1:86]^2)
  ParDraws.temp <- mvtnorm::rmvnorm(n = 5000, mean =  MuFINAL.temp, sigma = diag(VarFINAL.temp))
  W.Draws.temp <- matrix(0,nrow=5000,ncol=(8*4))
  Y.Draws.temp <- matrix(0,nrow=5000,ncol=(8*4))
  for(l in 1:5000){
    
    u1.draw.MAT[,l] <- SkewT_cdf(Y.input[(t.train.end-lag.p+1):t.train.end,1], xi = ParDraws.temp[l,1], log_omega = ParDraws.temp[l,2], alpha = ParDraws.temp[l,3], log_nu = ParDraws.temp[l,4])
    u2.draw.MAT[,l] <- SkewT_cdf(Y.input[(t.train.end-lag.p+1):t.train.end,2], xi = ParDraws.temp[l,5], log_omega = ParDraws.temp[l,6], alpha = ParDraws.temp[l,7], log_nu = ParDraws.temp[l,8])
    u3.draw.MAT[,l] <- TruncatedSkewT_cdf(Y.input[(t.train.end-lag.p+1):t.train.end,3], xi = ParDraws.temp[l,9], log_omega = ParDraws.temp[l,10], alpha = ParDraws.temp[l,11], log_nu = ParDraws.temp[l,12])
    u4.draw.MAT[,l] <- TruncatedSkewT_cdf(Y.input[(t.train.end-lag.p+1):t.train.end,4], xi = ParDraws.temp[l,13], log_omega = ParDraws.temp[l,14], alpha = ParDraws.temp[l,15], log_nu = ParDraws.temp[l,16])
    w.temp <- qnorm(c(rbind(u1.draw.MAT[,l],u2.draw.MAT[,l],u3.draw.MAT[,l],u4.draw.MAT[,l])))
    PhiArray.temp <- MakePhiArrayRcpp(ParDraws.temp[l,17:86])
    CorMAT.temp <- ConvertPCorToCorcpp(PhiArray = PhiArray.temp, m = 4)
    Gaussian.Cov.Temp <- CorMAT.temp[17:20,17:20] - CorMAT.temp[17:20,1:16]%*%solve(CorMAT.temp[1:16,1:16])%*%CorMAT.temp[1:16,17:20]
    Historical.temp <- w.temp
    for(h in 1:8){
      
      Gaussian.Mean.Temp <- c(CorMAT.temp[17:20,1:16]%*%solve(CorMAT.temp[1:16,1:16])%*%Historical.temp)
      W.Draws.temp[l,((h-1)*4+1):(4*h)] <- mvtnorm::rmvnorm(n=1, mean = Gaussian.Mean.Temp, sigma = Gaussian.Cov.Temp)
      Historical.temp <- c(Historical.temp[-(1:4)],W.Draws.temp[l,((h-1)*4+1):(4*h)])
    }
    U1.predict <- pnorm(W.Draws.temp[l,seq(1,32,4)])
    U2.predict <- pnorm(W.Draws.temp[l,seq(2,32,4)])
    U3.predict <-pnorm(W.Draws.temp[l,seq(3,32,4)])
    U4.predict <- pnorm(W.Draws.temp[l,seq(4,32,4)])
    
    U1.predict <- pmin(pmax( rep(0.001,8), U1.predict ), rep(0.999,8))
    U2.predict <- pmin(pmax( rep(0.001,8), U2.predict ), rep(0.999,8))
    Y.Draws.temp[l,seq(1,32,4)] <- qst(U1.predict, xi = ParDraws.temp[l,1], omega = exp(ParDraws.temp[l,2]), alpha = ParDraws.temp[l,3], nu = exp(ParDraws.temp[l,4]))
    Y.Draws.temp[l,seq(2,32,4)] <- qst(U2.predict, xi = ParDraws.temp[l,5], omega = exp(ParDraws.temp[l,6]), alpha = ParDraws.temp[l,7], nu = exp(ParDraws.temp[l,8]))
    normConst3.temp <- pst(0,xi = ParDraws.temp[l,9], omega = exp(ParDraws.temp[l,10]), alpha = ParDraws.temp[l,11], nu = exp(ParDraws.temp[l,12]))
    U3prime.temp <- U3.predict*(1 - normConst3.temp) + normConst3.temp
    U3prime.temp <- pmin(pmax( rep(0.01,8), U3prime.temp ), rep(0.85,8))
    Y.Draws.temp[l,seq(3,32,4)] <- qst(U3prime.temp , xi = ParDraws.temp[l,9], omega = exp(ParDraws.temp[l,10]), alpha = ParDraws.temp[l,11], nu = exp(ParDraws.temp[l,12]))
    normConst4.temp <- pst(0,xi = ParDraws.temp[l,13], omega = exp(ParDraws.temp[l,14]), alpha = ParDraws.temp[l,15], nu = exp(ParDraws.temp[l,16]))
    U4prime.temp <- U4.predict*(1 - normConst4.temp) + normConst4.temp
    U4prime.temp <- pmin(pmax( rep(0.01,8), U4prime.temp ), rep(0.85,8))
    Y.Draws.temp[l,seq(4,32,4)] <- qst(U4prime.temp , xi = ParDraws.temp[l,13], omega = exp(ParDraws.temp[l,14]), alpha = ParDraws.temp[l,15], nu = exp(ParDraws.temp[l,16]))
    
  }
  

  for(j in seq(1,32,4)){
    
    PredictiveDensityValue[t.train.end-83,j] <- log(mean(dnorm( (Y.input[t.train.end+Forward.t.fix[j],1]- Y.Draws.temp[,j])/BANDWIDTH.fix[j] ), na.rm=TRUE)/BANDWIDTH.fix[j])
  
  }
  

  for(j in seq(2,32,4)){
    
    PredictiveDensityValue[t.train.end-83,j] <- log(mean(dnorm( (Y.input[t.train.end+Forward.t.fix[j],2]- Y.Draws.temp[,j])/BANDWIDTH.fix[j] ), na.rm=TRUE)/BANDWIDTH.fix[j])
    
  }
  
  for(j in seq(3,32,4)){
    
    PredictiveDensityValue[t.train.end-83,j] <- log(mean(dnorm( (Y.input[t.train.end+Forward.t.fix[j],3]- Y.Draws.temp[,j])/BANDWIDTH.fix[j] ), na.rm=TRUE)/BANDWIDTH.fix[j])
    
  }
  
  for(j in seq(4,32,4)){
    
    PredictiveDensityValue[t.train.end-83,j] <- log(mean(dnorm( (Y.input[t.train.end+Forward.t.fix[j],4]- Y.Draws.temp[,j])/BANDWIDTH.fix[j] ), na.rm=TRUE)/BANDWIDTH.fix[j])
    
  }
  save.image("USmacroUncut.RData")
  
}



