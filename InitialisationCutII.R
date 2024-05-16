VBmarginal1.obj <- NewVB.cut.ST.USMacro(y = y1, mu.init = MarginalPar[1:4], C.init = diag(0.1,4), maxRuns = 20000)
VBmarginal2.obj <- NewVB.cut.ST.USMacro(y = y2, mu.init = MarginalPar[5:8], C.init = diag(0.1,4), maxRuns = 20000)
VBmarginal3.obj <- NewVB.cut.TruncST.USMacro(y = y3, mu.init = MarginalPar[9:12], C.init = diag(0.1,4), maxRuns = 20000)
VBmarginal4.obj <- NewVB.cut.TruncST.USMacro(y = y4, mu.init = MarginalPar[13:16], C.init = diag(0.1,4), maxRuns = 20000)


C.init.margin <- as.matrix(bdiag(VBmarginal1.obj$Cstore[,,20000], VBmarginal2.obj$Cstore[,,20000], VBmarginal3.obj$Cstore[,,20000], VBmarginal4.obj$Cstore[,,20000]))
margin1.est <- c(VBmarginal1.obj$mustore[20000,1], exp(VBmarginal1.obj$mustore[20000,2] + 0.5*VBmarginal1.obj$Cstore[2,2,20000]^2), VBmarginal1.obj$mustore[20000,3],  exp(VBmarginal1.obj$mustore[20000,4] + 0.5*VBmarginal1.obj$Cstore[4,4,20000]^2))
margin2.est <- c(VBmarginal2.obj$mustore[20000,1], exp(VBmarginal2.obj$mustore[20000,2] + 0.5*VBmarginal2.obj$Cstore[2,2,20000]^2), VBmarginal2.obj$mustore[20000,3],  exp(VBmarginal2.obj$mustore[20000,4] + 0.5*VBmarginal2.obj$Cstore[4,4,20000]^2))
margin3.est <- c(VBmarginal3.obj$mustore[20000,1], exp(VBmarginal3.obj$mustore[20000,2] + 0.5*VBmarginal3.obj$Cstore[2,2,20000]^2), VBmarginal3.obj$mustore[20000,3],  exp(VBmarginal3.obj$mustore[20000,4] + 0.5*VBmarginal3.obj$Cstore[4,4,20000]^2))
margin4.est <- c(VBmarginal4.obj$mustore[20000,1], exp(VBmarginal4.obj$mustore[20000,2] + 0.5*VBmarginal4.obj$Cstore[2,2,20000]^2), VBmarginal4.obj$mustore[20000,3],  exp(VBmarginal4.obj$mustore[20000,4] + 0.5*VBmarginal4.obj$Cstore[4,4,20000]^2))


u1.margin <- pst(y1, xi = margin1.est[1], omega = margin1.est[2], alpha = margin1.est[3], nu = margin1.est[4])
u2.margin <- pst(y2, xi = margin2.est[1], omega = margin2.est[2], alpha = margin2.est[3], nu = margin2.est[4])
u3.margin <- (pst(y3, xi = margin3.est[1], omega = margin3.est[2], alpha = margin3.est[3], nu = margin3.est[4]) - pst(0, xi = margin3.est[1], omega = margin3.est[2], alpha = margin3.est[3], nu = margin3.est[4]))/(1 - pst(0, xi = margin3.est[1], omega = margin3.est[2], alpha = margin3.est[3], nu = margin3.est[4]))
u4.margin <- (pst(y4, xi = margin4.est[1], omega = margin4.est[2], alpha = margin4.est[3], nu = margin4.est[4]) - pst(0, xi = margin4.est[1], omega = margin4.est[2], alpha = margin4.est[3], nu = margin4.est[4]))/(1 - pst(0, xi = margin4.est[1], omega = margin4.est[2], alpha = margin4.est[3], nu = margin4.est[4]))

OptimTry <- optim(par = c(expandedPhi,rep(-2,5)), fn = NewCopulaObjectiveFunctionCut, method = "L-BFGS-B", lower = c(rep(-6,16),rep(-6,16)), upper = c(rep(6,16),rep(4,16)), control = list(fnscale=-1), hessian = TRUE, U1=u1.margin, U2=u2.margin, U3 = u3.margin, U4 = u4.margin )
V <- diag(c(rep(0,16),sqrt(abs(-diag(solve(OptimTry$hessian))))))
V[1:16,1:16] <- C.init.margin
VB.obj <- NewCopulaVBCut(y1 = y1, y2 = y2, y3 = y3, y4 = y4, mu.init = c(VBmarginal1.obj$mustore[20000,],VBmarginal2.obj$mustore[20000,],VBmarginal3.obj$mustore[20000,], VBmarginal4.obj$mustore[20000,],OptimTry$par), C.init = V, maxRuns = 10000)