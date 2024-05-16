#define ARMA_USE_OPENMP

#include<cmath>
#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/expint.hpp>


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]


using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat ConvertMarginalCorCpp(const arma::mat& Omegamat){
  
  const int d = Omegamat.n_rows;
  const arma::vec a = Omegamat.diag(-1); 
  arma::mat R = eye(d,d);
  R.diag(1) = a;
  R.diag(-1) = a;
  arma::vec r1,r3;
  arma::mat R2;
  arma::rowvec r1primeR2Inv, r3primeR2Inv;
  double pieVar, r1primeR2Invr1, r1primeR2Invr3, r3primeR2Invr3;
  for(int k = 2; k<(d+1); k++){
    for(int j = 1; j<(d-k+1);j++){
      
      r1 = R.submat(j,(j-1),(j+k-2),(j-1));
      r3 = (R.submat((j+k-1),j,(j+k-1),(j+k-2))).t();
      R2 = R.submat( j, j, (j+k-2), (j+k-2) );
      pieVar = Omegamat((j+k-1),(j-1));
      r1primeR2Inv = r1.t()*R2.i();
      r3primeR2Inv = r3.t()*R2.i();
      r1primeR2Invr1 = (r1primeR2Inv*r1).eval()(0,0);
      r1primeR2Invr3 = (r1primeR2Inv*r3).eval()(0,0);
      r3primeR2Invr3 = (r3primeR2Inv*r3).eval()(0,0);
      R((j+k-1),(j-1)) = r1primeR2Invr3 + pieVar*sqrt( (1.0 - r1primeR2Invr1)*(1.0 - r3primeR2Invr3) );
      R((j-1),(j+k-1)) = R((j+k-1),(j-1));

      
    }
    
    
  }
  
  return R;
  
  
}


// [[Rcpp::export]]
arma::mat ConvertSemiPartialCorCpp(const arma::mat& Rmat){
  
  const int P = Rmat.n_rows;
  arma::mat SPcorMat(P,P,fill::zeros);
  arma::vec r1;
  arma::vec r3;
  arma::mat R2inv;
  double temp1;
  double temp2;
  double temp3;
  for(int j =0; j < (P-1); j++){
    
    SPcorMat(j, (j+1)) = Rmat(j,(j+1));
    SPcorMat((j+1), j) = SPcorMat(j, (j+1));
    SPcorMat(j,j) = 1.0;
    
  }
  for(int k = 2; k <= (P-1); k++){
    for(int j = 1; j <= (P-k); j++){
      
      //cout << "j= " << j <<  " k = " << k << endl;
      r1 = (Rmat.submat(j-1,j,j-1,j+k-2)).t();
      r3 = (Rmat.submat(j+k-1,j,j+k-1,j+k-2)).t();
      R2inv = (Rmat.submat(j, j, j+k-2, j+k-2)).i();
      temp1 = (r1.t() * R2inv * r3).eval()(0,0);
      temp2 = (r1.t() * R2inv * r1).eval()(0,0);
      temp3 = (r3.t() * R2inv * r3).eval()(0,0);
      SPcorMat(j-1,j+k-1) = (Rmat(j-1,j+k-1) - temp1)/sqrt( (1.0 - temp2)*( 1.0 - temp3 )  );
      SPcorMat(j+k-1, j-1) = SPcorMat(j-1, j+k-1);
      if((j==1)&(k==6)){
        
        cout << temp3  << endl;
        
      }
      
    }
    
    
  }
  SPcorMat(P-1,P-1) = 1.0;
  
  return SPcorMat;
  
}

// [[Rcpp::export]]
arma::mat PartialU(const arma::mat& BigU, const int& p, const arma::cube& PhiArray){
  
  const int m = BigU.n_cols;
  const int bigT = BigU.n_rows;
  const int N = m*bigT;
  int j,s,t,k,l1,l2;
  double num1, num2, utilde, vtilde;
  
  arma::mat OutputU(N,N,fill::zeros);
  int iprime = 0;
  for(int row_count = 0; row_count < BigU.n_rows; row_count++){
    for(int col_count = 0;  col_count< m; col_count++){
      
      OutputU(iprime,iprime) = max(min(BigU(row_count,col_count),0.999999999), 0.000000001);
      //OutputU(iprime,iprime) = BigU(row_count,col_count);
      iprime++;
      
    }
    
  }
  
  for(int r = 1; r < ((p+1)*m); r++){
    for(int i = (r+1); i <= N; i++){
      
      j = i - r;
      s = (int)(ceil( ( j/(double)m) ));
      t = (int)(ceil( ( i/(double)m) ));
      k = t - s;
      l1 = i - m*(t-1);
      l2 = j - m*(s-1);
      
      if(k>p){
        
        OutputU((i-1),(j-1)) = OutputU((i-1),j);
        OutputU((j-1),(i-1)) = OutputU((j-1),(i-2));
        
      }
      else{
        
        //utilde = R::qnorm( OutputU( (i-1),j), 0, 1, 1, 0);
        utilde = max(min( (R::qnorm( OutputU( (i-1),j), 0, 1, 1, 0)  ) , 6.000), -6.000);
        //vtilde = R::qnorm( OutputU( (j-1), (i-2) ), 0, 1, 1, 0);
        vtilde = max(min( (R::qnorm( OutputU( (j-1), (i-2) ), 0, 1, 1, 0)  ) , 6.000), -6.000);
        num1 = (utilde -  PhiArray((l1-1),(l2-1),k)*vtilde )/sqrt(1 - pow( PhiArray((l1-1),(l2-1),k), 2));
        num2 = (vtilde -  PhiArray((l1-1),(l2-1),k)*utilde )/sqrt(1 - pow( PhiArray((l1-1),(l2-1),k), 2));
        num1 = max(min(num1, 6.00),-6.00);
        num2 = max(min(num2, 6.00),-6.00);
        OutputU((i-1),(j-1)) = R::pnorm(num1, 0, 1, 1, 0);
        OutputU((j-1),(i-1)) = R::pnorm(num2 , 0, 1, 1, 0);
        
      }
      
    }
    
    
    
  }

  
  
  return OutputU;
  
}


// [[Rcpp::export]]
double Ktt(const double& t_ind, const arma::mat& Ucond, const arma::mat& Phi0, const int& m){
  
  const int at = m*(t_ind-1) + 1;
  const int bt = m*t_ind;
  double utilde, vtilde, rhotemp;
  int l1,l2;
  double ans_curr = 0.0;
    
  for(int i = at; i <= bt; i++){
    for(int j = at; j < i; j++){
      
      l1 = i - m*(t_ind-1);
      l2 = j - m*(t_ind-1);
      //rhotemp = Phi0(l1-1,l2-1);
      rhotemp = max(min(Phi0(l1-1,l2-1), 0.999999999), -0.999999999);
      utilde = R::qnorm( Ucond(i-1,j), 0, 1, 1, 0);
      vtilde = R::qnorm( Ucond(j-1,i-2), 0, 1, 1, 0);
      ans_curr += -0.5*log(1-rhotemp*rhotemp) - 0.5*( rhotemp*rhotemp*( utilde*utilde + vtilde*vtilde ) - 2.0*rhotemp*utilde*vtilde )/(1.0- rhotemp*rhotemp); 
      
    }
    
    
  }
  
  return ans_curr;
  
}

// [[Rcpp::export]]
double Kts(const double& t_ind, const double& s_ind, const arma::mat& Ucond, const arma::cube& PhiArray, const int& m){
  
  const int at = m*(t_ind-1) + 1;
  const int bt = m*t_ind;
  const int as = m*(s_ind-1) + 1;
  const int bs = m*s_ind;
  const int k = t_ind - s_ind;
  double utilde, vtilde, rhotemp;
  int l1,l2;
  double ans_curr = 0.0;
  
  for(int i = at; i <= bt; i++){
    for(int j = as; j <= bs; j++){
      
      l1 = i - m*(t_ind-1);
      l2 = j - m*(s_ind-1);
      //rhotemp = (PhiArray.slice(k))(l1-1,l2-1);
      rhotemp = max(min((PhiArray.slice(k))(l1-1,l2-1), 0.999999999), -0.999999999);
      utilde = R::qnorm( Ucond(i-1,j), 0, 1, 1, 0);
      vtilde = R::qnorm( Ucond(j-1,i-2), 0, 1, 1, 0);
      ans_curr += -0.5*log(1-rhotemp*rhotemp) - 0.5*( rhotemp*rhotemp*( utilde*utilde + vtilde*vtilde ) - 2.0*rhotemp*utilde*vtilde )/(1.0- rhotemp*rhotemp); 

      
    }
    
    
  }
  
  
  return ans_curr;
  
}

// [[Rcpp::export]]
double VineCopLikelihood(const arma::mat& Umat, const arma::cube& PhiArray, const int& p){
  
  const int m = Umat.n_cols;
  const int bigT = Umat.n_rows;
  const arma::mat& UcondMat = PartialU(Umat, p, PhiArray);
  double ans_curr = 0.0;
  const arma::mat Phi0 = PhiArray.slice(0);
  
  for(int tind = 1; tind <= bigT; tind++){
    for(int sind = max(1,tind-4); sind <= tind; sind++){
      
      if(tind==sind){
        
        ans_curr += Ktt(tind, UcondMat, Phi0, m);
        
      }
      else{
        
        ans_curr += Kts(tind, sind, UcondMat, PhiArray, m);
        
      }

      
    }
  }
  
  return ans_curr;
  
  
}

// [[Rcpp::export]]
arma::vec LogScoreMetricUPart(const arma::mat& Umat, const arma::cube& PhiArray, const int& p){
  
  const int bigT = Umat.n_rows;
  double ans_curr;
  const int m = Umat.n_cols;
  const arma::mat& UcondMat = PartialU(Umat, p, PhiArray);
  const arma::mat Phi0 = PhiArray.slice(0);
  arma::vec ansVec((bigT-p), fill::zeros);

  for(int tind = 5; tind <= bigT; tind++){
    
    ans_curr = 0.0;
    for(int sind = (tind-p); sind <= tind; sind++){
      
      if(tind==sind){
        
        ans_curr += Ktt(tind, UcondMat, Phi0, m);
        
      }
      else{
        
        ans_curr += Kts(tind, sind, UcondMat, PhiArray, m);
        
      }
      
      
    }
    ansVec((tind-p-1)) = ans_curr;
  }
  
  return ansVec;
  
  
}

// [[Rcpp::export]]
arma::rowvec GenerateWOneStepAhead(const arma::vec& Whistory, const arma::mat& Omega0, const arma::mat& Omega1, const arma::mat& Omega2, const arma::mat& Omega3, const arma::mat& Omega4){
  
  const arma::mat SigmaC = join_rows(Omega1.t(), Omega2.t(), Omega3.t(),  Omega4.t());
  const arma::mat SigmaB = join_cols(join_rows(Omega0, Omega1.t(), Omega2.t(), Omega3.t()), join_rows(Omega1, Omega0, Omega1.t(), Omega2.t()), join_rows(Omega2, Omega1, Omega0, Omega1.t()), join_rows(Omega3, Omega2, Omega1, Omega0) );
  const arma::mat CondCBinv = SigmaC*SigmaB.i();
  const arma::vec CondmeanVec = CondCBinv*Whistory;
  const arma::mat CondVar = Omega0 - CondCBinv*SigmaC.t();
  arma::arma_rng::set_seed_random();
  const arma::vec Zrand = randn(4);
  const arma::rowvec ZC = Zrand.t()*chol(CondVar);
  const arma::rowvec ansVec = (CondmeanVec.t() + ZC);
  
  return ansVec;
  
}

// [[Rcpp::export]]
arma::vec GenerateWvStepAhead(const arma::vec& WhistoryInitial, const arma::mat& Omega0, const arma::mat& Omega1, const arma::mat& Omega2, const arma::mat& Omega3, const arma::mat& Omega4, const int& v){
  
  const arma::mat SigmaC = join_rows(Omega1.t(), Omega2.t(), Omega3.t(),  Omega4.t());
  const arma::mat SigmaB = join_cols(join_rows(Omega0, Omega1.t(), Omega2.t(), Omega3.t()), join_rows(Omega1, Omega0, Omega1.t(), Omega2.t()), join_rows(Omega2, Omega1, Omega0, Omega1.t()), join_rows(Omega3, Omega2, Omega1, Omega0) );
  const arma::mat CondCBinv = SigmaC*SigmaB.i();
  const arma::mat CholCondVar = chol(Omega0 - CondCBinv*SigmaC.t());
  arma::vec CondmeanVec(4);
  arma::arma_rng::set_seed_random();
  const arma::mat Zrand(v,4,fill::randn);
  const int WhistoryLength = (v+4)*4;
  arma::vec Whistory(WhistoryLength);
  Whistory.rows(4*(v+4)-16,4*(v+4)-1) = WhistoryInitial;
  
  for(int g = 0; g < v; g++){
    
    CondmeanVec = CondCBinv*Whistory.rows(4*(v+4-g)-16,4*(v+4-g)-1);
    Whistory.rows(4*(v+4-g)-20,4*(v+4-g)-17) = (CondmeanVec.t() + Zrand.row(g)*CholCondVar).t();
    
  }
  
  return Whistory.rows(0,4*(v+4)-17);
  
}

// [[Rcpp::export]]
arma::vec GenerateWvStepAheadPredictiveDensity(const arma::vec& WhistoryInitial, const arma::mat& Omega0, const arma::mat& Omega1, const arma::mat& Omega2, const arma::mat& Omega3, const arma::mat& Omega4, const int& v){
  
  const arma::mat SigmaC = join_rows(Omega1.t(), Omega2.t(), Omega3.t(),  Omega4.t());
  //cout << "Inside1" << endl;
  const arma::mat SigmaB = join_cols(join_rows(Omega0, Omega1.t(), Omega2.t(), Omega3.t()), join_rows(Omega1, Omega0, Omega1.t(), Omega2.t()), join_rows(Omega2, Omega1, Omega0, Omega1.t()), join_rows(Omega3, Omega2, Omega1, Omega0) );
  //cout << "Inside2" << endl;
  const arma::mat CondCBinv = SigmaC*SigmaB.i();
  //cout << "Inside3" << endl;
  const arma::mat CholCondVar = chol(Omega0 - CondCBinv*SigmaC.t());
  //cout << "Inside4" << endl;
  arma::vec CondmeanVec(4);
  arma::arma_rng::set_seed_random();
  const arma::mat Zrand(v,4,fill::randn);
  const int WhistoryLength = (v+4)*4;
  arma::vec WhistoryVEC(WhistoryLength,fill::zeros);
  //Whistory.rows(4*(v+4)-16,4*(v+4)-1) = WhistoryInitial;
  
  
    for(int g = 0; g < v; g++){
      
      WhistoryVEC.rows(4*(v+4)-16,4*(v+4)-1) = WhistoryInitial;
      CondmeanVec = CondCBinv*WhistoryVEC.rows(4*(v+4-g)-16,4*(v+4-g)-1);
      WhistoryVEC.rows(4*(v+4-g)-20,4*(v+4-g)-17) = (CondmeanVec.t() + Zrand.row(g)*CholCondVar).t();
      
    }

  
  
  return WhistoryVEC.rows(0,4*(v+4)-17);
  
}

// [[Rcpp::export]]
arma::cube MakePhiArrayRcpp(const arma::vec& x){
  
  const arma::vec xprime = 2*normcdf(x) - ones(x.n_elem);
  arma::cube outputCube(4,4,5,fill::zeros);
  //c(1,xprime[1],xprime[2],xprime[3],xprime[1],1,xprime[4],xprime[5],xprime[2],xprime[4],1,xprime[6],xprime[3],xprime[5],xprime[6],1)
  outputCube.slice(0) = { { 1.0,xprime(0),xprime(1),xprime(2) }, { xprime(0), 1.0, xprime(3), xprime(4) }, {  xprime(1), xprime(3), 1.0, xprime(5) }, {  xprime(2), xprime(4), xprime(5), 1.0 } };
  outputCube.slice(0) = (outputCube.slice(0)).t();
  outputCube.slice(1) = { { xprime(6), xprime(7), xprime(8), xprime(9) }, { xprime(10), xprime(11), xprime(12), xprime(13) }, { xprime(14), xprime(15), xprime(16), xprime(17) }, { xprime(18), xprime(19), xprime(20), xprime(21) } };
  outputCube.slice(1) = (outputCube.slice(1)).t();
  outputCube.slice(2) = { { xprime(22), xprime(23), xprime(24), xprime(25) }, { xprime(26), xprime(27), xprime(28), xprime(29) }, { xprime(30), xprime(31), xprime(32), xprime(33) }, { xprime(34), xprime(35), xprime(36), xprime(37) } };
  outputCube.slice(2) = (outputCube.slice(2)).t();
  outputCube.slice(3) = { { xprime(38), xprime(39), xprime(40), xprime(41) }, { xprime(42), xprime(43), xprime(44), xprime(45) }, { xprime(46), xprime(47), xprime(48), xprime(49) }, { xprime(50), xprime(51), xprime(52), xprime(53) } };
  outputCube.slice(3) = (outputCube.slice(3)).t();
  outputCube.slice(4) = { { xprime(54), xprime(55), xprime(56), xprime(57) }, { xprime(58), xprime(59), xprime(60), xprime(61) }, { xprime(62), xprime(63), xprime(64), xprime(65) }, { xprime(66), xprime(67), xprime(68), xprime(69) } };
  outputCube.slice(4) = (outputCube.slice(4)).t();
  
  return outputCube;
  
}

// [[Rcpp::export]]
arma::mat FormPhiBlock(const arma::cube& PhiArray){
  
  arma::mat PhiDrawBlock(20,20);
  PhiDrawBlock.submat(0,4,3,7) = PhiArray.slice(1);
  PhiDrawBlock.submat(4,8,7,11) = PhiArray.slice(1);
  PhiDrawBlock.submat(8,12,11,15) = PhiArray.slice(1);
  PhiDrawBlock.submat(12,16,15,19) = PhiArray.slice(1);
  PhiDrawBlock.submat(4,0,7,3) = (PhiArray.slice(1)).t();
  PhiDrawBlock.submat(8,4,11,7) = (PhiArray.slice(1)).t();
  PhiDrawBlock.submat(12,8,15,11) = (PhiArray.slice(1)).t();
  PhiDrawBlock.submat(16,12,19,15) = (PhiArray.slice(1)).t();
  PhiDrawBlock.submat(0,8,3,11) = PhiArray.slice(2);
  PhiDrawBlock.submat(4,12,7,15) = PhiArray.slice(2);
  PhiDrawBlock.submat(8,16,11,19) = PhiArray.slice(2);
  PhiDrawBlock.submat(8,0,11,3) = (PhiArray.slice(2)).t();
  PhiDrawBlock.submat(12,4,15,7) = (PhiArray.slice(2)).t();
  PhiDrawBlock.submat(16,8,19,11) = (PhiArray.slice(2)).t();
  PhiDrawBlock.submat(0,12,3,15) = PhiArray.slice(3);
  PhiDrawBlock.submat(4,16,7,19) = PhiArray.slice(3);
  PhiDrawBlock.submat(12,0,15,3) = (PhiArray.slice(3)).t();
  PhiDrawBlock.submat(16,4,19,7) = (PhiArray.slice(3)).t();
  PhiDrawBlock.submat(0,16,3,19) = PhiArray.slice(4);
  PhiDrawBlock.submat(16,0,19,3) = (PhiArray.slice(4)).t();
  PhiDrawBlock.submat(0,0,3,3) = PhiArray.slice(0);
  PhiDrawBlock.submat(4,4,7,7) = PhiArray.slice(0);
  PhiDrawBlock.submat(8,8,11,11) = PhiArray.slice(0);
  PhiDrawBlock.submat(12,12,15,15) = PhiArray.slice(0);
  PhiDrawBlock.submat(16,16,19,19) = PhiArray.slice(0);
  
  
  return PhiDrawBlock;
  
}

// [[Rcpp::export]]
arma::mat PredictiveDrawRun(const arma::mat& ParDraws, const arma::mat& WhistoryMAT){
  
  arma::mat DrawsMAT(32,WhistoryMAT.n_rows,fill::zeros);
  arma::cube tempPhiArray(4,4,5);
  arma::vec ParCurr(70);
  arma::mat tempPhiBlock(20,20);
  arma::mat tempCorBlock(20,20);
  //cout << "here1" << endl;
  arma::vec tempWHistory(16);
  for(int l = 0; l < WhistoryMAT.n_rows; l++){
    
    ParCurr = (ParDraws.row(l)).t();
    tempPhiArray = MakePhiArrayRcpp(ParCurr);
    //cout << "here2" << endl;
    tempPhiBlock = FormPhiBlock(tempPhiArray);
    //cout << "here3" << endl;
    tempCorBlock = ConvertMarginalCorCpp(tempPhiBlock);
    //cout << "here4" << endl;
    tempWHistory = (WhistoryMAT.row(l)).t();
    //tempWHistory.rows(32,47) takes entries 33 to 48 on temp history (the most recent four quarters of observed W's at any given tprime) to forecast the next 8 quarters
    DrawsMAT.col(l) = GenerateWvStepAheadPredictiveDensity(tempWHistory.rows(32,47), tempCorBlock.submat(0,0,3,3),  tempCorBlock.submat(0,4,3,7), tempCorBlock.submat(0,8,3,11),tempCorBlock.submat(0,12,3,15), tempCorBlock.submat(0,16,3,19), 8);
    
  }
  return DrawsMAT;
  
}

//vcopulaPar = (vphi, log tau0^2, ...., log tau4^2)^top
// [[Rcpp::export]]
arma::rowvec GradLogQ(const arma::mat& UMAT, const arma::mat& Amat, const arma::mat& Bmat, const arma::vec& vcopulaPar, const arma::vec& vmu, const arma::mat& mC, const arma::mat& veta, const arma::mat& vwSq){
  

  const int bigT = UMAT.n_rows;
  const int J = UMAT.n_cols;
  const int d = mC.n_rows;
  const int dReduced = d*(d+1)/2;
  const arma::mat mCmC = mC*mC.t();
  const arma::rowvec gradMu = (vcopulaPar - vmu).t() * mCmC.i();
  double tempZ;
  arma::rowvec gradEta((bigT*J));
  arma::rowvec gradw((bigT*J));
  int counterPoint = 0;
  //const arma::mat Zmat = stats::qnorm((UMAT - Amat)/(Bmat - Amat), 0.0, 1.0 );
  for(int t = 0; t < bigT; t++){
    for(int j = 0; j < J; j++){
      
      tempZ = R::qnorm( ((UMAT(t,j) - Amat(t,j))/(Bmat(t,j) - Amat(t,j))),0.0,1.0, true, false);
      gradEta(counterPoint) = ( tempZ - veta(t,j) )*exp(-vwSq(t,j));
      gradw(counterPoint) = 0.5*pow( (tempZ - veta(t,j)), 2.0)*exp(-vwSq(t,j)) - 0.5;
      counterPoint++;
    }
  }
  
  arma::mat DG(mC.n_rows,mC.n_cols,fill::zeros);
  DG.diag(0) = 1.0/mC.diag(0);
  arma::rowvec DGvec = vectorise(DG, 1);
  arma::rowvec gradvecC = (vcopulaPar - vmu).t()*(mC.i()).t()* kron((vcopulaPar - vmu).t(), eye( d, d ))*(  kron((mC.i()).t(), mC.i()) ) - DGvec;
  arma::rowvec gradvecCReduced(dReduced, fill::zeros);
  int tcounter = 0;
  int tReducedcounter = 0;
  for(int t1 = 0; t1 < d; t1++){
    for(int t2 = 0; t2 < d; t2++){
      
      if(t1 <= t2){
        
        gradvecCReduced(tReducedcounter) = gradvecC(tcounter);
        tReducedcounter++;
        
      }
      tcounter++;
      
    }
  }
  arma::rowvec finalGrad = join_rows(gradMu, gradvecCReduced, gradEta, gradw);
  
  //Rcpp::List::create(Rcpp::Named("m_mu1") = m_mu1, Rcpp::Named("m_mu0") = m_mu0, Rcpp::Named("m_mu") = m_mu, Rcpp::Named("Sigma_Mu1") = Sigma_Mu1, Rcpp::Named("Sigma_Mu0") = Sigma_Mu0, Rcpp::Named("Sigma_Mu") = Sigma_Mu); 
  //return Rcpp::List::create(Rcpp::Named("finalGrad") = finalGrad, Rcpp::Named("gradvecC") = gradvecC);
  return finalGrad;
  
  
}

//vcopulaPar = (vphi, log tau0^2, ...., log tau4^2)^top
// [[Rcpp::export]]
double LogQEval(const arma::mat& UMAT, const arma::mat& Amat, const arma::mat& Bmat, const arma::vec& vcopulaPar, const arma::vec& vmu, const arma::mat& mC, const arma::mat& veta, const arma::mat& vwSq){
  
  
  const int bigT = UMAT.n_rows;
  const int J = UMAT.n_cols;
  const arma::mat mCmC = mC*mC.t();

  double tempZ;
  //const arma::mat Zmat = stats::qnorm((UMAT - Amat)/(Bmat - Amat), 0.0, 1.0 );
  double ansFinal = 0.0;
  for(int t = 0; t < bigT; t++){
    for(int j = 0; j < J; j++){
      
      tempZ = R::qnorm( ((UMAT(t,j) - Amat(t,j))/(Bmat(t,j) - Amat(t,j))),0.0,1.0, true, false);
      //ansFinal = ansFinal + log_normpdf( (tempZ - veta(t,j))/sqrt( exp(vwSq(t,j)) ) ) - log(Bmat(t,j) - Amat(t,j)) - log_normpdf(tempZ);
      ansFinal = ansFinal - 0.5*pow((tempZ - veta(t,j)),2.0)*exp(-vwSq(t,j)) - 0.5*vwSq(t,j) - log(Bmat(t,j) - Amat(t,j)) +0.5*tempZ*tempZ;
    }
  }

  arma::mat quadterm = ( vcopulaPar - vmu ).t()*mCmC.i()*( vcopulaPar - vmu );
  double val = max(det(mC),0.0000000001);

  ansFinal = ansFinal - log(val) - 0.5*vmu.n_elem*log(2*3.1416) - 0.5*quadterm(0,0);
  //return Rcpp::List::create(Rcpp::Named("ansFinal") = ansFinal, Rcpp::Named("quadterm") = quadterm, Rcpp::Named("val") = val, Rcpp::Named("tempZ") = tempZ); 
  return ansFinal;
  
  
}

// [[Rcpp::export]]
double NewCopulaObjectiveFunctionReverseCutCpp(const arma::vec& x, const arma::vec& U1, const arma::vec& U2, const arma::vec& U3, const arma::vec& U4){
  
  const arma::mat& UmatLocal = join_rows(U1, U2, U3, U4);
  const arma::cube& PhiArrayLocal = MakePhiArrayRcpp(x.rows(0,69));
  const double& LikVal = VineCopLikelihood(UmatLocal, PhiArrayLocal, 4);
  const double& PriorVal = -0.5*6.0*log(2*3.1416*exp(x(70))) - 0.5*accu(pow(x.rows(0,5),2.0))/exp(x(70)) -0.5*16.0*log(2*3.1416*exp(x(71))) - 0.5*accu(pow(x.rows(6,21),2.0))/exp(x(71)) -0.5*16.0*log(2*3.1416*exp(x(72))) - 0.5*accu(pow(x.rows(22,37),2.0))/exp(x(72)) -0.5*16.0*log(2*3.1416*exp(x(73))) - 0.5*accu(pow(x.rows(38,53),2.0))/exp(x(73)) - 0.5*16.0*log(2*3.1416*exp(x(74))) - 0.5*accu(pow(x.rows(54,69),2.0))/exp(x(74));
  const double& HyperpriorVal = 5.0*log(2.0/3.1416) - 5.0*0.5*log(1) - accu(x.rows(70,74)) - accu(log(1.0 + exp(2.0*x.rows(70,74))));
  double answer = LikVal+PriorVal+HyperpriorVal;
  return answer;
  
}

// [[Rcpp::export]]
arma::mat DrawQuCpp(const arma::mat& AMAT, const arma::mat& BMAT, const arma::mat& etamat, const arma::mat& wmat){
  
  const int TSlength = AMAT.n_rows;
  arma::arma_rng::set_seed_random();
  arma::mat zDrawMat(TSlength, 4, fill::randn);
  arma::mat uDrawMat(TSlength, 4);
  double anstemp;
  for(int tind = 0; tind < TSlength; tind++){
    for(int j = 0; j < 4; j++){
      
      anstemp = normcdf( (zDrawMat(tind, j)*exp(0.5*wmat(tind,j)) + etamat(tind,j)) );
      anstemp = max(min(anstemp,0.99),0.01);
      uDrawMat(tind,j) = anstemp*(BMAT(tind,j) - AMAT(tind,j)) + AMAT(tind,j);
      
    }
  }
  
  return uDrawMat;
  
}

// [[Rcpp::export]]
arma::vec DrawQPhiLogTauCpp(const arma::vec& vmu, const arma::mat& mC){
  
  const int drawLength = vmu.n_elem;
  arma::arma_rng::set_seed_random();
  const arma::vec drawsVec(drawLength,fill::randn);
  const arma::vec finaldraw = vmu + mC*drawsVec;
  
  return finaldraw;
  
}


// [[Rcpp::export]]
arma::mat invvechCpp(const arma::vec& vx){
  
  const int R = vx.n_elem;
  const int d = (sqrt(1+8*R)-1)/2;
  int vectorCounter = 0;
  arma::mat ansMAT(d,d,fill::zeros);
  for(int ind1 = 0; ind1 < d; ind1++){
    for(int ind2 = ind1; ind2 < d; ind2++){
    
      ansMAT(ind2,ind1) = vx(vectorCounter);
      vectorCounter++;
    
    }
  }
  
  return ansMAT;
}

// [[Rcpp::export]]
arma::vec EstimateXiCpp(const arma::mat& AMAT, const arma::mat& BMAT, const arma::mat& etamat, const arma::mat& wmat, const arma::vec& vmu, const arma::mat& mC, const int numDraws){

  const int d = vmu.n_elem;
  const int TSlength = AMAT.n_rows;
  const int numPars = d + d*(d+1)/2 + 2*etamat.n_elem;
  arma::vec E1(numPars, fill::zeros);
  arma::vec E2(numPars, fill::zeros);
  arma::vec E3(numPars, fill::zeros);
  arma::vec E4(numPars, fill::zeros);
  arma::mat UdrawCurr(AMAT.n_rows, AMAT.n_cols);
  arma::vec DrawTauPhiCurr(d);
  arma::arma_rng::set_seed_random();
  arma::mat Cz(vmu.n_elem ,numDraws, fill::randn);
  arma::mat zDrawMat(TSlength, 4, fill::randn);
  Cz = mC*Cz;
  double logHcurr, logQcurr;
  arma::vec GradLogQVec(numPars);
  arma::vec GradLogQVecSQ(numPars);
  arma::vec E1est(numPars,fill::zeros);
  arma::vec E2est(numPars,fill::zeros);
  arma::vec E3est(numPars,fill::zeros);
  arma::vec E4est(numPars,fill::zeros);
  double anstemp;

  for(int b = 0; b < numDraws; b++){

    //UdrawCurr = DrawQuCpp(AMAT, BMAT, etamat, wmat);
    DrawTauPhiCurr = vmu + Cz.col(b);
    anstemp = 0.0;
    for(int tind = 0; tind < TSlength; tind++){
      for(int j = 0; j < 4; j++){

        anstemp = normcdf( (zDrawMat(tind, j)*exp(0.5*wmat(tind,j)) + etamat(tind,j)) );
        anstemp = max(min(anstemp,0.99),0.01);
        UdrawCurr(tind,j) = anstemp*(BMAT(tind,j) - AMAT(tind,j)) + AMAT(tind,j);

      }
    }

    logHcurr = NewCopulaObjectiveFunctionReverseCutCpp(DrawTauPhiCurr, UdrawCurr.col(0), UdrawCurr.col(1), UdrawCurr.col(2), UdrawCurr.col(3));
    logQcurr = LogQEval(UdrawCurr, AMAT, BMAT, DrawTauPhiCurr, vmu, mC, etamat, wmat);
    GradLogQVec = (GradLogQ(UdrawCurr, AMAT, BMAT, DrawTauPhiCurr, vmu, mC, etamat, wmat)).t();
    GradLogQVecSQ = pow(GradLogQVec,2.0);
    E1est = E1est + (logHcurr - logQcurr)*GradLogQVecSQ;
    E2est = E2est + (logHcurr - logQcurr)*GradLogQVec;
    E3est = E3est + GradLogQVec;
    E4est = E4est + GradLogQVecSQ;

  }

  arma::vec denom = (E4est/numDraws) - pow( (E3est/numDraws), 2.0);
  arma::vec numer = (E1est/numDraws) - (E2est/numDraws)%(E3est/numDraws);
  int sgntemp;
  for(int iddenom = 0; iddenom < numPars; iddenom++){

    sgntemp = sign(denom(iddenom));
    if(sgntemp==0){ sgntemp = 1; }
    denom(iddenom) = sgntemp*max(abs(denom(iddenom)),0.0000001);

  }
  arma::vec vxi = numer/denom;

  return(vxi);

}

// [[Rcpp::export]]
arma::mat RearrangeVecMat(const arma::vec& vv, const int& numCols){
  
  const int numRows = vv.n_elem/numCols;
  arma::mat RetMAT(numRows,numCols);
  for(int j = 0; j < numRows; j++){
    
    RetMAT.row(j) = (vv.rows( numCols*j, (numCols*j+(numCols-1)) )).t();
    
  }
  return RetMAT;
  
}

//function(AMAT, BMAT, mu.init, C.init, etamat.init, wmat.init, maxRuns = 50000)
// [[Rcpp::export]]
List NewCopulaVBReverseCut(const arma::mat& AMAT, const arma::mat& BMAT, const arma::vec& muInit, const arma::mat& CInit, const arma::mat& EtaMatInit, const arma::mat& WMatInit, const int maxRuns){

  const int d = muInit.n_elem;
  const int S = 100;
  const int LengthTS = EtaMatInit.n_rows;
  const int numPar = d + d*(d+1)/2 + 8*LengthTS;

  arma::vec muCurr = muInit;
  arma::mat CCurr = CInit;
  const double rho = 0.85;
  const double eps = 1.0/1000000.0;

  arma::vec AccugradMu(d,fill::zeros);
  arma::vec AccuchangeMu(d,fill::zeros);
  arma::mat AccugradC(d,d,fill::zeros);
  arma::mat AccuchangeC(d,d,fill::zeros);
  arma::mat AccugradEta(LengthTS,4,fill::zeros);
  arma::mat AccuchangeEta(LengthTS,4,fill::zeros);
  arma::mat AccugradW(LengthTS,4,fill::zeros);
  arma::mat AccuchangeW(LengthTS,4,fill::zeros);

  arma::mat MuStore(maxRuns,d);
  arma::cube CStore(d,d,maxRuns);
  arma::cube EtaStore(LengthTS,4,maxRuns);
  arma::cube WStore(LengthTS,4,maxRuns);

  arma::mat EtaMatCurr = EtaMatInit;
  arma::mat WMatCurr = WMatInit;
  arma::vec XiEstCurr(numPar);
  arma::vec GradLogQVec(numPar);
  arma::vec gEstCurr(numPar, fill::zeros);

  arma::mat UdrawCurr(AMAT.n_rows, AMAT.n_cols);
  arma::vec DrawTauPhiCurr(d);

  arma::mat Cz(d ,S);
  arma::mat zDrawMat(LengthTS, 4);
  double anstemp, logHcurr, logQcurr;

  arma::vec gradMu(d);
  arma::vec ChangeMu(d);
  arma::vec muNew(d);

  arma::mat gradC(d,d);
  arma::mat ChangeC(d,d);
  arma::mat CNew(d,d);

  arma::mat gradEta(LengthTS,4);
  arma::mat ChangeEta(LengthTS,4);
  arma::mat EtaNew(LengthTS,4);

  arma::mat gradW(LengthTS,4);
  arma::mat ChangeW(LengthTS,4);
  arma::mat WNew(LengthTS,4);

  //const arma::mat EpsMATVech =  eps*invvechCpp(ones( (d + d*(d+1)/2) ));
  //const arma::vec EpsVec =  eps*ones(d);
  //const arma::mat EpsMATLong =  eps*ones(LengthTS,4);
  //cout << "Checkpoint1 cleared" << endl;
  for(int t = 0; t < maxRuns; t++){

    XiEstCurr = EstimateXiCpp(AMAT, BMAT, EtaMatCurr, WMatCurr, muCurr, CCurr, 100);
    arma::arma_rng::set_seed_random();
    Cz.randn();
    zDrawMat.randn();
    gEstCurr.zeros();
    //cout << "Checkpoint2 cleared" << endl;
    for(int s = 0; s < S; s++){

      DrawTauPhiCurr = muCurr + Cz.col(s);
      anstemp = 0.0;
      for(int tind = 0; tind < LengthTS; tind++){
        for(int j = 0; j < 4; j++){

          anstemp = normcdf( (zDrawMat(tind, j)*exp(0.5*WMatCurr(tind,j)) + EtaMatCurr(tind,j)) );
          anstemp = max(min(anstemp,0.99),0.01);
          UdrawCurr(tind,j) = anstemp*(BMAT(tind,j) - AMAT(tind,j)) + AMAT(tind,j);

        }
      }
      logHcurr = NewCopulaObjectiveFunctionReverseCutCpp(DrawTauPhiCurr, UdrawCurr.col(0), UdrawCurr.col(1), UdrawCurr.col(2), UdrawCurr.col(3));
      logQcurr = LogQEval(UdrawCurr, AMAT, BMAT, DrawTauPhiCurr, muCurr, CCurr, EtaMatCurr, WMatCurr);
      GradLogQVec = (GradLogQ(UdrawCurr, AMAT, BMAT, DrawTauPhiCurr, muCurr, CCurr, EtaMatCurr, WMatCurr)).t();
      gEstCurr = gEstCurr + (logHcurr - logQcurr - XiEstCurr)%GradLogQVec;

    }
    gEstCurr = gEstCurr/S;
    //cout << "Checkpoint3 cleared" << endl;
    gradMu = gEstCurr.rows(0,d-1);
    AccugradMu = rho*AccugradMu + (1-rho)*pow(gradMu, 2);
    ChangeMu = (sqrt( AccuchangeMu + eps )/sqrt( AccugradMu + eps ))%gradMu;
    muNew = muCurr + ChangeMu;
    AccuchangeMu = rho*AccuchangeMu + (1-rho)*pow(ChangeMu,2);
    //cout << "Checkpoint4 cleared" << endl;
    gradC = invvechCpp(gEstCurr.rows(d,d+d*(d+1)/2-1));
    AccugradC = rho*AccugradC + (1-rho)*pow(gradC,2);
    ChangeC = (sqrt( AccuchangeC + eps )/sqrt( AccugradC + eps ))%gradC;
    CNew = CCurr + ChangeC;
    AccuchangeC = rho*AccuchangeC + (1-rho)*pow(ChangeC,2);
    //cout << "Checkpoint5 cleared" << endl;
    gradEta = RearrangeVecMat( gEstCurr.rows(d+d*(d+1)/2, (d+d*(d+1)/2 + 4*LengthTS-1)), 4);
    AccugradEta = rho*AccugradEta + (1-rho)*pow(gradEta, 2);
    ChangeEta = (sqrt( AccuchangeEta + eps )/sqrt( AccugradEta + eps ))%gradEta;
    EtaNew = EtaMatCurr + ChangeEta;
    AccuchangeEta = rho*AccuchangeEta + (1-rho)*pow(ChangeEta,2);
    //cout << "Checkpoint6 cleared" << endl;
    gradW = RearrangeVecMat( gEstCurr.rows((d+d*(d+1)/2 + 4*LengthTS), (d+d*(d+1)/2 + 8*LengthTS-1)), 4);
    AccugradW = rho*AccugradW + (1-rho)*pow(gradW, 2);
    ChangeW = (sqrt( AccuchangeW + eps )/sqrt( AccugradW + eps ))%gradW;
    WNew = WMatCurr + ChangeW;
    AccuchangeW = rho*AccuchangeW + (1-rho)*pow(ChangeW,2);
    //cout << "Checkpoint7 cleared" << endl;
    // for(int dind = 0; dind < d; dind ++){
    //
    //   gradMu(dind) = gEstCurr(dind);
    //   AccugradMu(dind) = rho*AccugradMu(dind) + (1-rho)*gradMu(dind)*gradMu(dind);
    //   ChangeMu(dind) = sqrt(AccuchangeMu(dind) + eps)/sqrt(AccugradMu(dind) + eps)*gradMu(dind);
    //   muNew(dind) = muCurr(dind) + ChangeMu(dind);
    //   AccuchangeMu(dind) = rho*AccuchangeMu(dind) + (1-rho)*ChangeMu(dind)*ChangeMu(dind);
    //
    //   gradC = invvechCpp(gEstCurr.rows(d,d+d*(d+1)/2));
    //   for(int dind2 = dind; dind2 < d; dind2++){
    //
    //     AccugradC(dind2, dind) = rho*AccugradC(dind2, dind) + (1-rho)*gradC(dind2, dind)*gradC(dind2, dind);
    //     ChangeC(dind2, dind) = sqrt( AccuchangeC(dind2, dind) + eps )/sqrt( AccugradC(dind2, dind) + eps)*gradC(dind2, dind);
    //     CNew(dind2, dind) =  CCurr(dind2, dind) + ChangeC(dind2, dind);
    //     AccuchangeC(dind2, dind) = rho*AccuchangeC(dind2, dind) + (1-rho)*ChangeC(dind2, dind)*ChangeC(dind2, dind);
    //
    //   }
    //
    // }


    muCurr = muNew;
    CCurr = CNew;
    EtaMatCurr = EtaNew;
    WMatCurr = WNew;
    MuStore.row(t) = muCurr.t();
    CStore.slice(t) = CCurr;
    EtaStore.slice(t) = EtaMatCurr;
    WStore.slice(t) = WMatCurr;
    if(isnan(accu(muCurr))){

      break;

    }
    //cout << "Checkpoint8 cleared" << endl;
    if(( (t+1)%10)==0){

      cout << "Completed iter " << t+1 << endl;

    }

  }

  //return return Rcpp::List::create(Rcpp::Named("finalGrad") = finalGrad, Rcpp::Named("gradvecC") = gradvecC);
  return Rcpp::List::create(Rcpp::Named("MuStore") = MuStore, Rcpp::Named("CStore") = CStore, Rcpp::Named("EtaStore") = EtaStore, Rcpp::Named("WStore") = WStore);

}