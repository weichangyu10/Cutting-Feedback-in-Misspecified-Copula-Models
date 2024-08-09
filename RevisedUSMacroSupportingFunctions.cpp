#include<cmath>
#include <RcppArmadillo.h>
#include <RcppNumerical.h>


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]


using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace Numer;

// [[Rcpp::export]]
double dLogExpHalfCauchyRcpp(const arma::vec& x, const double S_sq){
  
  const double ans = -arma::sum(log(0.5*M_PI*sqrt(S_sq)) + log(exp(-x) + exp(x)/S_sq));
  
  return ans;
  //returns the sum of the log-density values instead of a vector of log-densities
  
  
}

// [[Rcpp::export]]
double dLogNormRcpp(const arma::vec& x, const double S){
  
  const int xLong = x.n_elem;
  const double ans = -0.5*xLong*log(2*M_PI*S*S) - 0.5*pow(norm(x,2),2)/(S*S);
  
  return ans;
  //returns the sum of the log-density values instead of a vector of log-densities
  
  
}

// [[Rcpp::export]]
double dLogLaplaceRcpp(const arma::vec& x, const double S){
  
  const double ans = - arma::sum( log(2*S) + abs(x)/S );
  
  return ans;
  //returns the sum of the log-density values instead of a vector of log-densities
  
  
}

// [[Rcpp::export]]
double dLogExpExpRcpp(const arma::vec& x, const double Lambda_Rate){
  
  const double ans = arma::sum( log(Lambda_Rate) + x - Lambda_Rate*exp(x) );
  return ans;
}

// [[Rcpp::export]]
arma::vec ApproxQnormRcpp( const arma::vec& vp  ){
  
  const Rcpp::NumericVector x = Rcpp::wrap(vp);
  const NumericVector QnormOutput = Rcpp::qnorm(x,  0.0, 1.0, true, false);
  arma::vec QnormOutputARMA = as<arma::vec>(QnormOutput);
  const arma::uvec Extreme = arma::find(vp > 0.9999683 || vp < 0.00003167);
  for (arma::uword i = 0; i < Extreme.n_elem; ++i) {
    
    if(vp(Extreme(i)) > 0.9999683){
      
      QnormOutputARMA(Extreme(i)) = sqrt( -2*log(1.0 - vp(Extreme(i))) - log(-2*log(1.0 - vp(Extreme(i)))) - log(2*M_PI)  );
      
    }
    else{
      
      QnormOutputARMA(Extreme(i)) = -sqrt( -2*log(vp(Extreme(i))) - log(-2*log(vp(Extreme(i)))) - log(2*M_PI)  );
      
    }
    
  }
  
  return QnormOutputARMA;
  
}

// [[Rcpp::export]]
arma::vec ApproxPnormRcpp( const arma::vec& vq  ){
  
  arma::vec PnormOutputARMA = normcdf(vq);
  return PnormOutputARMA;
  
}

class SkewTCDFIntegrand: public Func
{
private:
  const double xi;
  const double log_omega;
  const double alpha;
  const double log_nu;
public:
  SkewTCDFIntegrand(double xi_, double log_omega_, double alpha_, double log_nu_) : xi(xi_), log_omega(log_omega_), alpha(alpha_), log_nu(log_nu_) {}
  
  double operator()(const double& x) const
  {
    //const double xi = theta(0);
    // const double omega = exp(theta(1));
    // const double alpha = theta(2);
    // const double nu = exp(theta(3));
    // 
    // const Rcpp::NumericVector Xvec = Rcpp::NumericVector( Rcpp::wrap( ((Yvec - xi)/omega)  ) );
    // const arma::vec ANSvec =  log(2) - log(omega) + Rcpp::dt( Xvec, nu,  true ) + Rcpp::pt( alpha*Xvec*sqrt( (nu+1.0)/(nu+pow(Xvec,2)) ) , (nu+1.0),  true, true );
    // 
    const double& xnorm = (x - xi)*exp(-log_omega);
    const double& nu = exp(log_nu);
    //const double omega = exp(theta(1));
    //const double alpha = theta(2);
    //const double nu = exp(theta(3));
    
    //const Rcpp::NumericVector Xvec = Rcpp::NumericVector( Rcpp::wrap( ((Yvec - xi)/omega)  ) );
    //const arma::vec ANSvec =  log(2) - log(omega) + Rcpp::dt( Xvec, nu,  true ) + Rcpp::pt( alpha*Xvec*sqrt( (nu+1.0)/(nu+pow(Xvec,2)) ) , (nu+1.0),  true, true );
    return exp(log(2) - log_omega + R::dt(  xnorm, nu, true) + R::pt( (alpha*xnorm*sqrt( (nu+1.0)/(nu+xnorm*xnorm ) )), nu+1.0, true, true ));
    //return std::exp(t * x) * R::dbeta(x, a, b, 0);
  }
};

// [[Rcpp::export]]
arma::vec SkewT_cdf(const arma::vec& EvalPoints, double xi, double log_omega, double alpha, double log_nu)
{
  SkewTCDFIntegrand f(xi, log_omega, alpha, log_nu);
  double err_est;
  int err_code;
  const int NumPoints = EvalPoints.n_elem;
  arma::vec ResultsVec(NumPoints);
  for(int i = 0; i < NumPoints; i++){
    
    ResultsVec(i) = integrate(f, -1000.0, EvalPoints(i), err_est, err_code);
    
  }
  
  
  return ResultsVec;
}

// [[Rcpp::export]]
arma::vec TruncatedSkewT_cdf(const arma::vec& EvalPoints, double xi, double log_omega, double alpha, double log_nu)
{
  SkewTCDFIntegrand f(xi, log_omega, alpha, log_nu);
  double err_est;
  int err_code;
  const int NumPoints = EvalPoints.n_elem;
  arma::vec ResultsVec(NumPoints);
  for(int i = 0; i < NumPoints; i++){
    
    ResultsVec(i) = integrate(f, 0.0, EvalPoints(i), err_est, err_code);
    
  }
  double ComplementNormConst = integrate(f, -2000, 0.0, err_est, err_code);
  ResultsVec = ResultsVec/(1.0 - ComplementNormConst);
  
  return ResultsVec;
}

// [[Rcpp::export]]
arma::mat ComputeConditionalUsRcpp( const int M, const int p, const arma::vec& U, const arma::cube& PhiArray ){
  
  const int N = U.n_elem;
  arma::vec QnormU = ApproxQnormRcpp(U);
  arma::mat ConditionalQnormUMAT = diagmat(QnormU);
  arma::mat ConditionalUMAT = diagmat(U);
  int j, s, t, k, l1, l2;
  double QNORM_u_i_jPlusOneTemp, QNORM_u_j_iMinusOneTemp, phiTemp;
  
  for(int r = 1; r <= ((p+1)*M-1); r++){
    for(int i = (r+1); i <= N; i++){
      
      j = i - r;
      s = (j+M-1)/M;
      t = (i+M-1)/M;
      k = t - s;
      l1 = i - M*(t-1);
      l2 = j - M*(s-1);
      
      
      QNORM_u_i_jPlusOneTemp = ConditionalQnormUMAT(i-1,j );
      QNORM_u_j_iMinusOneTemp = ConditionalQnormUMAT(j-1,i-2);
      
      phiTemp = PhiArray(l1-1,l2-1, (min(k,p)) );
      
      if( k > p ){
        
        ConditionalQnormUMAT(i-1,j-1) = ConditionalQnormUMAT(i-1,j);
        ConditionalQnormUMAT(j-1,i-1) = ConditionalQnormUMAT(j-1,i-2);
        
      }
      else{
        
        ConditionalQnormUMAT(i-1,j-1) = (QNORM_u_i_jPlusOneTemp - phiTemp*QNORM_u_j_iMinusOneTemp)/sqrt(1.0 - phiTemp*phiTemp);
        ConditionalQnormUMAT(j-1,i-1) = ( QNORM_u_j_iMinusOneTemp - phiTemp*QNORM_u_i_jPlusOneTemp  )/sqrt(1.0 - phiTemp*phiTemp);
        
        
      }
      ConditionalUMAT(i-1,j-1) = normcdf(ConditionalQnormUMAT(i-1,j-1));
      ConditionalUMAT(j-1,i-1) = normcdf(ConditionalQnormUMAT(j-1,i-1));
      
    }
    
    
  }
  
  
  
  return ConditionalQnormUMAT;
  // return Rcpp::List::create(
  //   Rcpp::Named("ConditionalQnormUMAT") = ConditionalQnormUMAT, Rcpp::Named("ConditionalUMAT") = ConditionalUMAT
  // );
  
  
}

// [[Rcpp::export]]
arma::vec SkewTlogPDFcpp(const arma::vec& theta, const arma::vec& Yvec){
  
  const double xi = theta(0);
  const double omega = exp(theta(1));
  const double alpha = theta(2);
  const double nu = exp(theta(3));
  
  const Rcpp::NumericVector Xvec = Rcpp::NumericVector( Rcpp::wrap( ((Yvec - xi)/omega)  ) );
  const arma::vec ANSvec =  log(2) - log(omega) + Rcpp::dt( Xvec, nu,  true ) + Rcpp::pt( alpha*Xvec*sqrt( (nu+1.0)/(nu+pow(Xvec,2)) ) , (nu+1.0),  true, true );
  
  
  return ANSvec;
}

// [[Rcpp::export]]
double SkewTlogLikcpp(const arma::vec& theta, const arma::vec& Yvec){
  
  const int datLength = Yvec.n_elem;
  const double xi = theta(0);
  const double omega = exp(theta(1));
  const double alpha = theta(2);
  const double nu = exp(theta(3));
  
  const Rcpp::NumericVector Xvec = Rcpp::NumericVector( Rcpp::wrap( ((Yvec - xi)/omega)  ) );
  const double ANSOUT = datLength*log(2) - datLength*log(omega) + sum(Rcpp::dt( Xvec, nu,  true )) + sum(Rcpp::pt( alpha*Xvec*sqrt( (nu+1.0)/(nu+pow(Xvec,2)) ) , (nu+1.0),  true, true ));
  
  
  return ANSOUT;
}




// [[Rcpp::export]]
arma::vec TruncatedSkewTlogPDFcpp(const arma::vec& theta, const arma::vec& Yvec){
  
  const double xi = theta(0);
  const double omega = exp(theta(1));
  const double alpha = theta(2);
  const double nu = exp(theta(3));
  
  const Rcpp::NumericVector Xvec = Rcpp::NumericVector( Rcpp::wrap( ((Yvec - xi)/omega)  ) );
  arma::vec RawProb = SkewT_cdf(zeros(1), theta(0), theta(1), theta(2), theta(3));
  const double normConst = 1.0 - RawProb(0);
  const arma::vec ANSvec =  log(2) - log(omega) - log(normConst) + Rcpp::dt( Xvec, nu,  true ) + Rcpp::pt( alpha*Xvec*sqrt( (nu+1.0)/(nu+pow(Xvec,2)) ) , (nu+1.0),  true, true );
  
  
  return ANSvec;
}

// [[Rcpp::export]]
arma::mat EstimateCorrelationsFromUcpp( const arma::mat& BigU, const int& p){
  
  const int BigT = BigU.n_rows;
  const int m = BigU.n_cols;
  const int N = m*BigT;
  arma::cube CovArray(m,m,p+1);
  CovArray.slice(0) = cov( BigU );
  const arma::vec DiagSDs = 1.0/sqrt(diagvec(CovArray.slice(0)));
  const arma::mat DiagSDmat = DiagSDs * DiagSDs.t();
  arma::cube CorArray(m,m,p+1);
  CorArray.slice(0) = cor( BigU );
  //arma::mat CovMAT(N,N,fill::zeros);
  arma::mat CorMAT(N,N,fill::zeros);
  arma::mat TempCov(2*m,2*m);
  for(int k = 1; k <= p; k++){
    
    TempCov = cov( join_rows( BigU.rows( 0, BigT - k - 1 ), BigU.rows( k, BigT - 1 )  ) );
    CorArray.slice(k) = TempCov.submat(0,m, m-1, 2*m-1) % DiagSDmat;
    
  }
  
  int kTemp, a_t, b_t, a_tprime, b_tprime;
  for(int t = 1; t <= BigT; t++){
    for(int tprime = t; tprime <= min(t+p, BigT); tprime++){
      
      kTemp = tprime - t;
      a_t = (t-1)*m+1;
      b_t = t*m;
      a_tprime = (tprime-1)*m+1;
      b_tprime = tprime*m;
      CorMAT.submat( a_t-1,a_tprime - 1, b_t-1, b_tprime-1 ) = CorArray.slice(kTemp);
      CorMAT.submat( a_tprime - 1, a_t-1, b_tprime-1, b_t - 1) = (CorArray.slice(kTemp)).t();
      
    }
    
  }
  
  
  return CorMAT;
  
}

// [[Rcpp::export]]
arma::cube ConvertCorToPCorcpp(const arma::mat& M, const double& BigT, const int& m, const int& p){
  
  const int D = M.n_cols;
  arma::mat PCorMAT = eye(D, D);
  PCorMAT.diag(1) = M.diag(1);
  PCorMAT.diag(-1) = M.diag(-1);
  arma::rowvec r1;
  arma::rowvec r3;
  arma::mat R2;
  double Denom;
  arma::vec Intermediate1(1);
  arma::vec Intermediate2(1);
  arma::vec Intermediate3(1);
  for(int k = 2; k <= (D-1); k++){
    for(int j = 1; j <= (D-k); j++){
      
      r1 = M.submat( j-1, j, j-1, j+k-2 );
      r3 = M.submat(j+k-1, j, j+k-1, j+k-2);
      R2 = M.submat(j,j,j+k-2,j+k-2);
      Intermediate1 = r1*R2.i()*r1.t();
      Intermediate2 = r3*R2.i()*r3.t();
      Denom = sqrt( (1.0 - Intermediate1(0) )*( 1.0 - Intermediate2(0) ) );
      Intermediate3 = r1*R2.i()*r3.t();
      PCorMAT(j-1, j+k-1) = ( M(j-1, j+k-1) - Intermediate3(0) )/Denom;
      PCorMAT(j+k-1, j-1) = PCorMAT(j-1, j+k-1);
      
    }
    
  }
  
  
  arma::cube PhiArray(m,m,p+1);
  for(int k = 1; k <= (p+1); k++){
    
    PhiArray.slice(k-1) = PCorMAT.submat( 0, (k-1)*m, m-1, k*m-1 );
    
  }
  
  return PhiArray;
  
}

// [[Rcpp::export]]
arma::mat ConvertPCorToCorcpp(const arma::cube& PhiArray, const int& m){
  
  const int p = PhiArray.n_slices-1;
  const int MDIM = (p+1)*m;
  arma::mat PhiMAT(MDIM, MDIM);
  int k;
  
  int a_t, b_t, a_tprime, b_tprime;
  for(int t = 1; t <= (p+1); t++){
    for(int tprime = t;  tprime <= (p+1); tprime++){
      
      a_t = (t-1)*m;
      b_t = t*m-1;
      a_tprime = (tprime-1)*m;
      b_tprime = tprime*m-1;
      k = tprime - t;
      PhiMAT.submat(a_t, a_tprime, b_t, b_tprime) = PhiArray.slice(k);
      PhiMAT.submat(a_tprime, a_t, b_tprime, b_t) = (PhiArray.slice(k)).t();
      
    }
  }
  
  arma::mat R = eye(MDIM,MDIM);
  const arma::vec AdjacentCor = PhiMAT.diag(-1);
  R.diag(1) = AdjacentCor;
  R.diag(-1) = AdjacentCor;
  arma::rowvec r1,r3;
  arma::mat R2;
  double Denom;
  double r1primeR2Invr3prime, r3primeR2Invr3prime, r1primeR2Invr1prime;
  for(int k = 2; k <= (MDIM-1); k++){
    for(int j = 1; j <= (MDIM-k); j++){
      
      r1 = R.submat( j-1, j, j-1, j+k-2 );
      r3 = R.submat(j+k-1, j, j+k-1, j+k-2);
      R2 = R.submat(j,j,j+k-2,j+k-2);
      r1primeR2Invr1prime = (r1*R2.i()*r1.t()).eval()(0,0);
      r3primeR2Invr3prime = (r3*R2.i()*r3.t()).eval()(0,0);
      r1primeR2Invr3prime = (r1*R2.i()*r3.t()).eval()(0,0);
      Denom = sqrt( (1.0 - r1primeR2Invr1prime )*( 1.0 - r3primeR2Invr3prime ) );
      R(j-1,j+k-1) = r1primeR2Invr3prime + Denom*PhiMAT(j-1,j+k-1);
      R(j+k-1, j-1) = R(j-1,j+k-1);

    }
  }
  
  return R;
  
}

// [[Rcpp::export]]
arma::cube MakePhiArrayRcpp(const arma::vec& x){
  
  const arma::vec xprime = 2.0*normcdf(x) - ones(x.n_elem);
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
double ComputeLogKfuncSameTcpp( const arma::mat& QNormUcondMatSegment, const arma::mat& phiMat){
  
  //QNormUcondMatSegment should be of the form
  //u_a(t)|a(t)....u_a(t)|b(t)
  //....
  //u_b(t)|a(t)....u_b(t)|b(t)
  
  //phiMat should be the first slice of phiArray (discounting the cpp index convention)
  
  arma::vec left_U(6);
  arma::vec right_U(6);
  arma::vec phiVec(6);
  left_U(0) = QNormUcondMatSegment(1,1);
  left_U(1) = QNormUcondMatSegment(2,1);
  left_U(2) = QNormUcondMatSegment(2,2);
  left_U(3) = QNormUcondMatSegment(3,1);
  left_U(4) = QNormUcondMatSegment(3,2);
  left_U(5) = QNormUcondMatSegment(3,3);
  right_U(0) = QNormUcondMatSegment(0,0);
  right_U(1) = QNormUcondMatSegment(0,1);
  right_U(2) = QNormUcondMatSegment(1,1);
  right_U(3) = QNormUcondMatSegment(0,2);
  right_U(4) = QNormUcondMatSegment(1,2);
  right_U(5) = QNormUcondMatSegment(2,2);
  phiVec(0) = phiMat(1,0);
  phiVec(1) = phiMat(2,0);
  phiVec(2) = phiMat(2,1);
  phiVec(3) = phiMat(3,0);
  phiVec(4) = phiMat(3,1);
  phiVec(5) = phiMat(3,2);
  
  const double AnsRet = sum( -0.5*log(1.0 - square(phiVec)) - 0.5*(square(phiVec)%( square(left_U) + square(right_U)  ) - 2*phiVec%left_U%right_U  )/(1.0 - square(phiVec)) );
  return AnsRet;
  
}


// [[Rcpp::export]]
double ComputeLogKfuncDifferentTcpp( const arma::mat& QNormUcondMatLeftSegment, const arma::mat& QNormUcondMatRightSegment, const arma::mat& phiMat){
  
  //QNormUcondMatLeftSegment should be a m by m matrix with
  //u_a(t)|a(s)+1.....u_a(t)|b(s)+1
  //........
  //u_b(t)|a(s)+1.....u_b(t)|b(s)+1
  
  //QNormUcondMatRightSegment should be a m by m matrix with
  //u_a(s)|a(t)-1.....u_a(s)|b(t)-1
  //........
  //u_b(s)|a(t)-1.....u_b(s)|b(t)-1
  
  //phiMat should be the (t-k) entry of the phiArray (not accounting for cpp starting index)
  
  const arma::vec& vecQNormUcondMatLeftSegment = vectorise( QNormUcondMatLeftSegment.t() );
  const arma::vec& vecQNormUcondMatRightSegment = vectorise( QNormUcondMatRightSegment );
  const arma::vec& vecphiMat = vectorise( phiMat.t() );
  
  const double AnsRet = sum( -0.5*log(1.0 - square(vecphiMat)) - 0.5*(square(vecphiMat)%( square(vecQNormUcondMatLeftSegment) + square(vecQNormUcondMatRightSegment)  ) - 2*vecphiMat%vecQNormUcondMatLeftSegment%vecQNormUcondMatRightSegment )/(1.0 - square(vecphiMat)) );
  return AnsRet;
  
}

// [[Rcpp::export]]
double ComputeLikelihoodCopulaPartcpp( const arma::mat& QNormCondUinput, const arma::cube& PhiArrayInput, const int& p, const int& m){
  
  const int maxT = QNormCondUinput.n_cols/m;
  double AnsLikCop = ComputeLogKfuncSameTcpp( QNormCondUinput.submat(0,0,m-1,m-1), PhiArrayInput.slice(0) );
  int a_tInd, b_tInd, a_tPrimeInd,b_tPrimeInd ;
  
  for(int tInd = 2; tInd <= maxT; tInd++){
    
    a_tInd = (tInd-1)*m + 1;
    b_tInd = m*tInd;
    AnsLikCop +=  ComputeLogKfuncSameTcpp( QNormCondUinput.submat(a_tInd-1,a_tInd-1,b_tInd-1,b_tInd-1), PhiArrayInput.slice(0) );
    
    for(int tPrimeInd = max(1, tInd - p); tPrimeInd < tInd; tPrimeInd++){
      
      a_tPrimeInd = (tPrimeInd-1)*m + 1;
      b_tPrimeInd = m*tPrimeInd;
      AnsLikCop += ComputeLogKfuncDifferentTcpp( QNormCondUinput.submat( a_tInd-1, a_tPrimeInd, b_tInd-1, b_tPrimeInd  ), QNormCondUinput.submat( a_tPrimeInd-1, a_tInd - 2, b_tPrimeInd-1, b_tInd - 2), PhiArrayInput.slice(tInd - tPrimeInd)    );
    }
    
  }
  return AnsLikCop;
  
}


// [[Rcpp::export]]
double NewLogJointPosteriorUSMacrocpp( const arma::mat& Ymat, const arma::vec& PARsub, const int& p){

  const int m = Ymat.n_cols;
  const int numPars = PARsub.n_elem;
  
  const arma::vec Xivec = PARsub.rows(0,m-1);
  const arma::vec LogOmegavec = PARsub.rows(m,2*m-1);
  const arma::vec Alphavec = PARsub.rows(2*m, 3*m-1);
  const arma::vec LogNuvec = PARsub.rows(3*m, 4*m-1);

  arma::cube PhiCube(m,m,p+1);
  PhiCube(0,0,0) = 999999999;
  PhiCube(1,1,0) = 999999999;
  PhiCube(2,2,0) = 999999999;
  PhiCube(3,3,0) = 999999999;
  PhiCube(1,0,0) = PARsub(4*m);
  PhiCube(2,0,0) = PARsub(4*m+1);
  PhiCube(3,0,0) = PARsub(4*m+2);
  PhiCube(2,1,0) = PARsub(4*m+3);
  PhiCube(3,1,0) = PARsub(4*m+4);
  PhiCube(3,2,0) = PARsub(4*m+5);
  PhiCube.slice(0) = PhiCube.slice(0) + (PhiCube.slice(0)).t();
  arma::cube PnormPhiCube(m,m,p+1);
  PnormPhiCube.slice(0) = 2*normcdf(PhiCube.slice(0)) - 1.0;
  for(int mInd = 0; mInd < m; mInd++){
    for(int mIndPrime = 0; mIndPrime < m; mIndPrime++){
      
      PhiCube(mIndPrime, mInd, 1) = PARsub(4*m+6+mInd*m+mIndPrime);
      PhiCube(mIndPrime, mInd, 2) = PARsub(4*m+6+mInd*m+mIndPrime+m*m);
      PhiCube(mIndPrime, mInd, 3) = PARsub(4*m+6+mInd*m+mIndPrime+2*m*m);
      PhiCube(mIndPrime, mInd, 4) = PARsub(4*m+6+mInd*m+mIndPrime+3*m*m);
      
    }
   
    
  }
  for(int pInd = 1; pInd <= p; pInd++){
    
    PnormPhiCube.slice(pInd) = 2*normcdf( PhiCube.slice(pInd) ) - 1.0;
    
  }
  const arma::vec LogTauvec = PARsub.rows(numPars-p-1,numPars-1);
  const arma::vec Marg1Par = {PARsub(0), PARsub(4), PARsub(8), PARsub(12)};
  const arma::vec Marg2Par = {PARsub(1), PARsub(5), PARsub(9), PARsub(13)};
  const arma::vec Marg3Par = {PARsub(2), PARsub(6), PARsub(10), PARsub(14)};
  const arma::vec Marg4Par = {PARsub(3), PARsub(7), PARsub(11), PARsub(15)};
  
  const arma::vec U1sub = SkewT_cdf( Ymat.col(0), PARsub(0), PARsub(4), PARsub(8), PARsub(12) );
  const arma::vec U2sub = SkewT_cdf( Ymat.col(1), PARsub(1), PARsub(5), PARsub(9), PARsub(13) );
  const arma::vec U3sub = TruncatedSkewT_cdf( Ymat.col(2), PARsub(2), PARsub(6), PARsub(10), PARsub(14) );
  const arma::vec U4sub = TruncatedSkewT_cdf( Ymat.col(3), PARsub(3), PARsub(7), PARsub(11), PARsub(15) );
  
  
  const arma::mat UMAT = join_rows( U1sub, U2sub, U3sub, U4sub );
  arma::vec UMATvec = vectorise(UMAT.t());
  UMATvec.elem(arma::find(UMATvec>0.9999)).fill(0.9999);
  UMATvec.elem(arma::find(UMATvec<0.0001)).fill(0.0001);
  const arma::mat ConditionalUObjsub = ComputeConditionalUsRcpp(m, p, UMATvec, PnormPhiCube);
  const double LikContri = ComputeLikelihoodCopulaPartcpp(ConditionalUObjsub, PnormPhiCube, p, m);
  const double LikMarg1 = SkewTlogLikcpp(Marg1Par, Ymat.col(0));
  const double LikMarg2 = SkewTlogLikcpp(Marg2Par, Ymat.col(1));
  const double LikMarg3 = sum(TruncatedSkewTlogPDFcpp(Marg3Par, Ymat.col(2)));
  const double LikMarg4 = sum(TruncatedSkewTlogPDFcpp(Marg4Par, Ymat.col(3)));

  
  //Probably need to decrease scale for alpha and increase rate for ExpExp
  double PriorContri = dLogNormRcpp(Xivec, 100) + dLogExpHalfCauchyRcpp( LogOmegavec, 5.0 ) + dLogLaplaceRcpp( Alphavec, 10 ) + dLogExpExpRcpp( LogNuvec, 10 ) + dLogNormRcpp( PARsub.rows(4*m, 4*m+5), exp(0.5*LogTauvec(0)) ) + dLogNormRcpp( PARsub.rows(4*m+6, 4*m+5+m*m), exp(0.5*LogTauvec(1)) ) + dLogNormRcpp( PARsub.rows(4*m+6+m*m, 4*m+5+2*m*m), exp(0.5*LogTauvec(2)) ) + dLogNormRcpp( PARsub.rows(4*m+6+2*m*m, 4*m+5+3*m*m), exp(0.5*LogTauvec(3)) ) + dLogNormRcpp( PARsub.rows(4*m+6+3*m*m, 4*m+5+4*m*m), exp(0.5*LogTauvec(4)) ) + dLogExpHalfCauchyRcpp( LogTauvec, 10.0 );
  const double FinalAns = LikContri + LikMarg1 + LikMarg2 + LikMarg3 + LikMarg4 + PriorContri;
  
  return FinalAns;
}

// [[Rcpp::export]]
double NewLogJointPosteriorUSMacroReverseStage1cpp( const arma::mat& RankMAT, const arma::vec& PARsub, const int& p){
  
  const int Tlength = RankMAT.n_rows;
  const int m = RankMAT.n_cols;
  const int numPars = PARsub.n_elem;
 
  
  arma::cube PhiCube(m,m,p+1);
  PhiCube(0,0,0) = 999999999;
  PhiCube(1,1,0) = 999999999;
  PhiCube(2,2,0) = 999999999;
  PhiCube(3,3,0) = 999999999;
  PhiCube(1,0,0) = PARsub(0);
  PhiCube(2,0,0) = PARsub(1);
  PhiCube(3,0,0) = PARsub(2);
  PhiCube(2,1,0) = PARsub(3);
  PhiCube(3,1,0) = PARsub(4);
  PhiCube(3,2,0) = PARsub(5);
  PhiCube.slice(0) = PhiCube.slice(0) + (PhiCube.slice(0)).t();
  arma::cube PnormPhiCube(m,m,p+1);
  PnormPhiCube.slice(0) = 2*normcdf(PhiCube.slice(0)) - 1.0;
  const int NumPhiPar = m*(m-1)/2 + p*m*m;
  const arma::mat AMAT = (RankMAT - 1.0)/(Tlength+1.0);
  const arma::mat BMAT = RankMAT/(Tlength+1.0);
  const arma::vec Zvec = PARsub.rows(NumPhiPar+p+1, numPars-1);
  arma::mat UMAT(Tlength, m);
  
  for(int mInd = 0; mInd < m; mInd++){
    for(int mIndPrime = 0; mIndPrime < m; mIndPrime++){
      
      PhiCube(mIndPrime, mInd, 1) = PARsub(6+mInd*m+mIndPrime);
      PhiCube(mIndPrime, mInd, 2) = PARsub(6+mInd*m+mIndPrime+m*m);
      PhiCube(mIndPrime, mInd, 3) = PARsub(6+mInd*m+mIndPrime+2*m*m);
      PhiCube(mIndPrime, mInd, 4) = PARsub(6+mInd*m+mIndPrime+3*m*m);
      
    }
    
    UMAT.col(mInd) = ( BMAT.col(mInd) - AMAT.col(mInd)  )%normcdf(Zvec.rows( mInd*Tlength, (mInd+1)*Tlength-1 )) + AMAT.col(mInd);
    
  }
  for(int pInd = 1; pInd <= p; pInd++){
    
    PnormPhiCube.slice(pInd) = 2*normcdf( PhiCube.slice(pInd) ) - 1.0;
    
  }
  const arma::vec LogTauvec = PARsub.rows(NumPhiPar,NumPhiPar+p);
  const arma::vec UMATvec = vectorise(UMAT.t());
  const arma::mat ConditionalUObjsub = ComputeConditionalUsRcpp(m, p, UMATvec, PnormPhiCube);
  const double LikContri = ComputeLikelihoodCopulaPartcpp(ConditionalUObjsub, PnormPhiCube, p, m);
  
  //Probably need to decrease scale for alpha and increase rate for ExpExp
  double PriorContri = dLogNormRcpp( PARsub.rows(0, 5), exp(0.5*LogTauvec(0)) ) + dLogNormRcpp( PARsub.rows(6, 5+m*m), exp(0.5*LogTauvec(1)) ) + dLogNormRcpp( PARsub.rows(6+m*m, 5+2*m*m), exp(0.5*LogTauvec(2)) ) + dLogNormRcpp( PARsub.rows(6+2*m*m, 5+3*m*m), exp(0.5*LogTauvec(3)) ) + dLogNormRcpp( PARsub.rows(6+3*m*m, 5+4*m*m), exp(0.5*LogTauvec(4)) ) + dLogExpHalfCauchyRcpp( LogTauvec, 10.0 ) + dLogNormRcpp(Zvec, 1.0);
  const double FinalAns = LikContri + PriorContri;
  
  return FinalAns;
}

// [[Rcpp::export]]
double NewLogJointPosteriorSkewTUSMacroCutStage1cpp(const arma::vec& Yvec, const arma::vec& PARsub){
  
  const arma::vec Xi = PARsub.rows(0,0);
  const arma::vec LogOmega = PARsub.rows(1,1);
  const arma::vec Alpha = PARsub.rows(2,2);
  const arma::vec LogNu = PARsub.rows(3,3);
  
  const double LikPart = arma::sum(SkewTlogPDFcpp(PARsub, Yvec));
  const double PriorPart = dLogNormRcpp(Xi, 100) + dLogExpHalfCauchyRcpp( LogOmega, 5.0 ) + dLogLaplaceRcpp( Alpha, 10 ) + dLogExpExpRcpp( LogNu, 10 );
  
  return LikPart+PriorPart;
}

// [[Rcpp::export]]
double NewLogJointPosteriorTruncatedSkewTUSMacroCutStage1cpp(const arma::vec& Yvec, const arma::vec& PARsub){
  
  const arma::vec Xi = PARsub.rows(0,0);
  const arma::vec LogOmega = PARsub.rows(1,1);
  const arma::vec Alpha = PARsub.rows(2,2);
  const arma::vec LogNu = PARsub.rows(3,3);
  
  const double LikPart = arma::sum(TruncatedSkewTlogPDFcpp(PARsub, Yvec));
  const double PriorPart = dLogNormRcpp(Xi, 100) + dLogExpHalfCauchyRcpp( LogOmega, 5.0 ) + dLogLaplaceRcpp( Alpha, 10 ) + dLogExpExpRcpp( LogNu, 10 );
  
  return LikPart+PriorPart;
}
