#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <cmath>
#include <Rcpp.h>
// #include "subfun.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// // [[Rcpp::export]]
// arma::vec wpow(arma::vec x, double gamma0) {
//   //uvec ind = find(x==0);
//   vec y = pow(abs(x),-1.0*gamma0);
//   y.elem(find(x==0)).zeros();
//   y.elem(find(x==0)).fill(10000*y.max());
//   return y;
// }

// [[Rcpp::export]]
int nzcount(arma::vec x) {
  vec y = nonzeros(x) ;
  return y.n_elem;
}



void SampleNoReplace( arma::uvec &index, int nOrig, int size) {
  int ii, jj;
  arma::uvec sub(nOrig);
  for (ii = 0; ii < nOrig; ii++) {
    sub(ii) = ii;
  }
  for (ii = 0; ii < size; ii++) {
    jj = nOrig * unif_rand();
    index(ii) = sub(jj);
    // replace sampled element with last, decrement
    sub(jj) = sub(--nOrig);
  }
}


// // [[Rcpp::export]]
// arma::uvec mySdiff(arma::uvec x, arma::uvec y){
//   // cout<<zeros<mat>(4,5);
//   for (int j = 0; j < (int) y.n_elem; j++)
//     x = x.elem(find(x != y(j)));
//   return(x);
// }

// // [[Rcpp::export]]
// arma::mat xx1(arma::mat x, arma::vec y) {
//   mat sfac = stddev(x);
//   for(int i=0; i < x.n_cols; i++){
//     if(sfac(0,i)!=0) sfac(0,i) = 1/sfac(0,i); else sfac(0,i) = 1;
//   }
//   x  = x.each_row()%sfac;
//   cout << sqrt(accu(square(x.t()*y)))<< std::endl;
//   return x.t()*y;
// }

// // [[Rcpp::export]]
// arma::mat xx2(arma::mat x, arma::vec y) {
//   cout << y(1)+ x.col(1)%y<< std::endl;
//   return y;
// }


// // [[Rcpp::export]]
// arma::mat xx3(arma::mat x) {
//   x.insert_cols(0, ones<vec>(5));
//   return x;
// }

// [[Rcpp::export]]
double softThres(const double x, const double lambda) {
  return((x > lambda) ? x - lambda :
           (x < -lambda) ? x + lambda : 0.);
}


// [[Rcpp::export]]
double hardThres(const double x, const double lambda) {
  return((x > lambda) ? x  :
           (x < -lambda) ? x : 0.);
}


// [[Rcpp::export]]
arma::mat clsq(arma::mat X, arma::vec y, arma::mat C) {
  int p = X.n_cols, k = C.n_rows;
  arma::mat xx = X.t()*X; 
  xx.insert_rows(p,C);
  C = C.t();
  C.insert_rows(p,zeros<mat>(k,k));
  xx.insert_cols(p,C);
  arma::vec xy = X.t()*y;
  xy.insert_rows(p,zeros<vec>(k));
  xy = pinv(xx + 0.0*diagmat(ones<vec>(xx.n_cols)))*xy;
  return xy.subvec(0,p-1);
}


class thres {
  double x,lambda;
  double hThres(){
    return((x > lambda) ? x  :
             (x < -lambda) ? x : 0.);
  }
  double sThres(){
    return((x > lambda) ? x - lambda :
             (x < -lambda) ? x + lambda : 0.);
  }
public:
  double calculate (int operatr, double operand1, double operand2){
    double (thres::*opPtr)() = NULL;
    x = operand1;
    lambda = operand2;
    if (operatr == 1) opPtr = &thres::sThres;
    if (operatr == 0) opPtr = &thres::hThres;
    return (this->*opPtr)();
  }
};


// [[Rcpp::export]]
double thresC(int operatr, double x, double lambda) {
  // x = operand1;
  // lambda = operand2;
  if (operatr == 1) {
    return((x > lambda) ? x - lambda :
             (x < -lambda) ? x + lambda : 0.);
  }
  if (operatr == 0) {
    return((x > lambda) ? x  :
             (x < -lambda) ? x : 0.);
  }    else { return 0; }
}





// //' Engine function
// //'
////' @export
// [[Rcpp::export]]
Rcpp::List classoshe(arma::mat Xt, arma::vec y, arma::mat C, arma::vec we, double lam0,
                     Rcpp::List control){
  
  int i,n = Xt.n_rows, p = Xt.n_cols, maxiter = control["maxiter"], kk = C.n_rows;
  Rcpp::List out;
  mat sfac = stddev(Xt);
  double nu = control["nu"],mu,tol = control["tol"];
  for(i=0; i < p; i++){
    if(sfac(0,i)!=0) sfac(0,i) = 1/sfac(0,i); else sfac(0,i) = 1;
  }
  Xt  = Xt.each_row()%sfac;
  C  = C.each_row()%sfac;
  mat xx = Xt.t()*Xt/n , cc = C.t()*C;
  vec xy = Xt.t()*y/n;  
  vec lam, betap = zeros<vec>(p) , betau = zeros<vec>(p), psi;
  vec xxdiag = xx.diag(), ccdiag = cc.diag();
  double ssq2, ssq = sqrt(accu(square(y-Xt*betap))/n),err,du;
  we = we*lam0;
  int k,j;
  xx.diag().zeros();  cc.diag().zeros();
  for(i=0; i < maxiter; i++){
    ssq2=ssq;
    lam = we*ssq2;
    mu =  control["mu"]; psi = zeros<vec>(kk);
    for(k=0; k < maxiter; k++){
      for(j=0; j < p; j++){
        du = (xy(j) - sum(xx.row(j)*betau)); //// check 
        du = du-mu*(sum(cc.row(j)*betau) + sum(C.col(j)%psi));
        betau(j) = softThres(du,lam(j))/(xxdiag(j) + mu*ccdiag(j));
      }
      mu = mu*nu;
      psi = psi+ C*betau;
      err = sqrt(sum(square(betap-betau))); 
      betap = betau;
      if(err < tol && k >0) break;
    }
    ssq = sqrt(accu(square(y-Xt*betap))/n);
    if(abs(ssq-ssq2) < tol && i >0) break;
  }
  out["beta"] = betap%vectorise(sfac);
  return(out);
}




// [[Rcpp::export]]
Rcpp::List classol2(arma::mat Xt, arma::vec y, arma::mat C, arma::vec we, double lam0,
                    Rcpp::List control){
  int i,n = Xt.n_rows, p = Xt.n_cols, maxiter = control["maxiter"], kk = C.n_rows;
  Rcpp::List out;
  mat sfac = stddev(Xt);
  double tol = control["tol"]; //mu, nu = control["nu"],
  for(i=0; i < p; i++){
    if(sfac(0,i)!=0) sfac(0,i) = 1/sfac(0,i); else sfac(0,i) = 1;
  }
  Xt  = Xt.each_row()%sfac;
  C  = C.each_row()%sfac;
  mat xx = Xt.t()*Xt/n , cc = C.t()*C;
  mat U4,t5,V4;
  vec s4;
  svd(U4,s4,V4,xx + cc);
  
  vec xy = Xt.t()*y/n, err2;  
  vec lam, betap = zeros<vec>(p) , betau = zeros<vec>(p), psi;
  vec xxdiag = xx.diag(), ccdiag = cc.diag();
  double ssq2, ssq = sqrt(accu(square(y-Xt*betap))/n); //,du
  we = we*lam0;
  int k; //,j
  xx.diag().zeros();  cc.diag().zeros();
  for(i=0; i < maxiter; i++){
    ssq2=ssq;
    // cout << ssq << std::endl;
    lam = we*ssq2; lam = s4 + lam;
    psi = zeros<vec>(kk); // mu =  control["mu"];
    t5 = U4*diagmat(1/lam)*U4.t();
    // cout << lam << std::endl;
    // cout << sum(U4) << std::endl;
    // cout << sum(t5) << std::endl;
    for(k=0; k < maxiter; k++){
      betau = t5*(xy - C.t()*psi);
      err2 = C*betau;
      psi = psi+ err2;
      if(sum(abs(err2)) < tol && k >0) break;
    }
    
    ssq = sqrt(accu(square(y-Xt*betap))/n);
    if(abs(ssq-ssq2) < tol && i >0) break;
  }
  
  out["beta"] = betau%vectorise(sfac);
  return(out);
}





// //' Engine function 
// //' 
// @export
// [[Rcpp::export]]
Rcpp::List robregcc_nsp5(arma::mat X, arma::vec y, arma::mat C, int intercept, arma::vec gammawt,
                        arma::vec lampath, arma::vec shwt, Rcpp::List control, int ptype){
  // c('adaptive','soft','hard')[3]
  // clsq(arma::mat X, arma::vec y, arma::mat C)
  int n = X.n_rows, p = X.n_cols, outMiter = control["outMiter"];
  int nlam = lampath.n_rows,j,pind =1,kk,jj,fpath = control["fullpath"];
  double spn= control["sp"]; // ,lminf = control["lminfac"], gamma0 = control["gamma"],
  uvec tm;
  thres pt;
  Rcpp::List out;
  if(intercept==1){
    X.insert_cols(0, ones<vec>(n));
    p = p+1;
  }
  mat U, V, h2,ginvx;
  vec d,re,lam,tm2;
  // control parameter
  svd_econ(U,d,V,X.t()*X);
  mat x11 = diagmat(ones<vec>(p)) - C.t()*pinv(C*C.t())*C;
  x11 = X*x11;
  // h2 = X*pinv(X.t()*X)*X.t(); // U*U.t();
  h2 = x11*pinv(x11.t()*x11)*x11.t(); // U*U.t();
  re = y - h2*y;
  // d = pow(d,2);
  tm = find(d>1.490116e-08);
  ginvx  = V.cols(tm)*diagmat(1/d(tm))*V.cols(tm).t();
  // preparing function and weights
  // shwt = ones<vec>(n);
  if(ptype==3) pind = 0;
  arma::mat betapath =  zeros<mat>(p,nlam);
  arma::mat betapathR =  zeros<mat>(p,nlam);
  arma::mat gammapath =  zeros<mat>(n,nlam);
  vec gammap=zeros<vec>(n), gammau = zeros<vec>(n);
  spn = n*spn;
  double err,tol = control["tol"];
  // for loop of the body 
  for(j=0; j < nlam; j++){
    lam = shwt*lampath(j);
    if(ptype==3) gammap = gammawt;
    // gammap = gammawt;
    err = 1e6; kk=0;
    while ( (kk < outMiter)){
      tm2 = h2*gammap + re;
      for(jj=0; jj<n; jj++)
        gammau(jj) = pt.calculate(pind,tm2(jj),lam(jj));
      err = max(abs(gammap-gammau));//accu(abs(gammap)); 
      gammap = gammau;
      kk = kk+1;
      if((kk >2) && ( err < tol) ) { break;}
    }
    // cout << err << std::endl;
    // cout << sum(gammau) << std::endl;
    gammapath.col(j)=gammau;
    betapath.col(j) = clsq(X, y-gammau,C);
    // betapath.col(j) = ginvx*X.t()*(y-gammau);
    // tm = find(abs(gammau)>0);
    // betapathR.col(j) = pinv(X.rows(tm).t()*X.rows(tm))*X.rows(tm).t()*y(tm);
    // betapathR.col(j) = clsq(X.rows(tm), y(tm),C);
    // betapathR.col(j) = c_ridge(X.rows(tm), y(tm), C, 5,10);
    // cout << j << std::endl;
    
    if(((nzcount(gammau)+1) > spn) && (fpath == 0)) {break;}
  }
  j=std::min(j,nlam-1);
  betapath = betapath.cols(0,j);
  // betapathR = betapathR.cols(0,j);
  gammapath = gammapath.cols(0,j);
  lampath = lampath.subvec(0,j);
  out["Method"] = ptype;
  out["betapath"] = betapath;
  // out["betapathR"] = betapathR;
  out["gammapath"] = gammapath;
  out["lampath"] = lampath;
  // out["X"] = X;
  // out["y"] = y;
  // out["U"] = U;
  return(out);
}








// // [[Rcpp::export]]
// double xx5(double a, double b) {
//   thres t;
//   cout << t.calculate (1, a, b);
//   cout << t.calculate (0, a, b);
//   return t.calculate (0, a, b);
// }

// // [[Rcpp::export]]
// double max1(double a, double b) {
//   if(a>b) return a; else return b;
// }








// [[Rcpp::export]]
Rcpp::List robregcc_sp5(arma::mat X, arma::vec y, arma::mat C, 
                        arma::vec paramin, Rcpp::List control,
                        arma::vec shwt, arma::vec lampath,
                        int ptype, double k0, double alpha){
  // Mean centered data:
  int i,j,n = X.n_rows, p = X.n_cols, outMiter = control["outMiter"];
  int nlam = lampath.n_rows,pind =1, kk=0,jj=0 ,k=C.n_rows, np = n+p;
  double mu = control["mu"],nu = control["nu"],spb = control["spb"],spn= control["spy"];
  double k01 = 1-(1/k0); // n=1
  
  if(ptype==3) pind = 0;
  lampath = lampath/k0;
  
  // thres pt;
  std::function<double(double , double )> thr;
  if(pind==1)
    thr = softThres;
  if(pind==0)
    thr = hardThres;
  
  //  check for uninitialized variable 
  // store some value 
  vec xy1 = X.t()*y, y1(n);
  xy1 = xy1/k0;y1 = sqrt(1)*y/k0; // n=1
  mat C2 = (mu/k0)*C, xx = (X.t()*X)/k0, cc = (mu/k0)*C.t()*C;
  mat xx11 = diagmat(ones<vec>(p))- (xx+ cc); 
  mat x1 = sqrt(1)*X/k0, x1t = x1.t();
  
  // output structure for storing  
  arma::mat betapath =  zeros<mat>(p,nlam);
  arma::mat gammapath =  zeros<mat>(n,nlam);
  vec gammau=zeros<vec>(n), betau = zeros<vec>(p),elfac,lamgamma(n),lambeta(p);
  vec psi1,psi,xy2,betap,gammap,tempx,tempy;
  spn = n*spn; spb = p*spb;
  
  //  variables defined for the loop 
  double errout,tol = control["tol"],spgamma,spbeta, outtol = control["out.tol"];
  double errin,erbet,ergm;
  int inMiter = control["inMiter"];
  uvec tm2(n),tm3(p);
  psi = zeros<vec>(k); psi1 = zeros<vec>(k);
  int fpath = control["fullpath"];
  // for loop of the body 
  for(i=0; i < nlam; i++){
    // lambeta = alpha*lampath(i)*shwt.subvec(0,p-1); //
    // elfac = 1/(1+(1-alpha)*shwt.subvec(0,p-1)*lampath(i));// 
    // lamgamma = shwt.subvec(p,np-1)*lampath(i); // 
    elfac = 1/(1+(1-alpha)*lampath(i));// 
    lambeta = alpha*shwt.subvec(0,p-1)*lampath(i)*elfac; //
    // elfac = 1/(1+(1-alpha)*lampath(i));// 
    lamgamma = alpha*shwt.subvec(p,np-1)*lampath(i)*elfac; // 
    if(i > 0) {
      betau = betapath.col(i-1);
      gammau = gammapath.col(i-1);
    }
    
    if(ptype==3) {
      betau = paramin.subvec(0,p-1);
      gammau = paramin.subvec(p,np-1);
    }
    
    
    errout = 1e6;kk = 0;
    mu = control["mu"]; psi = psi*0; psi1 = psi1*0;
    spgamma = 1;spbeta = 1;
    //  while loop 
    while ((errout > outtol) && (kk < outMiter) && (0 < spgamma) && (spgamma < spn) && (0 < spbeta) && (spbeta < spb) ){
      errin = 1e6; j = 0;
      xy2 = xy1 - C2.t()*psi;
      
      betap = betau ; gammap = gammau; erbet = 10; ergm = 10;
      // inner core while loop 
      while ( (j < inMiter) && (erbet>0) && (ergm>0)){  //(errin > tol) &&
        tm2 = find(gammau!=0); 
        tm3 = find(betau!=0); 
        tempx = xy2 +  xx11.cols(tm3)*betau(tm3) - x1t.cols(tm2)*gammau(tm2); 
        tempy = y1 - x1.cols(tm3)*betau(tm3) + k01*gammau;
        
        tempy =tempy*elfac;
        for(jj=0; jj<n; jj++)
          gammau(jj) = thr(tempy(jj),lamgamma(jj));
        // gammau(jj) = pt.calculate(pind,tempy(jj),lamgamma(jj));
        
        
        tempx =tempx*elfac;
        for(jj=0; jj<p; jj++)
          betau(jj) = thr(tempx(jj),lambeta(jj));
        // betau(jj) = pt.calculate(pind,tempx(jj),lambeta(jj));
        
        
        erbet =max(abs(betap - betau));
        ergm = max(abs(gammap - gammau));
        // errin = erbet/accu(abs(betap)) + ergm/accu(abs(gammap));
        errin = erbet + ergm;
        gammap = gammau;
        betap = betau;
        j = j + 1;
        if((j>2) && ( errin < tol) ) { break;}
      }
      
      spgamma = nzcount(gammau); mu = mu*nu;
      tm3 = find(betau!=0); spbeta = tm3.n_elem;
      psi1 = C.cols(tm3)*betau(tm3);
      psi = psi + psi1;
      // errout = max(abs(psi - psi1));//accu(abs(psi1));
      errout = max(abs(psi1));
      // psi1 = psi;
      kk = kk+1;
    }
    // cout << j << " a " << kk << std::endl;
    
    gammapath.col(i)=gammau;
    betapath.col(i) = betau;
    if(((spgamma > spn) || (spbeta > spb)) && (fpath == 0)) {break;}
  }
  i=std::min(i,nlam-1);
  betapath = betapath.cols(0,i);
  gammapath = sqrt(1)*gammapath.cols(0,i); // n=1
  lampath = lampath.subvec(0,i);
  // betapath = betapath.each_col()%vectorise(sfac);
  // 
  return List::create(Named("betapath") = betapath,
                      Named("Method") = ptype,
                      Named("gammapath") = gammapath,
                      // Named("lampath") = lampath,
                      Named("nlambda") = i+1
  );
}









//'  Ridge regression with compositional covariates
//'
//' @param Xt CLR transformed predictor matrix. 
//' @param y model response vector
//' @param C Subcomposition matrix
//' @param nfold number of folds for crossvalidation
//' @param nlam number of lambda to generate solution path
//' @param control controling parameter for the model 
//' @return model parameter estimate 
//// @export
// [[Rcpp::export]]
Rcpp::List c_ridge2(arma::mat Xt, arma::vec y, arma::mat C, int nfold,
                    int nlam, Rcpp::List control){
  int i,n = Xt.n_rows, p = Xt.n_cols, maxiter = control["maxiter"], kk = C.n_rows;
  Rcpp::List out;
  mat sfac = stddev(Xt);
  double tol = control["tol"], lmfac = control["lminfac"]; //nu = control["nu"],mu,
  for(i=0; i < p; i++){
    if(sfac(0,i)!=0) sfac(0,i) = 1/sfac(0,i); else sfac(0,i) = 1;
  }
  Xt  = Xt.each_row()%sfac;
  C  = C.each_row()%sfac;
  uvec tm(n),tm1,tm2;
  
  arma::vec fid;
  fid = repmat(linspace<vec>(1,nfold,nfold), (n/nfold)+1, 1);
  fid = fid.subvec(0,n-1);
  SampleNoReplace(tm,n,n);
  fid = fid(tm);
  
  double lmax = max(abs(Xt.t()*y));
  arma::vec lampath =  exp(linspace<vec>(log(lmax),  log(lmax*lmfac), nlam));
  
  mat xx(p,p),xxtr(p,p) , cc = C.t()*C,U4,V4, dev = zeros<mat>(nfold,nlam);
  mat C2 = C.t();
  
  vec betau = zeros<vec>(p), psi,s4,xytr(p),lam(p),err(p);
  int k, ne,j; //,j
  // xx.diag().zeros();  cc.diag().zeros();
  for(i=1; i<=nfold; i++){
    tm2 = find(fid!=i);
    tm1 = find(fid==i);
    ne = tm2.n_elem;
    xytr = (Xt.rows(tm2).t()*y(tm2))/ne;
    svd(U4,s4,V4,(Xt.rows(tm2).t()*Xt.rows(tm2))/ne + cc,"std");
    for(j=0; j<nlam;j++){
      lam = 1/(s4+lampath(j));
      xx = U4*diagmat(lam)*U4.t();
      psi = zeros<vec>(kk);
      // convergence of the inner loop
      for(k=0; k < maxiter; k++){
        betau = xx*(xytr - C2*psi);
        err = C*betau; //sqrt(sum(square(betap-betau))); 
        psi = psi+ err;
        if(sum(abs(err)) < tol && k >0) break;
      }
      // cout << k << std::endl;
      ne = tm1.n_elem;
      dev(i-1,j) = accu(pow(y(tm1) - Xt.rows(tm1)*betau, 2))/ne;
    }
  }
  vec devsd = vectorise(stddev(dev,0,0))/as_scalar(sqrt(nfold));
  vec mndev = vectorise(mean(dev,0));
  uword gind = index_min(mndev);
  
  lmax = as_scalar( max(lampath(find(mndev <= (mndev(gind) + 0*devsd(gind))) ) ));
  
  svd(U4,s4,V4,(Xt.t()*Xt)/n + cc,"std");
  lam = 1/(s4+lmax);
  xx = U4*diagmat(lam)*U4.t();
  psi = zeros<vec>(kk);
  // convergence of the inner loop
  for(k=0; k < maxiter; k++){
    betau = xx*(xytr - C2*psi);
    err = C*betau; //sqrt(sum(square(betap-betau))); 
    psi = psi+ err;
    if(sum(abs(err)) < tol && k >0) break;
  }
  
  out["dev"] = dev;
  out["beta"] = betau%vectorise(sfac);
  return(out);
}

