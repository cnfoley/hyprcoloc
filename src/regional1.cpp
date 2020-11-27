//Includes/namespaces
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/Core>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace RcppEigen;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Rcpp::as;
using Eigen::LLT;
// [[Rcpp::export]]

List regional1(NumericMatrix zTemp1, NumericMatrix wTemp1, NumericVector traitsCo, NumericMatrix trait_cor, NumericVector EPS){

 // map Z-scores and prior ratios W by co-loc (1) and no co-loc (2)
 const Map<MatrixXd> matZ(as<Map<MatrixXd> >(zTemp1));
 const Map<MatrixXd> matW(as<Map<MatrixXd> >(wTemp1));

 // map n.cvs

 // map trait correlation matrix 
 const Map<MatrixXd> TRAIT_COR(as<Map<MatrixXd> >(trait_cor));

 // map correlation matrix for traits under consideration
 const Map<VectorXd> trtCo(as<Map<VectorXd> >(traitsCo));

 // map SNPs by co-loc (1) and no co-loc (2)
 const Map<VectorXd> eps(as<Map<VectorXd> >(EPS));

 // Set-up ordered Z and W score matrices
 int m(matZ.cols());

 MatrixXd sigZtmp(TRAIT_COR.rows(), m);
 MatrixXd sigZ(m, m);

 int Q(matZ.rows());
 MatrixXd Wnew = 0*matW;

 MatrixXd Ep = eps(0)*(VectorXd::Ones(m)).asDiagonal();
 MatrixXd I = (VectorXd::Ones(m)).asDiagonal();

 for (int i = 0; i < trtCo.size(); i++){
  Wnew.col(trtCo(i)-1) = matW.col(trtCo(i)-1);
 }

 for (int i = 0; i < m; i++){
  sigZtmp.col(i) = TRAIT_COR.col(2*i);
 }
 for (int j = 0; j < m; j++){
  sigZ.row(j) = sigZtmp.row(2*j);
 }

 sigZ = (sigZ + Ep).eval();

 MatrixXd Z = matZ;
 MatrixXd W(Q, m);
 MatrixXd Wtmp(1, m);

 // Define the diagonal prior ratio matrix "Wdiag"
 MatrixXd Wdiag(m, m);

 // Define the adjusted prior matrix "adjZ"
 MatrixXd adjZ(m, m);

 // Define the inverse sigma matrix "SIGINV"
 MatrixXd SIGINV(m, m);

 // Define the output matrix "out" containing a row of 
 // determinant prefactor constants and BF exponents
 VectorXd out(Q);


 for (int i = 0; i < Q; i++){

  Wtmp.row(0) = Wnew.row(i);
  VectorXd Wco(Map<VectorXd>(Wtmp.data(), m));
  W.row(i) = Wco.transpose();

  Wdiag = W.row(i).asDiagonal();
  adjZ = Wdiag*sigZ*Wdiag;

  MatrixXd ZL(sigZ.llt().matrixL().solve(I));
  MatrixXd ZAL((sigZ + adjZ).llt().matrixL().solve(I));

  SIGINV = 0.5*(ZAL.transpose()*ZAL)*adjZ*(ZL.transpose()*ZL);
  out(i) = -log(ZL.diagonal().prod()/ZAL.diagonal().prod()) + ((Z.row(i))*(SIGINV)*Z.row(i).transpose());

 }

 VectorXd out1_co(3);
 ArrayXd out_exp(Q);

 MatrixXd::Index maxIndex;

 out1_co(0) = out.maxCoeff(&maxIndex);
 out_exp = (out - out1_co(0)*VectorXd::Ones(Q)).array().exp();
 out1_co(1) = out_exp.sum();
 out1_co(2) = maxIndex + 1;

 if (trtCo.size() == m){
  return List::create(Named("out_co_1") = out1_co, Named("out_snp_scores") = out_exp/out1_co(1));
 }
 else {
  return wrap(out1_co);
 }
}
