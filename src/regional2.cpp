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
// [[Rcpp::export]]

List regional2(NumericMatrix zTemp1, NumericMatrix wTemp1, NumericMatrix sTemp1, NumericVector traitsCo, NumericMatrix rho, NumericMatrix trait_cor, NumericVector EPS){

 // map Z-scores and prior ratios W by co-loc (1) and no co-loc (2)
 const Map<MatrixXd> matZ(as<Map<MatrixXd> >(zTemp1));
 const Map<MatrixXd> matW(as<Map<MatrixXd> >(wTemp1));

 // map SNPs by co-loc (1) and no co-loc (2)
 const Map<MatrixXd> snps(as<Map<MatrixXd> >(sTemp1));

 // map SNP correlation matrix 
 const Map<MatrixXd> RHO(as<Map<MatrixXd> >(rho));

 // map trait correlation matrix 
 const Map<MatrixXd> TRAIT_COR(as<Map<MatrixXd> >(trait_cor));

 // map correlation matrix for traits under consideration
 const Map<VectorXd> trtCo(as<Map<VectorXd> >(traitsCo));

 // map SNPs by co-loc (1) and no co-loc (2)
 const Map<VectorXd> eps(as<Map<VectorXd> >(EPS));

 int n_cv(snps.rows()), ncol(snps.cols());
 int m(matZ.cols());

 MatrixXd sigZtmp(TRAIT_COR.rows(), m);
 MatrixXd sigZ(m, m);
 MatrixXd sigZ2 = TRAIT_COR;

 int Q(matZ.rows()), tmp(1);
 MatrixXd Wnew = 0*matW;

 MatrixXd Ep = eps(0)*(VectorXd::Ones(m)).asDiagonal();
 MatrixXd Ep2 = eps(0)*(VectorXd::Ones(n_cv*m)).asDiagonal();

 MatrixXd I = (VectorXd::Ones(m)).asDiagonal();
 MatrixXd I2 = (VectorXd::Ones(n_cv*m)).asDiagonal();

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
 sigZ2 = (sigZ2 + Ep2).eval();

 MatrixXd Z = matZ;
 MatrixXd W(Q, m);
 MatrixXd Wtmp(n_cv-1, m);

 MatrixXd Z2(ncol, n_cv*m);
 MatrixXd Ztmp2(n_cv, m);
 MatrixXd W2(ncol, n_cv*m);
 MatrixXd Wtmp2(n_cv, m);

 // Define the diagonal prior ratio matrix "Wdiag"
 MatrixXd Wdiag(m, m);
 MatrixXd Wdiag2(n_cv*m, n_cv*m);

 // Define the adjusted prior matrix "adjZ"
 MatrixXd adjZ(m, m);
 MatrixXd adjZ2(n_cv*m, n_cv*m);

 // Define the inverse sigma matrix "SIGINV"
 MatrixXd SIGINV(m, m);
 MatrixXd SIGINV2(n_cv*m, n_cv*m);

 // Define the output matrix "out" containing a row of 
 // determinant prefactor constants and BF exponents
 VectorXd out(Q);
 VectorXd out2(ncol);

 VectorXd tmp_rho(ncol);
 VectorXd tmp_rho_2(n_cv*m);
 VectorXd tmp_rho_21 = (VectorXd::Ones(n_cv*m));
 VectorXd tmp_rho_22 = (VectorXd::Ones(n_cv*m));

 MatrixXd sigZ2_new = sigZ2;

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

  for (int j = 0; j < n_cv; j++){
   tmp = snps(j,i)-1;
   Ztmp2.row(j) = matZ.row(tmp);
   Wtmp2.row(j) = Wnew.row(tmp);
  }

  VectorXd Zco2(Map<VectorXd>(Ztmp2.data(), n_cv*m));
  VectorXd Wco2(Map<VectorXd>(Wtmp2.data(), n_cv*m));

  Z2.row(i) = Zco2.transpose();
  W2.row(i) = Wco2.transpose();

  // compute the subseted correlation matrix for co and no-co Z-scores
  for (int k = 0; k < m; k++){
   tmp_rho(i) = RHO(snps(0 , i)-1, snps(1 , i)-1);
   tmp_rho_21(2*k+1) = tmp_rho(i);
   tmp_rho_22(2*k) = tmp_rho(i);
  }

  for (int k = 0; k < m; k++){
   tmp_rho_2 = sigZ2.col(2*k);
   sigZ2.col(2*k) = tmp_rho_2.cwiseProduct(tmp_rho_21);
   tmp_rho_2 = sigZ2.col(2*k + 1);
   sigZ2.col(2*k + 1) = tmp_rho_2.cwiseProduct(tmp_rho_22);
  }

  Wdiag2 = W2.row(i).asDiagonal();
  adjZ2 = Wdiag2*sigZ2*Wdiag2;

  MatrixXd ZL2(sigZ2.llt().matrixL().solve(I2));
  MatrixXd ZAL2((sigZ2 + adjZ2).llt().matrixL().solve(I2));

  SIGINV2 = 0.5*(ZAL2.transpose()*ZAL2)*adjZ2*(ZL2.transpose()*ZL2);
  out2(i) = -log(ZL2.diagonal().prod()/ZAL2.diagonal().prod()) + ((Z2.row(i))*(SIGINV2)*Z2.row(i).transpose());

  sigZ2 = sigZ2_new;

 }

for (int i = Q; i < ncol; i++){

  for (int j = 0; j < n_cv; j++){
   tmp = snps(j,i)-1;
   Ztmp2.row(j) = matZ.row(tmp);
   Wtmp2.row(j) = Wnew.row(tmp);
  }

  VectorXd Zco2(Map<VectorXd>(Ztmp2.data(), n_cv*m));
  VectorXd Wco2(Map<VectorXd>(Wtmp2.data(), n_cv*m));

  Z2.row(i) = Zco2.transpose();
  W2.row(i) = Wco2.transpose();

  // compute the subseted correlation matrix for coloc Z-scores
  for (int k = 0; k < m; k++){
   tmp_rho(i) = RHO(snps(0 , i)-1, snps(1 , i)-1);
   tmp_rho_21(2*k+1) = tmp_rho(i);
   tmp_rho_22(2*k) = tmp_rho(i);
  }

  for (int k = 0; k < m; k++){
   tmp_rho_2 = sigZ2.col(2*k);
   sigZ2.col(2*k) = tmp_rho_2.cwiseProduct(tmp_rho_21);
   tmp_rho_2 = sigZ2.col(2*k + 1);
   sigZ2.col(2*k + 1) = tmp_rho_2.cwiseProduct(tmp_rho_22);
  }

  Wdiag2 = W2.row(i).asDiagonal();
  adjZ2 = Wdiag2*sigZ2*Wdiag2;

  MatrixXd ZL2(sigZ2.llt().matrixL().solve(I2));
  MatrixXd ZAL2((sigZ2 + adjZ2).llt().matrixL().solve(I2));
  SIGINV2 = 0.5*(ZAL2.transpose()*ZAL2)*adjZ2*(ZL2.transpose()*ZL2);
  out2(i) = -log(ZL2.diagonal().prod()/ZAL2.diagonal().prod()) + ((Z2.row(i))*(SIGINV2)*Z2.row(i).transpose());
  sigZ2 = sigZ2_new;

 }

 VectorXd out1_co(3);
 VectorXd out2_co(3);
 ArrayXd out_exp1(Q);
 ArrayXd out_exp2(ncol);

 MatrixXd::Index maxIndex;

 out1_co(0) = out.maxCoeff(&maxIndex);
 out_exp1 = (out - out1_co(0)*VectorXd::Ones(Q)).array().exp();
 out1_co(1) = out_exp1.sum();
 out1_co(2) = maxIndex + 1;

 out2_co(0) = out2.maxCoeff(&maxIndex);
 out_exp2 = (out2 - (out2_co(0)*VectorXd::Ones(ncol))).array().exp();
 out2_co(1) = out_exp2.sum();
 out2_co(2) = maxIndex + 1;

 if (trtCo.size() == m){
  return List::create(Named("out_co_1") = out1_co, Named("out_co_2") = out2_co, Named("out_snp_scores_1") = out_exp1/out1_co(1), Named("out_snp_scores_2") = out_exp2/out2_co(1));
 }
 else {
  return List::create(Named("out_co_1") = out1_co, Named("out_co_2") = out2_co);
 }
}
