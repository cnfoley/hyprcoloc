//Includes/namespaces
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <Eigen/Core>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace RcppEigen;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::as;
// [[Rcpp::export]]

List align12(NumericMatrix Zmatrix, NumericMatrix Wmatrix, NumericMatrix sTemp1, NumericVector traitsCo, NumericVector traitNo, NumericMatrix rho, NumericMatrix trait_cor, NumericVector EPS){

 // map Z-scores and prior ratios W by co-loc (1) and no co-loc (2)
 const Map<MatrixXd> matZ(as<Map<MatrixXd> >(Zmatrix));
 const Map<MatrixXd> matW(as<Map<MatrixXd> >(Wmatrix));

 // map SNPs by co-loc (1) and no co-loc (2)
 const Map<MatrixXd> snps1(as<Map<MatrixXd> >(sTemp1));

 // map correlation matrix 
 const Map<MatrixXd> RHO(as<Map<MatrixXd> >(rho));

 // map trait correlation matrix 
 const Map<MatrixXd> TRAIT_COR(as<Map<MatrixXd> >(trait_cor));

 // map correlation matrix for traits under consideration
 const Map<VectorXd> trtCo(as<Map<VectorXd> >(traitsCo));
 const Map<VectorXd> trtNo(as<Map<VectorXd> >(traitNo));

 // increase diagonal elements by epsilon to avoid large condition number issues with covariance matrix sigZ 
 const Map<VectorXd> eps(as<Map<VectorXd> >(EPS));

 // Generate the constants that will define the looping parameters across Z and W. 
 // NOTE column space of snps1 == column space of snps2
 // NOTE row space of SNP
 int n_cv(snps1.rows()), ncol(snps1.cols());

 // Set-up ordered Z and W score matrices
 int sigZdim(trtCo.size()+n_cv);

 MatrixXd sigZtmp1(TRAIT_COR.rows(),sigZdim-1);
 MatrixXd sigZ1(sigZdim-1, sigZdim-1);
 MatrixXd sigZtmp2(TRAIT_COR.rows(), sigZdim);
 MatrixXd sigZ2(sigZdim, sigZdim);

 int Q(matZ.rows()), tmp(1);
 MatrixXd matZco(Q, trtCo.size());
 MatrixXd matWco(Q, trtCo.size());
 VectorXd vecZno(Q);
 VectorXd vecWno(Q);
 VectorXd tmp_tet(TRAIT_COR.rows());
 VectorXd tmp_tet1(sigZdim-1);
 VectorXd tmp_tet2(sigZdim);

 MatrixXd Ep1 = eps(0)*(VectorXd::Ones(sigZdim-1)).asDiagonal();
 MatrixXd Ep2 = eps(0)*(VectorXd::Ones(sigZdim)).asDiagonal();

 MatrixXd I1 = (VectorXd::Ones(sigZdim-1)).asDiagonal();
 MatrixXd I2 = (VectorXd::Ones(sigZdim)).asDiagonal();

 for (int i = 0; i < trtCo.size(); i++){
  matZco.col(i) = matZ.col(trtCo(i)-1);
  matWco.col(i) = matW.col(trtCo(i)-1);
 }

 vecZno = matZ.col(trtNo(0)-1);
 vecWno = matW.col(trtNo(0)-1);

 for (int i = 0; i < trtCo.size(); i++){
  tmp_tet = TRAIT_COR.col(2*trtCo(i)-1);
  sigZtmp1.col(i) = tmp_tet;
  sigZtmp2.col(i) = tmp_tet;
 }

 for (int k = 0; k < n_cv; k++){
  sigZtmp2.col(trtCo.size()+k) = TRAIT_COR.col(2*trtNo(0)-1);
 }	

 sigZtmp1.col(trtCo.size()) = TRAIT_COR.col(2*trtNo(0)-1);

 for (int j = 0; j < trtCo.size(); j++){
  tmp_tet1 = sigZtmp1.row(2*trtCo(j)-1);
  sigZ1.row(j) = tmp_tet1;
  tmp_tet2 = sigZtmp2.row(2*trtCo(j)-1);
  sigZ2.row(j) = tmp_tet2;
 }

 for (int l = 0; l < n_cv; l++){
  sigZ2.row(trtCo.size()+l) = sigZtmp2.row(2*trtNo(0)-1);
 }

 sigZ1.row(trtCo.size()) = sigZtmp1.row(2*trtNo(0)-1);
 sigZ1 = (sigZ1 + Ep1).eval();
 sigZ2 = (sigZ2 + Ep2).eval();

 int ncolZco(trtCo.size());
 int sigDim(sigZdim);
 int iters(ncol*Q);

 MatrixXd Z1(n_cv-1, ncolZco);
 VectorXd Zno_1(n_cv-1);
 MatrixXd Z_1cv(Q*(Q-1), sigDim-1);
 VectorXd Zno_2(n_cv);
 MatrixXd Z_2cv(iters, sigDim);

 MatrixXd W1(n_cv-1, ncolZco);
 VectorXd Wno_1(n_cv-1);
 MatrixXd W_1cv(Q*(Q-1), sigDim-1);
 VectorXd Wno_2(n_cv);
 MatrixXd W_2cv(iters, sigDim);

 VectorXd snp_cor_1cv(sigZdim - 1);
 VectorXd tet_1cv = sigZ1.row(sigZdim-n_cv);

 VectorXd snp_cor_2cv(sigZdim);
 VectorXd tet_2cv = sigZ2.row(sigZdim-1);

 int snp1(0), snp2(0); 

 // Define the diagonal prior ratio matrix "Wdiag"
 MatrixXd Wdiag_1(sigDim-1, sigDim-1);
 MatrixXd Wdiag_2(sigDim, sigDim);

 // Define the adjusted prior matrix "adjZ"
 MatrixXd adjZ_1(sigDim-1, sigDim-1);
 MatrixXd adjZ_2(sigDim, sigDim);

 // Define the inverse sigma matrix "SIGINV"
 MatrixXd SIGINV_1(sigDim-1, sigDim-1);
 MatrixXd SIGINV_2(sigDim, sigDim);

 // Define the output matrix "out" containing a row of 
 // determinant prefactor constants and BF exponents
 int n_no_1(Q*(Q - 1));
 VectorXd out_no_1(n_no_1);
 int n_no_2( Q*(ncol - Q + 1)), n_co_2( 2*ncol);
 VectorXd out_no_2(n_no_2);
 MatrixXd out_co_2(2, n_co_2);

 int iter_no_2(0), iter_co_2(0); 
 int iter_1(0), iter_2(0);

 MatrixXd sigZ1_new = sigZ1;
 MatrixXd sigZ2_new = sigZ2;

 for (int i = 0; i < Q; i++){

  tmp = i;
  Z1.row(0) = matZco.row(tmp);
  W1.row(0) = matWco.row(tmp);
  VectorXd Zco(Map<VectorXd>(Z1.data(), ncolZco));
  VectorXd Wco(Map<VectorXd>(W1.data(), ncolZco));

  for (int j2 = 0; j2 < ncol; j2++){
   if (j2 < Q){
    if (j2 != i){

     tmp = j2;
     Zno_1(0) = vecZno(tmp);
     Wno_1(0) = vecWno(tmp);

     Z_1cv.row(iter_1) << Zco.transpose(), Zno_1.transpose();
     W_1cv.row(iter_1) << Wco.transpose(), Wno_1.transpose();

     // compute the subseted correlation matrix for co and no-co Z-scores
     Wdiag_1 = W_1cv.row(iter_1).asDiagonal();
     snp_cor_1cv = RHO(i, j2)*tet_1cv;
     snp_cor_1cv(sigZdim - n_cv) = sigZ1_new(sigZdim - n_cv, sigZdim - n_cv);
     sigZ1.row(sigZdim - n_cv) = snp_cor_1cv;
     sigZ1.col(sigZdim - n_cv) = snp_cor_1cv;

     adjZ_1 = Wdiag_1*sigZ1*Wdiag_1;

     MatrixXd ZL1(sigZ1.llt().matrixL().solve(I1));
     MatrixXd ZAL1((sigZ1 + adjZ_1).llt().matrixL().solve(I1));

     SIGINV_1 = 0.5*(ZAL1.transpose()*ZAL1)*adjZ_1*(ZL1.transpose()*ZL1);
     out_no_1(iter_1) = -log(ZL1.diagonal().prod()/ZAL1.diagonal().prod()) + ((Z_1cv.row(iter_1))*(SIGINV_1)*Z_1cv.row(iter_1).transpose());

     iter_1 = iter_1 + 1;

    }
   }

   // 1 coloc CVs and 2 no-coloc CVs
   for (int k2 = 0; k2 < n_cv; k2++){
    tmp = snps1(k2,j2);
    Zno_2(k2) = vecZno(tmp-1);
    Wno_2(k2) = vecWno(tmp-1);
   }

   Z_2cv.row(iter_2) << Zco.transpose(), Zno_2.transpose();
   W_2cv.row(iter_2) << Wco.transpose(), Wno_2.transpose();

   // compute the subseted correlation matrix for co and no-co Z-scores
   Wdiag_2 = W_2cv.row(iter_2).asDiagonal();

   snp1 = snps1(0 , j2);
   snp2 = snps1(1 , j2);

   snp_cor_2cv = RHO(i, snp1 - 1)*tet_2cv;
   snp_cor_2cv(sigZdim - n_cv) = sigZ2_new(sigZdim - n_cv, sigZdim - n_cv);
   sigZ2.row(sigZdim - n_cv) = snp_cor_2cv;
   sigZ2.col(sigZdim - n_cv) = snp_cor_2cv;

   snp_cor_2cv = RHO(i, snp2 - 1)*tet_2cv;
   snp_cor_2cv(sigZdim - n_cv + 1) = sigZ2_new(sigZdim - n_cv + 1, sigZdim - n_cv + 1);
   sigZ2.row(sigZdim - n_cv + 1) = snp_cor_2cv;
   sigZ2.col(sigZdim - n_cv + 1) = snp_cor_2cv;

   sigZ2(ncolZco+1, ncolZco) = RHO(snp1 - 1, snp2 - 1);
   sigZ2(ncolZco, ncolZco+1) = RHO(snp1 - 1, snp2 - 1);

   adjZ_2 = Wdiag_2*sigZ2*Wdiag_2;

   MatrixXd ZL2(sigZ2.llt().matrixL().solve(I2));
   MatrixXd ZAL2((sigZ2 + adjZ_2).llt().matrixL().solve(I2));

   SIGINV_2 = 0.5*(ZAL2.transpose()*ZAL2)*adjZ_2*(ZL2.transpose()*ZL2);

   if (i+1 == snp1 || i+1 == snp2){
    out_co_2(0,iter_co_2) = -log(ZL2.diagonal().prod()/ZAL2.diagonal().prod()) + ((Z_2cv.row(iter_2))*(SIGINV_2)*Z_2cv.row(iter_2).transpose());
    out_co_2(1,iter_co_2) = i+1;
    iter_co_2 = iter_co_2 + 1;
   }
   else {
    out_no_2(iter_no_2) = -log(ZL2.diagonal().prod()/ZAL2.diagonal().prod()) + ((Z_2cv.row(iter_2))*(SIGINV_2)*Z_2cv.row(iter_2).transpose());
    iter_no_2 = iter_no_2 + 1;
   }
   iter_2 = iter_2 +1;

  }
 }

 VectorXd out1_no(2);
 VectorXd out2_no(2);
 VectorXd out2_co(3);

 MatrixXd::Index maxIndex;

 out1_no(0) = out_no_1.maxCoeff();
 out1_no(1) = (out_no_1 - out1_no(0)*VectorXd::Ones(iter_1)).array().exp().sum();
 out2_no(0) = out_no_2.maxCoeff();
 out2_no(1) = (out_no_2 - out2_no(0)*VectorXd::Ones(iter_no_2)).array().exp().sum();
 out2_co(0) = out_co_2.row(0).maxCoeff(&maxIndex);
 out2_co(1) = ((out_co_2.row(0).transpose()) - (out2_co(0)*VectorXd::Ones(iter_co_2))).array().exp().sum();
 out2_co(2) = out_co_2(1,maxIndex);

 return List::create(Named("out_no_1") = out1_no, Named("out_no_2") = out2_no, Named("out_co_2") = out2_co);
}
