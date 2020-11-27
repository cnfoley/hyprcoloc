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

List align2(NumericMatrix Zmatrix, NumericMatrix Wmatrix, NumericMatrix sTemp1, NumericVector traitsCo, NumericVector traitNo, NumericMatrix rho, NumericMatrix trait_cor, NumericVector EPS){

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
 int sigZdim(n_cv*(trtCo.size()+1));

 MatrixXd sigZtmp1(TRAIT_COR.rows(),sigZdim-1);
 MatrixXd sigZ1(sigZdim-1, sigZdim-1);
 MatrixXd sigZtmp2(TRAIT_COR.rows(), sigZdim);
 MatrixXd sigZ2(sigZdim, sigZdim);

 int Q(matZ.rows()), tmp(1);
 VectorXd tmp_rho(ncol);
 VectorXd tmp_rho_1(sigZdim-1);
 VectorXd tmp_rho_11 = (VectorXd::Ones(sigZdim-1));
 VectorXd tmp_rho_12 = (VectorXd::Ones(sigZdim-1));

 VectorXd tmp_rho_2(sigZdim);
 VectorXd tmp_rho_21 = (VectorXd::Ones(sigZdim));
 VectorXd tmp_rho_22 = (VectorXd::Ones(sigZdim));

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
  sigZtmp1.col(2*i) = tmp_tet;
  sigZtmp1.col(2*i+1) = tmp_tet;
  sigZtmp2.col(2*i) = tmp_tet;
  sigZtmp2.col(2*i+1) = tmp_tet; 
 }

 for (int k = 0; k < n_cv; k++){
  sigZtmp2.col(2*trtCo.size()+k) = TRAIT_COR.col(2*trtNo(0)-1);
 }	

 sigZtmp1.col(2*trtCo.size()) = TRAIT_COR.col(2*trtNo(0)-1);

 for (int j = 0; j < trtCo.size(); j++){
  tmp_tet1 = sigZtmp1.row(2*trtCo(j)-1);
  sigZ1.row(2*j) = tmp_tet1;
  sigZ1.row(2*j+1) = tmp_tet1;
  tmp_tet2 = sigZtmp2.row(2*trtCo(j)-1);
  sigZ2.row(2*j) = tmp_tet2;
  sigZ2.row(2*j+1) = tmp_tet2;
 }

 for (int l = 0; l < n_cv; l++){
  sigZ2.row(2*trtCo.size()+l) = sigZtmp2.row(2*trtNo(0)-1);
 }

 sigZ1.row(2*trtCo.size()) = sigZtmp1.row(2*trtNo(0)-1);
 sigZ1 = (sigZ1 + Ep1).eval();
 sigZ2 = (sigZ2 + Ep2).eval();

 int ncolZco(trtCo.size());
 int sigDim(sigZdim);
 int iters(ncol*(ncol-1));

 MatrixXd Z1(n_cv, ncolZco);
 VectorXd Zno_1(n_cv-1);
 MatrixXd Z_1cv(ncol*Q, sigDim-1);
 VectorXd Zno_2(n_cv);
 MatrixXd Z_2cv(iters, sigDim);

 MatrixXd W1(n_cv, ncolZco);
 VectorXd Wno_1(n_cv-1);
 MatrixXd W_1cv(ncol*Q, sigDim-1);
 VectorXd Wno_2(n_cv);
 MatrixXd W_2cv(iters, sigDim);

 VectorXd snp_cor_1cv(sigZdim - 1);
 VectorXd tet_1cv = sigZ1.row(sigZdim-n_cv);

 VectorXd snp_cor_2cv(sigZdim);
 VectorXd snp_cor_2cv_1 = snp_cor_2cv;
 VectorXd tet_2cv = sigZ2.row(sigZdim-1);

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
 int n_no_1( ncol*( Q - 2)), n_co_1( 2*ncol);
 VectorXd out_no_1(n_no_1);
 MatrixXd out_co_1(3, n_co_1);

 int n_no_2( ncol*( ncol - (2*Q - 3))), n_co_2( ncol*(2*(Q-2)));
 VectorXd out_no_2(n_no_2);
 MatrixXd out_co_2(3, n_co_2);

 int iter_no_1(0), iter_co_1(0), iter_no_2(0), iter_co_2(0); 
 int iter_1(0), iter_2(0);

 MatrixXd sigZ1_new = sigZ1;
 MatrixXd sigZ2_new = sigZ2;

 for (int i = 0; i < ncol; i++){

  sigZ1 = sigZ1_new;
  sigZ2 = sigZ2_new;

  for (int k1 = 0; k1 < n_cv; k1++){
   tmp = snps1(k1,i);
   Z1.row(k1) = matZco.row(tmp-1);
   W1.row(k1) = matWco.row(tmp-1);
  }

  VectorXd Zco(Map<VectorXd>(Z1.data(), n_cv*ncolZco));
  VectorXd Wco(Map<VectorXd>(W1.data(), n_cv*ncolZco));

  for (int m = 0; m < ncolZco; m++){
   tmp_rho(i) = RHO(snps1(1 , i)-1, snps1(0 , i)-1);
   tmp_rho_11(2*m+1) = tmp_rho(i);
   tmp_rho_12(2*m) = tmp_rho(i);
   tmp_rho_21(2*m+1) = tmp_rho(i);
   tmp_rho_22(2*m) = tmp_rho(i);
  }

  for (int m = 0; m < ncolZco; m++){
   tmp_rho_1 = sigZ1.col(2*m);
   sigZ1.col(2*m) = tmp_rho_1.cwiseProduct(tmp_rho_11);
   tmp_rho_1 = sigZ1.col(2*m + 1);
   sigZ1.col(2*m + 1) = tmp_rho_1.cwiseProduct(tmp_rho_12);
   tmp_rho_2 = sigZ2.col(2*m);
   sigZ2.col(2*m) = tmp_rho_2.cwiseProduct(tmp_rho_21);
   tmp_rho_2 = sigZ2.col(2*m + 1);
   sigZ2.col(2*m + 1) = tmp_rho_2.cwiseProduct(tmp_rho_22);
  }

  for (int j2 = 0; j2 < ncol; j2++){

   // 2 coloc CVs and 1 no-coloc CV
   if (j2 < Q){

    tmp = j2;
    Zno_1(0) = vecZno(tmp);
    Wno_1(0) = vecWno(tmp);
    Z_1cv.row(iter_1) << Zco.transpose(), Zno_1.transpose();
    W_1cv.row(iter_1) << Wco.transpose(), Wno_1.transpose();

    // compute the subseted correlation matrix for co and no-co Z-scores
    Wdiag_1 = W_1cv.row(iter_1).asDiagonal();

    for (int m = 0; m < ncolZco; m++){
     snp_cor_1cv(2*m) = RHO(snps1(0,i)-1, j2)*tet_1cv(2*m);
     snp_cor_1cv(2*m+1) = RHO(snps1(1,i)-1, j2)*tet_1cv(2*m+1);
    }

    snp_cor_1cv(sigZdim - n_cv) = sigZ1_new(sigZdim - n_cv,sigZdim - n_cv);
    sigZ1.row(sigZdim - n_cv) = snp_cor_1cv;
    sigZ1.col(sigZdim - n_cv) = snp_cor_1cv;

    adjZ_1 = Wdiag_1*sigZ1*Wdiag_1;

    MatrixXd ZL1(sigZ1.llt().matrixL().solve(I1));
    MatrixXd ZAL1((sigZ1 + adjZ_1).llt().matrixL().solve(I1));

    SIGINV_1 = 0.5*(ZAL1.transpose()*ZAL1)*adjZ_1*(ZL1.transpose()*ZL1);

    if (snps1.col(i)(0) == tmp+1){
     out_co_1(0,iter_co_1) = -log(ZL1.diagonal().prod()/ZAL1.diagonal().prod()) + ((Z_1cv.row(iter_1))*(SIGINV_1)*Z_1cv.row(iter_1).transpose());
     out_co_1(1, iter_co_1) = 1;
     out_co_1(2,iter_co_1) = i+1;
     iter_co_1 = iter_co_1 + 1;
    }
    else if (snps1.col(i)(1) == tmp+1){
     out_co_1(0,iter_co_1) = -log(ZL1.diagonal().prod()/ZAL1.diagonal().prod()) + ((Z_1cv.row(iter_1))*(SIGINV_1)*Z_1cv.row(iter_1).transpose());
     out_co_1(1, iter_co_1) = 2;
     out_co_1(2,iter_co_1) = i+1;
     iter_co_1 = iter_co_1 + 1;
    }
    else {
     out_no_1(iter_no_1) = -log(ZL1.diagonal().prod()/ZAL1.diagonal().prod()) + ((Z_1cv.row(iter_1))*(SIGINV_1)*Z_1cv.row(iter_1).transpose());
     iter_no_1 = iter_no_1 + 1;
    }

    iter_1 = iter_1 + 1;

   }

   // 2 coloc CVs and 2 no-coloc CVs
   if (j2 != i){

    for (int k2 = 0; k2 < n_cv; k2++){
     tmp = snps1(k2,j2);
     Zno_2(k2) = vecZno(tmp-1);
     Wno_2(k2) = vecWno(tmp-1);
    }

    Z_2cv.row(iter_2) << Zco.transpose(), Zno_2.transpose();
    W_2cv.row(iter_2) << Wco.transpose(), Wno_2.transpose();

    // compute the subseted correlation matrix for co and no-co Z-scores
    Wdiag_2 = W_2cv.row(iter_2).asDiagonal();

    for (int m = 0; m < ncolZco; m++){
     snp_cor_2cv(2*m) = RHO(snps1(0,i)-1, snps1(0,j2)-1)*tet_2cv(2*m);
     snp_cor_2cv(2*m+1) = RHO(snps1(1,i)-1, snps1(0,j2)-1)*tet_2cv(2*m+1);
     snp_cor_2cv_1(2*m) = RHO(snps1(0,i)-1, snps1(1,j2)-1)*tet_2cv(2*m);
     snp_cor_2cv_1(2*m+1) = RHO(snps1(1,i)-1, snps1(1,j2)-1)*tet_2cv(2*m+1);
    }

    snp_cor_2cv(sigZdim - n_cv) = sigZ2_new(sigZdim - n_cv, sigZdim - n_cv);
    snp_cor_2cv_1(sigZdim - 1) = sigZ2_new(sigZdim - 1, sigZdim - 1);
    sigZ2.row(sigZdim - n_cv) = snp_cor_2cv;
    sigZ2.col(sigZdim - n_cv) = snp_cor_2cv;
    sigZ2.row(sigZdim - 1) = snp_cor_2cv_1;
    sigZ2.col(sigZdim - 1) = snp_cor_2cv_1;
    sigZ2(2*ncolZco+1, 2*ncolZco) = RHO(snps1(1 , j2)-1, snps1(0 , j2)-1);
    sigZ2(2*ncolZco, 2*ncolZco+1) = RHO(snps1(1 , j2)-1, snps1(0 , j2)-1);

    adjZ_2 = Wdiag_2*sigZ2*Wdiag_2;

    MatrixXd ZL2(sigZ2.llt().matrixL().solve(I2));
    MatrixXd ZAL2((sigZ2 + adjZ_2).llt().matrixL().solve(I2));

    SIGINV_2 = 0.5*(ZAL2.transpose()*ZAL2)*adjZ_2*(ZL2.transpose()*ZL2);

    if (snps1.col(i)(0) == snps1.col(j2)(0) || snps1.col(i)(0) == snps1.col(j2)(1)){
     out_co_2(0,iter_co_2) = -log(ZL2.diagonal().prod()/ZAL2.diagonal().prod()) + ((Z_2cv.row(iter_2))*(SIGINV_2)*Z_2cv.row(iter_2).transpose());
     out_co_2(1, iter_co_2) = 1;
     out_co_2(2,iter_co_2) = i+1;
     iter_co_2 = iter_co_2 + 1;
    }
    else if (snps1.col(i)(1) == snps1.col(j2)(0) || snps1.col(i)(1) == snps1.col(j2)(1)){
     out_co_2(0,iter_co_2) = -log(ZL2.diagonal().prod()/ZAL2.diagonal().prod()) + ((Z_2cv.row(iter_2))*(SIGINV_2)*Z_2cv.row(iter_2).transpose());
     out_co_2(1, iter_co_2) = 2;
     out_co_2(2,iter_co_2) = i+1;
     iter_co_2 = iter_co_2 + 1;
    }
    else {
     out_no_2(iter_no_2) = -log(ZL2.diagonal().prod()/ZAL2.diagonal().prod()) + ((Z_2cv.row(iter_2))*(SIGINV_2)*Z_2cv.row(iter_2).transpose());
     iter_no_2 = iter_no_2 + 1;
    }

    iter_2 = iter_2 +1;

   }

  }

 } 

 VectorXd out1_no(2);
 VectorXd out1_co(4);
 VectorXd out2_no(2);
 VectorXd out2_co(4);

 MatrixXd::Index maxIndex;

 out1_no(0) = out_no_1.maxCoeff();
 out1_no(1) = (out_no_1 - out1_no(0)*VectorXd::Ones(iter_no_1)).array().exp().sum();
 out1_co(0) = out_co_1.row(0).maxCoeff(&maxIndex);
 out1_co(1) = ((out_co_1.row(0).transpose()) - out1_co(0)*VectorXd::Ones(iter_co_1)).array().exp().sum();
 out1_co(2) = out_co_1(1,maxIndex);
 out1_co(3) = out_co_1(2,maxIndex);

 out2_no(0) = out_no_2.maxCoeff();
 out2_no(1) = (out_no_2 - out2_no(0)*VectorXd::Ones(iter_no_2)).array().exp().sum();
 out2_co(0) = out_co_2.row(0).maxCoeff(&maxIndex);
 out2_co(1) = ((out_co_2.row(0).transpose()) - (out2_co(0)*VectorXd::Ones(iter_co_2))).array().exp().sum();
 out2_co(2) = out_co_2(1,maxIndex);
 out2_co(3) = out_co_2(2,maxIndex);

 return List::create(Named("out_no_1") = out1_no, Named("out_no_2") = out2_no, Named("out_co_1") = out1_co, Named("out_co_2") = out2_co);

}
