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

List align1(NumericMatrix Zmatrix, NumericMatrix Wmatrix, NumericVector N_CV, NumericVector traitsCo, NumericVector traitNo, NumericMatrix trait_cor, NumericMatrix rho, NumericVector EPS){

 // map Z-scores and prior ratios W by co-loc (1) and no co-loc (2)
 const Map<MatrixXd> matZ(as<Map<MatrixXd> >(Zmatrix));
 const Map<MatrixXd> matW(as<Map<MatrixXd> >(Wmatrix));

 // map SNPs by co-loc (1) and no co-loc (2)
 const Map<VectorXd> ncv(as<Map<VectorXd> >(N_CV));

 // map trait correlation matrix 
 const Map<MatrixXd> TRAIT_COR(as<Map<MatrixXd> >(trait_cor));

 // map SNP correlation matrix 
 const Map<MatrixXd> RHO(as<Map<MatrixXd> >(rho));

 // map correlation matrix for traits under consideration
 const Map<VectorXd> trtCo(as<Map<VectorXd> >(traitsCo));
 const Map<VectorXd> trtNo(as<Map<VectorXd> >(traitNo));
 const Map<VectorXd> eps(as<Map<VectorXd> >(EPS));

 // Generate the constants that will define the looping parameters across Z and W. 
 // NOTE row space of SNP
 int n_cv(ncv(0));

 // Set-up ordered Z and W score matrices
 int sigZdim(trtCo.size()+n_cv);
 MatrixXd sigZtmp(TRAIT_COR.rows(),sigZdim);
 MatrixXd sigZ(sigZdim, sigZdim);

 int Q(matZ.rows());
 MatrixXd matZco(Q, trtCo.size());
 MatrixXd matWco(Q, trtCo.size());
 VectorXd vecZno(Q);
 VectorXd vecWno(Q);

 int iters(Q*(Q-1));

 MatrixXd Ep = eps(0)*(VectorXd::Ones(sigZdim)).asDiagonal();
 MatrixXd I = (VectorXd::Ones(sigZdim)).asDiagonal();

 for (int i = 0; i < trtCo.size(); i++){
  matZco.col(i) = matZ.col(trtCo(i)-1);
  matWco.col(i) = matW.col(trtCo(i)-1);
 }

 vecZno = matZ.col(trtNo(0)-1);
 vecWno = matW.col(trtNo(0)-1);

 for (int i = 0; i < trtCo.size(); i++){
  sigZtmp.col(i) = TRAIT_COR.col(trtCo(i)-1); 
 }

 sigZtmp.col(trtCo.size()) = TRAIT_COR.col(trtNo(0)-1); 

 for (int i = 0; i < trtCo.size(); i++){
  sigZ.row(i) = sigZtmp.row(trtCo(i)-1);
 }

 sigZ.row(trtCo.size()) = sigZtmp.row(trtNo(0)-1);
 sigZ = (sigZ + Ep).eval();

 int ncolZco(matZco.cols());
 int sigDim(sigZdim);

 MatrixXd Z1(n_cv, ncolZco);
 VectorXd Zno(n_cv);
 MatrixXd Z(iters, sigDim);

 MatrixXd W1(n_cv, ncolZco);
 VectorXd Wno(n_cv);
 MatrixXd W(iters, sigDim);

 VectorXd snp_cor(sigZdim);
 VectorXd tet = sigZ.row(sigZdim-1);

 // Define the diagonal prior ratio matrix "Wdiag"
 MatrixXd Wdiag(sigDim, sigDim);

 // Define the adjusted prior matrix "adjZ"
 MatrixXd adjZ(sigDim, sigDim);

 // Define the inverse sigma matrix "SIGINV"
 MatrixXd SIGINV(sigDim, sigDim);

 // Define the output matrix "out" containing a row of 
 // determinant prefactor constants and BF exponents
 VectorXd out(Q*(Q-1));

 int iter_1(0);

 for (int i = 0; i < Q; i++){

  Z1.row(0) = matZco.row(i);
  W1.row(0) = matWco.row(i);

  VectorXd Zco(Map<VectorXd>(Z1.data(), n_cv*ncolZco));
  VectorXd Wco(Map<VectorXd>(W1.data(), n_cv*ncolZco));

  // 1 coloc CVs and 1 no-coloc CV
  for (int j1 = 0; j1 < Q; j1++){
   if (j1 != i){

    Zno(0) = vecZno(j1);
    Wno(0) = vecWno(j1);

    Z.row(iter_1) << Zco.transpose(), Zno.transpose();
    W.row(iter_1) << Wco.transpose(), Wno.transpose();

    // compute the subseted correlation matrix for co and no-co Z-scores
    Wdiag = W.row(iter_1).asDiagonal();
    snp_cor = RHO(i, j1)*tet;
    snp_cor(sigZdim-1) = 1;
    sigZ.row(sigZdim-1) = snp_cor;
    sigZ.col(sigZdim-1) = snp_cor;
    adjZ = Wdiag*sigZ*Wdiag;

    MatrixXd ZL(sigZ.llt().matrixL().solve(I));
    MatrixXd ZAL((sigZ + adjZ).llt().matrixL().solve(I));

    SIGINV = 0.5*(ZAL.transpose()*ZAL)*adjZ*(ZL.transpose()*ZL);
    out(iter_1) = -log(ZL.diagonal().prod()/ZAL.diagonal().prod()) + ((Z.row(iter_1))*(SIGINV)*Z.row(iter_1).transpose());

    iter_1 = iter_1 + 1;

   }
  }
 }

 VectorXd out1_no(2);

 out1_no(0) = out.maxCoeff();
 out1_no(1) = (out - out1_no(0)*VectorXd::Ones(Q*(Q-1))).array().exp().sum();

 return wrap(out1_no);
}
