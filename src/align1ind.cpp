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

List align1ind(NumericMatrix Zmatrix, NumericMatrix Wmatrix, NumericVector traitsCo, NumericVector traitNo){

 // map Z-scores and prior ratios W by co-loc (1) and no co-loc (2)
 const Map<MatrixXd> matZsq(as<Map<MatrixXd> >(Zmatrix));
 const Map<MatrixXd> matWsq(as<Map<MatrixXd> >(Wmatrix));

 // map correlation matrix for traits under consideration
 const Map<VectorXd> trtCo(as<Map<VectorXd> >(traitsCo));
 const Map<VectorXd> trtNo(as<Map<VectorXd> >(traitNo));

 // Set-up ordered Z and W score matrices
 // dimension integers
 int n_trts(trtCo.size()), m(matZsq.cols());
 int Q(matZsq.rows());

 MatrixXd Zsq(Q,n_trts);
 MatrixXd Wsq(Q,n_trts);
 VectorXd I = VectorXd::Ones(n_trts);
 VectorXd Zsq_no(Q); 
 Zsq_no = matZsq.col(trtNo(0)-1);
 VectorXd Wsq_no(Q);
 Wsq_no = matWsq.col(trtNo(0)-1);

 for (int i = 0; i < n_trts; i++){
  Zsq.col(i) = matZsq.col(trtCo(i)-1);
  Wsq.col(i) = matWsq.col(trtCo(i)-1);
 }

 VectorXd tmpZ(n_trts);
 VectorXd tmpW(n_trts);
 VectorXd tmp(n_trts);

 int iters(Q*(Q-1));

 MatrixXd Z(iters, m);
 MatrixXd W(iters, m);
 VectorXd Zno(1);
 VectorXd Wno(1);
 VectorXd W2no(1);

 // Define the output matrix "out" containing a row of 
 // determinant prefactor constants and BF exponents
 VectorXd out(iters);
 MatrixXd tmp_all(iters, m);

 int iter_1(0);

 for (int i = 0; i < Q; i++){

  tmpZ = Zsq.row(i);
  tmpW = I.transpose() - Wsq.row(i);
  tmp = Wsq.row(i);

  // 1 coloc CVs and 1 no-coloc CV
  for (int j1 = 0; j1 < Q; j1++){
   if (j1 != i){

    Zno(0) = Zsq_no(j1);
    Wno(0) = 1 - Wsq_no(j1);
    W2no(0) = Wsq_no(j1);

    Z.row(iter_1) << tmpZ.transpose(), Zno.transpose();
    W.row(iter_1) << tmpW.transpose(), Wno.transpose();
    tmp_all.row(iter_1) << tmp.transpose(), W2no.transpose();

    // compute the subseted correlation matrix for co and no-co Z-scores
    out(iter_1) = 0.5*log(tmp_all.row(iter_1).prod()) + 0.5*(Z.row(iter_1).dot(W.row(iter_1)));

    iter_1 = iter_1 + 1;

   }
  }
 }

 VectorXd out1_no(2);

 out1_no(0) = out.maxCoeff();
 out1_no(1) = (out - out1_no(0)*VectorXd::Ones(Q*(Q-1))).array().exp().sum();

 return wrap(out1_no);

}
