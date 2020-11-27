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

List regional1ind(NumericMatrix zTemp1, NumericMatrix wTemp1, NumericVector traitsCo){

 // map Z-scores and prior ratios W by co-loc (1) and no co-loc (2)
 const Map<MatrixXd> matZsq(as<Map<MatrixXd> >(zTemp1));
 const Map<MatrixXd> matWsq(as<Map<MatrixXd> >(wTemp1));

 // co-localised traits
 const Map<VectorXd> trtCo(as<Map<VectorXd> >(traitsCo));

 // dimension integers
 int n_trts(trtCo.size()), m(matZsq.cols());
 int Q(matZsq.rows());

 MatrixXd Zsq(Q,n_trts);
 MatrixXd Wsq(Q,n_trts);
 VectorXd I = VectorXd::Ones(n_trts);

 for (int i = 0; i < n_trts; i++){
  Zsq.col(i) = matZsq.col(trtCo(i)-1);
  Wsq.col(i) = matWsq.col(trtCo(i)-1);
 }

 // Define the output matrix "out" containing a row of 
 // determinant prefactor constants and BF exponents
 VectorXd out(Q);

 VectorXd tmpZ(n_trts);
 VectorXd tmpW(n_trts);
 MatrixXd tmp(1, n_trts);

 for (int i = 0; i < Q; i++){
  tmpZ = Zsq.row(i);
  tmpW = I.transpose() - Wsq.row(i);
  tmp.row(0) = Wsq.row(i);
  out(i) = 0.5*log(tmp.row(0).prod()) + 0.5*(tmpZ.dot(tmpW));
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
