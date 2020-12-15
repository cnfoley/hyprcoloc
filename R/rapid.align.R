##########################################################
##### Alignment (rapid) #####
##########################################################

#' rapid.align
#'
#' @param Zsq matrix of Z-scores
#' @param Wsq ratio matrix of prior standard deviation and observed standard errors squared
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.2 1 - prior probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param unifrom.priors uniform priors
#' @export
rapid.align <- function(Zsq, Wsq, prior.1, prior.2, uniform.priors){
  
  m = dim(Zsq)[2];
  Q = dim(Zsq)[1];
  
  if(uniform.priors==T){
    I.unif = 1;
  }else{
    I.unif = 0;
  }
  
  p.1 = prior.1;
  p.1.m = p.1/Q;
  p.m.1 = p.1/(m*Q*(Q-1));
  if(m==2){
    p.m.1 = 2*p.m.1;
    cnst = 0.25;	
  }else{
    cnst = 1
  }
  
  prior.all = I.unif*p.1.m + (1-I.unif)*prior(prior.1, prior.2, k = m);
  prior.align= I.unif*p.m.1 + (1-I.unif)*prior.1*prior(prior.1, prior.2, k = m-1);
  
  labf = 0.5*(log(Wsq) + (Zsq)*(1- Wsq));
  sum.labf = rowSums(labf);
  mx.labf = which.max(labf);
  clc.max = which.max(sum.labf);
  
  chi = sum.labf - sum.labf[clc.max];
  sum.exp.chi = sum(exp(chi));
  col.max.labf = colMax(labf);
  exp.labf.std = exp(t(t(labf) - col.max.labf));
  labf.tmp = t(colSums(exp.labf.std) - t(exp.labf.std));
  sum.tmp = colSums(exp(t(t(chi - labf)+col.max.labf))*labf.tmp);
  rm.trt = which.max(sum.tmp);
  
  align.prob = sum.exp.chi/(sum.exp.chi + cnst*(prior.align/prior.all)*sum(sum.tmp));
  
  return(list(align.prob, rm.trt))
  
}
