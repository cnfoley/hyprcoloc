##########################################################
##### Alignment #####
##########################################################

#' align.ABF.1
#'
#' @param Z matrix of Z-scores
#' @param W ratio matrix of prior standard deviation and observed standard errors
#' @param trait.cor correlation matrix between traits
#' @param sample.overlap matrix of sample overlap between traits
#' @param ld.matrix LD matrix
#' @param epsilon tolerance parameter
#' @param reg.res regional result
#' @param align.thresh alignment probability threshold
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.2 1 - prior probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param cor.adj.priors correlation adjusted priors
#' @param uniform.priors uniform priors
#' @param Zsq matrix of Z-scores squared
#' @param Wsq matrix of W squared
#' @param ind.traits are the traits independent or to be treated as independent
#' @export
align.ABF.1 <- function(Z, W, trait.cor, sample.overlap, ld.matrix,  epsilon, reg.res, align.thresh, prior.1, prior.2, cor.adj.priors, uniform.priors, Zsq, Wsq, ind.traits){
  
  m = dim(Z)[2];
  Q = dim(Z)[1];
  traits = c(1:m);
  trait.cor =trait.cor*sample.overlap;
  
  if(uniform.priors == T){
    I.unif = 1;
  }else{
    I.unif = 0;
  }
  
  p.1 = prior.1;
  p.m.1 = p.1/(m*Q*(Q-1));
  if(m==2){
    p.m.1 = 2*p.m.1;
  }
  
  prior.3 = prior.1;
  kappa = 2;
  if(cor.adj.priors==T){
    ave.cor = trait.cor[lower.tri(trait.cor)];
    ave.cor = mean(abs(ave.cor));	
    prior.2 = min(prior.2, 1-ave.cor^2);
    prior.3 = prior.1*(10^(-kappa*ave.cor^2));
    p.m.1 = p.m.1*(1-ave.cor)^kappa;
  }
  
  max.ABF.overall = sum.max.ABF.co = mpfr(0,120);
  log.sum.ABF.all.combs = sum.ABF = sum.ABF.co = mpfr(0,120);
  NaNs = 0;  
  
  df = data.frame(matrix(vector(), 0,3, dimnames=list(c(), c("log.sum.ABF.all.combs", "Trait.no.clc", "NaNs"))), stringsAsFactors=F)
  j = 1;
  align.prob = 1;
  trt.combn = sample(c(1:m));
  prior.algn= I.unif*p.m.1 + (1-I.unif)*prior.3*prior(prior.1, prior.2, k = m-1);
  
  if(m>2){
    while(j <=m){
      trt.no.clc = trt.combn[j] + 0.0;
      trt.clc = traits[-trt.no.clc] + 0.0;
      if(ind.traits==T){
        output = align1ind(Zsq, Wsq, trt.clc, trt.no.clc);
      }else{
        output = align1(Z, W, 1, trt.clc, trt.no.clc, trait.cor, ld.matrix, epsilon);
      }
      max.ABF = prior.algn*exp(mpfr(output[[1]][1], 120));
      sum.ABF.tmp = max.ABF*output[[2]][1];
      if(!is.nan(sum.ABF.tmp)){
        sum.ABF= sum.ABF + sum.ABF.tmp;
        if(max.ABF.overall<sum.ABF.tmp){
          max.ABF.overall = sum.ABF.tmp;  
          df[1,]$Trait.no.clc = toString(trt.no.clc); 
        }
      }else{NaNs = NaNs + 1}
      align.prob = (reg.res)/(reg.res + sum.ABF);
      j = j+1;
    }
  }else{
    trt.no.clc = 1.0;
    trt.clc = traits[-trt.no.clc] + 0.0;
    if(ind.traits == TRUE){
      output = align1ind(Zsq, Wsq, trt.clc, trt.no.clc);
    }else{
      output = align1(Z, W, 1, trt.clc, trt.no.clc, trait.cor, ld.matrix, epsilon);
    }
    max.ABF = prior.algn*exp(mpfr(output[[1]][1], 120));
    sum.ABF.tmp = max.ABF*output[[2]][1];
    if(!is.nan(sum.ABF.tmp)){
      sum.ABF= (sum.ABF + sum.ABF.tmp)/2;
      df[1,]$Trait.no.clc = toString(trt.no.clc); 
    }else{NaNs = NaNs + 1}
  }
  
  log.sum.ABF.all.combs=asNumeric(log(sum.ABF)); 
  df[1,]$log.sum.ABF.all.combs=log.sum.ABF.all.combs;
  df[1,]$NaNs = NaNs;
  
  return(df)
  
}
