##########################################################
##### Alignment 1CV and 2CV #####
########################################################## 

#' align.ABF.2
#'
#' @param Z matrix of Z-scores
#' @param W ratio matrix of prior standard deviation and observed standard errors
#' @param snps.clc SNPs colocalisation 
#' @param trait.cor correlation matrix between traits
#' @param sample.overlap matrix of sample overlap between traits
#' @param ld.matrix LD matrix
#' @param epsilon tolerance parameter
#' @param reg.res regional result
#' @param align.thresh alignment probability threshold
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.2 1 - prior probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param prior.3 prior probability that a trait contains a second causal variant given it contains one already
#' @param prior.4 1 - prior probability that trait two co-localises with trait one given traits one and two already share a causal variant and trait one contains a second causal variant
#' @param cor.adj.priors correlation adjusted priors
#' @param unifrom.priors uniform priors
#' @export
align.ABF.2 <- function(Z, W, snps.clc, trait.cor, sample.overlap, ld.matrix, epsilon, reg.res, align.thresh, prior.1, prior.2, prior.3, prior.4, cor.adj.priors, uniform.priors){
  
  m = dim(Z)[2];
  Q = dim(Z)[1];
  traits = c(1:m);
  n.cv = dim(snps.clc)[1];
  trait.cor =trait.cor*sample.overlap;
  
  if(uniform.priors == T){
    I.unif = 1;
  }else{
    I.unif = 0;
  }
  
  p.1 = prior.1;
  p.m.1 = p.1/(m*Q*(Q-1));
  p.1.m.2 = 2*p.1/(m*Q*(Q-1)*(Q-2));
  p.2.m.1 = 2*p.1/(m*Q*(Q-1)*(Q-2));
  p.2.m.2 = 4*p.1/(m*Q*(Q-1)*(Q-2)*(Q-3));
  p.co.m.1.2 = p.1/(m*Q*(Q-1));
  p.co.m.2.1 = p.1/(m*Q*(Q-1));
  p.co.m.2.2 = p.1/(m*Q*(Q-1)*(Q-2));
  
  if(m==2){
    p.m.1 = 2*p.m.1;
    p.2.m.2 = 2*p.2.m.2
    p.co.m.1.2 = p.co.m.1.2/2;
    p.co.m.2.1 = p.co.m.2.1/2;
    p.co.m.2.2 = 2*p.co.m.2.2
  }
  
  prior.5 = prior.1
  kappa = 10;
  
  if(cor.adj.priors==T){
    ave.cor = trait.cor[lower.tri(trait.cor)];
    ave.cor = mean(abs(ave.cor));	
    prior.2 = min(prior.2, 1-ave.cor^2);
    prior.3 = prior.3*(10^(-kappa*ave.cor^2));
    prior.4 = min(prior.4, 1-ave.cor^4);
    prior.5 = prior.1*(10^(-kappa*ave.cor^2));
    p.m.1 = p.m.1*(1-ave.cor)^2;
    p.1.m.2 = p.1.m.2*(1-ave.cor)^4;
    p.2.m.1 = p.2.m.1/(1-ave.cor)^2;
    p.2.m.2 = p.2.m.2*(1-ave.cor)^4;
    p.co.m.1.2 = p.co.m.1.2*(1-ave.cor)^2 
    p.co.m.2.1 = p.co.m.2.1*(1-ave.cor)^2
    p.co.m.2.2 = p.co.m.2.2*(1-ave.cor)^2  
  }
  
  ones=matrix(1,nrow=dim(snps.clc)[1], ncol=dim(snps.clc)[1]);
  trait.cor = as.matrix(kronecker(trait.cor, ones));
  
  max.ABF.overall = sum.max.ABF.co = mpfr(0,120);
  log.sum.ABF.align = sum.ABF = sum.ABF.co = mpfr(0,120);
  NAs = 0;
  
  df = data.frame(matrix(vector(), 0,7, dimnames=list(c(), c("log.sum.ABF.co.loc", "log.sum.ABF.align", "log.max.SNP.ABF.full", "Trait.no.clc", "Traits.clc", "SNPs.clc", "NaNs"))), stringsAsFactors=F)
  j = 1;
  align.prob = 1;
  trt.combn = sample(c(1:m));
  prior.no.11 = I.unif*p.m.1 + (1-I.unif)*prior.5*prior( prior.1, prior.2, k = m-1);
  prior.no.12 = I.unif*p.1.m.2 + (1-I.unif)*(prior.5*prior.3)*prior(prior.1, prior.2, k = m-1);
  prior.co.12 = I.unif*p.co.m.1.2 + (1-I.unif)*prior.3*prior(prior.1, prior.2, k = m);
  prior.no.21 = I.unif*p.2.m.1 + (1-I.unif)*prior.5*prior( prior.1, prior.2, k = m-1)*prior(prior.3, prior.4, k = m-1);
  prior.no.22 = I.unif*p.2.m.2 + (1-I.unif)*(prior.5*prior.3)*prior( prior.1, prior.2, k = m-1)*prior(prior.3, prior.4, k = m-1);
  prior.co.21 = I.unif*p.co.m.2.1 + (1-I.unif)*prior(prior.1, prior.2, k = m)*prior(prior.3, prior.4, k = m-1);
  prior.co.22 = I.unif*p.co.m.2.2 + (1-I.unif)*prior.3*prior(prior.1, prior.2, k = m)*prior(prior.3, prior.4, k = m-1);
  
  while(align.prob > align.thresh & j <=m){
    trt.no.clc = trt.combn[j] + 0.0;
    trt.clc = traits[-trt.no.clc] + 0.0;
    output_1 = align12(Z, W, snps.clc, trt.clc, trt.no.clc, ld.matrix, trait.cor, epsilon);
    max.ABF = c(prior.no.11*exp(mpfr(output_1[[1]][1], 120)),prior.no.12*exp(mpfr(output_1[[2]][1], 120)));
    sum.ABF.no = sum(max.ABF*c(output_1[[1]][2], output_1[[2]][2]));
    if(!is.nan(sum.ABF.no)){
      sum.ABF= sum.ABF + sum.ABF.no;
      if(max.ABF.overall<sum.ABF.no){
        max.ABF.overall = sum.ABF.no;  
        df[1,]$Trait.no.clc = toString(trt.no.clc); 
      }
    }else{NAs = NAs + 1}
    max.ABF.1 = prior.co.12*exp(mpfr(output_1[[3]][1], 120));
    sum.ABF.1 = max.ABF.1*output_1[[3]][2];
    clc.max.cvs.1 = output_1[[3]][3];
    if(!is.na(sum.ABF.1)){
      sum.ABF.co = sum.ABF.co + sum.ABF.1; 
      if(sum.max.ABF.co<sum.ABF.1){
        sum.max.ABF.co = sum.ABF.1; 
        log.max.ABF.co = log(max(max.ABF.1))
        df[1,]$Traits.clc = toString(trt.clc); 
        df[1,]$SNPs.clc = toString(clc.max.cvs.1)
      }
    }else{NAs = NAs + 1;}
    output_2 = align2(Z, W, snps.clc, trt.clc, trt.no.clc, ld.matrix, trait.cor, epsilon);
    max.ABF = c(prior.no.21*exp(mpfr(output_2[[1]][1], 120)),prior.no.22*exp(mpfr(output_2[[2]][1], 120)));
    sum.ABF.no = sum(max.ABF*c(output_2[[1]][2], output_2[[2]][2]));
    if(!is.nan(sum.ABF.no)){
      sum.ABF= sum.ABF + sum.ABF.no;
      if(max.ABF.overall<sum.ABF.no){
        max.ABF.overall = sum.ABF.no;  
        df[1,]$Trait.no.clc = toString(trt.no.clc); 
      }
    }else{NAs = NAs + 1}
    max.ABF.2 = c(prior.co.21*exp(mpfr(output_2[[3]][1], 120)), prior.co.22*exp(mpfr(output_2[[4]][1], 120)));
    sum.ABF.2 = sum(max.ABF.2*c(output_2[[3]][2], output_2[[4]][2]));
    clc.max.cvs.2 = snps.clc[output_2[[2+which(max.ABF.2==max(max.ABF.2))]][3],output_2[[2+which(max.ABF.2==max(max.ABF.2))]][4]];
    if(!is.na(sum.ABF.2)){
      sum.ABF.co = sum.ABF.co + sum.ABF.2; 
      if(sum.max.ABF.co<sum.ABF.2){
        sum.max.ABF.co = sum.ABF.2; 
        log.max.ABF.co = log(max(max.ABF.2));
        df[1,]$Traits.clc = toString(trt.clc); 
        df[1,]$SNPs.clc = toString(clc.max.cvs.2)
      }
    }else{NAs = NAs + 1;}
    align.prob = (reg.res + sum.ABF.co)/(reg.res + sum.ABF.co + sum.ABF);
    j = j+1;
  }
  
  log.sum.ABF.align=asNumeric(log(sum.ABF)); 
  df[1,]$log.sum.ABF.align=log.sum.ABF.align;
  df[1,]$NaNs=NAs;
  log.sum.ABF.co.loc = asNumeric(log(sum.ABF.co));
  df[1,]$log.sum.ABF.co.loc=log.sum.ABF.co.loc;
  df[1,]$log.max.SNP.ABF.full=asNumeric(log.max.ABF.co);
  
  return(df)
  
}
