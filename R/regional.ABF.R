##########################################################
##### Regional colocalisation #####
##########################################################

#' regional.ABF
#'
#' @param Z matrix of Z-scores
#' @param W ratio matrix of prior standard deviation and observed standard errors
#' @param snps.clc SNPs colocalisation
#' @param rho LD matrix
#' @param trait.cor correlation matrix between traits
#' @param sample.overlap matrix of sample overlap between traits
#' @param epsilon tolerance parameter
#' @param reg.thresh regional probability threshold
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.2 1 - prior probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param prior.3 prior probability that a trait contains a second causal variant given it contains one already
#' @param prior.4 1 - prior probability that trait two co-localises with trait one given traits one and two already share a causal variant and trait one contains a second causal variant
#' @param flag flag variable
#' @param test.2 test for 2CV
#' @param reg.steps regional step paramter
#' @param cor.adj.priors correlation adjusted priors
#' @param unifrom.priors uniform priors
#' @param branch.jump branch jump
#' @param Zsq matrix of Z-scores squared
#' @param Wsq matrix of W squared
#' @param ind.traits are the traits independent or to be treated as independent
#' @export
regional.ABF <- function(Z, W, snps.clc, rho, trait.cor, sample.overlap, epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, test.2, reg.steps, cor.adj.priors, uniform.priors, branch.jump, Zsq, Wsq, ind.traits){
  
  m = dim(Z)[2];
  Q = dim(Z)[1];
  trait.cor =trait.cor*sample.overlap;
  if(reg.steps>m){reg.steps = m;}
  
  if(uniform.priors == T){
    I.unif = 1;
  }else{
    I.unif = 0;
  }
  
  p.1 = prior.1;
  p.1.m = p.1/Q;
  p.2.m = p.1*2/(Q*(Q-1));
  
  if(cor.adj.priors==T){
    ave.cor = trait.cor[lower.tri(trait.cor)];
    ave.cor = mean(abs(ave.cor));	
    prior.2 = min(prior.2, 1-ave.cor^2);
    prior.4 = min(prior.4, 1-ave.cor^4);
  }
  
  if(branch.jump == T & reg.steps > 2){
    jmp = 1;
  }else{
    jmp = 0;
  }
  
  ones=matrix(1,nrow=dim(snps.clc)[1],ncol=dim(snps.clc)[1]);
  trait.cor = as.matrix(kronecker(trait.cor, ones));
  
  log.sum.max.ABF.1 = sum.max.ABF.1 = log.sum.max.ABF.2 = sum.max.ABF.2 = mpfr(0,120);
  log.sum.ABF.1 = sum.ABF.1 = sum.ABF.all.1 = log.sum.ABF.2 = sum.ABF.2 = sum.ABF.all.2 = mpfr(0,120);
  log.all.traits.ABF.1 = log.all.traits.ABF.2 = mpfr(0,120);
  log.max.ABF.1 = log.max.ABF.2 = mpfr(0,120)
  mx.ABF.k = mpfr(0,120);
  NAs = 0;
  
  if(test.2 ==0){
    df = data.frame(matrix(vector(), 0,10, dimnames=list(c(), c("log.sum.ABF.full", "log.sum.ABF.all.combs", "log.max.sum.ABF.all.traits", "log.max.SNP.ABF.full", "Post.ratio.2.vs.1and2",  "Traits", "SNPs", "NaNs", "reg_bb_alg" , "reg_only_trts" ))), stringsAsFactors=F);
    df[1,]$reg_bb_alg = FALSE;
    reg.prob = 1;
    i = m;
    count = 0;
    while(reg.prob > (1-jmp)*reg.thresh & i > 0){
      clc.trt=combn(m,i);
      count = count + dim(clc.trt)[2];
      prior1 = I.unif*p.1.m + (1-I.unif)*prior(prior.1,prior.2, k = i);
      prior2 = I.unif*p.2.m + (1-I.unif)*prior(prior.1,prior.2, k = i)*prior(prior.3, prior.4, k = i);
      for(j in 1:dim(clc.trt)[2]){
        trt.clc = clc.trt[,j] + 0.0;
        if(flag==0){
          if(ind.traits == T){
            output = regional1ind(Zsq, Wsq, trt.clc);
          }else{
            output = regional1(Z, W, trt.clc, trait.cor, epsilon);
          }
          if(i == m){
            max.ABF.1 = prior1*exp(mpfr(output[[1]][1], 120));
            sum.ABF.1 = max.ABF.1*output[[1]][2];
            clc.max.cvs.1 = output[[1]][3];
          }else{
            max.ABF.1 = prior1*exp(mpfr(output[[1]][1], 120));
            sum.ABF.1 = max.ABF.1*output[[2]][1];
            clc.max.cvs.1 = output[[3]][1];
            if(mx.ABF.k < sum.ABF.1){
              mx.ABF.k =  sum.ABF.1;
              df[1,]$reg_only_trts = toString(clc.trt[,j]); 
            }
          }
          if(!is.na(sum.ABF.1)){
            sum.ABF.all.1 = sum.ABF.all.1 + sum.ABF.1; 
            if(sum.max.ABF.1<sum.ABF.1){
              sum.max.ABF.1 = sum.ABF.1; 
              df[1,]$Traits = toString(clc.trt[,j]); 
              df[1,]$SNPs = toString(clc.max.cvs.1)
            }
          }else{NAs = NAs + 1;}
        }else{
          output = regional2(Z, W, snps.clc, trt.clc, rho, trait.cor, epsilon);
          max.ABF.1 = prior1*exp(mpfr(output[[1]][1], 120));
          sum.ABF.1 = max.ABF.1*output[[1]][2];
          clc.max.cvs.1 = output[[1]][3];
          if(!is.na(sum.ABF.1)){
            sum.ABF.all.1 = sum.ABF.all.1 + sum.ABF.1; 
            if(sum.max.ABF.1<sum.ABF.1){
              sum.max.ABF.1 = sum.ABF.1; 
              df[1,]$Traits = toString(clc.trt[,j]); 
              df[1,]$SNPs = toString(clc.max.cvs.1)
            }
          }else{NAs = NAs + 1;}
          max.ABF.2 = prior2*exp(mpfr(output[[2]][1], 120));
          sum.ABF.2 = max.ABF.2*output[[2]][2];
          clc.max.cvs.2 = snps.clc[,output[[2]][3]];
          if(!is.na(sum.ABF.2)){
            sum.ABF.all.2 = sum.ABF.all.2 + sum.ABF.2; 
            if(sum.max.ABF.2<sum.ABF.2){
              sum.max.ABF.2 = sum.ABF.2; 
              df[2,]$Traits = toString(clc.trt[,j]); 
              df[2,]$SNPs = toString(clc.max.cvs.2)
            }
          }else{NAs = NAs + 1;}
        }
      }
      if(i==m){
        log.max.ABF.1 = log(max.ABF.1)
        sum.ABF.full.1 = sum.ABF.all.1;
        log.sum.ABF.1=asNumeric(log(sum.ABF.1)); 
        df[1,]$log.sum.ABF.full=log.sum.ABF.1;
        df[1,]$log.max.SNP.ABF.full=asNumeric(log.max.ABF.1);
        snp.scores.1 = output[[2]];
        if(flag==1){
          log.max.ABF.2 = log(max.ABF.2);
          sum.ABF.full.2 = sum.ABF.all.2;
          log.sum.ABF.2=asNumeric(log(sum.ABF.2)); 
          df[2,]$log.sum.ABF.full=log.sum.ABF.2;
          df[2,]$log.max.SNP.ABF.full=asNumeric(log.max.ABF.2);
          snp.scores.1 = output[[3]];
          snp.scores.2.tmp = output[[4]];
          snp.scores.2  = ind.snp.score(Q, snp.scores.2.tmp);
        }
      }
      if(flag==0){
        reg.prob = sum.ABF.full.1/(sum.ABF.all.1);
      }else{
        reg.prob = (sum.ABF.full.1+sum.ABF.full.2)/(sum.ABF.all.1+sum.ABF.all.2);
      }
      i = i - 1;
      if(reg.steps!=0 & i == (m-reg.steps-1)){i=0;}
    }
    if(jmp == 1 | sum.ABF.full.1/mx.ABF.k >= (1/min(9,m))*count){df[1,]$reg_bb_alg = TRUE;}  
    log.sum.max.ABF.1=asNumeric(log(sum.max.ABF.1)); 
    log.all.traits.ABF.1=asNumeric(log(sum.ABF.all.1));
    df[1,]$log.sum.ABF.all.combs=log.all.traits.ABF.1; 
    df[1,]$log.max.sum.ABF.all.traits=log.sum.max.ABF.1;
    df[1,]$NaNs=NAs;
    out.list = list(df, snp.scores.1);
    if(flag==1){
      log.sum.max.ABF.2=asNumeric(log(sum.max.ABF.2)); 
      log.all.traits.ABF.2=asNumeric(log(sum.ABF.all.2));
      df[2,]$log.sum.ABF.all.combs=log.all.traits.ABF.2;
      df[2,]$log.max.sum.ABF.all.traits=log.sum.max.ABF.2;
      df[2,]$NaNs=NAs;
      df[2,]$Post.ratio.2.vs.1and2 = asNumeric(sum.ABF.all.2/(sum.ABF.all.1 + sum.ABF.all.2));
      out.list = list(df, snp.scores.1, snp.scores.2);
    }
  }else{
    prior1 = I.unif*p.1.m + (1-I.unif)*prior(prior.1,prior.2, k = m);
    prior2 = I.unif*p.2.m + (1-I.unif)*prior(prior.1,prior.2, k = m)*prior(prior.3, prior.4, k = m); 
    output = regional2(Z, W, snps.clc, c(1:m)+0.0, rho, trait.cor, epsilon);
    max.ABF.1 = prior1*exp(mpfr(output[[1]][1], 120));
    sum.ABF.1 = max.ABF.1*output[[1]][2];
    max.ABF.2 = prior2*exp(mpfr(output[[2]][1], 120));
    sum.ABF.2 = max.ABF.2*output[[2]][2];
    out.list = list(sum.ABF.2/(sum.ABF.2 + sum.ABF.1));
  }
  
  return(out.list)
  
} 

