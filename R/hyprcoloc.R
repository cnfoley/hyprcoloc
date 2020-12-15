##########################################################
##### HyPrColoc #####
##########################################################

#' HyPrColoc
#'
#' hyprcoloc is a function which allows the user to perform multi-trait colocalisation analyses in genomic regions
#' @param effect.est matrix of beta values
#' @param effect.se matrix of se values
#' @param binary.outcomes a binary vector depicting binary traits
#' @param trait.subset vector of traits from the full trait list for trageted coloclaisation analysis
#' @param trait.names vector of trait names corresponding to the columns in the effect.est matrix
#' @param snp.id vector of SNP IDs
#' @param ld.matrix LD matrix
#' @param trait.cor correlation matrix between traits
#' @param sample.overlap matrix of sample overlap between traits
#' @param n.cvs number of causal variants
#' @param bb.alg branch and bound algorithm
#' @param bb.selection branch and bound algorithm type
#' @param reg.steps regional step paramter
#' @param window.size size of window for 2CV testing
#' @param sentinel sentinel variant
#' @param epsilon tolerance parameter
#' @param reg.thres regional probability threshold
#' @param align.thresh alignment probability threshold
#' @param reg.tol regional tolerance parameter
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.c conditional colocalization prior: prior probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param prior.12 COLOC users: prior probability that a SNP is associated with two traits
#' @param prior.3 prior probability that a trait contains a second causal variant given it contains one already
#' @param prior.4 1 - prior probability that trait two co-localises with trait one given traits one and two already share a causal variant and trait one contains a second causal variant
#' @param sensitivity perform senstivity analysis
#' @param sens.1 first sensitivity analysis
#' @param sens.2 second sensitivity analysis
#' @param cor.adj.priors correlation adjusted priors
#' @param unifrom.priors uniform priors
#' @param branch.jump branch jump
#' @param ind.traits are the traits independent or to be treated as independent
#' @param snpscores output estimated posterior probability explained each SNP
#' @return data.frame of HyPrColoc colocalization results: rows denote clusters of colocalized traits
#' @return snpscores: a list of estimated posterior probabilities explained by the SNPs; for the BB algorithm there is a set of SNP probabilities for each iteration
#' @import compiler Rmpfr iterpc Matrix
#' @importFrom Rcpp evalCpp
#' @author Christopher N Foley <chris.neal.foley@gmail.com> and James R Staley <jrstaley95@gmail.com>
#' @examples
#' # Regression coefficients and standard errors from ten GWAS studies (Traits 1-5, 6-8 & 9-10 colocalize)
#' betas <- hyprcoloc::test.betas
#' head(betas)
#' ses <- hyprcoloc::test.ses
#' head(ses)
#'   
#' # Trait names and SNP IDs
#' traits <- paste0("T", 1:10)
#' rsid <- rownames(betas)
#' 
#' # Colocalisation analyses
#' results <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
#' @useDynLib hyprcoloc, .registration=TRUE
hyprcoloc <- function(effect.est, effect.se, binary.outcomes = rep(0, dim(effect.est)[2]),
                      trait.subset = c(1:dim(effect.est)[2]), trait.names = c(1:dim(effect.est)[2]),
                      snp.id = c(1:dim(effect.est)[1]), ld.matrix = diag(1, dim(effect.est)[1], dim(effect.est)[1]),
                      trait.cor = diag(1, dim(effect.est)[2], dim(effect.est)[2]), 
                      sample.overlap = matrix(rep(1,dim(effect.est)[2]^2) , nrow = dim(effect.est)[2]), 
                      bb.alg = TRUE, bb.selection = "regional", reg.steps = 1, 
                      reg.thresh = "default", align.thresh = "default", 
                      prior.1 = 1e-4, prior.c = 0.02, prior.12 = NULL,
                      sensitivity = FALSE, sense.1 = 1, sense.2 = 2, 
                      uniform.priors = FALSE, ind.traits = FALSE, snpscores=FALSE){

  if(any(is.na(effect.est))) stop("there are missing values in effect.est")
  if(any(is.na(effect.se))) stop("there are missing values in effect.se")
  if(any(effect.se==0)) stop("there are zero values in effect.se")
  if(any(is.na(binary.outcomes))) stop("there are missing values in binary.outcomes")
  if(any(!(binary.outcomes %in% c(0,1)))) stop("there are missing values in binary.outcomes")
  if(any(is.na(trait.subset))) stop("there are missing values in trait.subset")
  if(any(is.na(trait.names))) stop("there are missing values in trait.names")
  if(any(is.na(snp.id))) stop("there are missing values in snp.id")
  if(any(is.na(ld.matrix))) stop("there are missing values in ld.matrix")
  if(any(is.na(trait.cor))) stop("there are missing values in trait.cor")
  if(any(is.na(sample.overlap))) stop("there are missing values in sample.overlap")
  
  if(!is.null(prior.12)){
    prior.c = prior.12/(prior.1+prior.12);
  }
  prior.2 = 1 - prior.c;
  n.cvs = 1;
  test.2 = F;
  sentinel = 0;
  epsilon = 0;
  window.size = dim(effect.est)[1];
  prior.3 = 1e-3;
  prior.4 = 0.995;
  reg.tol = 0.699;
  cor.adj.priors = F;
  branch.jump = F;

  Z = effect.est/effect.se;    
  m = dim(Z)[2];
  Q = dim(Z)[1];
  
  if(uniform.priors==F & (reg.thresh == "default" | align.thresh == "default")){
    if(reg.thresh == "default"){reg.thresh = 0.5; reg.tol = 0.499;}
    if(align.thresh == "default"){align.thresh = 0.5;}
  }else if(uniform.priors==T & (reg.thresh == "default" | align.thresh == "default")){
    if(reg.thresh == "default"){reg.thresh =0.7;}
    if(align.thresh == "default"){align.thresh = 0.7;}
  }

  if(reg.tol >= reg.thresh){         
    reg.tol = reg.thresh-0.001;
  } 

  if(bb.alg == T & uniform.priors==F & (reg.thresh < 0.5 | align.thresh < 0.5)){stop("Do not run branch and bound algorithm with reg.thresh or align.thresh < 0.5 when using non-uniform priors")};
  if(bb.alg == T & uniform.priors==T & (reg.thresh < 0.6 | align.thresh < 0.6)){stop("Do not run branch and bound algorithm with reg.thresh or align.thresh < 0.6 when using uniform priors")};
  if(! n.cvs %in% c(1,2)){stop("Current version of HyPrMTMC is limited to assessment of a maximum of 2 CVs in a region: set n.cvs to either 1 or 2")}; 
  if(n.cvs == 2 & window.size >= 4E2){print("Computation likely to take a long time. Consider utilising the window.size and sentinel variables to focus assessment within a specified window of a sentinel SNP, default is the lead SNP across all traits")};
  if(sensitivity == T){bb.alg = F; print("Performing a regional probability sensitivity assessment: posterior probability of co-localisation will not be computed")};
  if(! reg.steps %in% 0:m){stop("reg.steps parameter must be an integer between 0 and number of traits m")};
  if(! sense.1 %in% 0:m | ! sense.2 %in% 0:m ){stop("Sensitivity parameters must be intergers between 0 and number of traits m")};
  if(sense.1 >= sense.2){stop("Sensitivity parameter sense.2 must be greater than sensitivity parameter sense.1")};
  if(window.size > Q){stop("Window size cannot be larger than the number of SNPs Q in the region")}
  if((length(trait.subset)<m & typeof(trait.subset)!="integer") & length(trait.names) < m){stop("When using character names in 'trait.subset' the variable 'trait.names' must be of length m, containing names for all traits considered")}
  if(bb.selection == "reg.only" & reg.thresh < 0.9){warning("Posterior evaluation and the BB algorithm are based on values of the regional statistic only.\n Consider setting a larger value for the regional threshold, e.g. > 0.9, to avoid difficulties in interpreting any regional association signals.");} 

  if(class(trait.subset)=="character" & class(trait.names)=="character"){
    tmp.traits = c(1:m);
    traits = tmp.traits[tolower(trait.names) %in% tolower(trait.subset)];
    trait.names = trait.names[traits];
  }else if(class(trait.subset)=="integer"){
    traits = trait.subset;
    trait.names = trait.names[traits];
  }else{stop("trait.subset must be character or intger valued")}
  
  if(n.cvs==1 & test.2==F){
    W = 0.15/effect.se;
    W[, which(binary.outcomes==1)] = 0.2/(effect.se[, which(binary.outcomes==1)]);
  }else if(n.cvs==1 & test.2==T){
    W = 0.15/effect.se;
    W[, which(binary.outcomes==1)] = 0.2/(effect.se[, which(binary.outcomes==1)]);
    W.nudge = matrix(0.15 + runif(m*Q, -0.015, 0.015), nrow=Q)/effect.se;
    W.nudge[, which(binary.outcomes==1)] = matrix(0.2 + runif(Q*length(which(binary.outcomes==1)), -0.02, 0.02), nrow=Q)/(effect.se[, which(binary.outcomes==1)]);
  }else{
    if(epsilon == 0){epsilon = 0.01;}
    W = matrix(0.15 + runif(m*Q, -0.015, 0.015), nrow=Q)/effect.se;
    W[, which(binary.outcomes==1)] = matrix(0.2 + runif(Q*length(which(binary.outcomes==1)), -0.02, 0.02), nrow=Q)/(effect.se[, which(binary.outcomes==1)]);
  }
  
  if(n.cvs==1 & (all(trait.cor[lower.tri(trait.cor)] == 0) | ind.traits == TRUE)){
    ind.traits = TRUE;
    Zsq = Z^2;
    Wsq = 1/(1+ W^2);
    if(reg.steps == 1 & test.2 == FALSE){rapid = TRUE}else{rapid = FALSE};
  }else{
    Zsq = Wsq = sparseMatrix(i=1, j=1, dims=c(Q,m));
    rapid = FALSE;
  }
  
  if(isTRUE(test.2)){
    test.2 = 1;
  }else{
    test.2 = 0;
  }

  snp.combin = function(x, y, vec){I = iterpc(x, y, labels = vec);return(getall(I)+0.0)};
  snp.scores = NA;
  
  if(class(sentinel)=="character"){sentinel = which(snp.id == sentinel);}
  tmp.vars = samp.reduc(Z, W, ld.matrix, trait.cor, sample.overlap, traits, window.size, sentinel, Zsq = Zsq, Wsq = Wsq);
  Z = tmp.vars[[1]];
  W = tmp.vars[[2]];
  Zsq = tmp.vars[[7]];
  Wsq = tmp.vars[[8]]
  ld.matrix = tmp.vars[[3]];
  trait.cor = tmp.vars[[4]];
  sample.overlap = tmp.vars[[5]];
  m = dim(Z)[2];
  Q = dim(Z)[1];
  snp.id = snp.id[tmp.vars[[6]]];
  traits  = 1:length(traits);
  
  cor.trts = identical(trait.cor, diag(m));
  if(branch.jump == T & cor.trts == F){branch.jump = F};
  
  snp.scores = list();
  bb.selection = tolower(bb.selection);

  if(n.cvs==1){
    if(rapid == F){
      snps = c(1:Q);
      snps.clc = as.matrix(t(snp.combin(Q, 2, snps)));
      flag = 0;        
    }
    if(bb.alg==T){    
      df = data.frame(matrix(vector(), 0,9, dimnames=list(c(), c("Algorithm_iteration", "Putatively_co_localised_traits", "HyPr_MTC_posterior", "Regional_probability", "Putatively_causal_SNP",  "Posterior_explained_by_SNP", "Posterior_ratio_2CVs_vs_1_or_2CVs", "Dropped_trait", "NaNs"))), stringsAsFactors=F)
      i=m;
      count = 1;
      reg.prob = align.prob = 0;
      while(i>1){
        trts=traits;
        while(reg.prob <=  reg.thresh | align.prob <= align.thresh){
          reg.prob = reg.only = 0;
          if(bb.selection == "align"){
            tmp.vars = samp.reduc(Z, W, ld.matrix, trait.cor, sample.overlap, trts, Zsq = Zsq, Wsq = Wsq);
            if(rapid == F){
              reg.result = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, 0, reg.steps, cor.adj.priors, uniform.priors, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
              reg.res = reg.result[[1]];
              reg.full = exp(mpfr(reg.res$log.sum.ABF.full,120));
              reg.prob = reg.full/(1 + exp(mpfr(reg.res$log.sum.ABF.all.combs,120)));
            }else{
              reg.res = rapid.reg(tmp.vars[[7]], tmp.vars[[8]], prior.1, prior.2, uniform.priors)
              reg.prob = reg.res[[1]];
            }  
          }else{
            while(reg.prob <= reg.thresh){
              tmp.vars = samp.reduc(Z, W, ld.matrix, trait.cor, sample.overlap, trts, Zsq = Zsq, Wsq = Wsq);
              if(rapid == F){
                reg.result = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, 0, reg.steps, cor.adj.priors, uniform.priors, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
                reg.res = reg.result[[1]];
                trts.tmp = as.numeric(unlist(strsplit(reg.res$Traits, split=", ")));
                reg.full = exp(mpfr(reg.res$log.sum.ABF.full,120));
                reg.prob = reg.full/(1 + exp(mpfr(reg.res$log.sum.ABF.all.combs,120)));
                if(reg.prob <= reg.thresh & length(trts) > length(trts.tmp)){
                  trts = trts[trts.tmp];
                }else if(reg.prob <= reg.thresh & length(trts) == length(trts.tmp) & reg.prob < reg.tol){
                  reg.only.trts = as.numeric(unlist(strsplit(reg.res$reg_only_trts, split=", ")));
                  trts = trts[reg.only.trts]; 
                }else if(reg.prob <= reg.thresh & length(trts) == length(trts.tmp) & reg.prob >= reg.tol){
                  break; 
                }else{
                  if(length(trts) != length(trts.tmp)){reg.prob = 0;trts = trts[trts.tmp];}
                }
              }else{
                reg.res = rapid.reg(tmp.vars[[7]], tmp.vars[[8]], prior.1, prior.2, uniform.priors)
                reg.prob = reg.res[[1]];
                if(reg.prob <= reg.thresh){trts = trts[-reg.res[[3]]]} 
              }
              if(length(trts)<=1){reg.only = 1; break;}
            }
          }  
          if(length(trts)<=1){break;}
          if(bb.selection != "reg.only"){
            if(ind.traits == F){
              align.res = align.ABF.1(tmp.vars[[1]], tmp.vars[[2]], tmp.vars[[4]], tmp.vars[[5]], tmp.vars[[3]], epsilon, reg.full, align.thresh, prior.1, prior.2, cor.adj.priors, uniform.priors, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
              align.prob = asNumeric(reg.full/(reg.full + exp(mpfr(align.res$log.sum.ABF.all.combs, 120))));
              if(align.prob <= align.thresh | reg.prob <=  reg.thresh){
                drop.trait = as.numeric(unlist(strsplit(align.res$Trait.no.clc, split=", ")))
                trts = trts[-drop.trait];
              }
            }else{
              align.res = rapid.align(tmp.vars[[7]], tmp.vars[[8]], prior.1, prior.2, uniform.priors);
              align.prob = align.res[[1]];
              if(align.prob <= align.thresh | reg.prob <=  reg.thresh){trts = trts[-align.res[[2]]];}
            }
          }else{
            reg.only = 1; break;
          }
          if(length(trts)<=1)break
        }
        df[count,]$Regional_probability = asNumeric(reg.prob);
        if(length(trts)>=2){
          df[count,]$Algorithm_iteration = count;          
          df[count,]$Putatively_co_localised_traits = toString(trait.names[trts]);          
          df[count,]$HyPr_MTC_posterior = asNumeric(reg.prob*align.prob);
          if(test.2 == 1){
            tmp.vars = samp.reduc(Z, W.nudge, ld.matrix, trait.cor, sample.overlap, trts, Zsq = Zsq, Wsq = Wsq);
            reg.tmp = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, 1, reg.steps, cor.adj.priors, uniform.priors=FALSE, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
            df[count,]$Posterior_ratio_2CVs_vs_1_or_2CVs = asNumeric(reg.tmp[[1]]); 
          }
          if(rapid == F){
            if(bb.selection != "reg.only" & ind.traits == F){
              df[count,]$NaNs = reg.res$NaNs + align.res$NaNs;
            }else{df[count,]$NaNs = reg.res$NaNs;}   
            df[count,]$Putatively_causal_SNP = toString(snp.id[as.numeric(reg.res$SNPs)]);
            df[count,]$Posterior_explained_by_SNP = exp(reg.res$log.max.SNP.ABF.full-reg.res$log.sum.ABF.full);          
            snp.scores[[count]] = reg.result[[2]];
          }else{
            df[count,]$Putatively_causal_SNP = toString(snp.id[reg.res[[2]]]);
            df[count,]$Posterior_explained_by_SNP = reg.res[[4]][reg.res[[2]]];
            snp.scores[[count]] = reg.res[[4]];
          }
        }else{
          df[count,]$Algorithm_iteration = count;
          df[count,]$Putatively_co_localised_traits = "None";   
          df[count,]$Dropped_trait = toString(trait.names[trts]);
          if(reg.only == 1 & rapid == F){df[count,]$NaNs = reg.res$NaNs;
          }else if(rapid == F){df[count,]$NaNs = reg.res$NaNs + align.res$NaNs;}
        }
        traits = traits[! traits %in% trts];
        reg.prob = align.prob = reg.only = 0;
        i = length(traits);
        count = count + 1;
      }
    }else{
      if(sensitivity == F){
        df = data.frame(matrix(vector(), 0,7, dimnames=list(c(), c("Traits", "HyPr_MTC_posterior", "Regional_probability", "Putatively_causal_SNP",  "Posterior_explained_by_SNP", "Posterior_ratio_2CVs_vs_1_or_2CVs", "NaNs"))), stringsAsFactors=F)
        df[1,]$Traits = toString(trait.names);
        tmp.vars = samp.reduc(Z, W, ld.matrix, trait.cor, sample.overlap, traits, Zsq = Zsq, Wsq = Wsq);
        if(rapid == F){
          reg.result = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, 0, reg.steps, cor.adj.priors, uniform.priors, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
          reg.res = reg.result[[1]];
          reg.full = exp(mpfr(reg.res$log.sum.ABF.full,120));
          reg.prob = asNumeric(reg.full/(1 + exp(mpfr(reg.res$log.sum.ABF.all.combs,120))));
          df[1,]$Regional_probability = reg.prob;
          if(reg.prob < reg.thresh*align.thresh){
            df[1,]$NaNs = reg.res$NaNs;
          }else{
            if(ind.traits == F){
              align.res = align.ABF.1(tmp.vars[[1]], tmp.vars[[2]], tmp.vars[[4]], tmp.vars[[5]], tmp.vars[[3]], epsilon, reg.full, align.thresh, prior.1, prior.2, cor.adj.priors, uniform.priors, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
              align.prob = asNumeric(reg.full/(reg.full + exp(mpfr(align.res$log.sum.ABF.all.combs, 120))));
              df[1,]$NaNs = reg.res$NaNs + align.res$NaNs;      
            }else{
              align.res = rapid.align(tmp.vars[[7]], tmp.vars[[8]], prior.1, prior.2, uniform.priors);
              align.prob = align.res[[1]];
              df[1,]$NaNs = reg.res$NaNs;
            }
            hypr_post = reg.prob*align.prob;
            df[1,]$HyPr_MTC_posterior = hypr_post;   
            if(hypr_post>=(reg.thresh*align.thresh)){
              df[1,]$Putatively_causal_SNP = toString(snp.id[as.numeric(reg.res$SNPs)])
              df[1,]$Posterior_explained_by_SNP = exp(reg.res$log.max.SNP.ABF.full-reg.res$log.sum.ABF.full); 
              if(test.2 == 1){
                reg.tmp = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, 1, reg.steps, cor.adj.priors, uniform.priors = FALSE, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
                df[1,]$Posterior_ratio_2CVs_vs_1_or_2CVs = asNumeric(reg.tmp[[1]]); 
              }
              snp.scores = reg.result[[2]];
            }
          }
        }else{
          rapid.result = rapid.hyprcoloc(tmp.vars[[7]], tmp.vars[[8]], prior.1, prior.2, uniform.priors);
          df[1,]$Regional_probability = rapid.result[[1]];   
          df[1,]$HyPr_MTC_posterior = rapid.result[[2]];   
          df[1,]$Putatively_causal_SNP = toString(snp.id[rapid.result[[3]]]);
          df[1,]$Posterior_explained_by_SNP = rapid.result[[4]][rapid.result[[3]]];
          snp.scores = rapid.result[[4]];
        }
      }else{
        df = data.frame(matrix(vector(), 0,4, dimnames=list(c(), c("Traits", "Regional_Pr_1", "Regional_Pr_2", "Relative_diff"))), stringsAsFactors=F)
        df[1,]$Traits = toString(trait.names);          
        tmp.vars = samp.reduc(Z, W, ld.matrix, trait.cor, sample.overlap, traits, Zsq = Zsq, Wsq = Wsq);
        reg.steps = sense.1 ;
        reg.result = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, 0, reg.steps, cor.adj.priors, uniform.priors, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
        reg.res = reg.result[[1]];
        reg.full = exp(mpfr(reg.res$log.sum.ABF.full,120));
        reg.prob.1 = asNumeric(reg.full/(1 + exp(mpfr(reg.res$log.sum.ABF.all.combs,120))));
        df[1,]$Regional_Pr_1 = reg.prob.1;
        reg.steps = sense.2;
        reg.result = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, 0, reg.steps, cor.adj.priors, uniform.priors, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
        reg.res = reg.result[[1]];
        reg.full = exp(mpfr(reg.res$log.sum.ABF.full,120));
        reg.prob.2 = asNumeric(reg.full/(1 + exp(mpfr(reg.res$log.sum.ABF.all.combs,120))));
        df[1,]$Regional_Pr_2 = reg.prob.2;
        df[1,]$Relative_diff = abs(reg.prob.1-reg.prob.2)/reg.prob.2;
      }
    }
  }else{
    flag = 1;
    test.2 = 0;
    snps = c(1:Q);
    snps.clc = as.matrix(t(snp.combin(Q, 2, snps)));
    if(bb.alg==T){
      NaNs = 0;
      i=m;
      count = 1;
      reg.prob = align.prob = 0;  
      df = data.frame(matrix(vector(), 0,8, dimnames=list(c(), c("Algorithm_iteration", "Putatively_co_localised_traits", "HyPr_MTC_posterior", "Regional_probability", "Putatively_causal_SNPs", "Posterior_explained_by_SNPs", "Dropped_trait", "NaNs"))), stringsAsFactors=F)
      while(i>1){
        trts=traits;
        tmp.vars = samp.reduc(Z, W, ld.matrix, trait.cor, sample.overlap, trts, Zsq = Zsq, Wsq = Wsq)
        while(reg.prob <=  reg.thresh | align.prob <= align.thresh){
          reg.prob = 0;
          while(reg.prob <= reg.thresh){
            reg.full = reg.all = reg.only = 0;
            reg.result = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, test.2, reg.steps, cor.adj.priors, uniform.priors, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
            reg.res = reg.result[[1]];
            reg.full = sum(exp(mpfr(reg.res$log.sum.ABF.full,120)));
            reg.all = sum(exp(mpfr(reg.res$log.sum.ABF.all.combs,120)));
            trts.tmp = as.numeric(unlist(strsplit(reg.res$Traits[1], split=", ")));
            tmp.mx = which(reg.res$log.sum.ABF.full==max(reg.res$log.sum.ABF.full));
            clc.snp = reg.res$SNPs[tmp.mx];
            reg.max.snp = reg.res$log.max.SNP.ABF.full[tmp.mx];
            reg.max.full = reg.res$log.sum.ABF.full[tmp.mx];
            reg.prob = asNumeric(reg.full/(1+reg.all));
            if(reg.prob < reg.thresh){
              if(length(trts) == length(trts.tmp)){
                trts = sample(trts, (length(trts)-1));
              }else{trts = trts[trts.tmp];}
              if(length(trts)>1){
                tmp.vars = samp.reduc(Z, W, ld.matrix, trait.cor, sample.overlap, trts, Zsq = Zsq, Wsq = Wsq);
              }
            }else{
              if(length(trts) != length(trts.tmp)){
                reg.prob = 0;
                trts = trts[trts.tmp];
                tmp.vars = samp.reduc(Z, W, ld.matrix, trait.cor, sample.overlap, trts, Zsq = Zsq, Wsq = Wsq);
              }             
            }
            if(length(trts)<=1){reg.only = 1; break;}
          }
          if(reg.only == 1){break;}
          align.res = align.ABF.2(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[4]], tmp.vars[[5]], tmp.vars[[3]], epsilon, reg.full, align.thresh, prior.1, prior.2, prior.3, prior.4, cor.adj.priors, uniform.priors);
          align.tot = exp(mpfr(align.res$log.sum.ABF.align,120));
          reg.full.align = sum(exp(mpfr(align.res$log.sum.ABF.co.loc,120)));
          if(reg.full.align > reg.full){
            clc.snp = align.res$SNPs.clc; 
            reg.max.snp = exp(mpfr(align.res$log.max.SNP.ABF.full, 120));
            reg.max.full = reg.full.align;
            tmp.mx = 1;
          }
          reg.full = reg.full + reg.full.align;
          drop.trait = as.numeric(unlist(strsplit(align.res$Trait.no.clc, split=", ")));
          align.prob = asNumeric(reg.full/(reg.full + align.tot));
          if(align.prob < align.thresh){
            trts = trts[-drop.trait];
          }
          if(length(trts)<=1)break
        }
        if(reg.only !=1){reg.prob = asNumeric(reg.full/(reg.full.align + reg.all));}
        df[count,]$Regional_probability = reg.prob;
        if(length(trts)>=2){
          NaNs = sum(reg.res$NaNs) + align.res$NaNs
          hypr_post = reg.prob*align.prob;
          df[count,]$Algorithm_iteration = count;          
          df[count,]$Putatively_co_localised_traits = toString(trait.names[trts]);          
          df[count,]$HyPr_MTC_posterior = hypr_post;   
          df[count,]$Putatively_causal_SNPs = toString(snp.id[as.numeric(unlist(strsplit(clc.snp, split=", ")))]);
          df[count,]$Posterior_explained_by_SNPs = asNumeric(reg.max.snp/reg.max.full);          
          df[count,]$NaNs = NaNs; 
          snp.scores[[count]] = reg.result[[tmp.mx+1]];
        }else{
          df[count,]$Algorithm_iteration = count;          
          df[count,]$Dropped_trait =toString(trait.names[trts]);
          if(reg.only == 1){df[count,]$NaNs = sum(reg.res$NaNs);
          }else{df[count,]$NaNs = NaNs;}
        }
        traits = traits[! traits %in% trts];
        reg.prob = align.prob = reg.only = 0;
        i = length(traits);
        count = count + 1;
      }
    }else{
      if(sensitivity == F){
        df = data.frame(matrix(vector(), 0,6, dimnames=list(c(), c("Traits", "HyPr_MTC_posterior", "Regional_probability", "Putatively_causal_SNPs",  "Posterior_explained_by_SNPs", "NaNs"))), stringsAsFactors=F)
        reg.full = reg.all = NaNs = 0;
        align.tot = drop.trait = vector("numeric", window.size+2);
        df[1,]$Traits = toString(trait.names); 
        reg.result = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, test.2, reg.steps, cor.adj.priors, uniform.priors, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
        reg.res = reg.result[[1]];
        reg.full = sum(exp(mpfr(reg.res$log.sum.ABF.full,120)));
        reg.all = sum(exp(mpfr(reg.res$log.sum.ABF.all.combs,120)));
        reg.prob = asNumeric(reg.full/(1 + reg.all));
        df[1,]$Regional_probability = reg.prob;
        NaNs = sum(reg.res$NaNs)
        if(reg.prob < reg.thresh*align.thresh){
          df[1,]$NaNs = NaNs;      
        }else{
          trts.tmp = as.numeric(unlist(strsplit(reg.res$Traits[1], split=", ")));
          tmp.mx = which(reg.res$log.sum.ABF.full==max(reg.res$log.sum.ABF.full));
          clc.snp = reg.res$SNPs[tmp.mx];
          reg.max.snp = reg.res$log.max.SNP.ABF.full[tmp.mx];
          reg.max.full = reg.res$log.sum.ABF.full[tmp.mx]; 
          align.res = align.ABF.2(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[4]], tmp.vars[[5]], tmp.vars[[3]], epsilon, reg.full, align.thresh, prior.1, prior.2, prior.3, prior.4, cor.adj.priors, uniform.priors);
          align.tot = exp(mpfr(align.res$log.sum.ABF.align,120));
          reg.full.align = sum(exp(mpfr(align.res$log.sum.ABF.co.loc,120)));
          if(reg.full.align > reg.full){
            clc.snp = align.res$SNPs.clc; 
            reg.max.full = reg.full.align;
            tmp.mx = 1;
          }
          reg.full = reg.full + reg.full.align;
          drop.trait = as.numeric(unlist(strsplit(align.res$Trait.no.clc, split=", ")));
          NaNs = NaNs + align.res$NaNs;
          reg.prob = asNumeric(reg.full/(1 + reg.full.align + reg.all));
          align.prob = asNumeric(reg.full/(1 + reg.full + align.tot));
          hypr_post = reg.prob*align.prob;
          df[1,]$Regional_probability = reg.prob;   
          df[1,]$HyPr_MTC_posterior = hypr_post;   
          df[1,]$NaNs = NaNs;      
          if(hypr_post>=(reg.thresh*align.thresh)){
            df[1,]$Putatively_causal_SNPs = toString(snp.id[as.numeric(unlist(strsplit(clc.snp, split=", ")))]);
            df[1,]$Posterior_explained_by_SNPs = asNumeric(reg.max.snp/reg.max.full); 
            snp.scores = reg.result[[tmp.mx + 1]];
          }
        }
      }else{
        df = data.frame(matrix(vector(), 0,4, dimnames=list(c(), c("Traits", "Regional_Pr_1", "Regional_Pr_2", "Relative_diff"))), stringsAsFactors=F)
        df[1,]$Traits = toString(trait.names);          
        tmp.vars = samp.reduc(Z, W, ld.matrix, trait.cor, sample.overlap, traits, Zsq = Zsq, Wsq = Wsq);
        reg.steps = sense.1;
        reg.result = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, test.2, reg.steps, cor.adj.priors, uniform.priors, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
        reg.res = reg.result[[1]];
        reg.full = sum(exp(mpfr(reg.res$log.sum.ABF.full,120)));
        reg.all = sum(exp(mpfr(reg.res$log.sum.ABF.all.combs,120)));
        reg.prob.1 = asNumeric(reg.full/(1 + reg.all));
        df[1,]$Regional_Pr_1 = reg.prob.1;
        reg.steps = sense.2;
        reg.result = regional.ABF(tmp.vars[[1]], tmp.vars[[2]], snps.clc, tmp.vars[[3]], tmp.vars[[4]], tmp.vars[[5]], epsilon, reg.thresh, prior.1, prior.2, prior.3, prior.4, flag, test.2, reg.steps, cor.adj.priors, uniform.priors, branch.jump, tmp.vars[[7]], tmp.vars[[8]], ind.traits);
        reg.res = reg.result[[1]];
        reg.full = sum(exp(mpfr(reg.res$log.sum.ABF.full,120)));
        reg.all = sum(exp(mpfr(reg.res$log.sum.ABF.all.combs,120)));
        reg.prob.2 = asNumeric(reg.full/(1 + reg.all));
        df[1,]$Regional_Pr_2 = reg.prob.2;
        df[1,]$Relative_diff = abs(reg.prob.1-reg.prob.2)/reg.prob.2;
      }
    }
  }
  
  if(length(which(df$NaNs!=0)) > 1){print("Results include NAs. Possible the Z score correlation matrix or adjusted prior matrix has a high condition number somewhere: consider increasing epsilon by a small amount");}
  names(df)[names(df)=="Algorithm_iteration"] <- "iteration"
  names(df)[names(df)=="Putatively_co_localised_traits"] <- "traits"
  names(df)[names(df)=="Traits"] <- "traits"
  names(df)[names(df)=="HyPr_MTC_posterior"] <- "posterior_prob"
  names(df)[names(df)=="Regional_probability"] <- "regional_prob"
  names(df)[names(df)=="Posterior_ratio_2CVs_vs_1_or_2CVs"] <- "posterior_ratio_2cvs_vs_1_or_2cvs"
  names(df)[names(df)=="Putatively_causal_SNP"] <- "candidate_snp"
  names(df)[names(df)=="Putatively_causal_SNPs"] <- "candidate_snp"
  names(df)[names(df)=="Posterior_explained_by_SNP"] <- "posterior_explained_by_snp"
  names(df)[names(df)=="Posterior_explained_by_SNPs"] <- "posterior_explained_by_snp"
  names(df)[names(df)=="Dropped_trait"] <- "dropped_trait"
  names(df)[names(df)=="Regional_Pr_1"] <- "regional_prob_1"
  names(df)[names(df)=="Regional_Pr_2"] <- "regional_prob_2"
  names(df)[names(df)=="Relative_diff"] <- "relative_diff"
  df$NaNs <- NULL
  rownames(df) <- NULL
  if(test.2==F){df$posterior_ratio_2cvs_vs_1_or_2cvs <- NULL}
  if(any(!is.na(df$posterior_prob))){df$posterior_prob[!is.na(df$posterior_prob)] <- round(as.numeric(df$posterior_prob[!is.na(df$posterior_prob)]),4); df$regional_prob[!is.na(df$regional_prob)] <- round(as.numeric(df$regional_prob[!is.na(df$regional_prob)]),4); df$posterior_explained_by_snp[!is.na(df$posterior_explained_by_snp)] <- round(as.numeric(df$posterior_explained_by_snp[!is.na(df$posterior_explained_by_snp)]),4)}
  if(any(!is.na(df$regional_prob))){df$regional_prob[!is.na(df$regional_prob)] <- round(as.numeric(df$regional_prob[!is.na(df$regional_prob)]),4); df$posterior_explained_by_snp[!is.na(df$posterior_explained_by_snp)] <- round(as.numeric(df$posterior_explained_by_snp[!is.na(df$posterior_explained_by_snp)]),4)}
  if(length(df$posterior_ratio_2cvs_vs_1_or_2cvs)>0){if(any(!is.na(df$posterior_ratio_2cvs_vs_1_or_2cvs))){df$posterior_ratio_2cvs_vs_1_or_2cvs[!is.na(df$posterior_ratio_2cvs_vs_1_or_2cvs)] <- round(as.numeric(df$posterior_ratio_2cvs_vs_1_or_2cvs[!is.na(df$posterior_ratio_2cvs_vs_1_or_2cvs)]),4)}}
  if(sum(is.na(snp.scores))==0 & length(snp.scores)>0 & snpscores==T){results = list(results=df, snpscores=snp.scores)}else{results = list(results=df)};
  class(results) <- "hyprcoloc"; 

  return(results)

}