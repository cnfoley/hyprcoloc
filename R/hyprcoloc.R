##########################################################
##### Prior #####
##########################################################  

#' prior
#'
#' @param p.1 probability of one trait having a causal variant in the genetic region
#' @param gamma step probability
#' @param k number of traits
#' @export
prior <- function(p.1, gamma, k){
  if(k==1){
    p.1;
  }else{
    i=c(2:k);G=(1-gamma^(i-1));p.1*prod(G);
  }
}

##########################################################
##### SNP scores #####
##########################################################

#' p.0
#'
#' @param j number of traits
#' @param Q number of SNPs
#' @export
p.0 <- function(j, Q){
  (j-1)*Q - j*(j-1)/2 + 1;
}

#' ind.snp.score
#'
#' @param Q the number of traits 
#' @param snp.scores vector of per snp contributions to the posterior probability of colocalization
#' @export
ind.snp.score <- function(Q, snp.scores){
  loc = vector("numeric", Q-1);
  res = vector("numeric", Q);
  for(j in 1:Q){
    i = 1;
    while(i < j){
      loc[i] = p.0(i, Q) + j - (i + 1);
      i = i+1;
    }
    if(j<=Q-1){
      loc[i : (Q-1)] = (p.0(j, Q)):(p.0(j+1, Q) - 1);
    }
    res[j] = sum(snp.scores[loc]);
  }
  return(res)
}


##########################################################
##### Reduce dimensions of datasets #####
########################################################## 

#' samp.reduc
#'
#' @param Z matrix of Z-scores
#' @param W ratio matrix of prior standard deviation and observed standard errors
#' @param ld.matrix LD matrix
#' @param trait.cor correlation matrix between traits
#' @param sample.overlap matrix of sample overlap between traits
#' @param trts vector of traits
#' @param window.size size of window for 2CV testing
#' @param sentinel sentinel variant 
#' @param Zsq matrix of Z-scores squared
#' @param Wsq matrix of W squared
#' @export
samp.reduc <- function(Z, W, ld.matrix, trait.cor, sample.overlap, trts, window.size = dim(Z)[1], sentinel = 0, Zsq, Wsq){
  
  q = dim(Z)[1];
  Z = Z[,trts];
  W = W[,trts];
  Zsq = Zsq[,trts];
  Wsq = Wsq[,trts];
  
  trait.cor = trait.cor[trts,trts];
  sample.overlap = sample.overlap[trts,trts];
  
  if(window.size < q){
    if(sentinel == 0){
      max.Z = which(abs(Z) == max(abs(Z)));
      if(length(max.Z)>1){max.Z = sample(max.Z,1)};
      row.Z = max.Z%%dim(Z)[1];
      if(row.Z == 0){row.Z = dim(Z)[1]};
    }else{ row.Z = sentinel;}
    
    Q = window.size;
    if(Q%%2 == 1){Q = Q-1;}
    
    if(max(row.Z - Q/2,1) == 1){
      row.col = 1:(Q+1);
      Z = Z[row.col,];
      W = W[row.col,];
      Zsq = Zsq[row.col,];
      Wsq = Wsq[row.col,];
      ld.matrix = ld.matrix[row.col,row.col];
    }else if(min(row.Z + Q/2,q)==q){
      row.col = (q - Q):q;
      Z = Z[row.col,];
      W = W[row.col,];
      Zsq = Zsq[row.col,];
      Wsq = Wsq[row.col,];
      ld.matrix = ld.matrix[row.col,row.col];
    }else{
      row.col = (row.Z-Q/2):(row.Z+Q/2);
      Z = Z[row.col,];
      W = W[row.col,];
      Zsq = Zsq[row.col,];
      Wsq = Wsq[row.col,];
      ld.matrix = ld.matrix[row.col,row.col];
    }
  }else{row.col = 1:q}
  
  return(list(Z, W, ld.matrix, trait.cor, sample.overlap, row.col, Zsq, Wsq))
  
}

#' colMax
#'
#' @param x data.frame or matrix
#' @export
colMax <- function(x){
  apply(x, MARGIN = c(2), max)
}


################################################################################
##### Compute credible sets of snps for each clusetr of colocalized traits #####
################################################################################

#' cred.sets
#'
#' @param res list of results from hyprcoloc when setting "snpscores = TRUE" 
#' @param value the sum of the probabilities of the snps included in the credible set  
#' @export

cred.sets = function(res, value = 0.95){
  clusters = which(! is.na(res[[1]]$traits));
  results = vector('list', length(clusters));
  if(length(clusters)>0){
    count = 1;
    for(i in clusters){
      snp.scores = res[[2]][[i]];
      snp.scores = snp.scores[order(-snp.scores)];
      tmp.rmv = which(snp.scores < (1- value)/length(snp.scores));
      tmp.lgt = length(snp.scores) - length(tmp.rmv);
      snp.scores = snp.scores[- tmp.rmv];
      if(tmp.lgt >1){
        row.col = combn(tmp.lgt,2)
        row = c(1:tmp.lgt, row.col[1,]);
        col = c(1:tmp.lgt, row.col[2,]);
        snp.cols = sparseMatrix(i = row, j = col, x = 1, dims = c(tmp.lgt,tmp.lgt));
        wch.snps = snp.cols[,min(which(colSums(snp.scores*snp.cols) >= value))]==1;
        results[[count]] = snp.scores[wch.snps];
      }else{ results[[count]] = snp.scores[tmp.lgt];}
      count = count + 1;
    }
  }
  return(results)
}


###########################################################################################################################################
##### Perform a sensitivity analysis by varying the algorithm (regional and alignment) thresholds and coloclaization prior (prior.c)  #####
###########################################################################################################################################

#' sensitivity.plot
#'
#' sensitivity.plot is a function which repeatedly calls the hyprcoloc function to compute a similarity matrix which illustrates how strongly clustered/colocalized pairs of traits are across different input thresholds and priors     
#' @param effect.est matrix of snp regression coefficients (i.e. regression beta values) in the genomic region 
#' @param effect.se matrix of standard errors associated with the beta values
#' @param binary.outcomes a binary vector of dimension the number of traits: 1 represents a binary trait 0 otherwise
#' @param trait.subset vector of trait names (or number) from the full trait list: used for trageted colocalization analysis in a region
#' @param trait.names vector of trait names corresponding to the columns in the effect.est matrix
#' @param snp.id vector of SNP IDs
#' @param ld.matrix LD matrix
#' @param trait.cor matrix of pairwise correlations between traits
#' @param sample.overlap matrix of pairwise sample overlap between traits
#' @param bb.alg branch and bound algorithm: TRUE, employ BB algorithm; FALSE, do not
#' @param bb.selection branch and bound algorithm type, e.g. regional or alignment selection
#' @param reg.steps regional step paramter
#' @param reg.thresh a vector of regional probability thresholds
#' @param align.thresh a vector of alignment probability thresholds
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.c vector of conditional colocalization priors: where prior.c is the probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param prior.12 COLOC prior p12: prior probability of a SNP being associated with any two traits
#' @param uniform.priors uniform priors
#' @param ind.traits are the traits independent or to be treated as independent
#' @param equal.thresholds fix the regional and alignment thresholds to be equal
#' @param similarity.matrix To be viewed as a similarity matrix. Default FALSE.
#' @import pheatmap RColorBrewer
#' @export
sensitivity.plot = function(effect.est, effect.se, binary.outcomes = rep(0, dim(effect.est)[2]), 
                            trait.subset = c(1:dim(effect.est)[2]), trait.names = c(1:dim(effect.est)[2]),
                            snp.id = c(1:dim(effect.est)[1]), ld.matrix = diag(1, dim(effect.est)[1], dim(effect.est)[1]),
                            trait.cor = diag(1, dim(effect.est)[2], dim(effect.est)[2]), 
                            sample.overlap = matrix(rep(1,dim(effect.est)[2]^2), nrow = dim(effect.est)[2]),
                            bb.alg = TRUE, bb.selection = "regional", 
                            reg.steps = 1, reg.thresh = c(0.6,0.7,0.8,0.9), align.thresh = c(0.6,0.7,0.8,0.9),
                            prior.1 = 1e-4, prior.c = c(0.02, 0.01, 0.005), prior.12 = NULL,
                            uniform.priors = FALSE,
                            ind.traits = TRUE, equal.thresholds = FALSE, similarity.matrix = FALSE){
  
  m = dim(effect.est)[2];                            
  snp.combin = function(x, y, vec){I = iterpc(x, y, labels = vec);return(getall(I)+0.0)};
  sim.mat = diag(0,m);
  
  if(!is.null(prior.12)){
    prior.c = prior.12/(prior.1 + prior.12);
  }
  
  for(i in reg.thresh){
      for(k in prior.c){
          if(equal.thresholds){
              j = i;
              tmp.mat = diag(1,m);                                     
              res = hyprcoloc(effect.est, effect.se, binary.outcomes = binary.outcomes, trait.subset = trait.subset, trait.names = trait.names,
                              snp.id = snp.id, ld.matrix = ld.matrix, trait.cor = trait.cor, sample.overlap = sample.overlap, bb.alg = bb.alg, bb.selection = bb.selection,
                              reg.steps = reg.steps, reg.thresh = i, align.thresh = j,
                              prior.1 = prior.1, prior.c = k, uniform.priors = uniform.priors, ind.traits = ind.traits);
              trt.clusts = res[[1]]$traits;
              for(its in 1:length(trt.clusts)){
                tmp.clust = unlist(strsplit(trt.clusts[its], split=", "));
                if(tmp.clust[1]!="None"){
                  tmp.vec = which(trait.names %in% tmp.clust);
                  # coloc.pairs = t(snp.combin(m, 2, tmp.vec));
                  coloc.pairs = t(snp.combin(length(tmp.vec), 2, tmp.vec));
                  tmp.mat[t(coloc.pairs)] = 1;
                }
              }
              sim.mat = sim.mat + tmp.mat + t(tmp.mat) - diag(1,m);
          }else{
              for(j in align.thresh){
              tmp.mat = diag(1,m);                                     
              res = hyprcoloc(effect.est, effect.se, binary.outcomes = binary.outcomes, trait.subset = trait.subset, trait.names = trait.names,
                              snp.id = snp.id, ld.matrix = ld.matrix, trait.cor = trait.cor, sample.overlap = sample.overlap, bb.alg = bb.alg, bb.selection = bb.selection,
                              reg.steps = reg.steps, reg.thresh = i, align.thresh = j,
                              prior.1 = prior.1, prior.c = k, uniform.priors = uniform.priors, ind.traits = ind.traits);
              trt.clusts = res[[1]]$traits;
              for(its in 1:length(trt.clusts)){
                tmp.clust = unlist(strsplit(trt.clusts[its], split=", "));
                if(tmp.clust[1]!="None"){
                  tmp.vec = which(trait.names %in% tmp.clust);
                  # coloc.pairs = t(snp.combin(m, 2, tmp.vec));
                  coloc.pairs = t(snp.combin(length(tmp.vec), 2, tmp.vec));
                  tmp.mat[t(coloc.pairs)] = 1;
                }
              }
              sim.mat = sim.mat + tmp.mat + t(tmp.mat) - diag(1,m);
              }
        }
    }
  }
  sim.mat = sim.mat/length(reg.thresh)/length(align.thresh)/length(prior.c);
  if(equal.thresholds){
  sim.mat = sim.mat*length(align.thresh);
  }
  rownames(sim.mat) = trait.names;
  colnames(sim.mat) = trait.names;
  
  #                                  annotation_row = data.frame(
  #                                  Clusters = factor(dta$cluster.class[match(rownames(smat), obs.names)])
  #                                )
  #                                rownames(annotation_row) <- rownames(smat)
  breaksList = seq(0,1,by=0.02);
  plot = pheatmap(
    mat               = sim.mat,
    color             = colorRampPalette((brewer.pal(n = 9, name = "Blues")))(length(breaksList)),
    breaks = breaksList,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    border_color      = NA,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    #annotation_col    = annotation_row,
    drop_levels       = TRUE,
    fontsize          = 6,
    main              = "Prior sensitivity heatmap"
  )
  
  if(!similarity.matrix){
    return(plot)}else{
      return(list(plot, sim.mat))
    }        
  
}



##########################################################
##### Regional colocalisation (rapid) #####
##########################################################

#' rapid.reg
#'
#' @param Zsq matrix of Z-scores
#' @param Wsq ratio matrix of prior standard deviation and observed standard errors squared
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.2 1 - prior probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param uniform.priors uniform priors
#' @export
rapid.reg <- function(Zsq, Wsq, prior.1, prior.2, uniform.priors){
  
  m = dim(Zsq)[2];
  Q = dim(Zsq)[1];
  
  if(uniform.priors==T){
    I.unif = 1;
  }else{
    I.unif = 0;
  }
  
  p.1 = prior.1;
  p.1.m = p.1/Q;
  
  prior.all = I.unif*p.1.m + (1-I.unif)*prior(prior.1, prior.2, k = m);
  prior.sub = I.unif*p.1.m + (1-I.unif)*prior(prior.1, prior.2, k = m-1);  
  
  labf = 0.5*(log(Wsq) + (Zsq)*(1- Wsq));
  sum.labf = rowSums(labf);
  mx.labf = which.max(labf);
  clc.max = which.max(sum.labf);
  
  chi = sum.labf - sum.labf[clc.max];
  exp.chi = exp(chi);
  sum.exp.chi = sum(exp.chi);
  sum.labf.subset = colSums(exp(chi - labf));
  rm.trt = which.max(sum.labf.subset);
  
  reg.prob = sum.exp.chi/(exp(-sum.labf[clc.max] - log(prior.all))  + sum.exp.chi + (prior.sub/prior.all)*sum(sum.labf.subset));

  return(list(reg.prob, clc.max, rm.trt, exp.chi/sum.exp.chi))

}

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
#' @param uniform.priors uniform priors
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


##########################################################
##### Alignment (rapid) #####
##########################################################

#' rapid.align
#'
#' @param Zsq matrix of Z-scores
#' @param Wsq ratio matrix of prior standard deviation and observed standard errors squared
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.2 1 - prior probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param uniform.priors uniform priors
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
#' @param uniform.priors uniform priors
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

##########################################################
##### HyPrColoc (rapid) #####
##########################################################

#' rapid.hyprcoloc
#'
#' @param Zsq matrix of Z-scores
#' @param Wsq ratio matrix of prior standard deviation and observed standard errors squared
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.2 1 - prior probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param uniform.priors uniform priors
#' @export
rapid.hyprcoloc <- function(Zsq, Wsq, prior.1, prior.2, uniform.priors){
  
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
  prior.sub = I.unif*p.1.m + (1-I.unif)*prior(prior.1, prior.2, k = m-1);  
  prior.align= I.unif*p.m.1 + (1-I.unif)*prior.1*prior(prior.1, prior.2, k = m-1);
  
  labf = 0.5*(log(Wsq) + (Zsq)*(1- Wsq));
  sum.labf = rowSums(labf);
  mx.labf = which.max(labf);
  clc.max = which.max(sum.labf);
  
  chi = sum.labf - sum.labf[clc.max];
  exp.chi = exp(chi);
  sum.exp.chi = sum(exp.chi);
  sum.labf.subset = colSums(exp(chi - labf));
  rm.trt = which.max(sum.labf.subset);
  
  reg.prob = sum.exp.chi/(exp(-sum.labf[clc.max] - log(prior.all))  + sum.exp.chi + (prior.sub/prior.all)*sum(sum.labf.subset));
 
  col.max.labf = colMax(labf);
  exp.labf.std = exp(t(t(labf) - col.max.labf));
  labf.tmp = t(colSums(exp.labf.std) - t(exp.labf.std));
  sum.tmp = colSums(exp(t(t(chi - labf)+col.max.labf))*labf.tmp);

  align.prob = sum.exp.chi/(sum.exp.chi + cnst*(prior.align/prior.all)*sum(sum.tmp));
  
  return(list(reg.prob, reg.prob*align.prob, clc.max, exp.chi/sum.exp.chi))

}

##########################################################
##### HyPrColoc #####
##########################################################

#' HyPrColoc
#'
#' hyprcoloc is a function used to identify clusters of colocalized traits and candidate causal SNPs in genomic regions
#' @param effect.est matrix of snp regression coefficients (i.e. regression beta values) in the genomic region 
#' @param effect.se matrix of standard errors associated with the beta values
#' @param binary.outcomes a binary vector of dimension the number of traits: 1 represents a binary trait 0 otherwise
#' @param trait.subset vector of trait names (or number) from the full trait list: used for trageted colocalization analysis in a region
#' @param trait.names vector of trait names corresponding to the columns in the effect.est matrix
#' @param snp.id vector of SNP IDs
#' @param ld.matrix LD matrix
#' @param trait.cor matrix of pairwise correlations between traits
#' @param sample.overlap matrix of pairwise sample overlap between traits
#' @param bb.alg branch and bound algorithm: TRUE, employ BB algorithm; FALSE, do not
#' @param bb.selection branch and bound algorithm type, e.g. regional or alignment selection
#' @param reg.steps regional step paramter
#' @param reg.thresh threshold probability beyond which traits are believed to share a regional association signal 
#' @param align.thresh threshold probability beyond which traits are believed to align at a single causal variant
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.c conditional colocalization prior: probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param prior.12 COLOC prior p12: prior probability of a SNP being associated with any two traits
#' @param sensitivity perform senstivity analysis
#' @param sense.1 first sensitivity analysis
#' @param sense.2 second sensitivity analysis
#' @param uniform.priors uniform priors
#' @param ind.traits are the traits independent or to be treated as independent
#' @param snpscores output estimated posterior probability explained each SNP
#' @return A data.frame of HyPrColoc results: each row is a cluster of colocalized traits or is coded NA (if no colocalization is identified)
#' @return If snpscores=TRUE: additionally returns a list of posterior probability explained by each SNPs and for each cluster of colocalized traits identified
#' @import compiler Rmpfr iterpc Matrix
#' @importFrom Rcpp evalCpp
#' @useDynLib hyprcoloc
#' @author Christopher N Foley <chris.neal.foley@gmail.com> and James R Staley <jrstaley95@gmail.com>
#' @examples
#' # Regression coefficients and standard errors from ten GWAS studies
#' # (Traits 1-5, 6-8 & 9-10 are the clusters of colocalized traits)
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
#' @export
hyprcoloc <- function(effect.est, effect.se, binary.outcomes = rep(0, dim(effect.est)[2]), 
                      trait.subset = c(1:dim(effect.est)[2]), trait.names = c(1:dim(effect.est)[2]), 
                      snp.id = c(1:dim(effect.est)[1]), ld.matrix = diag(1, dim(effect.est)[1], dim(effect.est)[1]),
                      trait.cor = diag(1, dim(effect.est)[2], dim(effect.est)[2]), 
                      sample.overlap = matrix(rep(1,dim(effect.est)[2]^2) , nrow = dim(effect.est)[2]), 
                      bb.alg = TRUE, bb.selection = "regional", 
                      reg.steps = 1, reg.thresh = "default", align.thresh = "default", 
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
    prior.c = prior.12/(prior.1 + prior.12);
  }
  prior.2 = 1-prior.c;
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

#' Print HyPrColoc
#'
#' print method for class "hyprcoloc"
#' @param x an object of class "hyprcoloc"
#' @author Christopher N Foley (University of Cambridge) <chris.neal.foley@gmail.com> and James R Staley (University of Bristol) <jrstaley95@gmail.com>
#' @param ... Other arguments passed to or from other methods
#'
#' @export print.hyprcoloc
#' @export
print.hyprcoloc <- function(x, ...){
  cat("\nCall: \nhyprcoloc")
  cat("\n\nResults:\n")
  print(x$results)
  cat("\n")
}

# prevents R check throwing errors about unassigned global variables
utils::globalVariables(c("colorRampPalette", "combn", "runif"))
