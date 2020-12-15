###########################################################################################################################################
##### Perform a sensitivity analysis by varying the algorithm (regional and alignment) thresholds and coloclaization prior (prior.2)  #####
###########################################################################################################################################

#' sensitivity.plot
#'
#' sensitivity.plot is a function which repeatedly calls the hyprcoloc function to compute a similarity matrix which illustrates how strongly clustered/colocalized pairs of traits are across different input thresholds and priors     
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
#' @param reg.thres a vector of regional probability thresholds
#' @param align.thresh a vector of alignment probability thresholds
#' @param reg.tol regional tolerance parameter
#' @param prior.1 prior probability of a SNP being associated with one trait
#' @param prior.c a vector of prior probabilities where: prior.c is the "conditional colocalization prior", i.e. the probability of a SNP being associated with an additional trait given that the SNP is associated with at least 1 other trait
#' @param prior.3 prior probability that a trait contains a second causal variant given it contains one already
#' @param prior.4 1 - prior probability that trait two co-localises with trait one given traits one and two already share a causal variant and trait one contains a second causal variant
#' @param unifrom.priors uniform priors
#' @param ind.traits are the traits independent or to be treated as independent
#' @param equal.thresholds fix the regional and alignment thresholds to be equal
#' @export
sensitivity.plot = function(effect.est, effect.se, binary.outcomes = rep(0, dim(effect.est)[2]), 
                            trait.subset = c(1:dim(effect.est)[2]), trait.names = c(1:dim(effect.est)[2]),
                            snp.id = c(1:dim(effect.est)[1]), ld.matrix = diag(1, dim(effect.est)[1], dim(effect.est)[1]),
                            trait.cor = diag(1, dim(effect.est)[2], dim(effect.est)[2]), sample.overlap = matrix(rep(1,dim(effect.est)[2]^2), nrow = dim(effect.est)[2]),
                            bb.alg = TRUE, bb.selection = "regional", reg.steps = 1, reg.thresh = c(0.6,0.7,0.8,0.9), align.thresh = c(0.6,0.7,0.8,0.9),
                            prior.1 = 1e-4, prior.c = c(0.02, 0.01, 0.005), 
                            prior.12 = NULL, uniform.priors = FALSE,
                            ind.traits = TRUE, equal.thresholds = FALSE, similarity.matrix = FALSE){
  
  m = dim(effect.est)[2];                            
  snp.combin = function(x, y, vec){I = iterpc(x, y, labels = vec);return(getall(I)+0.0)};
  sim.mat = diag(0,m);
  
  if(!is.null(prior.12)){
    prior.c = prior.12/(prior.1+prior.12);
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
    color             = colorRampPalette((brewer.pal(n = 9, name = "OrRd")))(length(breaksList)),
    breaks = breaksList,
    border_color      = NA,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    #annotation_col    = annotation_row,
    drop_levels       = TRUE,
    fontsize          = 6,
    main              = "Default Heatmap"
  )
  
  if(!similarity.matrix){
    return(plot)}else{
      return(list(plot, sim.mat))
    }        
  
}
