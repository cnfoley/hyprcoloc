##########################################################
##### Prior #####
##########################################################  

#' prior
#'
#' @param p.1 probability of one trait having a causal variant in the genetic region
#' @param gamma step probability, i.e. 1 - pc (conditional colocalization probability)
#' @param k number of traits
#' @export
prior <- function(p.1, gamma, k){
  if(k==1){
    p.1;
  }else{
    i=c(2:k);G=(1-gamma^(i-1));p.1*prod(G);
  }
}
