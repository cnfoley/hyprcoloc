##########################################################
##### SNP scores #####
##########################################################

p.0 <- function(j, Q){
  (j-1)*Q - j*(j-1)/2 + 1;
}

#' ind.snp.score
#'
#' @param Q the number of traits 
#' @param snp.scores per snp posterior explained values
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
