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
