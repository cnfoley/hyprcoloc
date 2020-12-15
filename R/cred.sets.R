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
