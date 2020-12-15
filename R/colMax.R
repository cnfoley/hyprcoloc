#' colMax
#'
#' @param x data.frame or matrix
#' @export
colMax <- function(x){
  apply(x, MARGIN = c(2), max)
}
