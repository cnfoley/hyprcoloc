#' Print HyPrColoc
#'
#' print method for class "hyprcoloc"
#' @param x an object of class "hyprcoloc"
#' @author Christopher N Foley <chris.neal.foley@gmail.com> and James R Staley <jrstaley95@gmail.com>
#' @export
print.hyprcoloc <- function(x, ...){
  cat("\nCall: \nhyprcoloc")
  cat("\n\nResults:\n")
  print(x$results)
  cat("\n")
}