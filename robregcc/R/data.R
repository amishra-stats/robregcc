#' Simulated date for testing functions in the robregcc package (sparse setting).
#' 
#' A list of response (y), predictors (X) and sub-cpmposition matrix (C).
#' 
#' Vector y, response with a certain percentage of observations as outliers.
#' 
#' Matrix X, Compositional predictors.
#'
#' @docType data
#'
#' @usage data(simulate_robregcc)
#'
#' @format A list with three components: 
#' \describe{
#'   \item{X}{Compositional predictors.}
#'   \item{y}{Outcome with outliers.}
#'   \item{C}{Sub-cmposition matrix.}
#' }
#'
#' @keywords datasets
#' @source Similated data
#' 
#' @examples
#' 
#' library(robregcc)
#' data(simulate_robregcc)
#' X <- simulate_robregcc$X;
#' y <- simulate_robregcc$y
#' C <- simulate_robregcc$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
#' 
"simulate_robregcc"
