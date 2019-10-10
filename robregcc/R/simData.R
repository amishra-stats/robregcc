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
#' @usage data(simulate_robregcc_sp)
#'
#' @format A list with three components: 
#' \describe{
#'   \item{X}{Compositional predictors.}
#'   \item{y}{Outcome with outliers.}
#'   \item{C}{Sub-cmposition matrix.}
#' }
#'
#' @keywords datasets

#' @examples
#' data(simulate_robregcc_sp)
#' X <- simulate_robregcc_sp$X;
#' y <- simulate_robregcc_sp$y
#' C <- simulate_robregcc_sp$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
"simulate_robregcc_sp"



#' Simulated date for testing functions in the robregcc package (non-sparse setting).
#' 
#' A list of response (y), predictors (X) and sub-cpmposition matrix (C).
#' 
#' Vector y, response with a certain percentage of observations as outliers.
#' 
#' Matrix X, Compositional predictors.
#'
#' @docType data
#'
#' @usage data(simulate_robregcc_nsp)
#'
#' @format A list with three components: 
#' \describe{
#'   \item{X}{Compositional predictors.}
#'   \item{y}{Outcome with outliers.}
#'   \item{C}{Sub-cmposition matrix.}
#' }
#'
#' @keywords datasets

#' @examples
#' data(simulate_robregcc_nsp)
#' X <- simulate_robregcc_nsp$X;
#' y <- simulate_robregcc_nsp$y
#' C <- simulate_robregcc_nsp$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
"simulate_robregcc_nsp"