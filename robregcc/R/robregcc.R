if (getRversion() >= "3.5.0") utils::globalVariables(c("."))


#' Simulation data
#'
#' Simulate data for the robust regression with compositional covariates
#'
#' @param n sample size
#' @param betacc model parameter satisfying compositional covariates
#' @param beta0 intercept 
#' @param O number of outlier
#' @param Sigma covariance matrix of simulated predictors
#' @param levg 1/0 whether to include leveraged observation or not
#' @param snr noise to signal ratio
#' @param shft multiplying factor to model variance for creating outlier
#' @param m test sample size
#' @param C subcompositional matrix
#' @param out list for obtaining output with simulated data structure
#' @return a list containing simulated output.
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
#' @import magrittr
#' @importFrom stats sd
#' @export
#' @examples  
#' ## Simulation example:
#' library(robregcc)
#' library(magrittr)
#' 
#' ## n: sample size 
#' ## p: number of predictors
#' ## o: fraction of observations as outliers
#' ## L: {0,1} => leveraged {no, yes}, 
#' ## s: multiplicative factor 
#' ## ngrp: number of subgroup in the model 
#' ## snr: noise to signal ratio for computing true std_err
#' 
#' ## Define parameters to simulate example
#' p <- 80                # number of predictors  
#' n <- 300               # number of sample   
#' O <- 0.10*n            # number of outlier
#' L <- 1                             
#' s <- 8                          
#' ngrp <- 4              # number of sub-composition
#' snr <- 3               # Signal to noise ratio
#' example_seed <- 2*p+1  # example seed
#' set.seed(example_seed) 

#' # Simulate subcomposition matrix
#' C1 <- matrix(0,ngrp,23)
#' tind <- c(0,10,16,20,23)
#' for (ii in 1:ngrp)
#'   C1[ii,(tind[ii] + 1):tind[ii + 1]] <- 1
#' C <- matrix(0,ngrp,p)
#' C[,1:ncol(C1)] <- C1            
#' # model parameter beta
#' beta0 <- 0.5
#' beta <- c(1, -0.8, 0.4, 0, 0, -0.6, 0, 0, 0, 0, -1.5, 0, 1.2, 0, 0, 0.3)
#' beta <- c(beta,rep(0,p - length(beta)))
#' # Simulate response and predictor, i.e., X, y
#' Sigma  <- 1:p %>% outer(.,.,'-') %>% abs(); Sigma  <- 0.5^Sigma
#' data.case <- vector("list",1)
#' set.seed(example_seed)
#' data.case <- robregcc_sim(n,beta,beta0, O = O,
#'       Sigma,levg = L, snr,shft = s,0, C,out = data.case)
#' data.case$C <- C                         
#' # We have saved a copy of simulated data in the package 
#' # with name simulate_robregcc 
#' # simulate_robregcc = data.case;
#' # save(simulate_robregcc, file ='data/simulate_robregcc.rda')
#' 
#' X <- data.case$X                 # predictor matrix
#' y <- data.case$y                 # model response 
#' 
#' 
#' @references
#' Mishra, A., Mueller, C.,(2019) \emph{Robust regression with compositional covariates. In prepration.} arXiv:1909.04990.
robregcc_sim <- function(n, betacc, beta0, O, Sigma, levg, snr,
                         shft, m, C, out = list()) {
  out$p <- length(betacc)
  out$beta <- betacc
  out$n <- n
  out$L <- levg
  out$O <- O
  out$shft <- seq(shft + 1, shft, length.out = O)
  x <- exp(MASS::mvrnorm(n, c(
    rep(log(out$p / 2), 5),
    rep(0, out$p - 5)
  ), Sigma))
  if (out$L == 1) {
    x[1:round(out$O / 2), ] <- getLevObs(x, round(out$O / 2), C)
  }
  x <- x %>%
    apply(1, function(x) x / sum(x)) %>%
    t() %>%
    log()
  y <- beta0 + tcrossprod(x, t(betacc))
  sigma <- 1 # if(betafac)
  sigma <- norm(y, "2") / sqrt(n) / snr
  out$shift <- c(sign(y[1:out$O]) * out$shft * sigma, rep(0, n - out$O))
  y <- y + sigma * rnorm(n)
  out$yo <- y
  # out$shift <- c(rep(shft*sigma,out$O),rep(0,n-out$O))
  y <- y + out$shift
  out$sigma <- sigma
  out$SNR <- snr
  out$y <- y
  out$X <- x
  ## generate test  data set for the simulation:
  if (m > 0) {
    x <- exp(MASS::mvrnorm(m, c(
      rep(log(out$p / 2), 5),
      rep(0, out$p - 5)
    ), Sigma))
    out$Xte <- x %>%
      apply(1, function(x) x / sum(x)) %>%
      t() %>%
      log()
    out$Yte <- beta0 + tcrossprod(out$Xte, t(betacc)) + sigma * rnorm(m)
  }
  return(out)
}







#' Control  parameter for model estimation:
#'
#' The model approach use scaled lasoo approach for model selection.
#'
#' @param maxiter maximum number of iteration for convergence
#' @param tol tolerance value set for convergence
#' @param nlam number of lambda to be genrated to obtain solution path
#' @param out.tol tolernce value set for convergence of outer loop
#' @param lminfac a multiplier of determing lambda_min as a fraction of lambda_max
#' @param lmaxfac a multiplier of lambda_max
#' @param mu penalty parameter used in enforcing orthogonality
#' @param nu penalty parameter used in enforcing orthogonality (incremental rate of mu)
#' @param sp maximum proportion of nonzero elements in shift parameter
#' @param gamma adaptive penalty weight exponential factor
#' @param outMiter maximum number of outer loop iteration
#' @param inMiter maximum number of inner loop iteration
#' @param kmaxS maximum number of iteration for fast S estimator for convergence
#' @param tolS tolerance value set for convergence in case of fast S estimator
#' @param nlamx number of x lambda
#' @param nlamy number of y lambda
#' @param spb sparsity in beta
#' @param spy sparsity in shift gamma
#' @param lminfacX a multiplier of determing lambda_min as a fraction of lambda_max for sparsity in X
#' @param lminfacY a multiplier of determing lambda_min as a fraction of lambda_max for sparsity in shift parameter
#' @param kfold nummber of folds for crossvalidation
#' @param fullpath 1/0 to get full path yes/no
#' @param sigmafac multiplying factor for the range of standard deviation
#' @return a list of controling parameter.
#' @export
#' @examples  
#' # default options
#' library(robregcc)
#' control_default = robregcc_option()
#' # manual options
#' control_manual <- robregcc_option(maxiter=1000,tol = 1e-4,lminfac = 1e-7)
robregcc_option <- function(maxiter = 10000, tol = 1e-10, nlam = 100,
                            out.tol = 1e-8,
                            lminfac = 1e-8, lmaxfac = 10, mu = 1,
                            nu = 1.05, sp = 0.3, gamma = 2,
                            outMiter = 3000, inMiter = 500,
                            kmaxS = 500,
                            tolS = 1e-4, nlamx = 20,
                            nlamy = 20, spb = 0.3,
                            spy = 0.3,
                            lminfacX = 1e-6,
                            lminfacY = 1e-2, kfold = 10,
                            fullpath = 0, sigmafac = 2) {
  return(list(
    maxiter = maxiter, tol = tol, nlam = nlam, lminfac = lminfac,
    lmaxfac = lmaxfac,
    mu = mu, nu = nu, sp = sp, outMiter = outMiter, inMiter = inMiter,
    out.tol = out.tol, gamma = gamma, kmaxS = kmaxS, tolS = tolS,
    nlamx = nlamx, nlamy = nlamy, spb = spb, spy = spy,
    lminfacX = lminfacX, lminfacY = lminfacY, kfold = kfold,
    fullpath = fullpath, sigmafac = sigmafac
  ))
}



#' Compute solution path of constrained lasso.
#'
#' The model uses scaled lasoo approach for model selection.
#'
#' @param Xt CLR transformed predictor matrix.
#' @param y model response vector
#' @param C sub-compositional matrix
#' @param we specify weight of model parameter
#' @param control a list of internal parameters controlling the model fitting
#' @return
#'   \item{betapath}{solution path estimate}
#'   \item{beta}{model parameter estimate}
#' @export
#' @importFrom stats qnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib robregcc
#' @examples
#' library(robregcc)
#' library(magrittr)
#' 
#' data(simulate_robregcc)
#' X <- simulate_robregcc$X;
#' y <- simulate_robregcc$y
#' C <- simulate_robregcc$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
#' 
#' 
#' #
#' Xt <- cbind(1,X)                         # accounting for intercept in predictor
#' C <- cbind(0,C)                           # accounting for intercept in constraint
#' bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept 
#' 
#' # Non-robust regression
#' control <- robregcc_option(maxiter = 5000, tol = 1e-7, lminfac = 1e-12)
#' fit.path <- classo_path(Xt, y, C, we = bw, control = control)
classo_path <- function(Xt, y, C, we = NULL, control = list()) {
  ## centered Xt and centered y provided
  p <- ncol(Xt)
  if (is.null(control)) control <- robregcc::robregcc_option()
  if (is.null(we)) we <- rep(1, p)
  Xt <- crossprod(t(Xt),subtract(diag(1,p,p),svd(t(C))$u %>% tcrossprod()))
  fo <- classopath(Xt, y, C, we,control)
  return(fo)
}



#' Estimate parameters of linear regression model with compositional covariates using method suggested by Pixu shi.
#'
#' The model uses scaled lasoo approach for model selection.
#'
#' @param Xt CLR transformed predictor matrix.
#' @param y model response vector
#' @param C sub-compositional matrix
#' @param we specify weight of model parameter
#' @param type 1/2 for l1 / l2 loss in the model
#' @param control a list of internal parameters controlling the model fitting
#' @return
#'   \item{beta}{model parameter estimate}
#' @export
#' @importFrom stats qnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib robregcc
#' @examples
#' library(robregcc)
#' library(magrittr)
#' 
#' data(simulate_robregcc)
#' X <- simulate_robregcc$X;
#' y <- simulate_robregcc$y
#' C <- simulate_robregcc$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
#' 
#' # Predictor transformation due to compositional constraint:
#' Xt <- cbind(1,X)          # accounting for intercept in predictor
#' C <- cbind(0,C)            # accounting for intercept in constraint
#' bw <- c(0,rep(1,p))        # weight matrix to not penalize intercept 
#' 
#' # Non-robust regression, [Pixu Shi 2016]
#' control <- robregcc_option(maxiter = 5000, tol = 1e-7, lminfac = 1e-12)
#' fit.nr <- classo(Xt, y, C, we = bw, type = 1, control = control) 
#' @references
#' Shi, P., Zhang, A. and Li, H., 2016. \emph{Regression analysis for microbiome compositional data. The Annals of Applied Statistics}, 10(2), pp.1019-1040.
classo <- function(Xt, y, C, we = NULL,
                       type = 1, control = list()) {
  ## centered Xt and centered y provided
  p <- ncol(Xt)
  n <- nrow(Xt)
  k <- nrow(C)
  if (is.null(we)) we <- rep(1, p)
  cindex = which(we == 0)
  Xt <- crossprod(t(Xt),subtract(diag(1,p,p),svd(t(C))$u %>% tcrossprod()))
  xmean <- colMeans(Xt)
  xmean[cindex] <- 0
  sfac <- apply(Xt, 2, sd)
  sfac[sfac == 0] <- 1
  sfac[cindex] <- 1
  sfac <- drop(1 / sfac)
  
  ## scaling of X;
  Xt <- t(t(Xt) - xmean)  * rep(sfac, each = n)
  C <- C * rep(sfac, each = k)
  mu_y <- mean(y)
  y <- y - mu_y
  
  dtf <- cbind( seq(0, p, length.out = 2000),
    abs(fnk(seq(0, p, length.out = 2000), p)))
  
  kk <- dtf[which.min(dtf[, 2]), 1]
  lam0 <- sqrt(2 / n) * qnorm(1 - (kk / p))

  
  # return(c_ridge2( Xt, y, C, 5, 10, control))
  if (type == 1) {
    fo <-  classoshe(Xt, y, C, we, lam0, control)
  } else {
    fo <- classol2(Xt, y, C, we, lam0, control)
  }
  
  tem <- fo$beta[-1]*sfac[-1]
  fo$beta <- c(mu_y + fo$beta[1] - sum(tem*xmean[-1]),tem)
  
  class(fo) = "constrained lasso"
  return(fo)
}









#' Principal sensitivity component analysis with compositional covariates in sparse setting.
#'
#' Produce model and its residual estimate based on PCS analysis.
#'
#' @param X0 CLR transformed predictor matrix.
#' @param y0 model response vector
#' @param alp (0,0.5) fraction of data sample to be removed to generate subsample
#' @param cfac initial value of shift parameter for weight construction/initialization
#' @param b1 tukey bisquare function parameter producing desired breakdown point
#' @param cc1 tukey bisquare function parameter producing desired breakdown point
#' @param C sub-compositional matrix
#' @param we penalization index for model parameters beta
#' @param type 1/2 for l1 / l2 loss in the model
#' @param control a list of internal parameters controlling the model fitting
#' @return
#'   \item{betaf}{TModel parameter estimate}
#'   \item{residuals}{residual estimate}
#' @export
#' @import magrittr
#' @importFrom MASS ginv
#' @importFrom Rcpp evalCpp
#' @useDynLib robregcc
#' @examples  
#' 
#' library(robregcc)
#' library(magrittr)
#' 
#' data(simulate_robregcc)
#' X <- simulate_robregcc$X;
#' y <- simulate_robregcc$y
#' C <- simulate_robregcc$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
#' 
#' Xt <- cbind(1,X)  # include intercept in predictor
#' C <- cbind(0,C)    # include intercept in constraint
#' bw <- c(0,rep(1,p)) # weights not penalize intercept 
#' 
#' example_seed <- 2*p+1               
#' set.seed(example_seed) 
#' 
#' # Breakdown point for tukey Bisquare loss function 
#' b1 = 0.5                    # 50% breakdown point
#' cc1 =  1.567                # corresponding model parameter
#' b1 = 0.25; cc1 =  2.937   
#' 
#' # Initialization [PSC analysis for compositional data]
#' control <- robregcc_option(maxiter = 1000,
#'  tol = 1e-4,lminfac = 1e-7)
#' fit.init <- cpsc_sp(Xt, y,alp = 0.4, cfac = 2, b1 = b1,
#' cc1 = cc1,C,bw,1,control) 
#' 
#' @references
#' Mishra, A., Mueller, C.,(2019) \emph{Robust regression with compositional covariates. In prepration}. arXiv:1909.04990.
cpsc_sp <- function(X0, y0, alp = 0.4, cfac = 2, b1 = 0.25,
                    cc1 = 2.937, C = NULL,
                    we, type, control = list()) {
  # X0=Xt;y0 <- y;alp=0.4;cfac=2
  # H <- crossprod(X0) %>%
  #   MASS::ginv() %>%
  #   tcrossprod(X0, .) %>%
  #   tcrossprod(X0)
  # # ind <- 1-diag(H)<=(max((1-diag(H)))- sd(1-diag(H))/2)
  # ind <- order(1 - diag(H))[1:round(0.7 * nrow(X0))]
  X2 = X0; p <- ncol(X0); y2 <- y0; y0 <- y0/sd(y0)
  X0 <- crossprod(t(X0),subtract(diag(p),svd(t(C))$u %>% tcrossprod()))
  ind <- getLevIndex(X0)
  X1 <- X0[ind, ]
  y1 <- y0[ind]
  fo <- getscsfun.sp(X1, y1, alp, b1, cc1, C, we, type, control)
  X1 <- X0
  y1 <- y0
  # print(fo$scale)
  n1 <- nrow(X1)
  res0 <- y1 - X1 %*% fo$beta
  ind <- abs(res0) < cfac * fo$scale
  i <- 0
  err <- 1e6
  sc0 <- fo$scale
  while (i < 10 && sum(ind) > (n1 / 2) && err > 1e-2) {
    # print(i);print(dim(X1[ind,]))
    fo <- getscsfun.sp(X1[ind, ], y1[ind], alp, b1, cc1, C, we, type, control)
    res0 <- y1 - X1 %*% fo$beta
    ind <- abs(res0) < cfac * fo$scale
    i <- i + 1
    err <- abs(fo$scale - sc0)
    sc0 <- fo$scale
    # print(sum(ind))
    # print(sc0)
  }
  res0 <- y1 - X1 %*% fo$beta
  fo$inlier <- ind <- abs(res0) < cfac * fo$scale
  fo$scale <- fo$beta <- NULL
  
  fo$betaf <- sd(y2)*classo(X2[ind, ], y1[ind], C, we, type, control)$beta
  fo$residuals <- y2 - X2 %*% fo$betaf
  
  control$maxiter <- 2000
  control$tol <- 1e-10
  control$lminfac <- 1e-8
  # fo$betaR <- classo(X2[ind, ], y1[ind], C, we, 2, control)$beta
  fo$betaR <- sd(y2)*c_ridge2( X2[ind, ], y1[ind], C, 10, 50, control)$beta
  fo$residualR <- y2 - X2 %*% fo$betaR
  
  class(fo) = "initialization"
  return(fo)
}






#' Robust model estimation approach for regression with compositional covariates.
#'
#'
#' Fit regression model with compositional covariates for a  range of tuning parameter lambda. Model parameters is assumed to be sparse.
#'
#' @param X predictor matrix 
#' @param y phenotype/response vector
#' @param C conformable sub-compositional matrix
#' @param beta.init initial value of model parameter beta
#' @param gamma.init inital value of shift parameter gamma
#' @param cindex index of control (not penalized) variable in the model 
#' @param control a list of internal parameters controlling the model fitting
#' @param penalty.index a vector of length 2 specifying type of penalty for model parameter and shift parameter respectively. 1, 2, 3 corresponding to adaptive, soft and hard penalty
#' @param alpha elastic net penalty
#' @param verbose TRUE/FALSE for showing progress of the cross validation
#' @return
#'   \item{Method}{Type of penalty used}
#'   \item{betapath}{model parameter estimate along solution path}
#'   \item{gammapath}{shift parameter estimate along solution path}
#'   \item{lampath}{sequence of fitted lambda)}
#'   \item{k0}{scaling factor}
#'   \item{cver}{error from k fold cross validation }
#'   \item{selInd}{selected index from minimum and 1se rule cross validation error}
#'   \item{beta0}{beta estimate corresponding to selected index}
#'   \item{gamma0}{mean shift estimate corresponding to selected index}
#'   \item{residual0}{residual estimate corresponding to selected index}
#'   \item{inlier0}{inlier index corresponding to selected index}
#'   \item{betaE}{Post selection estimate corresponding to selected index}
#'   \item{residualE}{post selection residual corresponding to selected index}
#'   \item{inlierE}{post selection inlier index corresponding to selected index}
#' @export
#' @importFrom utils modifyList
#' @importFrom Rcpp evalCpp
#' @useDynLib robregcc
#' @examples
#' library(magrittr)
#' library(robregcc)
#' 
#' data(simulate_robregcc)
#' X <- simulate_robregcc$X;
#' y <- simulate_robregcc$y
#' C <- simulate_robregcc$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
#' 
#' Xt <- cbind(1,X)                         # accounting for intercept in predictor
#' C <- cbind(0,C)                           # accounting for intercept in constraint
#' bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept 
#' 
#' example_seed <- 2*p+1               
#' set.seed(example_seed) 
#' \donttest{
#' # Breakdown point for tukey Bisquare loss function 
#' b1 = 0.5                    # 50% breakdown point
#' cc1 =  1.567                # corresponding model parameter
#' b1 = 0.25; cc1 =  2.937   
#' 
#' # Initialization [PSC analysis for compositional data]
#' control <- robregcc_option(maxiter=1000,tol = 1e-4,lminfac = 1e-7)
#' fit.init <- cpsc_sp(Xt, y,alp = 0.4, cfac = 2, b1 = b1, 
#' cc1 = cc1,C,bw,1,control)  
#' 
#' ## Robust model fitting
#' 
#' # control parameters
#' control <- robregcc_option()
#' beta.wt <- fit.init$betaR    # Set weight for model parameter beta
#' beta.wt[1] <- 0
#' control$gamma = 1            # gamma for constructing  weighted penalty
#' control$spb = 40/p           # fraction of maximum non-zero model parameter beta
#' control$outMiter = 1000      # Outer loop iteration
#' control$inMiter = 3000       # Inner loop iteration
#' control$nlam = 50            # Number of tuning parameter lambda to be explored
#' control$lmaxfac = 1          # Parameter for constructing sequence of lambda
#' control$lminfac = 1e-8       # Parameter for constructing sequence of lambda 
#' control$tol = 1e-20;         # tolrence parameter for converging [inner  loop]
#' control$out.tol = 1e-16      # tolerence parameter for convergence [outer loop]
#' control$kfold = 10           # number of fold of crossvalidation
#' control$sigmafac = 2#1.345
#' # Robust regression using adaptive lasso penalty
#' fit.ada <- robregcc_sp(Xt,y,C,
#'                        beta.init = beta.wt,  cindex = 1, 
#'                        gamma.init = fit.init$residuals,
#'                        control = control, 
#'                        penalty.index = 1, alpha = 0.95)
#' 
#' # Robust regression using lasso penalty [Huber equivalent]   
#' fit.soft <- robregcc_sp(Xt,y,C, cindex = 1, 
#'                         control = control, penalty.index = 2, 
#'                         alpha = 0.95)
#' 
#' 
#' # Robust regression using hard thresholding penalty
#' control$lmaxfac = 1e2               # Parameter for constructing sequence of lambda
#' control$lminfac = 1e-3              # Parameter for constructing sequence of lambda
#' control$sigmafac = 2#1.345
#' fit.hard <- robregcc_sp(Xt,y,C, beta.init = fit.init$betaf, 
#'                         gamma.init = fit.init$residuals,
#'                         cindex = 1, 
#'                         control = control, penalty.index = 3, 
#'                         alpha = 0.95)
#'  }                       
#' @references
#' Mishra, A., Mueller, C.,(2019) \emph{Robust regression with compositional covariates. In prepration.} arXiv:1909.04990.
robregcc_sp <- function(X, y, C, beta.init = NULL, gamma.init = NULL,
                        cindex = 1,
                        control = list(), penalty.index = 3, alpha = 1,
                        verbose = TRUE) {
  truncLoss1 <- function(r, f) {
    r. <- r - stats::median(r)
    sigmae <- 1.483 * stats::median(abs(r.))  ## MAD : stats::median
    # sigmae <- (1.345 * sigmae)^2
    sigmae <- (f * sigmae)^2
    r2 <- r.^2
    gind <- r2 < sigmae
    sd(r[gind])
  }
  
  ## robust regression for mean centered data X and y
  p <- ncol(X); n <- nrow(X); k <- nrow(C)
  if (ncol(C) != p ) stop('Mismatch in the dimension of X and C')
  out <- list(X0 = X, y0 = y, C0 = C)
  X <- crossprod(t(X),subtract(diag(1,p,p),svd(t(C))$u %>% tcrossprod()))
  
  xmean <- colMeans(X)
  xmean[cindex] <- 0
  sfac <- apply(X, 2, sd)
  sfac[sfac == 0] <- 1
  sfac[cindex] <- 1
  sfac <- drop(1 / sfac)
  out$sfac <- sfac
  out$xmean <- xmean
  
  ## scaling of X;
  X <- t(t(X) - xmean)  * rep(sfac, each = n)
  C <- C * rep(sfac, each = k)
  mu_y <- mean(y)
  out$sdy <- sd(y)
  y <- (y - mu_y)/out$sdy
  out$X <- X
  out$y <- y
  out$C <- C
  out$ymean <- mu_y

  if (penalty.index == 1) {
    gamma.wt <- gamma.init
    beta.wt <- beta.init
    if (is.null(gamma.wt) || is.null(beta.wt)) {
      stop("Error: provide initial value of both shift and model parameter for weight 
           construction.")
    }
    if (length(beta.wt) != p) {
      stop("Weights not correctly specified:")
    }
    beta.wt[cindex] <- 0 
    beta.wt <- beta.wt / sfac
    ab <- shift4resid(gamma.wt)
    shwt <- abs(c(beta.wt, ab[,2] / sqrt(n)))/out$sdy
    shwt[shwt != 0] <- (shwt[shwt != 0])^(-control$gamma)
    gamma.init <- rep(1, n)
    beta.init <- rep(1, p)
    param.ini <- c( beta.init, gamma.init)
    # shwt <- shwt * sqrt(1 + log(n))
    }
  
  if (penalty.index == 2) {
    gamma.init <- rep(1, n)
    beta.init <- rep(1, p)
    param.ini <- c( beta.init, gamma.init)
    gamma.wt <- rep(1, n)
    beta.wt <- rep(1, p)
    beta.wt[cindex] <- 0
    shwt <- c(beta.wt, gamma.wt)
    # shwt <- shwt * sqrt(1 + log(n))
  }
  
  
  if (penalty.index == 3) {
    if (is.null(gamma.init) || is.null(beta.init)) {
      stop("Error: provide initial value of both shift and model parameter for weight initialization.")
    }
    # beta.init[1] <- beta.init[1] + mu_y - 
    # sum(xmean[-1]*beta.init[-1])
    
    beta.init[cindex] <- 0 
    beta.init <- beta.init / sfac
    ab <- shift4resid(gamma.init)
    # param.ini <- c(beta.init, gamma.init/ sqrt(n))/out$sdy
    param.ini <- c(beta.init, 0*ab[,1] + (ab[,2]/sqrt(n)) )/out$sdy
    gamma.wt <- rep(1, n)
    beta.wt <- rep(1, p)
    beta.wt[cindex] <- 0
    shwt <- c(beta.wt, gamma.wt)
    # shwt <- sqrt(shwt * sqrt(1 + log(n)))
    }
  shwt <- c(shwt[1:p]*sqrt(1 + log(p)), 
            shwt[p + (1:n)]*sqrt(1 + log(n)))
  
  
  tm <- abs(shwt) > 0
  lmax <- max(abs(c(crossprod(X, y)/n, y/sqrt(n))[tm] / shwt[tm]))
  ## scaling factor calculation
  k0 <- (svd(cbind(  rbind(X/sqrt(n), sqrt(control$mu) * C),
                     rbind( sqrt(1) * diag(n),
                            matrix(0, k, n)) ))$d[1])^2
  
  lampath <- exp(seq(log(lmax * control$lmaxfac),
                     log(lmax * control$lminfac),
                     length.out = control$nlam ))
  out$shwt <- shwt
  out$k0 <- k0
  out$lampath <- lampath / k0
  # Rcpp::sourceCpp('src/robregcc.cpp')
  out <- modifyList(out, robregcc_sp5(
    X, y, C, param.ini,
    control, shwt, lampath,
    penalty.index, k0, alpha ))
  out$lampathcv <- lampath[1:out$nlambda]
  control$fullpath <- 1
  
  
  ##  cross validation check of the model:
  kfold <- control$kfold
  out$sind <- sind <- sample(rep.int(1:kfold, 
                                     times = ceiling(n / kfold)), n)
  nspcv <- vector("list", kfold)
  bind <- rep(T, p)
  for (kind in 1:kfold) { # kind =2
    if (verbose) cat("Cross validation stage: ", kind, "\n")
    samInd <- sind != kind
    param.ini2 <- param.ini[c(bind, samInd)]
    nspcv[[kind]] <- robregcc_sp5(
      X[samInd, ], y[samInd], C, param.ini2,
      control, shwt[c(bind, samInd)], out$lampathcv,
      penalty.index, k0, alpha )
    nspcv[[kind]]$tre <- drop(y[samInd]) -
      (X[samInd, ] %*% nspcv[[kind]]$betapath +
         1 * nspcv[[kind]]$gammapath)
    nspcv[[kind]]$tresd <- drop(apply(nspcv[[kind]]$tre, 2, function(x) {
      lnc <- round(1 * length(x))
      sd((x[order(x)])[1:lnc])
    }))
    nspcv[[kind]]$te <- y[!samInd] - X[!samInd, ] %*%
      nspcv[[kind]]$betapath
    nspcv[[kind]]$scale.te <- nspcv[[kind]]$te /
      rep(nspcv[[kind]]$tresd, each = sum(!samInd))
  }
  out$cv <- nspcv
  
  #  Cross-validation implementation:
  ermat2 <- array(NA, dim = c(kfold, length(out$lampathcv)))
  for (i in 1:kfold) {
    nte <- dim(out$cv[[i]]$te)
    for (j in 1:nte[2]) {
      x <- out$cv[[i]]$te[, j] / out$cv[[i]]$tresd[j]
      ermat2[i, j] <- truncLoss1(x,control$sigmafac)
    }
    ermat2[i, ] <- abs(ermat2[i, ] - 1)
  }
  out$cver <- ermat2
  
  
  out$gind <- ((apply(abs(out$gammapath) > 0, 2, sum) != 0) *
                 (apply(abs(out$betapath) > 0, 2, sum) != 0)) != 0
  avger <- apply(ermat2, 2, mean, na.rm = T)
  avger <- avger * out$gind
  avger[avger == 0] <- NA
  
  sderr <- apply(ermat2, 2, sd, na.rm = T) / sqrt(kfold)
  sderr <- sderr * out$gind
  sderr[sderr == 0] <- NA
  
  dm1 <- max(which(avger == (min(avger, na.rm = T)), arr.ind = T))
  dm2 <- min(which(avger <= avger[dm1] + 1 * sderr[dm1]))
  out$selInd <- c(dm1, dm2)
  
  # Retrieval of the output 
  out$beta0 <- out$betapath[, out$selInd]
  out$gamma0 <- out$gammapath[, out$selInd]
  out$residuals0 <- drop(y) - X %*% out$beta0
  out$inlier0 <- abs(out$gamma0) < 0.01
  
  
  ## sanity check
  control <- robregcc_option(maxiter = 500,tol = 1e-6,
                             lminfac = 1e-7)
  control <- robregcc_option()
  bw2 <- rep(0, p)
  bw2[beta.wt != 0] <- 1
  out$beta1 <- with(out, cbind(
    classo(X0[inlier0[, 1], ], y0[inlier0[, 1]],
           C0,
           we = bw2, 1, control
    )$beta,
    classo(X0[inlier0[, 2], ], y0[inlier0[, 2]],
           C0,
           we = bw2, 1, control
    )$beta
  ))
  out$re1 <- with(out, drop(y0) - X0 %*% beta1)
  out$sd1 <- apply(out$re1 * out$inlier0, 2, function(x) sd(x[x != 0]))
  out$inlierE <- cbind(
    abs(out$re1[, 1]) < 3 * out$sd1[1],
    abs(out$re1[, 2]) < 3 * out$sd1[2]
  )
  
  
  ## final estimate
  out$betaE <- with(out, cbind(
    classo(X0[inlierE[, 1], ], y0[inlierE[, 1]], C0,
           we = bw2, 1, control
    )$beta,
    classo(X0[inlierE[, 2], ], y0[inlierE[, 2]], C0,
           we = bw2, 1, control
    )$beta
  ))
  out$residualsE <- drop(out$y0) - out$X0 %*% out$betaE
  out$sdE <- apply(out$residualsE * out$inlierE, 2, 
                   function(x) sd(x[x != 0]))
  
  ## rescaling of beta in the model 
  updatebeta = function(ymean, beta0, sfac, xmean, sdy){
    tem <- sdy*beta0[-1,]*sfac[-1]
    rbind(as.vector(ymean + sdy*beta0[1,drop = F] - 
                      crossprod(tem, xmean[-1])), tem)}
  
  out$beta0 <- with(out, updatebeta(ymean, beta0, sfac, xmean, sdy))
  out$residuals0 <- drop(out$y0) - out$X0 %*% out$beta0
  # out$beta1 <- with(out, updatebeta(ymean, beta1, sfac, xmean))
  # out$betaE <- with(out, updatebeta(ymean, betaE, sfac, xmean))
  out$X <- NULL; out$y <- NULL; out$C <- NULL; 
  names(out)[names(out) == "X0"] <- "X"
  names(out)[names(out) == "y0"] <- "y"
  names(out)[names(out) == "C0"] <- "C"
  out$lampath <- out$lampathcv
  out$gammapath <- out$sdy*out$gammapath
  out$gamma0 <- out$gammapath[, out$selInd]
  out$betapath <- with(out, updatebeta(ymean, betapath, sfac, xmean, sdy))
  
  out$sfac = NULL; out$sdy = NULL; out$ymean = NULL
  out$shwt = NULL; out$k0 = NULL; # out$Method = NULL
  out$cv = NULL; out$sind = NULL; out$gind = NULL
  out$xmean <- NULL; out$lampathcv <- NULL;
  out$beta1 = NULL; out$re1 = NULL; out$sd1 = NULL
  
  class(out) <- "robregcc"
  return(out)
  }






