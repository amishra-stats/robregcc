if (getRversion() >= "3.5.0") utils::globalVariables(c("."))


#' Simulation data
#'
#' Simulate data for the robust regression with compositional covariates
#'
#' @param n sample size
#' @param betacc model parameter satisfying compositional covariates
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
#' 
#' library(robregcc)
#' library(magrittr)
#' 
#' ## n: sample size 
#' ## p: number of predictors
#' ## o: fraction of observations as outliers
#' ## L: {0,1} => leveraged {no, yes}, indicator variable for outlier type
#' ## shFac: multiplicative factor of true standard deviation by which O, 
#' ##         i.e., outliers fraction of observations are shifted. 
#' ## ngrp: number of subgroup in the model 
#' ## snr: noise to signal ratio for computing true standard deviation of error 
#' 
#' p <- 80                            
#' n <- 300                           
#' o <- 0.10                            
#' L <- 1                              
#' shFac <- 6       # shFac = {6,8} corresponds to {moderate, high} outlier 
#' ngrp <- 4                         
#' snr <- 3   
#' sp_beta <- 1
#' 
#' # Set seed for reproducibility 
#' example_seed <- 2*p+1               
#' set.seed(example_seed) 
#' 
#' ## 1. coefficient and subcomposition matrix C
#' if(sp_beta == 1){         ## sparse model coefficient matrix 
#'   #' subcomposition matrix C
#'   C1 <- matrix(0,ngrp,23)
#'   tind <- c(0,10,16,20,23)
#'   for(ii in 1:ngrp)
#'     C1[ii,(tind[ii]+1):tind[ii+1]] <- 1
#'   C <- matrix(0,ngrp,p)
#'   C[,1:ncol(C1)] <- C1            
#'   
#'   
#'   # model coefficient beta; Follow examples from [Pixu Shi 2016]
#'   beta <- c(1, - 0.8, 0.4, 0, 0, - 0.6, 0, 0, 0, 0, -1.5, 
#'             0, 1.2, 0, 0, 0.3)
#'   beta <- c(beta,rep(0,p-length(beta)))
#'   tcrossprod(C,t(beta)) ##' sanity check
#' }  else if(sp_beta == 0) { ## non sparse model coefficient matrix 
#'   # subcomposition matrix C
#'   j <- 1; C <- matrix(0,ngrp,p)
#'   for(ii in 1:ngrp){
#'     tv <-  min(c(round(ii*p/ngrp),p))
#'     C[ii,j:tv] <- 1
#'     j <- tv+1
#'   }
#'   
#'   # model coefficient beta;
#'   beta <- sample(c(1,-1),p,replace = T)*runif(p,.3,.4)
#'   beta <- svd(t(C))$u %>% tcrossprod() %>% 
#'     subtract(diag(p),.) %>% 
#'     tcrossprod(.,t(beta))
#'   tcrossprod(C,t(beta)) ## sanity check
#' }
#' # number of outliers
#' O <- o*n  
#' 
#' ## 2. simulate response and predictor matrix, i.e., X, y
#' Sigma  <- 1:p %>% outer(.,.,'-') %>% abs(); Sigma  <- 0.5^Sigma
#' data.case <- vector("list",1)
#' data.case <- robregcc_sim(n,beta,O = O,Sigma,levg = L, snr,shft = shFac,0,
#'                           C,out=data.case)
#' 
#' # We have saved a copy of simulated data in the package 
#' # with name simulate_robregcc_sp and simulate_robregcc_nsp
#' 
#' X <- data.case$X                          # predictor matrix
#' y <- data.case$y                          # model response 
#' 
#' 
#' @references
#' Mishra, A., Mueller, C.,(2019) \emph{Robust regression with compositional covariates. In prepration.} arXiv:1909.04990.
robregcc_sim <- function(n, betacc, O, Sigma, levg, snr,
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
  y <- tcrossprod(x, t(betacc))
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
    out$Yte <- tcrossprod(out$Xte, t(betacc)) + sigma * rnorm(m)
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
#' @return a list of controling parameter.
#' @export
#' @examples  
#' library(robregcc)
#' # default options
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
                            lminfacY = 1e-2, kfold = 5,
                            fullpath = 0) {
  return(list(
    maxiter = maxiter, tol = tol, nlam = nlam, lminfac = lminfac,
    lmaxfac = lmaxfac,
    mu = mu, nu = nu, sp = sp, outMiter = outMiter, inMiter = inMiter,
    out.tol = out.tol, gamma = gamma, kmaxS = kmaxS, tolS = tolS,
    nlamx = nlamx, nlamy = nlamy, spb = spb, spy = spy,
    lminfacX = lminfacX, lminfacY = lminfacY, kfold = kfold,
    fullpath = fullpath
  ))
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
#' 
#' library(robregcc)
#' library(magrittr)
#' 
#' data(simulate_robregcc_sp)
#' X <- simulate_robregcc_sp$X;
#' y <- simulate_robregcc_sp$y
#' C <- simulate_robregcc_sp$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
#' 
#' # Predictor transformation due to compositional constraint:
#' # Equivalent to performing centered log-ratio transform 
#' Xt <- svd(t(C))$u %>% tcrossprod() %>% subtract(diag(p),.) %>% crossprod(t(X),.)
#' #
#' Xm <- colMeans(Xt)
#' Xt <- scale(Xt,Xm,FALSE)                  # centering of predictors 
#' mean.y <- mean(y)
#' y <- y - mean.y                           # centering of response 
#' Xt <- cbind(1,Xt)                         # accounting for intercept in predictor
#' C <- cbind(0,C)                           # accounting for intercept in constraint
#' bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept 
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

  dtf <- cbind( seq(0, p, length.out = 2000),
    abs(fnk(seq(0, p, length.out = 2000), p)))
  
  kk <- dtf[which.min(dtf[, 2]), 1]
  lam0 <- sqrt(2 / n) * qnorm(1 - (kk / p))

  if (is.null(we)) we <- rep(1, p)
  # return(c_ridge2( Xt, y, C, 5, 10, control))
  if (type == 1) {
    fo <-  classoshe(Xt, y, C, we, lam0, control)
  } else {
    fo <- classol2(Xt, y, C, we, lam0, control)
  }
  
  class(fo) = "constrained lasso"
  return(fo)
}












#' Robust model estimation approach for regression with compositional covariates.
#'
#'
#' Fit regression model with compositional covariates for a  range of tuning parameter lambda. Model parameters is assumed to be sparse.
#'
#' @param X CLR transformed predictor matrix. Centered X
#' @param y model response vector, centered y
#' @param C sub-compositional matrix
#' @param beta.init initial value of model parameter beta
#' @param gamma.init inital value of shift parameter gamma
#' @param beta.wt specify weight of model parameter
#' @param gamma.wt initial value of shift parameter for weight construction/initialization
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
#' @export
#' @importFrom utils modifyList
#' @importFrom Rcpp evalCpp
#' @useDynLib robregcc
#' @examples
#' 
#' library(robregcc)
#' library(magrittr)
#' 
#' data(simulate_robregcc_sp)
#' X <- simulate_robregcc_sp$X;
#' y <- simulate_robregcc_sp$y
#' C <- simulate_robregcc_sp$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
#' 
#' # Predictor transformation due to compositional constraint:
#' # Equivalent to performing centered log-ratio transform 
#' Xt <- svd(t(C))$u %>% tcrossprod() %>% subtract(diag(p),.) %>% crossprod(t(X),.)
#' #
#' Xm <- colMeans(Xt)
#' Xt <- scale(Xt,Xm,FALSE)                  # centering of predictors 
#' mean.y <- mean(y)
#' y <- y - mean.y                           # centering of response 
#' Xt <- cbind(1,Xt)                         # accounting for intercept in predictor
#' C <- cbind(0,C)                           # accounting for intercept in constraint
#' bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept 
#' 
#' example_seed <- 2*p+1               
#' set.seed(example_seed) 
#' 
#' # Breakdown point for tukey Bisquare loss function 
#' b1 = 0.5                    # 50% breakdown point
#' cc1 =  1.567                # corresponding model parameter
#' # b1 = 0.25; cc1 =  2.937   
#' 
#' \donttest{
#' # Initialization [PSC analysis for compositional data]
#' control <- robregcc_option(maxiter=1000,tol = 1e-4,lminfac = 1e-7)
#' fit.init <- cpsc_sp(Xt, y,alp=0.4, cfac=2, b1=b1,cc1=cc1,C,bw,1,control) 
#' 
#' # Robust procedure
#' # control parameters
#' control <- robregcc_option()
#' beta.wt <- fit.init$betaR           # Set weight for model parameter beta
#' beta.wt[1] <- 0
#' control$gamma = 2                   # gamma for constructing  weighted penalty
#' control$spb = 40/p                  # fraction of maximum non-zero model parameter beta
#' control$outMiter = 1000             # Outer loop iteration
#' control$inMiter = 3000              # Inner loop iteration
#' control$nlam = 50                   # Number of tuning parameter lambda to be explored
#' control$lmaxfac = 1                 # Parameter for constructing sequence of lambda
#' control$lminfac = 1e-8              # Parameter for constructing sequence of lambda 
#' control$tol = 1e-20;                # tolrence parameter for converging [inner  loop]
#' control$out.tol = 1e-16             # tolerence parameter for convergence [outer loop]
#' control$kfold = 5                   # number of fold of crossvalidation
#' 
#' 
#' # Robust regression using adaptive elastic net penalty [case III, Table 1]
#' fit.ada <- robregcc_sp(Xt,y,C, beta.init=fit.init$betaR, 
#'                        gamma.init = fit.init$residualR,
#'                        beta.wt=abs(beta.wt), 
#'                        gamma.wt = abs(fit.init$residualR),
#'                        control = control, 
#'                        penalty.index = 1, alpha = 0.95) 
#'                        
#' # Robust regression using lasso penalty [Huber equivalent]   [case II, Table 1]
#' fit.soft <- robregcc_sp(Xt,y,C, beta.init=NULL, gamma.init = NULL,
#'                         beta.wt=bw, gamma.wt = NULL,
#'                         control = control, penalty.index = 2, 
#'                         alpha = 0.95)
#' 
#' 
#' # Robust regression using hard thresholding penalty [case I, Table 1]
#' control$lmaxfac = 1e2        # Parameter for constructing sequence of lambda
#' control$lminfac = 1e-3       # Parameter for constructing sequence of lambda
#' fit.hard <- robregcc_sp(Xt,y,C, beta.init=fit.init$betaf, 
#'                         gamma.init = fit.init$residuals,
#'                         beta.wt=bw, gamma.wt = NULL,
#'                         control = control, penalty.index = 3, 
#'                         alpha = 0.95)
#'                         
#'                         
#'  }                       
#' @references
#' Mishra, A., Mueller, C.,(2019) \emph{Robust regression with compositional covariates. In prepration.} arXiv:1909.04990.
robregcc_sp <- function(X, y, C, beta.init = NULL, gamma.init = NULL,
                        beta.wt = NULL, gamma.wt = NULL,
                        control = list(), penalty.index = 3, alpha = 1,
                        verbose = TRUE) {
  ## robust regression for mean centered data X and y

  p <- ncol(X)
  n <- nrow(X)
  k <- nrow(C)
  if (penalty.index == 1) {
    if (is.null(gamma.wt) || is.null(beta.wt)) {
      stop("Error: provide initial value of both shift and model parameter for weight 
           construction.")
    }
    if (length(beta.wt) != p) {
      stop("Weights not correctly specified:")
    }
    if (is.null(gamma.init)) {
      gamma.init <- rep(1, n)
    }
    if (is.null(beta.init)) {
      beta.init <- rep(1, p)
    }
  }

  if (penalty.index == 2) {
    gamma.init <- rep(1, n)
    beta.init <- rep(1, p)
    if (is.null(gamma.wt)) {
      gamma.wt <- rep(1, n)
    }
    if (is.null(beta.wt)) {
      beta.wt <- rep(1, p)
    }
  }


  if (penalty.index == 3) {
    if (is.null(gamma.init) || is.null(beta.init)) {
      stop("Error: provide initial value of both shift and model parameter for weight 
           initialization.")
    }
    if (is.null(gamma.wt)) {
      gamma.wt <- rep(1, n)
    }
    if (is.null(beta.wt)) {
      beta.wt <- rep(1, p)
    }
  }


  if (penalty.index == 4) {
    if (is.null(gamma.wt) || is.null(beta.wt)) {
      stop("Error: provide initial value of both shift and model parameter for weight 
           construction.")
    }
    if (length(beta.wt) != p) {
      stop("Weights not correctly specified:")
    }
    if (is.null(gamma.init) || is.null(beta.init)) {
      stop("Error: provide initial value of both shift and model parameter for weight 
           initialization.")
    }
  }

  out <- list()
  out$X0 <- X
  out$y0 <- y
  out$C0 <- C

  # y <- y/sqrt(n); X <- X/sqrt(n); gamma.wt <- gamma.wt/sqrt(n)
  # X <- Xt
  sfac <- apply(X, 2, sd)
  sfac[sfac == 0] <- 1
  sfac <- drop(1 / sfac)
  out$sfac <- sfac

  ## scaling of X
  X <- X * rep(sfac, each = n)
  C <- C * rep(sfac, each = k)
  X <- X / sqrt(n)
  y <- y / sqrt(n) #
  out$X <- X
  out$y <- y
  out$C <- C

  # parameter initialization
  out$param.ini <- param.ini <- c(
    beta.init / sfac,
    gamma.init / sqrt(n)
  )

  # Weight construction
  h2 <- tcrossprod(tcrossprod(X, MASS::ginv(crossprod(X))), X)
  # re <- y - h2%*%y
  if (penalty.index == 1) {
    shwt <- abs(c(beta.wt / sfac, gamma.wt / sqrt(n)))
    shwt[shwt != 0] <- (shwt[shwt != 0])^(-control$gamma)
    # shwt[shwt==0] <- .Machine$integer.max
    shwt <- shwt * sqrt(1 + log(n)) #* sqrt(1-diag(h2))
  }

  if (penalty.index == 2) {
    # shwt <- c(beta.wt, gamma.wt * sqrt(1 - diag(h2)))
    shwt <- c(beta.wt, gamma.wt * sqrt(1 - 0))
    shwt <- shwt * sqrt(1 + log(n))
  }

  if (penalty.index == 3) {
    # shwt = c(rep(max(abs(param.ini[1:p])),p), 
    # rep(max(abs(param.ini[(p+1):np])),n) )
    # shwt = abs(param.ini)
    # shwt[shwt!=0] <- (shwt[shwt!=0])^(-control$gamma)
    # shwt[shwt==0] <- .Machine$integer.max
    # shwt <- c(beta.wt, gamma.wt * sqrt(1 - diag(h2))) #
    shwt <- c(beta.wt, gamma.wt * sqrt(1 - 0))
    shwt <- sqrt(shwt * sqrt(1 + log(n)))
  }

  if (penalty.index == 4) {
    shwt <- abs(c(beta.wt / sfac, gamma.wt / sqrt(n)))
    shwt[shwt != 0] <- (shwt[shwt != 0])^(-control$gamma)
    # shwt = c(rep(max(abs(param.ini[1:p])),p), 
    # rep(max(abs(param.ini[(p+1):np])),n) )
    # shwt = abs(param.ini)
    # shwt[shwt!=0] <- (shwt[shwt!=0])^(-control$gamma)
    # shwt[shwt==0] <- .Machine$integer.max
    # shwt <- c(beta.wt,gamma.wt*sqrt(1-diag(h2))) #
    shwt <- sqrt(shwt * sqrt(1 + log(n)))
  }


  ## lampath calculation
  k0 <- (svd(cbind(
    rbind(X, sqrt(control$mu) * C),
    rbind(
      sqrt(1) * diag(n),
      matrix(0, k, n)
    )
  ))$d[1])^2

  tm <- abs(shwt) > 0
  lmax <- max(abs(c(crossprod(X, y), sqrt(1) * y)[tm] / shwt[tm]))
  #save(list = ls(), file = "xxx.rda")
  lampath <- exp(seq(log(lmax * control$lmaxfac),
    log(lmax * control$lminfac),
    length.out = control$nlam
  ))

  out$k0 <- k0
  out$shwt <- shwt
  out$lampath <- lampath / k0


  out <- modifyList(out, robregcc_sp5(
    X, y, C, param.ini,
    control, shwt, lampath,
    penalty.index, k0, alpha
  ))
  out$betapath <- apply(out$betapath, 2, function(x) x * sfac)
  out$gammapath <- sqrt(n) * out$gammapath
  out$lampathcv <- lampath[1:out$nlambda]

  control$fullpath <- 1


  ##  cross validation check of the model:
  kfold <- control$kfold
  out$sind <- sind <- sample(rep.int(1:kfold,
    times = ceiling(n / kfold)
  ), n)
  nspcv <- vector("list", kfold)
  ermat2 <- matrix(nrow = kfold, ncol = length(out$lampathcv))
  bind <- rep(T, p)
  for (kind in 1:kfold) { # kind =2
    if (verbose) cat("Cross validation stage: ", kind, "\n")
    samInd <- sind != kind
    param.ini2 <- param.ini[c(bind, samInd)]
    nspcv[[kind]] <- robregcc_sp5(
      X[samInd, ], y[samInd], C, param.ini2,
      control, shwt[c(bind, samInd)], out$lampathcv,
      penalty.index, k0, alpha
    )
    nspcv[[kind]]$betapath <- apply(
      nspcv[[kind]]$betapath, 2,
      function(x) x * sfac
    )
    nspcv[[kind]]$gammapath <- sqrt(n) * nspcv[[kind]]$gammapath

    nspcv[[kind]]$tre <- tre <- drop(out$y0[samInd]) -
      (out$X0[samInd, ] %*% nspcv[[kind]]$betapath +
        1 * nspcv[[kind]]$gammapath)
    nspcv[[kind]]$tresd <- drop(apply(tre, 2, function(x) {
      lnc <- round(1 * length(x))
      sd((x[order(x)])[1:lnc])
    }))
    nspcv[[kind]]$te <- out$y0[!samInd] - out$X0[!samInd, ] %*%
      nspcv[[kind]]$betapath
    nspcv[[kind]]$scale.te <- nspcv[[kind]]$te /
      rep(nspcv[[kind]]$tresd, each = sum(!samInd))
  }
  out$cv <- nspcv




  #  Cross-validation implementation:
  nfold <- length(out$cv)
  ermat2 <- array(NA, dim = c(nfold, length(out$lampathcv)))
  for (i in 1:nfold) {
    nte <- dim(out$cv[[i]]$te)
    for (j in 1:nte[2]) {
      x <- out$cv[[i]]$te[, j] / out$cv[[i]]$tresd[j]
      ermat2[i, j] <- truncLoss1(x)
    }
    ermat2[i, ] <- abs(ermat2[i, ] - 1)
  }
  out$cver <- ermat2

  out$gind <- ((apply(abs(out$gammapath) > 0, 2, sum) != 0) *
    (apply(abs(out$betapath) > 0, 2, sum) != 0)) != 0
  avger <- apply(ermat2, 2, mean, na.rm = T)
  avger <- avger * out$gind
  avger[avger == 0] <- NA

  sderr <- apply(ermat2, 2, sd, na.rm = T) / sqrt(nfold)
  sderr <- sderr * out$gind
  sderr[sderr == 0] <- NA


  dm1 <- max(which(avger == (min(avger, na.rm = T)), arr.ind = T))
  dm2 <- max(which(avger <= avger[dm1] + 1 * sderr[dm1]))
  out$selInd <- c(dm1, dm2)
  out$beta0 <- out$betapath[, out$selInd]
  out$gamma0 <- out$gammapath[, out$selInd]
  out$residuals0 <- drop(out$y0) - out$X0 %*% out$beta0
  # dm1 <- which(avger == (min(avger,na.rm = T)),arr.ind = T)[1]
  # out$selInd <- max(which(avger <= (min(avger,na.rm = T)) + 0*sd(avger,na.rm = T)))


  out$tData <- abs(out$gammapath[, out$selInd]) < 0.05
  # save(list = ls(), file = "aditya.rda")
  # print(1234)
  ## sanity check
  control <- robregcc_option(
    maxiter = 1000, tol = 1e-2,
    lminfac = 1e-7
  )
  # control <- robregcc_option()
  bw2 <- rep(0, p)
  bw2[beta.wt != 0] <- 1
  out$beta1 <- with(out, cbind(
    classo(X0[tData[, 1], ], y0[tData[, 1]],
      C0,
      we = bw2, 1, control
    )$beta,
    classo(X0[tData[, 2], ], y0[tData[, 2]],
      C0,
      we = bw2, 1, control
    )$beta
  ))
  out$re1 <- drop(out$y0) - out$X0 %*% out$beta1
  sd1 <- apply(out$re1 * out$tData, 2, function(x) sd(x[x != 0]))
  out$trueDataInd <- cbind(
    abs(out$re1[, 1]) < 3.1 * sd1[1],
    abs(out$re1[, 2]) < 3.1 * sd1[2]
  )


  ## final estimate
  out$betaE <- with(out, cbind(
    classo(X0[trueDataInd[, 1], ], y0[trueDataInd[, 1]], C0,
      we = bw2, 1, control
    )$beta,
    classo(X0[trueDataInd[, 2], ], y0[trueDataInd[, 2]], C0,
      we = bw2, 1, control
    )$beta
  ))
  out$residualsE <- drop(out$y0) - out$X0 %*% out$betaE

  names(out)[names(out) == "X0"] <- "X"
  names(out)[names(out) == "y0"] <- "y"
  names(out)[names(out) == "C0"] <- "C"

  class(out) <- "robregcc"
  return(out)
}















#' Robust model estimation approach for regression with compositional covariates.
#'
#' Generate solution path for range of lambda in case of where model parameter beta is not assumed to be sparse(nsp).
#'
#' @param X CLR transformed predictor matrix.
#' @param y model response vector
#' @param C sub-compositional constraint matrix
#' @param intercept true/false to include intercept term in the model
#' @param gamma.wt initial value of shift parameter for weight construction/initialization
#' @param control a list of internal parameters controlling the model fitting
#' @param penalty.index 1, 2, 3 corresponding to adaptive, soft and hard penalty
#' @param verbose TRUE/FALSE for showing progress of the cross validation
#' @return
#'   \item{Method}{Type of penalty used}
#'   \item{betapath}{model parameter estimate along solution path}
#'   \item{gammapath}{shift parameter estimate along solution path}
#'   \item{lampath}{sequence of fitted lambda)}
#'   \item{X}{predictors}
#'   \item{y}{response}
#' @export
#' @importFrom utils modifyList
#' @importFrom Rcpp evalCpp
#' @importFrom MASS ginv
#' @useDynLib robregcc
#' @examples
#' 
#' library(robregcc)
#' library(magrittr)
#' 
#' data(simulate_robregcc_nsp)
#' X <- simulate_robregcc_nsp$X;
#' y <- simulate_robregcc_nsp$y
#' C <- simulate_robregcc_nsp$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
#' 
#' # Predictor transformation due to compositional constraint:
#' # Equivalent to performing centered log-ratio transform 
#' Xt <- svd(t(C))$u %>% tcrossprod() %>% subtract(diag(p),.) %>% crossprod(t(X),.)
#' #
#' Xm <- colMeans(Xt)
#' Xt <- scale(Xt,Xm,FALSE)                  # centering of predictors 
#' mean.y <- mean(y)
#' y <- y - mean.y                           # centering of response 
#' Xt <- cbind(1,Xt)                         # accounting for intercept in predictor
#' C <- cbind(0,C)                           # accounting for intercept in constraint
#' bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept 
#' 
#' example_seed <- 2*p+1               
#' set.seed(example_seed) 
#' 
#' # Breakdown point for tukey Bisquare loss function 
#' b1 = 0.5                    # 50% breakdown point
#' cc1 =  1.567                # corresponding model parameter
#' # b1 = 0.25; cc1 =  2.937   
#' 
#' \donttest{
#' # Initialization [PSC analysis for compositional data]
#' control <- robregcc_option(maxiter=3000,tol = 1e-6)
#' fit.init  <- cpsc_nsp(Xt, y,alp=0.4,cfac=2,b1 = b1, cc1 = cc1,C,control)
#' 
#' 
#' # Robust procedure
#' # control parameters
#' control <- robregcc_option()
#' control$tol <- 1e-30
#' control$nlam = 25; 
#' control$lminfac = 1e-5;
#' control$outMiter = 10000
#' control$gamma <- 2
#' # Robust regression using adaptive elastic net penalty [case III, Table 1]
#' fit.ada <- robregcc_nsp(Xt,y, C, intercept = FALSE,  
#'                         gamma.wt = fit.init$residuals,
#'                         control = control, penalty.index = 1)
#' 
#' 
#' # Robust regression using elastic net penalty [case II, Table 1]
#' control$lminfac = 1e-1;
#' fit.soft <- robregcc_nsp(Xt,y,C,intercept = FALSE, gamma.wt = NULL,
#'                          control = control, penalty.index = 2)
#' 
#' 
#' 
#' # Robust regression using hard-ridge penalty [case I, Table 1]
#' control$tol <- 1e-30
#' control$nlam = 25; 
#' control$lminfac = 1e-1; 
#' control$outMiter = 10000
#' fit.hard <- robregcc_nsp(Xt,y,C, intercept = FALSE, 
#'                          gamma.wt = fit.init$residuals,
#'                          control = control, penalty.index = 3) 
#' 
#' }
#' @references
#' Mishra, A., Mueller, C.,(2019) \emph{Robust regression with compositional covariates. In prepration}. arXiv:1909.04990.
robregcc_nsp <- function(X, y, C, intercept = FALSE, gamma.wt = NULL,
                         control = list(), penalty.index = 1, verbose = TRUE) {
  ## robust regression for mean centered data X and y

  n <- nrow(X)
  int0 <- 0
  if (intercept) int0 <- 1
  if (penalty.index != 2) {
    if (is.null(gamma.wt)) {
      stop("Error: provide initial value of shift parameter for weight 
           construction/initialization.")
    }
  } else {
    gamma.wt <- rep(1, n)
  }
  out <- list()
  out$X0 <- X
  out$y0 <- y
  out$C0 <- C
  k <- nrow(C)
  sfac <- apply(X, 2, sd)
  sfac[sfac == 0] <- 1
  sfac <- drop(1 / sfac)
  out$sfac <- sfac
  X <- X * rep(sfac, each = n)
  C <- C * rep(sfac, each = k)
  X <- X / sqrt(n)
  y <- y / sqrt(n) ##
  out$X <- X
  out$y <- y
  out$C <- C
  gamma.wt <- gamma.wt / sqrt(n)
  
  
  ## generate lampath for the case:
  h2 <- tcrossprod(tcrossprod(X, MASS::ginv(crossprod(X))), X)
  re <- y - h2 %*% y

  shwt <- rep(1, n)
  if (penalty.index == 1) {
    shwt <- abs(gamma.wt)
    shwt[shwt != 0] <- (shwt[shwt != 0])^(-control$gamma)
    shwt <- shwt * sqrt(1 + log(n)) #* sqrt(1-diag(h2));
  }

  if (penalty.index == 2) {
    shwt <- shwt * sqrt(1 + log(n)) #* sqrt(1-diag(h2));
  }

  if (penalty.index == 3) {
    shwt <- sqrt(shwt * sqrt(1 + log(n))) #* sqrt(1-diag(h2));
  }
  tm <- abs(shwt) > 0
  lmax <- max(abs(re[tm] / shwt[tm]))
  lampath <- exp(seq(log(lmax), log(lmax * control$lminfac),
    length.out = control$nlam
  ))

  out <- modifyList(out, robregcc_nsp5(
    X, y, C, int0, gamma.wt,
    lampath, shwt, control, penalty.index
  ))
  out$betapath <- apply(out$betapath, 2, function(x) x * sfac)
  out$gammapath <- sqrt(n) * out$gammapath

  lampath <- out$lampath
  control$fullpath <- 1
  out$shwt <- shwt
  # out$X <- X ; out$y <- y

  ##  cross validation check of the model:
  kfold <- control$kfold
  out$sind <- sind <- sample(rep.int(1:kfold,
    times = ceiling(n / kfold)
  ), n)
  nspcv <- vector("list", kfold)
  ermat2 <- matrix(nrow = kfold, ncol = length(lampath))
  for (kind in 1:kfold) { # kind =2
    if (verbose) cat("Cross validation stage: ", kind, "\n")
    samInd <- sind != kind
    nspcv[[kind]] <- robregcc_nsp5(
      X[samInd, ], y[samInd], C, int0, gamma.wt[samInd], lampath,
      shwt[samInd], control, penalty.index
    )
    nspcv[[kind]]$betapath <- apply(
      nspcv[[kind]]$betapath, 2,
      function(x) x * sfac
    )
    nspcv[[kind]]$gammapath <- sqrt(n) * nspcv[[kind]]$gammapath

    nspcv[[kind]]$tre <- tre <- drop(out$y0[samInd]) -
      (out$X0[samInd, ] %*% nspcv[[kind]]$betapath +
        1 * nspcv[[kind]]$gammapath)


    nspcv[[kind]]$tresd <- drop(apply(tre, 2, function(x) {
      lnc <- round(1 * length(x))
      sd((x[order(x)])[1:lnc])
    }))
    nspcv[[kind]]$te <- out$y0[!samInd] - out$X0[!samInd, ] %*%
      nspcv[[kind]]$betapath
    nspcv[[kind]]$scale.te <- nspcv[[kind]]$te /
      rep(nspcv[[kind]]$tresd, each = sum(!samInd))
  }
  out$cv <- nspcv
  out$lampathcv <- lampath

  # gind <- out$gammapath
  # gind[abs(gind) < 0.01] <- 0
  # gind <- gind==0
  # out$betapathSub <- 0*out$betapath
  # for(i in 1:ncol(out$betapath))
  #   out$betapathSub[,i] <- crossprod(ginv(crossprod(out$X[gind[,i],])),
  #                                    crossprod(out$X[gind[,i],],out$y[gind[,i]]))


  #  Cross-validation implementation:
  nfold <- length(out$cv)
  ermat2 <- array(NA, dim = c(nfold, length(out$lampathcv)))
  for (i in 1:nfold) {
    nte <- dim(out$cv[[i]]$te)
    for (j in 1:nte[2]) {
      x <- out$cv[[i]]$te[, j] / out$cv[[i]]$tresd[j]
      ermat2[i, j] <- truncLoss1(x)
    }
    ermat2[i, ] <- abs(ermat2[i, ] - 1)
  }
  out$cver <- ermat2

  out$gind <- ((apply(abs(out$gammapath) > 0, 2, sum) != 0) *
    (apply(abs(out$betapath) > 0, 2, sum) != 0)) != 0
  avger <- apply(ermat2, 2, mean, na.rm = T)
  avger <- avger * out$gind
  avger[avger == 0] <- NA

  sderr <- apply(ermat2, 2, sd, na.rm = T) / sqrt(nfold)
  sderr <- sderr * out$gind
  sderr[sderr == 0] <- NA

  dm1 <- max(which(avger == (min(avger, na.rm = T)), arr.ind = T))
  dm2 <- max(which(avger <= avger[dm1] + 1 * sderr[dm1]))
  out$selInd <- c(dm1, dm2)
  out$beta0 <- out$betapath[, out$selInd]
  out$gamma0 <- out$gammapath[, out$selInd]
  out$residuals0 <- drop(out$y0) - out$X0 %*% out$beta0
  out$tData <- abs(out$gammapath[, out$selInd]) == 0

  ## sanity check
  out$beta1 <- with(out, cbind(
    crossprod(
      MASS::ginv(crossprod(X0[tData[, 1], ])),
      crossprod(X0[tData[, 1], ], y0[tData[, 1]])
    ),
    crossprod(
      MASS::ginv(crossprod(X0[tData[, 2], ])),
      crossprod(X0[tData[, 2], ], y0[tData[, 2]])
    )
  ))
  out$re1 <- drop(out$y0) - out$X0 %*% out$beta1
  sd1 <- apply(out$re1 * out$tData, 2, function(x) sd(x[x != 0]))
  out$trueDataInd <- cbind(
    abs(out$re1[, 1]) < 3 * sd1[1],
    abs(out$re1[, 2]) < 3 * sd1[2]
  )

  ## final estimate
  out$betaE <- with(out, cbind(
    crossprod(
      MASS::ginv(crossprod(X0[trueDataInd[, 1], ])),
      crossprod(X0[trueDataInd[, 1], ], y0[trueDataInd[, 1]])
    ),
    crossprod(
      MASS::ginv(crossprod(X0[trueDataInd[, 2], ])),
      crossprod(X0[trueDataInd[, 2], ], y0[trueDataInd[, 2]])
    )
  ))
  out$residualsE <- drop(out$y0) - out$X0 %*% out$betaE

  names(out)[names(out) == "X0"] <- "X"
  names(out)[names(out) == "y0"] <- "y"
  names(out)[names(out) == "C0"] <- "C"

  class(out) <- "robregcc"
  return(out)
}





#' Principal sensitivity component analysis with compositional covariates in non-sparse setting.
#'
#' Produce model and its residual estimate based in PCS analysis.
#'
#' @param X0 CLR transformed predictor matrix.
#' @param y0 model response vector
#' @param alp (0,0.5) fraction of data sample to be removed to generate subsample
#' @param cfac initial value of shift parameter for weight construction/initialization
#' @param b1 tukey bisquare function parameter producing desired breakdown point
#' @param cc1 tukey bisquare function parameter producing desired breakdown point
#' @param C sub-compositional matrix
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
#' data(simulate_robregcc_nsp)
#' X <- simulate_robregcc_nsp$X;
#' y <- simulate_robregcc_nsp$y
#' C <- simulate_robregcc_nsp$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
#' 
#' # Predictor transformation due to compositional constraint:
#' # Equivalent to performing centered log-ratio transform 
#' Xt <- svd(t(C))$u %>% tcrossprod() %>% subtract(diag(p),.) %>% crossprod(t(X),.)
#' #
#' Xm <- colMeans(Xt)
#' Xt <- scale(Xt,Xm,FALSE)                  # centering of predictors 
#' mean.y <- mean(y)
#' y <- y - mean.y                           # centering of response 
#' Xt <- cbind(1,Xt)                         # accounting for intercept in predictor
#' C <- cbind(0,C)                           # accounting for intercept in constraint
#' bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept 
#' 
#' example_seed <- 2*p+1               
#' set.seed(example_seed) 
#' 
#' # Breakdown point for tukey Bisquare loss function 
#' b1 = 0.5                    # 50% breakdown point
#' cc1 =  1.567                # corresponding model parameter
#' # b1 = 0.25; cc1 =  2.937   
#' 
#' # Initialization [PSC analysis for compositional data]
#' control <- robregcc_option(maxiter=3000,tol = 1e-6)
#' fit.init  <- cpsc_nsp(Xt, y,alp=0.4,cfac=2,b1 = b1, cc1 = cc1,C,control)
#' 
#' 
#' @references
#' Mishra, A., Mueller, C.,(2019) \emph{Robust regression with compositional covariates. In prepration}. arXiv:1909.04990.
cpsc_nsp <- function(X0, y0, alp = 0.4, cfac = 2, b1 = 0.25, cc1 = 2.937,
                     C = NULL, control = list()) {
  # if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
  # X1=Xt;y1 <- y;alp=0.4;cfac=2
  H <- crossprod(X0) %>%
    MASS::ginv() %>%
    tcrossprod(X0, .) %>%
    tcrossprod(X0)
  ind <- order(1 - diag(H))[1:round(0.6 * nrow(X0))]
  # ind <- 1-diag(H)<=(max((1-diag(H)))- sd(1-diag(H))/2)
  # ind <- rep(T,nrow(X0))
  X1 <- X0[ind, ]
  y1 <- y0[ind]
  dm1 <- dim(X1)
  # if (dm1[1] <= dm1[2]) {
  #   fo <- getscsfun.sp(X1, y1, alp, b1, cc1, C, control)
  # } else {
  fo <- getscsfun(X1, y1, alp, b1, cc1)
  # }
  X1 <- X0
  y1 <- y0
  n1 <- nrow(X1)
  res0 <- y1 - X1 %*% fo$beta
  ind <- sind <- abs(res0) < cfac * fo$scale
  i <- 0
  err <- 1e6
  sc0 <- fo$scale
  # print(sc0)
  while (i < 200 && sum(ind) > (n1 / 2) && err > 1e-3) {
    # print(i);print(dim(X1[ind,]))
    dm1 <- dim(X1[ind, ])
    # if (dm1[1] <= dm1[2]) {
    #   fo <- getscsfun.sp(X1[ind, ], y1[ind], alp, b1, cc1, C, control)
    # } else {
    fo <- getscsfun(X1[ind, ], y1[ind], alp, b1, cc1)
    # }
    res0 <- y1 - X1 %*% fo$beta
    ind <- abs(res0) < cfac * fo$scale
    i <- i + 1
    err <- abs(fo$scale - sc0)
    sc0 <- fo$scale
    # print(sc0)
    # print(sum(ind))
  }
  res0 <- y1 - X1 %*% fo$beta
  ind <- abs(res0) < cfac * fo$scale
  fo$scale <- fo$beta <- NULL

  dm1 <- dim(X1[ind, ])
  # if (dm1[1] <= dm1[2]) {
  #   fo$betaf <- classo(X1[ind, ], y1[ind], C, we = NULL, control)$beta
  # } else {
  fo$betaf <- crossprod(
    MASS::ginv(crossprod(X1[ind, ])),
    crossprod(X1[ind, ], y1[ind])
  )
  # }

  fo$residuals <- y0 - X0 %*% fo$betaf
  
  class(fo) = "initialization"
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
#' data(simulate_robregcc_sp)
#' X <- simulate_robregcc_sp$X;
#' y <- simulate_robregcc_sp$y
#' C <- simulate_robregcc_sp$C
#' n <- nrow(X); p <- ncol(X); k <-  nrow(C)
#' 
#' # Predictor transformation due to compositional constraint:
#' # Equivalent to performing centered log-ratio transform 
#' Xt <- svd(t(C))$u %>% tcrossprod() %>% subtract(diag(p),.) %>% crossprod(t(X),.)
#' #
#' Xm <- colMeans(Xt)
#' Xt <- scale(Xt,Xm,FALSE)                  # centering of predictors 
#' mean.y <- mean(y)
#' y <- y - mean.y                           # centering of response 
#' Xt <- cbind(1,Xt)                         # accounting for intercept in predictor
#' C <- cbind(0,C)                           # accounting for intercept in constraint
#' bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept 
#' 
#' example_seed <- 2*p+1               
#' set.seed(example_seed) 
#' 
#' # Breakdown point for tukey Bisquare loss function 
#' b1 = 0.5                    # 50% breakdown point
#' cc1 =  1.567                # corresponding model parameter
#' # b1 = 0.25; cc1 =  2.937   
#' 
#' \donttest{
#' # Initialization [PSC analysis for compositional data]
#' control <- robregcc_option(maxiter=1000,tol = 1e-4,lminfac = 1e-7)
#' fit.init <- cpsc_sp(Xt, y,alp=0.4, cfac=2, b1=b1,cc1=cc1,C,bw,1,control)  
#' 
#' 
#' 
#' }
#' @references
#' Mishra, A., Mueller, C.,(2019) \emph{Robust regression with compositional covariates. In prepration}. arXiv:1909.04990.
cpsc_sp <- function(X0, y0, alp = 0.4, cfac = 2, b1 = 0.25,
                        cc1 = 2.937, C = NULL,
                        we, type, control = list()) {
  # X0=Xt;y0 <- y;alp=0.4;cfac=2
  H <- crossprod(X0) %>%
    MASS::ginv() %>%
    tcrossprod(X0, .) %>%
    tcrossprod(X0)
  # ind <- 1-diag(H)<=(max((1-diag(H)))- sd(1-diag(H))/2)
  ind <- order(1 - diag(H))[1:round(0.7 * nrow(X0))]
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
  while (i < 20 && sum(ind) > (n1 / 2) && err > 1e-2) {
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
  fo$ind <- ind <- abs(res0) < cfac * fo$scale
  fo$scale <- fo$beta <- NULL
  fo$betaf <- classo(
    X1[ind, ], y1[ind], C, we, type,
    control
  )$beta


  fo$residuals <- y0 - X0 %*% fo$betaf
  control$maxiter <- 10000
  control$tol <- 1e-20
  control$lminfac <- 1e-8
  fo$betaR <- c_ridge2(
    X1[ind, ], y1[ind], C, 10,
    control$nlam, control
  )$beta
  fo$residualR <- y0 - X0 %*% fo$betaR
  
  class(fo) = "initialization"
  return(fo)
}


