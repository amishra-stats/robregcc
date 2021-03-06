% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robregcc.R
\name{robregcc_sp}
\alias{robregcc_sp}
\title{Robust model estimation approach for regression with compositional covariates.}
\usage{
robregcc_sp(
  X,
  y,
  C,
  beta.init = NULL,
  gamma.init = NULL,
  cindex = 1,
  control = list(),
  penalty.index = 3,
  alpha = 1,
  verbose = TRUE
)
}
\arguments{
\item{X}{predictor matrix}

\item{y}{phenotype/response vector}

\item{C}{conformable sub-compositional matrix}

\item{beta.init}{initial value of model parameter beta}

\item{gamma.init}{inital value of shift parameter gamma}

\item{cindex}{index of control (not penalized) variable in the model}

\item{control}{a list of internal parameters controlling the model fitting}

\item{penalty.index}{a vector of length 2 specifying type of penalty for model parameter and shift parameter respectively. 1, 2, 3 corresponding to adaptive, soft and hard penalty}

\item{alpha}{elastic net penalty}

\item{verbose}{TRUE/FALSE for showing progress of the cross validation}
}
\value{
\item{Method}{Type of penalty used}
  \item{betapath}{model parameter estimate along solution path}
  \item{gammapath}{shift parameter estimate along solution path}
  \item{lampath}{sequence of fitted lambda)}
  \item{k0}{scaling factor}
  \item{cver}{error from k fold cross validation }
  \item{selInd}{selected index from minimum and 1se rule cross validation error}
  \item{beta0}{beta estimate corresponding to selected index}
  \item{gamma0}{mean shift estimate corresponding to selected index}
  \item{residual0}{residual estimate corresponding to selected index}
  \item{inlier0}{inlier index corresponding to selected index}
  \item{betaE}{Post selection estimate corresponding to selected index}
  \item{residualE}{post selection residual corresponding to selected index}
  \item{inlierE}{post selection inlier index corresponding to selected index}
}
\description{
Fit regression model with compositional covariates for a  range of tuning parameter lambda. Model parameters is assumed to be sparse.
}
\examples{
library(magrittr)
library(robregcc)

data(simulate_robregcc)
X <- simulate_robregcc$X;
y <- simulate_robregcc$y
C <- simulate_robregcc$C
n <- nrow(X); p <- ncol(X); k <-  nrow(C)

Xt <- cbind(1,X)                         # accounting for intercept in predictor
C <- cbind(0,C)                           # accounting for intercept in constraint
bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept 

example_seed <- 2*p+1               
set.seed(example_seed) 
\donttest{
# Breakdown point for tukey Bisquare loss function 
b1 = 0.5                    # 50\% breakdown point
cc1 =  1.567                # corresponding model parameter
b1 = 0.25; cc1 =  2.937   

# Initialization [PSC analysis for compositional data]
control <- robregcc_option(maxiter=1000,tol = 1e-4,lminfac = 1e-7)
fit.init <- cpsc_sp(Xt, y,alp = 0.4, cfac = 2, b1 = b1, 
cc1 = cc1,C,bw,1,control)  

## Robust model fitting

# control parameters
control <- robregcc_option()
beta.wt <- fit.init$betaR    # Set weight for model parameter beta
beta.wt[1] <- 0
control$gamma = 1            # gamma for constructing  weighted penalty
control$spb = 40/p           # fraction of maximum non-zero model parameter beta
control$outMiter = 1000      # Outer loop iteration
control$inMiter = 3000       # Inner loop iteration
control$nlam = 50            # Number of tuning parameter lambda to be explored
control$lmaxfac = 1          # Parameter for constructing sequence of lambda
control$lminfac = 1e-8       # Parameter for constructing sequence of lambda 
control$tol = 1e-20;         # tolrence parameter for converging [inner  loop]
control$out.tol = 1e-16      # tolerence parameter for convergence [outer loop]
control$kfold = 10           # number of fold of crossvalidation
control$sigmafac = 2#1.345
# Robust regression using adaptive lasso penalty
fit.ada <- robregcc_sp(Xt,y,C,
                       beta.init = beta.wt,  cindex = 1, 
                       gamma.init = fit.init$residuals,
                       control = control, 
                       penalty.index = 1, alpha = 0.95)

# Robust regression using lasso penalty [Huber equivalent]   
fit.soft <- robregcc_sp(Xt,y,C, cindex = 1, 
                        control = control, penalty.index = 2, 
                        alpha = 0.95)


# Robust regression using hard thresholding penalty
control$lmaxfac = 1e2               # Parameter for constructing sequence of lambda
control$lminfac = 1e-3              # Parameter for constructing sequence of lambda
control$sigmafac = 2#1.345
fit.hard <- robregcc_sp(Xt,y,C, beta.init = fit.init$betaf, 
                        gamma.init = fit.init$residuals,
                        cindex = 1, 
                        control = control, penalty.index = 3, 
                        alpha = 0.95)
 }                       
}
\references{
Mishra, A., Mueller, C.,(2019) \emph{Robust regression with compositional covariates. In prepration.} arXiv:1909.04990.
}
