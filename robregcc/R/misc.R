## Truncated loss function for calculatin standard deviations
truncLoss1 <- function(r) {
  r. <- r - stats::median(r)
  sigmae <- 1.483 * stats::median(abs(r.))  ## MAD : stats::median absolute deviation
  sigmae <- (1.345 * sigmae)^2
  r2 <- r.^2
  gind <- r2 < sigmae
  sd(r[gind])
}




## rho function for the calculation of the error function
rho <- function(u, cc = 1.56) {
  w <- abs(u) <= cc
  v <- (u^2 / (2) * (1 - (u^2 / (cc^2)) + (u^4 / (3 * cc^4)))) * w +
    (1 - w) * (cc^2 / 6)
  v <- v * 6 / cc^2
  return(v)
}




## Calculate M-scale estimate 
scale1 <- function(u, b, cc, initial.sc = stats::median(abs(u)) / .6745) {
  max.it <- 200
  sc <- initial.sc
  i <- 0
  eps <- 1e-8
  # magic number alert
  err <- 1
  while (((i <- i + 1) < max.it) && (err > eps)) {
    sc2 <- sqrt(sc^2 * mean(rho(u / sc, cc)) / b)
    err <- abs(sc2 / sc - 1)
    sc <- sc2
  }
  return(sc)
}




# Subfunction for principal sensitive component analysis:
# 
# Subfubction PCS non-sparse.
# 
# @param Xa CLR transformed predictor matrix.
# @param ya model response vector
# @param alp0 (0,0.5) fraction of data sample to be removed to generate subsample
# @param b1 tukey bisquare function parameter producing desired breakdown point
# @param cc1 tukey bisquare function parameter producing desired breakdown point
# @return
#   \item{beta}{Model parameter estimate}
#   \item{scale}{scale estimate}
# @references
# Mishra, A., Mueller, C.,(2018) \emph{Robust regression with compositional covariates. In prepration}.
#' @importFrom MASS ginv
#' @import magrittr
#' @useDynLib robregcc
getscsfun <- function(Xa, ya, alp0 = 0.4, b1 = 0.25, cc1 = 2.937) {
  # if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
  # ginv2 = function(x){
  #   svdx <- svd(x)
  #   gind <- svdx$d > 1e-6
  #   svdx$u[,gind]%*%diag(1/svdx$d[gind],nrow = sum(gind))%*%t(svdx$u[,gind])
  # }
  H <- crossprod(Xa) %>%
    MASS::ginv() %>%
    tcrossprod(Xa, .) %>%
    tcrossprod(Xa)
  # print(dim)
  e <- ya - crossprod(H, ya)
  W <- diag(drop(e / (1 - diag(H))))
  # print((1-diag(H)))
  P <- tcrossprod(crossprod(H, W^2), H)
  svdp <- svd(P)
  U <- svdp$u[, abs(svdp$d) > 1e-4]
  pu <- ncol(U)
  pa <- ncol(Xa)
  n1 <- nrow(Xa)
  alp <- round(alp0 * n1)
  betaOLSmat <- array(0, dim = c(pa, 3, pu))
  for (i in 1:pu) {
    ord <- vector("list", 3)
    tord <- order(U[, i])
    ord[[1]] <- tord[(alp + 1):n1]
    ord[[2]] <- tord[1:(n1 - alp)]
    tord <- order(abs(U[, i]))
    ord[[3]] <- tord[1:(n1 - alp)]
    for (j in 1:3) {
      # betaOLSmat[,j,i] <- c_ridge(Xt[ord[[j]],],y[ord[[j]]],C,5,40)
      betaOLSmat[, j, i] <- crossprod(
        MASS::ginv(crossprod(Xa[ord[[j]], ])),
        crossprod(
          Xa[ord[[j]], ],
          ya[ord[[j]]]
        )
      )
    }
  }
  scmat <- apply(betaOLSmat, 2:3, function(x) {
    res <- ya - Xa %*% x
    scale1(res,
           b = b1, cc = cc1,
           initial.sc = stats::median(abs(res)) / .6745
    )
    # mscale(res)
  })
  beta0 <- crossprod(MASS::ginv(crossprod(Xa)), crossprod(Xa, ya))
  res <- ya - Xa %*% beta0
  # scale0 <- mscale(res)
  scale0 <- scale1(res,
                   b = b1, cc = cc1, initial.sc =
                     stats::median(abs(res)) / .6745
  )
  ab <- which(scmat == min(scmat, scale0), arr.ind = T)
  if (dim(ab)[1] == 0) {
    out <- list(beta = beta0, scale = scale0)
  } else {
    out <- list(
      beta = betaOLSmat[, ab[1, 1], ab[1, 2]],
      scale = scmat[ab[1, 1], ab[1, 2]]
    )
  }
  out
}





# Subfunction for principal sensitive component analysis (sparsity):
#'
# Subfubction PCS sparse.
# 
# @param Xa CLR transformed predictor matrix.
# @param ya model response vector
# @param alp0 (0,0.5) fraction of data sample to be removed to generate subsample
# @param b1 tukey bisquare function parameter producing desired breakdown point
# @param cc1 tukey bisquare function parameter producing desired breakdown point
# @param C sub-compositional matrix
# @param we penalization index for model parameters beta
# @param type 1/2 for l1 / l2 loss in the model
# @param control a list of internal parameters controlling the model fitting
# @return
#   \item{beta}{Model parameter estimate}
#   \item{scale}{scale estimate}
# @references
# Mishra, A., Mueller, C.,(2018) \emph{Robust regression with compositional covariates. In prepration}.
#' @importFrom MASS ginv
#' @import magrittr
#' @useDynLib robregcc
getscsfun.sp <- function(Xa, ya, alp0 = 0.4, b1 = 0.25,
                         cc1 = 2.937, C = NULL,
                         we, type, control = list()) {
  
  # ginv2 = function(x){
  #   svdx <- svd(x)
  #   gind <- svdx$d > 1e-6
  #   svdx$u[,gind]%*%diag(1/svdx$d[gind],nrow =
  #   sum(gind))%*%t(svdx$u[,gind])
  # }
  beta0 <- classo(Xa, ya, C, we, type, control)$beta
  e <- ya - Xa %*% beta0
  # H <- crossprod(Xa)%>%ginv()%>%tcrossprod(Xa,.)%>%tcrossprod(Xa)
  # e <- ya-crossprod(H,ya)
  Xa2 <- Xa[, beta0 != 0]
  H <- crossprod(Xa2) %>%
    MASS::ginv() %>%
    tcrossprod(Xa2, .) %>%
    tcrossprod(Xa2)
  # e <- ya-crossprod(H,ya)
  W <- diag(drop(e / (1 - diag(H))))
  P <- tcrossprod(crossprod(H, W^2), H)
  svdp <- svd(P)
  U <- svdp$u[, abs(svdp$d) > 1e-6]
  pu <- ncol(U)
  pa <- ncol(Xa)
  n1 <- nrow(Xa)
  alp <- round(alp0 * n1)
  betaOLSmat <- array(0, dim = c(pa, 3, pu))
  for (i in 1:pu) {
    ord <- vector("list", 3)
    tord <- order(U[, i])
    ord[[1]] <- tord[(alp + 1):n1]
    ord[[2]] <- tord[1:(n1 - alp)]
    tord <- order(abs(U[, i]))
    ord[[3]] <- tord[1:(n1 - alp)]
    for (j in 1:3) {
      betaOLSmat[, j, i] <- classo(
        Xa[ord[[j]], ], ya[ord[[j]]], C, we, type,
        control
      )$beta
      # betaOLSmat[,j,i] <- crossprod(ginv(crossprod(Xa[ord[[j]],])),
      # crossprod(Xa[ord[[j]],],ya[ord[[j]]]))
    }
  }
  # save(list=ls(),file = "aditya.rda")
  
  
  scmat <- apply(betaOLSmat, 2:3, function(x) {
    x[is.na(x)] <- 0
    res <- ya - Xa %*% x
    scale1(res,
           b = b1, cc = cc1, initial.sc =
             stats::median(abs(res)) / .6745
    )
    # mscale(res)
  })
  # beta0 <- crossprod(ginv(crossprod(Xa)),crossprod(Xa,ya))
  # beta0 <- classo(Xa,ya,C,we=NULL,control)$beta
  # res <- ya-Xa%*%beta0
  # scale0 <- mscale(res)
  scale0 <- scale1(e,
                   b = b1, cc = cc1,
                   initial.sc = stats::median(abs(e)) / .6745
  )
  ab <- which(scmat == min(scmat, scale0), arr.ind = T)
  if (dim(ab)[1] == 0) {
    out <- list(beta = beta0, scale = scale0)
  } else {
    out <- list(
      beta = betaOLSmat[, ab[1, 1], ab[1, 2]], scale =
        scmat[ab[1, 1], ab[1, 2]]
    )
  }
  out
}




## Generate leveraged observations from x2;
## O2 is the index of observations in x2 which will be leveraged 
## C is the compositional costraint matrix; 
getLevObs <- function(x2, O2, C) {
  # x2 <- Xt
  n <- nrow(x2)
  levObs <- apply(x2, 2, order)
  for (i in 1:ncol(x2)) {
    x2[, i] <- x2[levObs[, i], i]
  }
  out2 <- x2[1:O2, ]
  x2 <- x2[n:1, ] + 4
  for (im in 1:nrow(C)) {
    gind <- C[im, ] == 1
    out2[, gind][, 1] <- x2[, gind][1:O2, 1]
  }
  out2
}



fnk <- function(kkkt, p) { ## Gradient of 'fr'
  (qnorm(1 - (kkkt / p)))^4 + 2 * ((qnorm(1 -
    (kkkt / p)))^2) - kkkt
}











