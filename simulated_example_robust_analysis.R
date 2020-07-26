###########################################################
## R code for simulation examples                        ##
## in "Robust regression with Compositional Covariates"  ##
###########################################################


rm(list = ls())
## Install and load the package "robregcc" from CRAN/GitHub
library(devtools)
# devtools::install_github("amishra-stats/robregcc/robregcc")
library(robregcc)

## Load required packages
library(MASS)
require(ggplot2)
require(reshape2)
library(magrittr)
library(graphics)


######################################
## R code for simulation examples   ##
######################################

## -------------------Simulation setting-------------------
## n: sample size 
## p: number of predictors
## o: fraction of observations as outliers
## L: {0,1} => leveraged {no, yes}, indicator variable for outlier type
## s: multiplicative factor of true standard deviation by which O, i.e., outliers fraction of observations are shifted. 
## ngrp: number of subgroup in the model 
## snr: noise to signal ratio for computing true standard deviation of error 
## nsim:  number of simulation run 


## For each setting, generate 100 simulation run
n <- 200
p <- c(100, 300)
s <- c(6, 8)
L <- c(0, 1)
o <- c(0.05, 0.10, 0.15, 0.20)
nsim <- 100
# Begin simulation run setting
j = 1
sim_setting <- NULL
for (tem_p in p) {
  for (i in 1:nsim) {
    sim_setting <- rbind(sim_setting, c(n, tem_p, 0, 8, 0, i, j))
    j <- j + 1
  }
}
for (tem_p in p) {
  for (tem_L in L) {
    for (tem_s in s) {
      for (tem_o in o) {
        for (i in 1:nsim) {
          sim_setting <- rbind(sim_setting, c(n, tem_p, tem_L, tem_s, 
                                              tem_o*n, i, j))
          j <- j + 1
        }
      }
    }
  }
}
sim_setting <- data.frame(sim_setting)
names(sim_setting) <- c('n','p','L','s','O','index','example_seed')
# End simulation run setting




## ---------------------- single replication ------------------------
## For a chosen setting_index execute following chunk of code;
## 
## 
setting_index <- 1600
setting_eval <- paste(paste(names(sim_setting), 
                            sim_setting[setting_index,] , 
                            sep = ' = '),collapse = ';')
# Evaluate expression to get simulation setting;
eval(parse(text = setting_eval))
# n = 200;p = 100;L = 1;s = 8;O = 0.1;index = 100;example_seed = 1600


#Additional parameter 
ngrp <- 4                           # number of sub-composition
snr <- 3                            # Signal to noise ratio

#--------------------------------------------------------------
## Simulate true model variables, i.e., y, X, C, beta
## Follow examples from [Pixu Shi 2016]

## -----------------------Generate data--------------------
## 1. coefficient and subcomposition matrix C
# Simulate subcomposition matrix
C1 <- matrix(0,ngrp,23)
tind <- c(0,10,16,20,23)
for (ii in 1:ngrp)
  C1[ii,(tind[ii] + 1):tind[ii + 1]] <- 1
C <- matrix(0,ngrp,p)
C[,1:ncol(C1)] <- C1            

# model parameter beta
beta0 <- 0.5
beta <- c(1, -0.8, 0.4, 0, 0, -0.6, 0, 0, 0, 0, -1.5, 0, 1.2, 0, 0, 0.3)
beta <- c(beta,rep(0,p - length(beta)))

# Simulate response and predictor, i.e., X, y
Sigma  <- 1:p %>% outer(.,.,'-') %>% abs(); Sigma  <- 0.5^Sigma
data.case <- vector("list",1)
set.seed(example_seed)
data.case <- robregcc_sim(n,beta,beta0, O = O,Sigma,levg = L, 
                          snr,shft = s,0, C,out = data.case)

#--------------------------------------------------------------
## Data preprocessing:
X <- data.case$X                          # predictor matrix
y <- data.case$y                          # model response 
# Account for intercept in the model
Xt <- cbind(1,X)                          # accounting for intercept in predictor
C <- cbind(0,C)                           # accounting for intercept in constraint
bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept 


set.seed(example_seed)        # unique seed
#--------------------------------------------------------------
# Non-robust regression, [Pixu Shi 2016]
control <- robregcc_option(maxiter = 5000, tol = 1e-8, lminfac = 1e-8)
fit.nr <- classo(Xt, y, C, we = bw, type = 1, control = control) 
if (O > 0) {
  fit.oracle <- classo(Xt[-(1:O),], data.case$yo[-(1:O)], C, 
                       we = bw, type = 1, control = control) 
} else {
  fit.oracle <- fit.nr
}



#--------------------------------------------------------------
## Initialization
# control parameter 
# Breakdown point for tukey Bisquare loss function 
# b1 = 0.5; cc1 =  1.567             # 50% breakdown point
b1 = 0.25; cc1 =  2.937               # initalization for scale parameter 
control <- robregcc_option(maxiter = 1000,tol = 1e-4,lminfac = 1e-7)
fit.init <- cpsc_sp(Xt, y,alp = 0.4, cfac = 2, b1 = b1,
                    cc1 = cc1,C,bw,1,control)  


#--------------------------------------------------------------
## Robust model fitting

# control parameters
control <- robregcc_option()
beta.wt <- fit.init$betaR          # Set weight for model parameter beta
beta.wt[1] <- 0
control$gamma = 1                   # gamma for constructing  weighted penalty
control$spb = 40/p                  # fraction of maximum non-zero model parameter beta
control$outMiter = 1000             # Outer loop iteration
control$inMiter = 3000              # Inner loop iteration
control$nlam = 50                   # Number of tuning parameter lambda to be explored
control$lmaxfac = 1                 # Parameter for constructing sequence of lambda
control$lminfac = 1e-8              # Parameter for constructing sequence of lambda 
control$tol = 1e-20;                # tolrence parameter for converging [inner  loop]
control$out.tol = 1e-16             # tolerence parameter for convergence [outer loop]
control$kfold = 10                   # number of fold of crossvalidation
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




## -----------------------Figure 2 in the manuscript ------------
fdir <-  getwd()
pdf(file = file.path(fdir,'cv_sim.pdf'), width = 15, height = 4)
par(mfrow = c(1,3), mar = 5*c(1,1,0.5,0.1))
plot_path(fit.ada)
plot_cv(fit.ada)
plot_resid(fit.ada, 1, s = 0)
dev.off()
# 



## -----------------------Model comparisons  ------------
## A function to compute the error measures for comparing different methods
## Table 2 (s = 8) in the manuscript and Table 1 (s = 6) in the supplementary material 
## Figure 3 compares the misspecification measure using HM 
## model_compare : provides metric to compare models 

# fit :  fit objects obtained using robregcc
# O: Number of outliers in the simulated data 
# beta: true parameter beta 
# s_index: 1/2 corresponding to minimum/1se rule cross validation error 
model_compare = function(fit, O, beta,s_index){
  outlier_identifier <- !fit$inlierE[,s_index]
  if (O) { 
    # O = 0 case
    FP2 <- sum(outlier_identifier[-c(1:O)])
    outlier_identifier <- abs(fit$gamma0[,s_index]) > 0.001
    
    FN <- sum(!outlier_identifier[1:O])
    FP1 <- sum(outlier_identifier[-c(1:O)])
    HM <- sum(!outlier_identifier[1:O]) + sum(outlier_identifier[-c(1:O)])
  } else { 
    # O = 0 case
    FP2 <- sum(outlier_identifier)
    outlier_identifier <- abs(fit$gamma0[,s_index]) > 0.001
    
    FN <- 0
    FP1 <- sum(outlier_identifier)
    HM <- FN + FP1
  }
  error_beta <- sqrt(sum((fit$betaE[,s_index] - beta)^2))/length(beta)
  return(c(error_beta, FN, FP1, HM, FP2))
}

beta_oracle <- c(beta0,beta)
error_beta_nr <- sqrt(sum(((beta_oracle - fit.nr$beta))^2))/length(beta_oracle)

s_index <- 1  # 1/2 corresponding to minimum/1se rule cross validation error 
compare.fit <- rbind(c(1,model_compare(fit.ada,O,beta_oracle,s_index)),
                     c(2,model_compare(fit.soft,O,beta_oracle,s_index)),
                     c(3,model_compare(fit.hard,O,beta_oracle,s_index)),
                     c(4,error_beta_nr,NA,NA,NA,NA))
colnames(compare.fit) <- c('Method', 'Er_beta', 'FN', 'FP_1','HM', 'FP_2')
compare.fit <- data.frame(compare.fit)
compare.fit$Method <- c('A','E','H','NR')[compare.fit$Method]
compare.fit <- cbind(n,p,L,s,O,index,s_index, compare.fit)
compare.fit








