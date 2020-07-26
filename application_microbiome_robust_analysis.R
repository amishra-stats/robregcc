###########################################################
## R code for  microbiome applications                   ##
## in "Robust regression with Compositional Covariates"  ##
###########################################################



#####################################################################

##########################################################
### Modeling inflammation matrker: HIV data analysis  ####
##########################################################


rm(list = ls())

# Load data
load("../HIV_data/sCD14wC.rda")         # Load data for estimation

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


## -----------------------------------
## Data analysis with:  C = 1_p : no-subcomposition  ####
## -----------------------------------
C <- matrix(1,nrow = 1, ncol = 60)

## Add pseudo count to impute zeros  in the observation 
W <- sCD14[,1:60]     # species abundance data 
# relative abindance data 
X <- t( apply(W, 1, function(x) {
  indx <- x != 0; x[!indx]  <- 0.5; x/sum(x) } ))

# response matrix 
y <- sCD14[,61]



# Data setting 
exp_seed <- 123
set.seed(exp_seed)
Xt <- cbind(1,log(X))
Ct <- cbind(matrix(0,nrow = nrow(C),ncol = 1),C)
p <- ncol(Xt)
bw <- c(0,rep(1,p - 1))                       # weight matrix to not penalize intercept 


# ------------ Fitting non-robust model [Pixy shi 2016]
control <- robregcc_option(maxiter = 5000, tol = 1e-12, lminfac = 1e-12)
fit.nr <- classo(Xt, y, Ct, we = bw, type = 1, control = control) 
resid_nr <- y - Xt %*% fit.nr$beta


# ------------ intialization -- [PSC analysis for compositional data]
# control parameter 
set.seed(exp_seed)
b1 = 0.25; cc1 =  2.937   # initalization for scale parameter
control <- robregcc_option(maxiter = 5000,tol = 1e-4,lminfac = 1e-7)
fit.init <- cpsc_sp(Xt, y,alp = 0.4, cfac = 2, b1 = b1, cc1 = cc1,Ct,bw,1,control)  




## Model fitting
set.seed(exp_seed)

# control parameters
control <- robregcc_option()
beta.wt <- fit.init$betaR          # Set weight for model parameter beta
beta.wt[1] <- 0
control$gamma = 1                   # gamma for constructing  weighted penalty
control$spb = 50/p                  # fraction of maximum non-zero model parameter beta
control$outMiter = 1000             # Outer loop iteration
control$inMiter = 3000              # Inner loop iteration
control$nlam = 100                   # Number of tuning parameter lambda to be explored
control$lmaxfac = 1                 # Parameter for constructing sequence of lambda
control$lminfac = 1e-8              # Parameter for constructing sequence of lambda 
control$tol = 1e-20;                # tolrence parameter for converging [inner  loop]
control$out.tol = 1e-16             # tolerence parameter for convergence [outer loop]
control$kfold = 10                   # number of fold of crossvalidation ----
control$sigmafac = 2

# Robust regression using adaptive lasso penalty
fit.ada <- robregcc_sp(Xt,y,Ct,
                       beta.init = beta.wt,  cindex = 1, 
                       gamma.init = fit.init$residuals,
                       control = control, 
                       penalty.index = 1, alpha = 0.95)




set.seed(exp_seed)
control$lmaxfac = 1                 # Parameter for constructing sequence of lambda
control$lminfac = 1e-8              # Parameter for constructing sequence of lambda 
# Robust regression using lasso penalty [Huber equivalent]   
fit.soft <- robregcc_sp(Xt,y,Ct, cindex = 1, 
                        control = control, penalty.index = 2, 
                        alpha = 0.95)




set.seed(exp_seed)
# Robust regression using hard thresholding penalty
control$lmaxfac = 1e1               # Parameter for constructing sequence of lambda
control$lminfac = 1e-2              # Parameter for constructing sequence of lambda
control$sigmafac = 2#1.345
fit.hard <- robregcc_sp(Xt,y,Ct, beta.init = fit.init$betaf, 
                        gamma.init = fit.init$residuals,
                        cindex = 1, 
                        control = control, penalty.index = 3, 
                        alpha = 0.95)


## save robust predictor estimate here as predictor_est_c1
s_index <- 1  # 1/2 corresponding to minimum/1se rule cross validation error 
beta_lam_cv <- cbind(fit.nr$beta,fit.ada$betaE[,s_index],
                     fit.soft$betaE[,s_index],fit.hard$betaE[,s_index])
predictor_name <- c('intercept',colnames(X))
ord_predictor <- order(predictor_name)
# save estimate of the predictor with C = 1_p
predictor_est_c1 <- data.frame(beta_lam_cv[ord_predictor,],
                               v.name = predictor_name[ord_predictor])


## we wil use the saved output to plot later 
save(predictor_est_c1, file = '../model_fit_plot1.rda')


## Compute R-square of the robust and the non-robust model 
shift_est_lam_cv <- cbind(0,fit.ada$gamma0[,s_index],
                          fit.soft$gamma0[,s_index],fit.hard$gamma0[,s_index])
# compute weight and calculate weighted r square 
residual_est_lam_cv <- drop(y) - Xt %*% beta_lam_cv - shift_est_lam_cv
R_square <-  1 - apply(residual_est_lam_cv, 2, function(x) sum(x^2) )/sum((y - mean(y))^2) 
names(R_square) <- c('NR','A','E','H')
R_square

# NR         A         E         H 
# 0.2742904 0.5714158 0.5330283 0.7169454 






# ---------------------------------------------
# HIV data analysis plots with C = 1_p
# Model fit plots are provided in the supplementary material 
## Plots for the supplementary material 
## Figure S1, S2, and S3 in the supplementary material 
## 
fdir <-  getwd()
pdf(file = file.path(fdir,'appl_hiv1.pdf'), width = 15, height = 4)

# [adaptive]
par(mfrow = c(1,3), mar = 5*c(1,1,0.5,0.1))
plot_path(fit.ada)
plot_cv(fit.ada)
plot_resid(fit.ada,1,0)


# [soft]
par(mfrow = c(1,3), mar = 5*c(1,1,0.5,0.1))
plot_path(fit.soft)
plot_cv(fit.soft)
plot_resid(fit.soft,1,0)


# [Hard]
par(mfrow = c(1,3), mar = 5*c(1,1,0.5,0.1))
plot_path(fit.hard)
plot_cv(fit.hard)
plot_resid(fit.hard,1,0)
par(mfrow = c(1,1))
dev.off()



## Predicted value:
Y.hat <- data.frame(Xt %*% beta_lam_cv,Index = 1:length(y),y = y)
nr_model <- lm(y~Y.hat[,1] -1)
obs.sel <- (!cbind(fit.ada$inlierE[,s_index],fit.soft$inlierE[,s_index],
                   fit.hard$inlierE[,s_index])) + 0
obs.sel2 <- cbind(0,obs.sel)     # observatio selected 0




## Figure S4 in the supplementary material 
pdf(file.path(fdir,'hiv_selbal_fit.pdf'),height = 5,width = 6)
## Adaptive vs non_robust approach 
par(mfrow=c(1,1), mar = 50*c(0.1,0.1,0.03,0.03))
selIndex <- 2
kindx <- obs.sel2[,selIndex] == 0
ctype <- obs.sel2[,selIndex] + 1
ptype <- obs.sel2[,selIndex] + 19
ptype[ptype == 20] <- 17
plot(y,Y.hat[,selIndex],col=ctype,ylab=expression(widehat(y)),
     pch=ptype,cex=1,cex.axis = 1.5,cex.lab=2, xlab = "y", font  = 2, 
     mgp = c(2.2,1,0), log = 'xy', ylim = c(5000,8500), 
     yaxt = "n" , xaxt = "n" )
abline(lm(y[kindx]~Y.hat[kindx,selIndex]-1), lty = 1, col = 'red', lwd = 4)
abline(nr_model, col = 'blue', lty = 3, lwd = 4)
legend("topleft", inset=.01, legend=c('[A]', "[NR]"),
       col=c("red", "blue"), lty=c(1,3), cex=1.4, lwd = 4,
       text.font = 2, bg="transparent")
axis(2, at = seq(5500, 8500, by = 1000), font = 2, cex.axis = 1.5)
axis(1, at = c(2000, 5000, 10000, 20000), font = 2, cex.axis = 1.5)
text(2300,6000, expression(R[A]^2*' = 0.57'), cex=1.9, lwd = 4, font = 2)
text(2300,5500, expression(R[NR]^2*' = 0.27'), cex=1.9, lwd = 4, font = 2)


## Soft vs non_robust approach 
par(mfrow=c(1,1), mar = 50*c(0.1,0.1,0.03,0.03))
selIndex <- 3
kindx <- obs.sel2[,selIndex] == 0
ctype <- obs.sel2[,selIndex] + 1
ptype <- obs.sel2[,selIndex] + 19
ptype[ptype == 20] <- 17
plot(y,Y.hat[,selIndex],col=ctype,ylab=expression(widehat(y)),
     pch=ptype,cex=0.8,cex.axis = 1.5,cex.lab=2, xlab = "y", font  = 2, 
     mgp = c(2.2,1,0), log = 'xy', ylim = c(5000,8500), 
     yaxt = "n" , xaxt = "n"  )
abline(lm(y[kindx]~Y.hat[kindx,selIndex]-1), lty = 1, col = 'red', lwd = 4)
abline(nr_model, col = 'blue', lty = 3, lwd = 4)
legend("topleft", inset=.01, legend=c('[E]', "[NR]"),
       col=c("red", "blue"), lty=c(1,3), cex=1.4, lwd = 4,
       text.font = 2, bg="transparent" )
axis(2, at = seq(5500, 8500, by = 1000), font = 2, cex.axis = 1.5)
axis(1, at = c(2000, 5000, 10000, 20000), font = 2,  cex.axis = 1.5)
text(2300,6000, expression(R[E]^2*' = 0.53'), cex=1.9, lwd = 4, font = 2)
text(2300,5500, expression(R[NR]^2*' = 0.27'), cex=1.9, lwd = 4, font = 2)


## Hard vs non_robust approach 
par(mfrow=c(1,1), mar = 50*c(0.1,0.1,0.03,0.03))
selIndex <- 4
kindx <- obs.sel2[,selIndex] == 0
ctype <- obs.sel2[,selIndex] + 1
ptype <- obs.sel2[,selIndex] + 19
ptype[ptype == 20] <- 17
plot(y,Y.hat[,selIndex],col=ctype,ylab=expression(widehat(y)),
     pch=ptype,cex=0.8,cex.axis = 1.5,cex.lab=2, xlab = "y", font  = 2, 
     mgp = c(2.2,1,0), log = 'xy' , 
     ylim = c(5000,8500), yaxt = "n" , xaxt = "n" )
abline(lm(y[kindx]~Y.hat[kindx,selIndex]-1), lty = 1, col = 'red', lwd = 4)
abline(nr_model, col = 'blue', lty = 3, lwd = 4)
legend("topleft", inset=.01, legend=c('[H]', "[NR]"),
       col=c("red", "blue"), lty=c(1,3), cex=1.4, lwd = 4,
       text.font = 2, bg="transparent")
axis(2, at = seq(5500, 8500, by = 1000), font = 2, cex.axis = 1.5)
axis(1, at = c(2000, 5000, 10000, 20000), font = 2, cex.axis = 1.5)
text(2300,6000, expression(R[H]^2*' = 0.71'), cex=1.9, lwd = 4, font = 2)
text(2300,5500, expression(R[NR]^2*' = 0.27'), cex=1.9, lwd = 4, font = 2)

dev.off()


















##########################################################
### Modeling inflammation matrker: HIV data analysis  ####
### Analysis with:  C =  Subcomposition coherence     ####
##########################################################
## -----------------------------------
## Data analysis with:  C =  Subcomposition coherence  

rm(list = ls())
# Load data
load("../HIV_data/sCD14wC.rda")         # Load data for estimation


## Add pseudo count to impute zeros  in the observation 
W <- sCD14[,1:60]     # species abundance data 
# relative abindance data 
X <- t( apply(W, 1, function(x) {
  indx <- x != 0; x[!indx]  <- 0.5; x/sum(x) } ))

# response matrix 
y <- sCD14[,61]



# Data setting 
exp_seed <- 1234
set.seed(exp_seed)
Xt <- cbind(1,log(X))
Ct <- cbind(matrix(0,nrow = nrow(C),ncol = 1),C)
p <- ncol(Xt)
bw <- c(0,rep(1,p - 1))                       # weight matrix to not penalize intercept 


# ------------ Fitting non-robust model [Pixy shi 2016]
control <- robregcc_option(maxiter = 5000, tol = 1e-12, lminfac = 1e-12)
fit.nr <- classo(Xt, y, Ct, we = bw, type = 1, control = control) 
resid_nr <- y - Xt %*% fit.nr$beta


# ------------ intialization -- [PSC analysis for compositional data]
# control parameter 
set.seed(exp_seed)
b1 = 0.25; cc1 =  2.937   # initalization for scale parameter
control <- robregcc_option(maxiter = 5000,tol = 1e-4,lminfac = 1e-7)
fit.init <- cpsc_sp(Xt, y,alp = 0.4, cfac = 2, b1 = b1, cc1 = cc1,Ct,bw,1,control)  




## Model fitting
set.seed(exp_seed)

# control parameters
control <- robregcc_option()
beta.wt <- fit.init$betaR          # Set weight for model parameter beta
beta.wt[1] <- 0
control$gamma = 1                   # gamma for constructing  weighted penalty
control$spb = 50/p                  # fraction of maximum non-zero model parameter beta
control$outMiter = 1000             # Outer loop iteration
control$inMiter = 3000              # Inner loop iteration
control$nlam = 100                   # Number of tuning parameter lambda to be explored
control$lmaxfac = 1                 # Parameter for constructing sequence of lambda
control$lminfac = 1e-8              # Parameter for constructing sequence of lambda 
control$tol = 1e-20;                # tolrence parameter for converging [inner  loop]
control$out.tol = 1e-16             # tolerence parameter for convergence [outer loop]
control$kfold = 10                   # number of fold of crossvalidation ----
control$sigmafac = 2

# Robust regression using adaptive lasso penalty
fit.ada <- robregcc_sp(Xt,y,Ct,
                       beta.init = beta.wt,  cindex = 1, 
                       gamma.init = fit.init$residuals,
                       control = control, 
                       penalty.index = 1, alpha = 0.95)




set.seed(exp_seed)
control$lmaxfac = 1                 # Parameter for constructing sequence of lambda
control$lminfac = 1e-8              # Parameter for constructing sequence of lambda 
# Robust regression using lasso penalty [Huber equivalent]   
fit.soft <- robregcc_sp(Xt,y,Ct, cindex = 1, 
                        control = control, penalty.index = 2, 
                        alpha = 0.95)




set.seed(exp_seed)
# Robust regression using hard thresholding penalty
control$lmaxfac = 1e1               # Parameter for constructing sequence of lambda
control$lminfac = 1e-2              # Parameter for constructing sequence of lambda
control$sigmafac = 2#1.345
fit.hard <- robregcc_sp(Xt,y,Ct, beta.init = fit.init$betaf, 
                        gamma.init = fit.init$residuals,
                        cindex = 1, 
                        control = control, penalty.index = 3, 
                        alpha = 0.95)


## save robust predictor estimate here as predictor_est_c1
s_index <- 1  # 1/2 corresponding to minimum/1se rule cross validation error 
beta_lam_cv <- cbind(fit.nr$beta,fit.ada$betaE[,s_index],
                     fit.soft$betaE[,s_index],fit.hard$betaE[,s_index])
predictor_name <- c('intercept',colnames(X))
ord_predictor <- order(predictor_name)
#
subcomp_name <- vector('character',length(predictor_name))
for (i in 1:nrow(Ct)) 
  subcomp_name[which(Ct[i,] == 1)] <- rownames(Ct)[i]
# save estimate of the predictor with C = subcompositional coherence 
predictor_est_ck <- data.frame(beta_lam_cv[ord_predictor,],
                               v.name = predictor_name[ord_predictor],
                               ordInfo = subcomp_name[ord_predictor])

## we wil use the saved output to plot later 
save(predictor_est_ck, file = '../model_fit_plotK.rda')



## Compute R-square of the robust and the non-robust model 
shift_est_lam_cv <- cbind(0,fit.ada$gamma0[,s_index],
                          fit.soft$gamma0[,s_index],fit.hard$gamma0[,s_index])
# compute weight and calculate weighted r square 
residual_est_lam_cv <- drop(y) - Xt %*% beta_lam_cv - shift_est_lam_cv
R_square <-  1 - apply(residual_est_lam_cv, 2, function(x) sum(x^2) )/sum((y - mean(y))^2) 
names(R_square) <- c('NR','A','E','H')
R_square <- signif(R_square,2)
R_square
# NR         A         E         H 
# 0.30 0.56 0.41 0.62 














# ---------------------------------------------
# HIV data analysis plots with C = Subcomposition coherence 
# Model fit plots are provided in the supplementary material 
## Plots for the supplementary material 
## Figure S9, S10, and S11 in the supplementary material 
## 



fdir <-  getwd()
pdf(file = file.path(fdir,'appl_hivk.pdf'), width = 15, height = 4)
## Plots for the supplementary material 
# [adaptive]
par(mfrow = c(1,3), mar = 5*c(1,1,0.5,0.1))
plot_path(fit.ada)
plot_cv(fit.ada)
plot_resid(fit.ada,1,0)


# [soft]
par(mfrow = c(1,3), mar = 5*c(1,1,0.5,0.1))
plot_path(fit.soft)
plot_cv(fit.soft)
plot_resid(fit.soft,1,0)


# [Hard]
par(mfrow = c(1,3), mar = 5*c(1,1,0.5,0.1))
plot_path(fit.hard)
plot_cv(fit.hard)
plot_resid(fit.hard,1,0)
par(mfrow = c(1,1))
dev.off()



## Predicted value:
Y.hat <- data.frame(Xt %*% beta_lam_cv,Index = 1:length(y),y = y)
nr_model <- lm(y~Y.hat[,1] -1)
obs.sel <- (!cbind(fit.ada$inlierE[,s_index],fit.soft$inlierE[,s_index],
                   fit.hard$inlierE[,s_index])) + 0
obs.sel2 <- cbind(0,obs.sel)     # observatio selected 0



## Figure S12 in the supplementary material 
pdf(file.path(fdir,'hiv_selbal_fit1.pdf'),height = 5,width = 6)
## Adaptive vs non_robust approach 
par(mfrow=c(1,1), mar = 50*c(0.1,0.1,0.03,0.03))
selIndex <- 2
kindx <- obs.sel2[,selIndex] == 0
ctype <- obs.sel2[,selIndex] + 1
ptype <- obs.sel2[,selIndex] + 19
ptype[ptype == 20] <- 17
plot(y,Y.hat[,selIndex],col=ctype,ylab=expression(widehat(y)),
     pch=ptype,cex=1,cex.axis = 1.5,cex.lab=2, xlab = "y", font  = 2, 
     mgp = c(2.2,1,0), log = 'xy', ylim = c(5000,8500), 
     yaxt = "n" , xaxt = "n" )
abline(lm(y[kindx]~Y.hat[kindx,selIndex]-1), lty = 1, col = 'red', lwd = 4)
abline(nr_model, col = 'blue', lty = 3, lwd = 4)
legend("topleft", inset=.01, legend=c('[A]', "[NR]"),
       col=c("red", "blue"), lty=c(1,3), cex=1.4, lwd = 4,
       text.font = 2, bg="transparent")
axis(2, at = seq(5500, 8500, by = 1000), font = 2, cex.axis = 1.5)
axis(1, at = c(2000, 5000, 10000, 20000), font = 2, cex.axis = 1.5)
text(2300,6000, bquote(R[A]^2 == .(R_square[2])), cex=1.9, lwd = 4, font = 2)
text(2300,5500, bquote(R[NR]^2 == .(R_square[1])), cex=1.9, lwd = 4, font = 2)



## Soft vs non_robust approach 
par(mfrow=c(1,1), mar = 50*c(0.1,0.1,0.03,0.03))
selIndex <- 3
kindx <- obs.sel2[,selIndex] == 0
ctype <- obs.sel2[,selIndex] + 1
ptype <- obs.sel2[,selIndex] + 19
ptype[ptype == 20] <- 17
plot(y,Y.hat[,selIndex],col=ctype,ylab=expression(widehat(y)),
     pch=ptype,cex=0.8,cex.axis = 1.5,cex.lab=2, xlab = "y", font  = 2, 
     mgp = c(2.2,1,0), log = 'xy', ylim = c(5000,8500), 
     yaxt = "n" , xaxt = "n"  )
abline(lm(y[kindx]~Y.hat[kindx,selIndex]-1), lty = 1, col = 'red', lwd = 4)
abline(nr_model, col = 'blue', lty = 3, lwd = 4)
legend("topleft", inset=.01, legend=c('[E]', "[NR]"),
       col=c("red", "blue"), lty=c(1,3), cex=1.4, lwd = 4,
       text.font = 2, bg="transparent" )
axis(2, at = seq(5500, 8500, by = 1000), font = 2, cex.axis = 1.5)
axis(1, at = c(2000, 5000, 10000, 20000), font = 2,  cex.axis = 1.5)
text(2300,6000, bquote(R[E]^2 == .(R_square[3])), cex=1.9, lwd = 4, font = 2)
text(2300,5500, bquote(R[NR]^2 == .(R_square[1])), cex=1.9, lwd = 4, font = 2)


## Hard vs non_robust approach 
par(mfrow=c(1,1), mar = 50*c(0.1,0.1,0.03,0.03))
selIndex <- 4
kindx <- obs.sel2[,selIndex] == 0
ctype <- obs.sel2[,selIndex] + 1
ptype <- obs.sel2[,selIndex] + 19
ptype[ptype == 20] <- 17
plot(y,Y.hat[,selIndex],col=ctype,ylab=expression(widehat(y)),
     pch=ptype,cex=0.8,cex.axis = 1.5,cex.lab=2, xlab = "y", font  = 2, 
     mgp = c(2.2,1,0), log = 'xy' , 
     ylim = c(5000,8500), yaxt = "n" , xaxt = "n" )
abline(lm(y[kindx]~Y.hat[kindx,selIndex]-1), lty = 1, col = 'red', lwd = 4)
abline(nr_model, col = 'blue', lty = 3, lwd = 4)
legend("topleft", inset=.01, legend=c('[H]', "[NR]"),
       col=c("red", "blue"), lty=c(1,3), cex=1.4, lwd = 4,
       text.font = 2, bg="transparent")
axis(2, at = seq(5500, 8500, by = 1000), font = 2, cex.axis = 1.5)
axis(1, at = c(2000, 5000, 10000, 20000), font = 2, cex.axis = 1.5)
text(2300,6000, bquote(R[H]^2 == .(R_square[4])), cex=1.9, lwd = 4, font = 2)
text(2300,5500, bquote(R[NR]^2 == .(R_square[1])), cex=1.9, lwd = 4, font = 2)


dev.off()







## -----------------------------------------------------
## Figure 4 and 5 of the main manuscript 
## 

rm(list = ls())
load(file = '../model_fit_plot1.rda')
tem_c1 <- predictor_est_c1
load(file = '../model_fit_plotK.rda')
tem_ck <- predictor_est_ck
tem_c1[,1:4] <- tem_c1[,1:4]*(abs(tem_c1[,1:4]) > 1e-5)
tem_ck[,1:4] <- tem_ck[,1:4]*(abs(tem_ck[,1:4]) > 1e-5)

tem_c1_sel <- tem_c1$v.name[rowSums(tem_c1[,1:4] != 0) > 0]
tem_ck_sel <- tem_ck$v.name[rowSums(tem_ck[,1:4] != 0) > 0]

selvar <- union(as.character(tem_c1_sel), as.character(tem_ck_sel))

dt_c1 <- tem_c1[match(selvar, tem_c1$v.name),]
dt_ck <- tem_ck[match(selvar, tem_ck$v.name),]
ordtem <- order(dt_ck$ordInfo)
dt_c1 <- dt_c1[ordtem, ]
dt_ck <- dt_ck[ordtem, ]
dt_c1$ordInfo <- dt_ck$ordInfo




## plots witth C = 1_p;  Figure 5 in the manuscript 
dtf <- dt_c1
names(dtf)[1:4] <- c('NR','A','E','H')
dtf[,1:4][abs(dtf[,1:4]) < 1e-3] <- 0
dtf2 <- dtf
##  ----- New modification to the plot function 
dtf3 <- dtf2
# dtf3 <- dtf3[rowSums(dtf3[,1:4] != 0) > 0,]
dtf3 <- dtf3[dtf3$v.name != 'intercept',]
# identify number of distinct beta in the model
beta0index <- data.frame(colSums(dtf3[,1:4] != 0))



# modify order names
library(stringr)
tvf <- as.character(dtf3$v.name)
tvf <- unlist(lapply(strsplit(tvf,'_'),paste, collapse = ' '))
tvf <- gsub('sensu stricto 1','sensu_stricto_1',tvf)
tvf <- gsub('f ','[f]',tvf)
tvf <- gsub('g ','[g]',tvf)
tvf <- gsub('o ','[o]',tvf)
tvf <- gsub(' S','_S',tvf)
dtf3$v.name <- str_wrap(tvf,width = 20)

# plot module 
dtf3 <- melt(dtf3,5:6)
# ii <- 1:nrow(dtf3)
dtf3$Name <- dtf3$v.name
names(dtf3)[names(dtf3) == 'variable'] <- 'M'

dtf3$M <- factor(as.character(dtf3$M),levels = c('NR',"A", "E","H"))
levels(dtf3$M) <- paste(paste(levels(dtf3$M), beta0index[levels(dtf3$M),],
                              sep = ' ['),']',sep = '')

ggplot(dtf3, aes(x = v.name, y = value)) +  
  geom_bar(aes(fill = M, linetype = M), stat="identity", position="dodge") +
  scale_fill_manual(values = c(4,1,3,2)) +
  theme_bw() + xlab('') + ylab('') +
  facet_grid(~ ordInfo, space = "free_x", scales = "free_x") + # , switch = "x"
  theme(strip.placement = "outside",
        legend.position = "top",
        strip.text  = element_text(size = 20,color = "transparent", 
                                   face = 'bold.italic',vjust = 1),
        strip.background = element_rect(fill = NA,colour = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 20,colour = 'black',
                                   face = 'bold',vjust = 0.4),
        axis.text.y = element_text(angle = 0, hjust = 1,size = 20,
                                   colour = 'black',face = 'bold'),
        panel.spacing = unit(.01,"cm"),
        legend.title = element_text(colour = "Black", size = 16, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20, face = "bold") )  

# Bacteroidales, Clostridiales, Erysipelotrichales, Selenomonadales, Uncategorized







## plots Figure 5
dtf <- dt_ck
names(dtf)[1:4] <- c('NR','A','E','H')
dtf[,1:4][abs(dtf[,1:4]) < 1e-3] <- 0
dtf2 <- dtf
##  ----- New modification to the plot function 
dtf3 <- dtf2
# dtf3 <- dtf3[rowSums(dtf3[,1:4] != 0) > 0,]
dtf3 <- dtf3[dtf3$v.name != 'intercept',]
# identify number of distinct beta in the model
beta0index <- data.frame(colSums(dtf3[,1:4] != 0))



# modify order names
library(stringr)
tvf <- as.character(dtf3$v.name)
tvf <- unlist(lapply(strsplit(tvf,'_'),paste, collapse = ' '))
tvf <- gsub('f ','[f]',tvf)
tvf <- gsub('g ','[g]',tvf)
tvf <- gsub('o ','[o]',tvf)
tvf <- gsub(' S','_S',tvf)
tvf <- gsub('sensu stricto 1','sensu_stricto_1',tvf)
dtf3$v.name <- str_wrap(tvf,width = 20)

# plot module 
dtf3 <- melt(dtf3,5:6)
# ii <- 1:nrow(dtf3)
dtf3$Name <- dtf3$v.name
names(dtf3)[names(dtf3) == 'variable'] <- 'M'

dtf3$M <- factor(as.character(dtf3$M),levels = c('NR',"A", "E","H"))
levels(dtf3$M) <- paste(paste(levels(dtf3$M), beta0index[levels(dtf3$M),],
                              sep = ' ['),']',sep = '')

ggplot(dtf3, aes(x = v.name, y = value)) +  
  geom_bar(aes(fill = M, linetype = M), stat="identity", position="dodge") +
  scale_fill_manual(values = c(4,1,3,2)) +
  theme_bw() + xlab('') + ylab('') +
  facet_grid(~ ordInfo, space = "free_x", scales = "free_x") + # , switch = "x"
  theme(strip.placement = "outside",
        legend.position = "top",
        strip.text  = element_text(size = 20,color = "transparent", 
                                   face = 'bold.italic',vjust = 1),
        strip.background = element_rect(fill = NA,colour = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 20,colour = 'black',
                                   face = 'bold',vjust = 0.4),
        axis.text.y = element_text(angle = 0, hjust = 1,size = 20,
                                   colour = 'black',face = 'bold'),
        panel.spacing = unit(.01,"cm"),
        legend.title = element_text(colour = "Black", size = 16, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20, face = "bold")  )  





