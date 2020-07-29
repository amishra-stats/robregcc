##########################################################
### Modeling inflammation matrker: HIV data analysis  ####
### Mislabeling experiement                           ####
##########################################################
## -----------------------------------
## Data analysis with:  C = 1_p : no-subcomposition  ####
## Data analysis with:  With or without corrupted observations
## -----------------------------------


rm(list = ls())

# # Load data
data_url = "https://github.com/amishra-stats/robregcc/raw/master/sCD14wC.rds"
download.file(data_url, destfile='temp.rds')
load("temp.rds") 
# load("../HIV_data/sCD14wC.rda")         # Load data for estimation

## Install and load the package "robregcc" from CRAN/GitHub
library(devtools)
# devtools::install_github("amishra-stats/robregcc/robregcc")
# install.packages("robregcc")
library(robregcc)

## Load required packages
library(MASS)
require(ggplot2)
require(reshape2)
library(magrittr)
library(graphics)



C <- matrix(1,nrow = 1, ncol = 60)

## Add pseudo count to impute zeros  in the observation 
W <- sCD14[,1:60]     # species abundance data 
# relative abindance data 
X <- t( apply(W, 1, function(x) {
  indx <- x != 0; x[!indx]  <- 0.5; x/sum(x) } ))

# response matrix 
y <- sCD14[,61]


# generate corrupted observation
O = round(0.06*length(y)) + 1     # number of outlier, e.g. 15% of observation
tindex <- order(y)
tindex1 <- tindex[1:(O/2)]
tindex2 <- tindex[ length(y) + 1 - (1:(O/2)) ]

# corrupted y dennoted by y_c
y_c <- y                        
y_c[tindex1] <- y[tindex2] + 2*0
y_c[tindex2] <- y[tindex1] - 2*0
tindex <- c(tindex1,tindex2)


# Data setting 
exp_seed <- 123
set.seed(exp_seed)
Xt <- cbind(1,log(X))
Ct <- cbind(matrix(0,nrow = nrow(C),ncol = 1),C)
p <- ncol(Xt)
bw <- c(0,rep(1,p - 1))                       # weight matrix to not penalize intercept 


# ------------ Fitting non-robust model [Pixy shi 2016]
# corrupted observation: NO 
set.seed(exp_seed)
control <- robregcc_option(maxiter = 5000, tol = 1e-12, lminfac = 1e-12)
fit.nr <- classo(Xt, y, Ct, we = bw, type = 1, control = control) 
resid_nr <- y - Xt %*% fit.nr$beta


# corrupted observation: Yes
set.seed(exp_seed)
fit.nr_c <- classo(Xt, y_c, Ct, we = bw, type = 1, control = control) 
resid_nr_c <- y_c - Xt %*% fit.nr_c$beta




# ------------ intialization -- [PSC analysis for compositional data]
# control parameter 
set.seed(exp_seed)
b1 = 0.25; cc1 =  2.937   # initalization for scale parameter
control <- robregcc_option(maxiter = 5000,tol = 1e-4,lminfac = 1e-7)
# corrupted observation: NO 
fit.init <- cpsc_sp(Xt, y,alp = 0.4, cfac = 2, b1 = b1, 
                    cc1 = cc1,Ct,bw,1,control)  


# corrupted observation: Yes
set.seed(exp_seed)     
fit.init_c <- cpsc_sp(Xt, y_c, alp = 0.4, cfac = 2, b1 = b1,
                      cc1 = cc1,Ct,bw,1,control) 



## Model fitting
set.seed(exp_seed)

# control parameters
control <- robregcc_option()
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
# corrupted observation: NO 
set.seed(exp_seed)
fit.ada <- robregcc_sp(Xt,y,Ct,
                       beta.init = fit.init$betaR,  cindex = 1, 
                       gamma.init = fit.init$residuals,
                       control = control, 
                       penalty.index = 1, alpha = 0.95)

# corrupted observation: Yes
set.seed(exp_seed)
fit.ada_c <- robregcc_sp(Xt,y_c,Ct,
                         beta.init = fit.init_c$betaR,  cindex = 1, 
                         gamma.init = fit.init_c$residuals,
                         control = control, 
                         penalty.index = 1, alpha = 0.95)



set.seed(exp_seed)
control$lmaxfac = 1                 # Parameter for constructing sequence of lambda
control$lminfac = 1e-8              # Parameter for constructing sequence of lambda 
# Robust regression using lasso penalty [Huber equivalent]   
# corrupted observation: NO
fit.soft <- robregcc_sp(Xt,y,Ct, cindex = 1, 
                        control = control, penalty.index = 2, 
                        alpha = 0.95)

# corrupted observation: Yes
set.seed(exp_seed)
fit.soft_c <- robregcc_sp(Xt,y_c,Ct, cindex = 1, 
                          control = control, penalty.index = 2, 
                          alpha = 0.95)



# Robust regression using hard thresholding penalty
control$lmaxfac = 1e1               # Parameter for constructing sequence of lambda
control$lminfac = 1e-2              # Parameter for constructing sequence of lambda

# corrupted observation: NO
set.seed(exp_seed)
fit.hard <- robregcc_sp(Xt,y,Ct, beta.init = fit.init$betaf, 
                        gamma.init = fit.init$residuals,
                        cindex = 1, 
                        control = control, penalty.index = 3, 
                        alpha = 0.95)

# corrupted observation: Yes
set.seed(exp_seed)
fit.hard_c <- robregcc_sp(Xt,y_c,Ct, beta.init = fit.init_c$betaf, 
                          gamma.init = fit.init_c$residuals,
                          cindex = 1, 
                          control = control, penalty.index = 3, 
                          alpha = 0.95)



## ----------------- Corrupted observation: No
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
# save(predictor_est_c1, file = 'hiv_c1.rda')


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
# 0.27 0.63 0.53 0.57 


## ----------------- Corrupted observation: Yes
beta_lam_cv_c <- cbind(fit.nr_c$beta,fit.ada_c$betaE[,s_index],
                       fit.soft_c$betaE[,s_index],fit.hard_c$betaE[,s_index])
## Compute R-square of the robust and the non-robust model 
shift_est_lam_cv_c <- cbind(0,fit.ada_c$gamma0[,s_index],
                            fit.soft_c$gamma0[,s_index],fit.hard_c$gamma0[,s_index])
# compute weight and calculate weighted r square 
residual_est_lam_cv_c <- drop(y_c) - Xt %*% beta_lam_cv_c - shift_est_lam_cv_c
R_square_c <-  1 - apply(residual_est_lam_cv_c, 2, function(x) sum(x^2) )/sum((y - mean(y))^2) 
names(R_square_c) <- c('NR','A','E','H')
R_square_c <- signif(R_square_c,2)
R_square_c
# NR     A     E     H 
# 0.065 0.550 0.550 0.690



save(list = ls(), file = '../mislabel_experiment.rda')




# ---------------------------------------------------
# Model fit plots are provided in the supplementary material 
## Figure S5, S6, and S7 in the supplementary material 

# Corrupted observation: Yes
fdir <-  getwd()
pdf(file = file.path(fdir,'appl_hiv_c.pdf'), width = 15, height = 4)
## Plots for the supplementary material 
# [adaptive]
par(mfrow = c(1,3), mar = 5*c(1,1,0.5,0.1))
plot_path(fit.ada_c)
plot_cv(fit.ada_c)
plot_resid(fit.ada_c,1,0)


# [soft]
par(mfrow = c(1,3), mar = 5*c(1,1,0.5,0.1))
plot_path(fit.soft_c)
plot_cv(fit.soft_c)
plot_resid(fit.soft_c,1,0)


# [Hard]
par(mfrow = c(1,3), mar = 5*c(1,1,0.5,0.1))
plot_path(fit.hard_c)
plot_cv(fit.hard_c)
plot_resid(fit.hard_c,1,0)
par(mfrow = c(1,1))
dev.off()




## ----------------------------  Predicted value:
# Corrupted observation: No
Y.hat <- data.frame(Xt %*% beta_lam_cv,Index = 1:length(y),y = y)
nr_model <- lm(y~Y.hat[,1] -1)
obs.sel <- (!cbind(fit.ada$inlierE[,s_index],fit.soft$inlierE[,s_index],
                   fit.hard$inlierE[,s_index])) + 0
obs.sel2 <- cbind(0,obs.sel)     # observatio selected 0

# Corrupted observation: Yes
Y.hat_c <- data.frame(Xt %*% beta_lam_cv_c,Index = 1:length(y),y = y_c)
nr_model <- lm(y_c~Y.hat_c[,1] -1)
obs.sel <- (!cbind(fit.ada_c$inlierE[,s_index],fit.soft_c$inlierE[,s_index],
                   fit.hard_c$inlierE[,s_index])) + 0
obs.sel2 <- cbind(0,obs.sel)     # observatio selected 0



## Figure S8 in the supplementary material 
pdf(file.path(fdir,'hiv_selbal_fit_c.pdf'),height = 5,width = 6)
## Adaptive vs non_robust approach 
par(mfrow=c(1,1), mar = 50*c(0.1,0.1,0.03,0.03))
selIndex <- 2
kindx <- obs.sel2[,selIndex] == 0
ctype <- obs.sel2[,selIndex] + 1
ptype <- obs.sel2[,selIndex] + 19
ptype[ptype == 20] <- 17
plot(y_c,Y.hat_c[,selIndex],col=ctype,ylab=expression(widehat(y)),
     pch=ptype,cex=1,cex.axis = 1.5,cex.lab=2, xlab = "y", font  = 2, 
     mgp = c(2.2,1,0), log = 'xy', ylim = c(5000,8500), 
     yaxt = "n" , xaxt = "n" )
abline(lm(y_c[kindx]~Y.hat_c[kindx,selIndex]-1), lty = 1, col = 'red', lwd = 4)
abline(nr_model, col = 'blue', lty = 3, lwd = 4)
legend("topleft", inset=.01, legend=c('[A]', "[NR]"),
       col=c("red", "blue"), lty=c(1,3), cex=1.4, lwd = 4,
       text.font = 2, bg="transparent")
axis(2, at = seq(5500, 8500, by = 1000), font = 2, cex.axis = 1.5)
axis(1, at = c(2000, 5000, 10000, 20000), font = 2, cex.axis = 1.5)
text(2300,6000, bquote(R[A]^2 == .(R_square_c[2])), cex=1.9, lwd = 4, font = 2)
text(2300,5500, bquote(R[NR]^2 == .(R_square_c[1])), cex=1.9, lwd = 4, font = 2)



## Soft vs non_robust approach 
par(mfrow=c(1,1), mar = 50*c(0.1,0.1,0.03,0.03))
selIndex <- 3
kindx <- obs.sel2[,selIndex] == 0
ctype <- obs.sel2[,selIndex] + 1
ptype <- obs.sel2[,selIndex] + 19
ptype[ptype == 20] <- 17
plot(y_c,Y.hat_c[,selIndex],col=ctype,ylab=expression(widehat(y)),
     pch=ptype,cex=0.8,cex.axis = 1.5,cex.lab=2, xlab = "y", font  = 2, 
     mgp = c(2.2,1,0), log = 'xy', ylim = c(5000,8500), 
     yaxt = "n" , xaxt = "n"  )
abline(lm(y_c[kindx]~Y.hat_c[kindx,selIndex]-1), lty = 1, col = 'red', lwd = 4)
abline(nr_model, col = 'blue', lty = 3, lwd = 4)
legend("topleft", inset=.01, legend=c('[E]', "[NR]"),
       col=c("red", "blue"), lty=c(1,3), cex=1.4, lwd = 4,
       text.font = 2, bg="transparent" )
axis(2, at = seq(5500, 8500, by = 1000), font = 2, cex.axis = 1.5)
axis(1, at = c(2000, 5000, 10000, 20000), font = 2,  cex.axis = 1.5)
text(2300,6000, bquote(R[E]^2 == .(R_square_c[3])), cex=1.9, lwd = 4, font = 2)
text(2300,5500, bquote(R[NR]^2 == .(R_square_c[1])), cex=1.9, lwd = 4, font = 2)


## Hard vs non_robust approach 
par(mfrow=c(1,1), mar = 50*c(0.1,0.1,0.03,0.03))
selIndex <- 4
kindx <- obs.sel2[,selIndex] == 0
ctype <- obs.sel2[,selIndex] + 1
ptype <- obs.sel2[,selIndex] + 19
ptype[ptype == 20] <- 17
plot(y_c,Y.hat_c[,selIndex],col=ctype,ylab=expression(widehat(y)),
     pch=ptype,cex=0.8,cex.axis = 1.5,cex.lab=2, xlab = "y", font  = 2, 
     mgp = c(2.2,1,0), log = 'xy' , 
     ylim = c(5000,8500), yaxt = "n" , xaxt = "n" )
abline(lm(y_c[kindx]~Y.hat_c[kindx,selIndex]-1), lty = 1, col = 'red', lwd = 4)
abline(nr_model, col = 'blue', lty = 3, lwd = 4)
legend("topleft", inset=.01, legend=c('[H]', "[NR]"),
       col=c("red", "blue"), lty=c(1,3), cex=1.4, lwd = 4,
       text.font = 2, bg="transparent")
axis(2, at = seq(5500, 8500, by = 1000), font = 2, cex.axis = 1.5)
axis(1, at = c(2000, 5000, 10000, 20000), font = 2, cex.axis = 1.5)
text(2300,6000, bquote(R[H]^2 == .(R_square_c[4])), cex=1.9, lwd = 4, font = 2)
text(2300,5500, bquote(R[NR]^2 == .(R_square_c[1])), cex=1.9, lwd = 4, font = 2)

dev.off()







## -----------------------------------------------------------------------------
## Comparison of corrupted model vs true model 
## Figure 6 in the main manuscript
## -----------------------------------------------------------------------------

### ------------------------ Figure 6 (b)
## error in model parameter beta estimate 
del_beta_cv <- sqrt(colSums((beta_lam_cv)^2)/nrow(beta_lam_cv))
del_beta <- sqrt(colSums((beta_lam_cv - beta_lam_cv_c)^2)/nrow(beta_lam_cv))
del_beta_fraction <- del_beta/del_beta_cv
# fraction of suppport mismatch 
del_beta_support <- 1 - colSums((beta_lam_cv!=0) * (beta_lam_cv_c!=0))/colSums((beta_lam_cv!=0))

# Hamming distance 
del_beta_support <- abs(((beta_lam_cv!=0) +0) - ((beta_lam_cv_c!=0) +0))
del_beta_support <- colSums(del_beta_support)


# use bar chart to represent support mismatch index; 
dtf <- data.frame(M = c('NR','A','E','H'), del_beta = 100*del_beta_fraction, 
                  supp_beta = del_beta_support)

# plot module 
dtf$M <- factor(as.character(dtf$M),levels = c('NR',"A", "E","H"))
df_plot1 <- ggplot(dtf, aes(x = M, y = del_beta)) +  
  geom_bar(aes(fill = M), stat="identity", position="dodge") +
  scale_fill_manual(values = c(4,1,3,2)) +
  theme_bw() + xlab('') + ylab('%') + 
  theme(#strip.placement = "outside",
    legend.position =  "none", #"top",
    axis.text.x = element_text( hjust = 0.5,size = 25,colour = 'black',
                                face = 'bold',vjust = 0.4),
    axis.text.y = element_text(angle = 0, hjust = 1,size = 25,
                               colour = 'black',face = 'bold'),
    axis.title = element_text(size = 25, face = "bold") )  


df_plot2 <- ggplot(dtf, aes(x = M, y = supp_beta)) +  
  geom_bar(aes(fill = M), stat="identity", position="dodge") +
  scale_fill_manual(values = c(4,1,3,2)) +
  theme_bw() + xlab('') + ylab('') + 
  theme(legend.position =  "none", #"top",
        axis.text.x = element_text( hjust = 0.5,size = 25,colour = 'black',
                                    face = 'bold',vjust = 0.4),
        axis.text.y = element_text(angle = 0, hjust = 1,size = 25,
                                   colour = 'black',face = 'bold'),
        axis.title = element_text(size = 25, face = "bold") )  


fdir <-  getwd()
pdf(file.path(fdir,'appl_hiv_delbeta.pdf'),onefile = TRUE,width = 4, height = 6)
print(df_plot1)
print(df_plot2)
dev.off()




### ------------------------ Figure 6 (c)
## Present outliers itentified by the model 
fdir <-  getwd()
pdf(file.path(fdir,'appl_hiv_corrupted_outlier.pdf'),
    onefile = TRUE,width = 8, height = 6)
n <- length(y)
obs <- data.frame(Index = 1:n, value = fit.ada_c$residualsE[,s_index], 
                  Outlier = fit.ada_c$inlierE[,s_index]+1)
obs$Outlier <- c('Y','N')[obs$Outlier]
ggplot(obs, aes(x = Index, y = value)) +  
  geom_point(aes(color = Outlier), size = 5) + 
  scale_color_manual(values=c(1, 2)) +
  theme_bw() + xlab('') + ylab(expression(bold( y - hat(y)) )) +
  theme(legend.position = "none",
        legend.box = "horizontal",
        axis.text=element_text(size=25),
        axis.title = element_text(color="black", size=30, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 25,colour = 'black',
                                   face = 'bold',vjust = 0.4),
        axis.text.y = element_text(angle = 0, hjust = 1,size = 25,
                                   colour = 'black',face = 'bold'),
        legend.title = element_text(colour = "Black", size = 25, face = "bold"),
        legend.text = element_text(size = 25, face = "bold"),
        legend.background = element_rect( linetype="solid", 
                                          colour ="black"))  
dev.off()



### ------------------------ Figure 6 (a)
# Plot corrupted observation
n <- length(y)
obs <- data.frame(Index = 1:n, y_hat = Y.hat[,2], y = y, Mislabel = 2)
# plot(obs$T, obs$Index)
obs$Mislabel[tindex1] <- 2+(1:length(tindex1))
obs$Mislabel[tindex2] <- 2+(1:length(tindex1))
obs$High <- obs$Mislabel
obs$Mislabel <- as.factor(obs$Mislabel)
library(reshape2)
library(ggplot2)
library(gghighlight)
p <- ggplot(obs, aes(x =  y , y = y_hat, shape = Mislabel,size = Mislabel)) +  # 
  geom_point(aes(color = Mislabel)) + 
  scale_shape_manual(values=c(16,8,10,15,17,18)) +
  scale_color_manual(values=c(1:6)) +
  scale_size_manual(values=c(4,6*rep(1,5)))+
  theme_bw() + xlab('y') + ylab(expression(bold(hat(y)))) +
  ggtitle('RobRegCC with adaptive penalty') + 
  geom_abline(intercept = 0, slope = 1.051, linetype = "dashed",size = 2,col='red') +
  theme(legend.position = 'None', 
        plot.title = element_text(hjust = 0.5,color="black", size=25, face="bold"),
        axis.text=element_text(size=25),
        axis.title = element_text(color="black", size=30, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 25,colour = 'black',
                                   face = 'bold',vjust = 0.4),
        axis.text.y = element_text(angle = 0, hjust = 1,size = 25,
                                   colour = 'black',face = 'bold')
  )   
p
fdir <-  getwd()
pdf(file.path(fdir,'appl_hiv_corupted.pdf'),onefile = TRUE,width = 8, height = 6)
print(p)
dev.off()

