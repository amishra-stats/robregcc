# Robust regression with compositional covariates

R code for simulation studies and gut microbiome data analysis in the paper ["Robust regression with compositional covariates" by Aditya Mishra and Christian L. M&uuml;ller](https://arxiv.org/abs/1909.04990).


## Manuscript simulation studies and microbiome application 
#### Simulation scripts
**simulated_example_robust_analysis.R** - generates all the simulation settings (each setting is replicated 100 times) and computes the robust and non-robust model parameters estimate for a single replication. In addition, it also generates Figure 2 and model performance statistics for Table 2 and Figure 3.

**application_microbiome_robust_analysis.R** - computes the robust and non-robust model parameters estimate in the case of C = 1_p and C = sub-compositional coherence. The manuscript generates Figures 4 and 5 of the main manuscript.


**mislabel_experiement_robust_analysis.R** - in the mislabel experiment, the script computes the robust and non-robust model parameters estimate in the case of C = 1_p. The manuscript generates Figure 6 of the main manuscript.



#### Supporting data and precomputed results are in the folder **Data**:

**model_fit_plot1.rda** - non-robust and robust model parameters estimate in case of C = 1_p. The output is used to generate Figure 4 and 5 in the manuscript.

**model_fit_plotK.rda** - non-robust and robust parameters estimate in case of a model with sub-compositional coherence. The output is used to generate Figure 4 and 5 of the manuscript.  

**mislabel_experiment.rda** - n the mislabel experiment, the script computes the robust and non-robust model parameters estimate in the case of C = 1_p. The output is used to generate Figure 6 of the manuscript.  






## Getting Started

We have implemented the algorithm performing the robust regression with compositional covariates in the R package **robregcc**. Currently, GitHub host a development version of the package. Here we demonstrate the usage of functions available for the robust model fitting and outlier detection.

### Package installation

##### Dependency

```
## For installation
install.packages("MASS")
install.packages("magrittr")

install.packages("Rcpp")
install.packages("RcppArmadillo")


## Package for visualization:
install.packages("ggplot2")
install.packages("reshape2")
```


#### Install

Please download the package from the host. Mac/Linux user should download source file "robregcc_1.0.tar.gz".  Windows user should download binary file "robregcc_1.0.tar.gz". 
In the instruction given below, specify "download_location" accordingly.

```
## Windows user:  Binary format of the package
install.packages("download_location/robregcc_1.0.tgz", repos = NULL, type = "source")

## Mac/linux user: Use source file to install the package
install.packages("download_location/robregcc_1.0.tar.gz", repos = NULL, type = "source")
```


## Simulation examples:

Here we demonstrate the effectiveness of the proposed procedure on the simulated data. Please execute the documented code below to understand the required robust analysis. 

See the output of the example code in "testRobRegCC.pdf".


```
## -----------------------------------------------------------------------------

## Load package
rm(list = ls())
library(robregcc)
library(magrittr)

#--------------------------------------------------------------
## Define parameters to simulate example
## 
p <- 80                             # number of predictors  
n <- 300                            # number of sample   
O <- 0.10*n                         # number of outlier, e.g. 15% of observation    
L <- 1                              # indicator variable for outlier type, 
                                    # L = {0,1} => leveraged {no, yes}
                                    
# generate outlier by shifting "O"observation by amount equals to shFac times 
# true error variance sigma.
# s = {6,8} corresponds to {moderate, high} outlier 
s <- 8                          
ngrp <- 4                           # number of sub-composition
snr <- 3                            # Signal to noise ratio
example_seed <- 2*p+1               # example seed
set.seed(example_seed) 

#--------------------------------------------------------------
## Simulate true model variables, i.e., y, X, C, beta
## Follow examples from [Pixu Shi 2016]

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
#
#
# Account for intercept in the model
Xt <- cbind(1,X)                          # accounting for intercept in predictor
C <- cbind(0,C)                           # accounting for intercept in constraint
bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept 



#--------------------------------------------------------------
## Initialization

# Breakdown point for tukey Bisquare loss function 
b1 = 0.5                    # 50% breakdown point
cc1 =  1.567                # corresponding model parameter
b1 = 0.25; cc1 =  2.937   # initalization for scale parameter

set.seed(example_seed)      # unique seed

# control parameter for intialization method
control <- robregcc_option(maxiter = 1000,tol = 1e-4,lminfac = 1e-7)

# intialization
fit.init <- cpsc_sp(Xt, y,alp = 0.4, cfac = 2, b1 = b1,
                    cc1 = cc1,C,bw,1,control)  


#--------------------------------------------------------------
## Model fitting

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



#--------------------------------------------------------------
## Extract fitted model parameters

# coefficient estimate: [adaptive] 
coef_cc(fit.ada, type = 0, s = 1)

# coefficient estimate: [lasso/Huber] 
coef_cc(fit.soft, type = 0, s = 1)

# coefficient estimate: [Hard] 
coef_cc(fit.hard, type = 0, s = 1)



# residual estimate: [adaptive] 
residuals(fit.ada)

# residual estimate: [lasso/Huber] 
residuals(fit.soft)

# residual estimate: [Hard] 
residuals(fit.hard)




#--------------------------------------------------------------
## Plot model output: solution path, cross validation error, estimated residual

# [adaptive]
par(mfrow=c(1,3))
plot_path(fit.ada)
plot_cv(fit.ada)
plot_resid(fit.ada,type = 0, s = 1)


# [soft]
par(mfrow=c(1,3))
plot_path(fit.soft)
plot_cv(fit.soft)
plot_resid(fit.soft,type = 0, s = 1)
#title(sub ='[Soft]: Solution path, Cross-validation error, residual')

# [Hard]
par(mfrow=c(1,3))
plot_path(fit.hard)
plot_cv(fit.hard)
plot_resid(fit.hard,type = 0, s = 1)
#title(sub ='[Hard]: Solution path, Cross-validation error, residual')
par(mfrow=c(1,1))






#--------------------------------------------------------------
## True/Estimated parameter comparison

library(reshape2)
library(ggplot2)

# [Adaptive]
tmp <- data.frame(c(beta0,beta),fit.ada$beta0[,1])
names(tmp) <- c('Simulated parameter','Estimated parameter')
tmp$Index <- 1:(p+1)
df <- melt(tmp,3)
names(df)[2] <- "Comparison"
ggplot(data=df, aes(x=Index, y=value, fill=Comparison)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_bw() + theme(legend.position="bottom") + ggtitle('Adaptive lasso') +
  theme(plot.title = element_text(hjust = 0.5))

# [Lasso/Huber] 
tmp <- data.frame(c(beta0,beta),fit.soft$beta0[,1])
names(tmp) <- c('Simulated parameter','Estimated parameter')
tmp$Index <- 1:(p+1)
df <- melt(tmp,3)
names(df)[2] <- "Comparison"
ggplot(data=df, aes(x=Index, y=value, fill=Comparison)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_bw() + theme(legend.position="bottom") + ggtitle('Lasso') +
  theme(plot.title = element_text(hjust = 0.5))


# [Hard] 
tmp <- data.frame(c(beta0,beta),fit.hard$beta0[,1])
names(tmp) <- c('Simulated','Estimated')
tmp$Index <- 1:(p+1)
df <- melt(tmp,3)
names(df)[2] <- "Comparison"
ggplot(data=df, aes(x=Index, y=value, fill=Comparison)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_bw() + theme(legend.position="bottom") + ggtitle('Hard penalty') +
  theme(plot.title = element_text(hjust = 0.5))

```


## Queries
Please contact authors and creators for any queries related to using the package **robregcc**. 



-   Aditya Mishra: [mailto](mailto:amishra@flatironinstitute.org)
-   Christian Mueller: [mailto](mailto:cmueller@flatironinstitute.org)
