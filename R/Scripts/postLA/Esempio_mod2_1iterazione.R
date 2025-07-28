# Introduction ------------------------------------------------------------
rm(list = ls())
setwd("C:/Users/colom/DynamicFeatureAllocation/R/Scripts/postLA")
Rcpp::sourceCpp("../../../src/RcppFunctions.cpp")
source("../../Rfunctions.R" )
source("../../genera_img_ilaria.R")

# What is here? --------------------------------------------------------

# In this example I show you that if I know "everything", (the centers and sig2X), then
# I am able to get the perfect solution by setting sig2A very small

# Immagini Base -----------------------------------------------------------

load("C:/Users/colom/DynamicFeatureAllocation/R/Scripts/save/img_base.Rdat")

# Data generation ---------------------------------------------------------
seed = 2927
set.seed(seed)

Ttot <- Ti <- 24    # Time length
di   = 8
D    = di*di # Size of the problem

sigma2_X <- sig2_X <- sig2X <- 0.01
sim_data = sim_images(sig2_X,Ti)

A = sim_data$A
Z = sim_data$Z
X = sim_data$X
zetas = rbind(img_base[1,],img_base[2,],img_base[3,],img_base[4,])

Nfeat_tot = nrow(A)
plot_multiple_imgs(A, plt_wnd = c(2,2))


par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  mean_t_vec = (Z%*%A)[t,]
  mean_t = matrix(data = mean_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img(mean_t, center_value = 0.5, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = paste0("Mean T = ",t))
}

par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  Xt_vec = X[t,]
  X_t = matrix(data = Xt_vec, nrow = di, ncol = di, byrow = F)
  plot_img(X_t, center_value = 0.5, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = paste0("Obs. T = ",t))
}



## Filtering  ---------------------------------------------------
N  = 500
M1 = 0.1
Psurv = 0.5
sigma2_X <- sigma2_X # vedi sopra
sigma2_A <- 1e-5
use_VS <- TRUE
H = 4
zetas_mcmc = lapply(1:H, function(h){zetas[h,]})


fit = CondSMC(X = X, N = N, D = D, Ttot = Ttot, 
              Bfix = rep(-1,Ttot), 
              Particle_fix = NULL,  
              M1 = M1, 
              Psurv = Psurv,
              sigma2_X = sigma2_X,
              sigma2_A = sigma2_A,
              zeta = zetas_mcmc, 
              proposal_Nnew_1T = NULL,
              use_VS = use_VS)

# Number of new features at each time step
temp = sapply(fit$Path_k[1:Ttot],function(x) x$Nnew)
temp

# Total number of features
nrow(fit$Amat_k)

# Filtered mean
plot_multiple_imgs(fit$mean_k, plt_wnd = c(3,4))

# Feature path
Afinal = fit$Amat_k; dim(Afinal)
Zfinal = fit$Zmat_k; dim(Zfinal)
for(i in 1:nrow(Afinal)){
  # plot figure
  par(mfrow = c(1,2), mar = c(2,2,2,1), bty = "l")
  res_t_vec = Afinal[i,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
  
  # plot survival path
  plot(x = 1:Ttot, y = Zfinal[,i], pch = 16, xlab = "", ylab = "")
  segments(x0 = 1:Ttot, x1 = 1:Ttot, y0 = rep(0,Ttot), y1 = Zfinal[,i])
}

## Filtering - Tests  ---------------------------------------------------

# Here I want to play with the hyperparameters to test different settings
# - works if N is smaller (N=10)
# - "fails" if N is small (N=10) and M1 is large (M1=1)
# - works if N is large (N=100;N=1000) and M1 is large (M1=1)
# - fails if sig2A = sig2X; works if sig2A = sig2X/10
# - works if I set H=6 with two centers equal to rep(0,D). 
# - total DISASTER if I set H=4, all equals to 0

N  = 500
M1 = 1
Psurv = 0.5
sigma2_X <- sigma2_X # vedi sopra
sigma2_A <- 1e-5
use_VS <- TRUE
H = 6
zetas_mcmc = vector("list", H)
zetas_mcmc[1:4] = lapply(1:4, function(h){zetas[h,]})
zetas_mcmc[5:6] = lapply(5:6, function(h){rep(0,D)})


fit = CondSMC(X = X, N = N, D = D, Ttot = Ttot, 
              Bfix = rep(-1,Ttot), 
              Particle_fix = NULL,  
              M1 = M1, 
              Psurv = Psurv,
              sigma2_X = sigma2_X,
              sigma2_A = sigma2_A,
              zeta = zetas_mcmc, 
              proposal_Nnew_1T = NULL,
              use_VS = use_VS)

# Number of new features at each time step
temp = sapply(fit$Path_k[1:Ttot],function(x) x$Nnew)
temp

# Total number of features
nrow(fit$Amat_k)

# Filtered mean
plot_multiple_imgs(fit$mean_k, plt_wnd = c(3,4))

# Feature path
Afinal = fit$Amat_k; dim(Afinal)
Zfinal = fit$Zmat_k; dim(Zfinal)
for(i in 1:nrow(Afinal)){
  # plot figure
  par(mfrow = c(1,2), mar = c(2,2,2,1), bty = "l")
  res_t_vec = Afinal[i,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
  
  # plot survival path
  plot(x = 1:Ttot, y = Zfinal[,i], pch = 16, xlab = "", ylab = "")
  segments(x0 = 1:Ttot, x1 = 1:Ttot, y0 = rep(0,Ttot), y1 = Zfinal[,i])
}

## Filtering - empty centers  ---------------------------------------------------

# If I set H=1 and zeta = 0, then I must be able to go back to model 1 --> yes

N  = 500
M1 = 0.5
Psurv = 0.5
sigma2_X <- sigma2_X # vedi sopra
sigma2_A <- 1
use_VS <- TRUE
H = 4
zetas_mcmc = vector("list", H)
zetas_mcmc = lapply(zetas_mcmc, function(x){rep(0,D)})


fit = CondSMC(X = X, N = N, D = D, Ttot = Ttot, 
              Bfix = rep(-1,Ttot), 
              Particle_fix = NULL,  
              M1 = M1, 
              Psurv = Psurv,
              sigma2_X = sigma2_X,
              sigma2_A = sigma2_A,
              zeta = zetas_mcmc, 
              proposal_Nnew_1T = NULL,
              use_VS = use_VS)

# Number of new features at each time step
temp = sapply(fit$Path_k[1:Ttot],function(x) x$Nnew)
temp

# Total number of features
nrow(fit$Amat_k)

# Filtered mean
plot_multiple_imgs(fit$mean_k, plt_wnd = c(3,4))

# Feature path
Afinal = fit$Amat_k; dim(Afinal)
Zfinal = fit$Zmat_k; dim(Zfinal)
for(i in 1:nrow(Afinal)){
  # plot figure
  par(mfrow = c(1,2), mar = c(2,2,2,1), bty = "l")
  res_t_vec = Afinal[i,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
  
  # plot survival path
  plot(x = 1:Ttot, y = Zfinal[,i], pch = 16, xlab = "", ylab = "")
  segments(x0 = 1:Ttot, x1 = 1:Ttot, y0 = rep(0,Ttot), y1 = Zfinal[,i])
}
