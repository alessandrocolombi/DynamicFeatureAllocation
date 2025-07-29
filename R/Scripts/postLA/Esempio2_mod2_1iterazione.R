# Introduction ------------------------------------------------------------
rm(list = ls())
setwd("C:/Users/colom/DynamicFeatureAllocation/R/Scripts/postLA")
Rcpp::sourceCpp("../../../src/RcppFunctions.cpp")
source("../../Rfunctions.R" )
source("../../genera_img_ilaria.R")

# What is here? --------------------------------------------------------

# In this example I show you that if I know "everything", (the centers and sig2X), then
# I am able to get the perfect solution by setting sig2A very small

# This time I generate data differently wrt Esempio_mod2_1iterazione.R

# Immagini Base -----------------------------------------------------------

load("C:/Users/colom/DynamicFeatureAllocation/R/Scripts/save/img_base.Rdat")

seed = 42
set.seed(seed)

Ttot = 24    # Time length
di   = 8
D    = di*di # Size of the problem

H = 4
zetas = rbind(img_base[1,],img_base[2,],img_base[3,],img_base[4,])
sigma2_X <- sig2_X <- sig2X <- 0.01
A = img_base[1:4, ]
Z = matrix(NA,Ttot,4)
Z[,1] = c(1,1,1,1,
          1,1,1,1,
          1,0,0,0,
          0,0,0,0,
          1,1,1,1,
          0,0,0,0) # cross
Z[,2] = c(0,0,0,0,
          1,1,1,1,
          1,1,1,1,
          1,1,1,1,
          1,1,1,0,
          0,0,0,0)
Z[,3] = c(0,0,0,0,
          0,0,0,0,
          1,1,1,1,
          1,1,1,1,
          0,0,0,0,
          1,1,1,1)
Z[,4] = c(0,0,0,0,
          0,0,0,0,
          0,0,0,0,
          0,0,1,1,
          1,1,1,1,
          1,1,1,1)
X = matrix(0,nrow = Ttot, ncol = D)
for(t in 1:Ttot){
  mean_t_vec = (Z%*%A)[t,]
  X[t,] = mvtnorm::rmvnorm(n = 1, mean = mean_t_vec, sigma = sig2X * diag(D))
}

par(mfrow = c(3,4), mar = c(1,1,1,1) , bty = "l")
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
  plot_img(X_t, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = paste0("Obs. T = ",t))
}


## Filtering  ---------------------------------------------------
N  = 1000
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

beepr::beep()


# Number of new features at each time step
temp = sapply(fit$Path_k[1:Ttot],function(x) x$Nnew)
temp

# Total number of features
nrow(fit$Amat_k) 


# Active labels (just a check)
for(t in 1:Ttot){
  cat("\n t = ",t,"; active: ",fit$Path_k[[t]]$label_actives)
}

plot_multiple_imgs(fit$mean_k, plt_wnd = c(3,4))

# Feature path
Afinal = fit$Amat_k
Zfinal = fit$Zmat_k
for(i in 1:Num_feat_tot){
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












# Check and break down ----------------------------------------------------

for(t in 1:Ttot){
  cat("\n t = ",t,"; active: ",fit$Path_k[[t]]$label_actives)
}
cat("\n")



fit$Path_k[[20]]$label_actives
ActFea = fit$Path_k[[20]]$active_values
plot_multiple_imgs(ActFea, plt_wnd = c(1,2))



## Important changes
# Custom Zfinal
Num_feat_tot = max(fit$cum_Ntot_k)
Zfinal = matrix(0,nrow = Ttot, ncol = Num_feat_tot)
for(t in 1:Ttot){
  Zfinal[t,fit$Path_k[[t]]$label_actives] = 1
}
# Custom Afinal - Filtered mean
Afinal = matrix(NA,nrow = Num_feat_tot, ncol = D)
FiltMean = matrix(NA,nrow = Ttot, ncol = D)
for(t in 1:Ttot){
  Afinal[ fit$Path_k[[t]]$label_actives , ] = fit$Path_k[[t]]$active_values
  FiltMean[t,] = colSums(fit$Path_k[[t]]$active_values)
}
# Filtered mean
plot_multiple_imgs(FiltMean, plt_wnd = c(3,4))
plot_multiple_imgs(fit$mean_k, plt_wnd = c(3,4))





# Number of new features at each time step
temp = sapply(fit$Path_k[1:Ttot],function(x) x$Nnew)
temp