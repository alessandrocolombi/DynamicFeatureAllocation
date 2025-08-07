# Introduction ------------------------------------------------------------
rm(list = ls())
setwd("C:/Users/colom/DynamicFeatureAllocation/R/Scripts/postLA")
Rcpp::sourceCpp("../../../src/RcppFunctions.cpp")
source("../../Rfunctions.R" )
source("../../genera_img_ilaria.R")

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

# K-means on X -----------------------------------------------------------------
### Fare Kmeans sulle immagini osservate non è una buona idea perché non riesce 
### a cogliere le features. Non ragiona in modo additivo

plot_multiple_imgs(X, plt_wnd = c(3,4))

H = 4
km0 = kmeans(x = X, centers = H, nstart = 100)
zeta00 = km0$centers
plot_multiple_imgs(zeta00, plt_wnd = c(2,2))


# Options (true sig2X) -----------------------------------------------------------------
D = di*di
Ttot = 24
N = 500
Bfix = rep(-1,Ttot); Particle_fix = NULL
M1 = 0.5
Psurv = 0.5
sig2A <- sig2_A <- sigma2_A <- 10
proposal_Nnew_1T = NULL
use_VS = TRUE

## Run Model 1  -------------------------------------------------------------

fit1 = Conditional_SeqMonteCarlo( X = X, N = N, D = D, Ttot = Ttot,
                                  Bfix = Bfix, Particle_fix = Particle_fix,  
                                  M1 = M1, Psurv = Psurv,
                                  sigma2_A = sigma2_A,
                                  sigma2_X = sigma2_X,
                                  proposal_Nnew_1T = proposal_Nnew_1T,
                                  use_VS = use_VS)
beepr::beep()

# Number of new features at each time step
temp = sapply(fit1$Path_k[1:Ttot],function(x) x$Nnew)
temp

# Total number of features
nrow(fit1$Amat_k) 
Num_feat_tot = nrow(fit1$Amat_k) 

# Features
plot_multiple_imgs(fit1$Amat_k, plt_wnd = c(3,3))

# Filtered mean
plot_multiple_imgs(fit1$mean_k, plt_wnd = c(3,4))

# Feature path
Afinal = fit1$Amat_k
Zfinal = fit1$Zmat_k
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


## K-means -----------------------------------------------------------------
plot_multiple_imgs(Afinal, plt_wnd = c(3,3))

H = 4
km = kmeans(x = Afinal, centers = H, nstart = 100)
zeta0 = km$centers
plot_multiple_imgs(zeta0, plt_wnd = c(2,2))



## Run Model 2 -------------------------------------------------------------
N  = N
M1 = 0.1
Psurv = 0.5
sigma2_X <- sigma2_X # vedi sopra
sigma2_A <- 1e-4
use_VS <- TRUE
H = H
zetas_mcmc = lapply(1:H, function(h){zeta0[h,]})


fit = CondSMC(X = X, N = N, D = D, Ttot = Ttot, 
              Bfix = Bfix,Particle_fix = Particle_fix,  
              M1 = M1, Psurv = Psurv,
              sigma2_X = sigma2_X,
              sigma2_A = sigma2_A,
              zeta = zetas_mcmc, 
              proposal_Nnew_1T = proposal_Nnew_1T,
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

# Filtered mean
plot_multiple_imgs(fit$mean_k, plt_wnd = c(3,4))


# Features
plot_multiple_imgs(fit$Amat_k, plt_wnd = c(2,2))


# Feature path
Afinal = fit$Amat_k
Zfinal = fit$Zmat_k
for(i in 1:nrow(fit$Amat_k) ){
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




## K-means 2 -----------------------------------------------------------------
plot_multiple_imgs(fit$Amat_k, plt_wnd = c(3,3))

H = 4
km2 = kmeans(x = fit$Amat_k, centers = H, nstart = 100)
zeta2 = km2$centers
plot_multiple_imgs(zeta2, plt_wnd = c(2,2))



















# Options (wrong sig2X) -----------------------------------------------------------------
D = di*di
Ttot = 24
N = 500
Bfix = rep(-1,Ttot); Particle_fix = NULL
M1 = 0.5
Psurv = 0.5
sig2A <- sig2_A <- sigma2_A <- 10
sig2X_fake <- 0.01 + 0
proposal_Nnew_1T = NULL
use_VS = TRUE

## Run Model 1  -------------------------------------------------------------

fit1 = Conditional_SeqMonteCarlo( X = X, N = N, D = D, Ttot = Ttot,
                                  Bfix = Bfix, Particle_fix = Particle_fix,  
                                  M1 = M1, Psurv = Psurv,
                                  sigma2_A = sigma2_A,
                                  sigma2_X = sig2X_fake,
                                  proposal_Nnew_1T = proposal_Nnew_1T,
                                  use_VS = use_VS )
beepr::beep()

# Number of new features at each time step
temp = sapply(fit1$Path_k[1:Ttot],function(x) x$Nnew)
temp

# Total number of features
nrow(fit1$Amat_k) 
Num_feat_tot = nrow(fit1$Amat_k) 

# Features
plot_multiple_imgs(fit1$Amat_k, plt_wnd = c(3,3))

# Filtered mean
plot_multiple_imgs(fit1$mean_k, plt_wnd = c(3,4))

# Feature path
Afinal = fit1$Amat_k
Zfinal = fit1$Zmat_k
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

