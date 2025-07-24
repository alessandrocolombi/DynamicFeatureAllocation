# Introduction ------------------------------------------------------------
rm(list = ls())
setwd("C:/Users/colom/DynamicFeatureAllocation/R/Scripts")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../Rfunctions.R")
source("../genera_img_ilaria.R")

# Custom functions --------------------------------------------------------


# Immagini Base -----------------------------------------------------------

di   = 8
D    = di*di # Size of the problem

par(mfrow = c(3,2), mar = c(2,2,2,1), bty = "l")
A1 = rep(0, D)
idx = c(2, 9:11, 18)
A1[idx] = rep(1, length(idx))
image((matrix(data = A1[1:D], nrow = di, ncol = di, byrow = F)), 
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 

A2 = rep(0, D)
idx = c(5:8, 13,16, 21,24, 29:32)
A2[idx] = rep(1, length(idx))
image((matrix(data = A2[1:D], nrow = di, ncol = di, byrow = F)),
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 


A3 = rep(0, D)
idx = c(33, 41:42, 49:51, 57:60)
A3[idx] = rep(1, length(idx))
image((matrix(data = A3[1:D], nrow = di, ncol = di, byrow = F)), 
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 

A4 = rep(0, D)
idx = c(38:40, 47, 55, 63)
A4[idx] = rep(1, length(idx))
image((matrix(data = A4[1:D], nrow = di, ncol = di, byrow = F)), 
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 

A5 = A2[D:1]
image((matrix(data = A5[1:D], nrow = di, ncol = di, byrow = F)), 
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 

A6 = c(rep(c(rep(0,3), 1,1,rep(0,3)), di))
image((matrix(data = A6[1:D], nrow = di, ncol = di, byrow = F)), 
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 



img_base = matrix(NA,nrow = 6, ncol = D)
img_base[1,] = A1; img_base[2,] = A2
img_base[3,] = A3; img_base[4,] = A4
img_base[5,] = A5; img_base[6,] = A6

# Sim. data from model 2  -------------------------------------------------

seed = 42
set.seed(seed)

Ttot = 24    # Time length
di   = 8
D    = di*di # Size of the problem

H = 4
zetas = rbind(A1,A2,A3,A4)
Psurv_true = 0.91
M1true = 1
c = 1/Psurv_true - 1
gamma = M1true * Psurv_true

sig2_A = 0.0001; sigma2_X = 0.01


# sim_data = sim_images_ale( img_base = zetas, Ti = Ttot, 
#                            sig2X = sigma2_X, sig2A = sigma2_A,
#                            Psurv = Psurv_true)
sim_data = sim_images_ale_mod2( zetas = zetas, Ti = Ttot, 
                                sig2X = 0.001, sig2A = 0.001,
                                Psurv = 0.89, M1 = 0.33 )

A = sim_data$A
Z = sim_data$Z
X = sim_data$X

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




# Independent fit ---------------------------------------------------------
# stop("Vecchio, non runnare")
# seed = 42
# set.seed(seed)
# 
# N = 10 # Number of particles
# c = 1/Psurv_true - 1
# gamma = M1true * Psurv_true
# sigma = 0
# sigma2_ker = 0.1 # Variance base measure
# fit_smc = lapply(1:H, function(h){
#   SeqMonteCarlo(X,N,D,Ttot,
#                 M1 = M1true,Psurv = Psurv_true,
#                 sigma2_A = sigma2_ker,sigma2_X, 
#                 mu0 = zetas[h,], 
#                 use_VS = FALSE)
# })
# 
# mean_all = lapply(fit_smc, function(x){x$Final_Est})
# for(h in 1:H){
#   par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
#   for(t in 1:Ttot){
#     res_t_vec = mean_all[[h]][t,]
#     res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
#     plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
#               horizontal = F, main = paste0("Gr. ",h,"; Est. T = ",t))
#   }
# }
# 
# 
# mean_all_sum = Reduce(`+`, mean_all)
# result = mean_all_sum
# par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
# for(t in 1:Ttot){
#   res_t_vec = result[t,]
#   res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
#   plot_img(res_t, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
#            horizontal = F, main = paste0("Est. T = ",t))
# }


# PG with independent chains and fix number of centers -------------------------------------------
# stop("Vecchio, non runnare")
# seed = 42
# set.seed(seed)
# 
# N  = 10
# G  = 4 + 1
# M1 = M1true
# use_VS = FALSE
# H  = 4 # number of centers
# sigma2_ker = 1e-5 # variance base measure
# 
# par_sig2X = set_par_invgamma(media = 0.01, var = 0.1)
# par_sig2A = set_par_invgamma(media = 0.1,   var = 2)
# par_gamma = set_par_gamma(media = 0.25,    var = 0.1)
# par_c     = set_par_gamma(media = 1.1,    var = 1)
# 
# # Gibbs sampler structures:
# 
# # Beta process (gamma, c, sigma)
# BetaProcess_params = matrix(NA, nrow = G, ncol = 3)
# colnames(BetaProcess_params) = c("gamma","c","sigma")
# 
# # Variance
# sigma2_X_mcmc = matrix(NA, nrow = G, ncol = 1)
# colnames(sigma2_X_mcmc) = c("sigma2_X")
# 
# # Centers 
# zetas_mcmc = lapply(1:G, function(g){
#   matrix(NA, nrow = H, ncol = D)
# })
# 
# # Auxiliaries
# var_prop <- M1_mcmc <- Psurv_mcmc <- rep(0,G)
# 
# # Latent allocation process
# G_th = lapply(1:H, function(h){
#   list("Particles_PG" = vector("list", G))
# })
# # G_th[[h]] is the Particle Gibbs for the h-th chain (h=1,...,H)
# Ktot = rep(NA,G)
# 
# 
# # Initialization
# var_prop[1] = 0.1 # initial adaptive variance
# BetaProcess_params[1,1] = gamma # gamma at first iteration
# BetaProcess_params[1,2] = c # c at first iteration
# BetaProcess_params[,3]  = 0 # set sigma = 0 for all iterations
# M1_mcmc[1] = computeM1(BetaProcess_params[1,])
# Psurv_mcmc[1] = computePsurv(BetaProcess_params[1,])
# 
# sigma2_X_mcmc[1,1] = 0.001 # initial value for sigma2_X
# 
# zetas_mcmc[[1]] = zetas
# 
# g = 1
# temp = lapply(1:H, function(h){
#   
#   Conditional_SeqMonteCarlo( X = X, N = N, D = D, Ttot = Ttot,
#                              mu0 = zetas_mcmc[[g]][h,],
#                              Bfix = rep(-1,Ttot), Particle_fix = NULL,  
#                              M1 = M1_mcmc[1],Psurv = Psurv_mcmc[1],
#                              sigma2_A = sigma2_ker, 
#                              sigma2_X = sigma2_X_mcmc[g,1],
#                              use_VS = use_VS)
# })
# for(h in 1:H){
#   # Save the sampled paths in the right spot
#   G_th[[h]]$Particles_PG[[g]] = temp[[h]]
# }
# 
# 
# 
# g = 2
# UpdateSig2X = FALSE
# Updatec = FALSE; Updategamma = FALSE
# Updatezetas = FALSE; # IMPLEMENTARE L'AGGIORNAMENTO DEI CENTRI!
# pb = txtProgressBar(min = 2, max = G, initial = 2) 
# for(g in 2:G){
#   
#   # cat("\n ++++++ g = ",g," ++++++ \n")
#   setTxtProgressBar(pb,g)
#   
#   ## Assemble G_th for h=1,...,H to define G_t (for each t=1,...,Ttot)
#   mean_vec_1T = Reduce(`+`, lapply(G_th, function(x){x$Particles_PG[[g-1]]$mean_k}))
#   Zmat        = Reduce(cbind, lapply(G_th, function(x){x$Particles_PG[[g-1]]$Zmat_k})) 
#   Amat        = Reduce(rbind, lapply(G_th, function(x){x$Particles_PG[[g-1]]$Amat_k})) 
#   Ktot[g-1]   = ifelse(is.null(nrow(Amat)), 0, nrow(Amat))
#   
#   if(Ktot[g-1] == 0)
#     cat("\n Ktot[",g-1,"] è 0, speriamo vada tutto bene \n")
#   
#   
#   # Update sigma2_X
#   if(UpdateSig2X){
#     suff_stat_X = sum( apply(X-mean_vec_1T, 1, function(x){t(x)%*%x}) )
#     sigma2_X_mcmc[g,1] = r_fullcond_sigma2_X(par_sig2X[1], par_sig2X[2], D, Ttot, suff_stat_X)
#     # sigma2_XA_mcmc[g,1] = sigma2_X # set true values
#   }else{
#     sigma2_X_mcmc[g,1] = sigma2_X_mcmc[g-1,1]
#   }
#   
#   
#   # Update_gamma
#   if(Updategamma){
#     BetaProcess_params[g,1] = r_fullcond_gamma(par_gamma[1], par_gamma[2], 
#                                                Ktot[g-1], Ttot, 
#                                                BetaProcess_params[g-1,2], 
#                                                BetaProcess_params[g-1,3])
#     # BetaProcess_params[g,1] = gamma # set true values
#   }else{
#     BetaProcess_params[g,1] = BetaProcess_params[g-1,1]
#   }
#   
#   # Update c
#   if(Updatec){
#     suff_stat_Z = stat_suff_Z(Zmat = Zmat)
#     BetaProcess_params[g,2] = r_fullcond_c(BetaProcess_params[g-1,2], 
#                                            par_c[1], par_c[2],
#                                            gamma = BetaProcess_params[g,1], 
#                                            sigma = BetaProcess_params[g,3],
#                                            NumTrailTot = suff_stat_Z[1], 
#                                            NumFailTot  = suff_stat_Z[3], 
#                                            NnewTot = Ktot[g-1], Ttot = Ttot,
#                                            var_prop[g-1])
#     # BetaProcess_params[g,2] = c # set true values
#     # Update var_prop
#     var_prop[g] = 0.01
#   }else{
#     BetaProcess_params[g,2] = BetaProcess_params[g-1,2]
#     var_prop[g] = 0.01
#   }
#   
#   # Update M1 and Psurv
#   M1_mcmc[g]    = computeM1(BetaProcess_params[g,])
#   Psurv_mcmc[g] = computePsurv(BetaProcess_params[g,])
#   
#   # cat("\n M1 = ",M1_mcmc[g],"; Psurv = ",Psurv_mcmc[g],"\n")
#   
#   
#   # Update zetas
#   if(Updatezetas){
#     stop("Non l'ho ancora fatto :(")
#   }else{
#     zetas_mcmc[[g]] = zetas_mcmc[[g-1]]
#   }
#   # Run Conditional SMC
#   temp = lapply(1:H, function(h){
#     
#     Conditional_SeqMonteCarlo( X = X, N = N, D = D, Ttot = Ttot,
#                                mu0 = zetas_mcmc[[g]][h,],
#                                Bfix = G_th[[h]]$Particles_PG[[g-1]]$B_k,
#                                Particle_fix = G_th[[h]]$Particles_PG[[g-1]]$Path_k,  
#                                M1 = M1_mcmc[g],Psurv = Psurv_mcmc[g],
#                                sigma2_A = sigma2_ker, 
#                                sigma2_X = sigma2_X_mcmc[g,1],
#                                use_VS = use_VS)
#     
#   })
#   for(h in 1:H){
#     # Save the sampled paths in the right spot
#     G_th[[h]]$Particles_PG[[g]] = temp[[h]]
#   }
# 
# }
# close(pb)
# 
# # Inference sigma2_X and sigma2_A
# quantiles   = apply(sigma2_X_mcmc, 2, quantile, probs = c(0.0025,0.5,0.9975))
# dens_thetas = apply(sigma2_X_mcmc, 2, function(x){density(x)$y})
# 
# par(mfrow = c(1,1), mar = c(4,4,2,2), bty = "l")
# plot(0,0,type = "n",  main = "Variances", xlab = "", ylab = "Dens.",
#      xlim = range(quantiles), ylim = c(0,max(dens_thetas)) )
# hist(sigma2_X_mcmc[,1], col = "darkgreen", nclass = "FD", freq = F, add = T)
# abline(v = apply(sigma2_X_mcmc, 2, median), col = c("darkgreen"), lwd = 2, lty = 3 )
# legend("topright", legend = colnames(sigma2_X_mcmc) ,col = c("darkgreen"), pch = 16)
# 
# 
# 
# # Gamma
# summary(BetaProcess_params[,1])
# plot(BetaProcess_params[,1], type = "l", main = "gamma")
# 
# # c
# summary(BetaProcess_params[,2])
# plot(BetaProcess_params[,2], type = "l", main = "c")
# 
# par(mfrow = c(1,3), bty = "l")
# plot(M1_mcmc, type = "l", main = "M1")
# plot(Psurv_mcmc, type = "l", main = "Prob. surv.");
# plot(Ktot, type = "l", main = "Num. features")
# 
# 
# 
# 
# Kt_mcmc_h = lapply(G_th, function(x){
#   t(sapply(x$Particles_PG, function(x){apply(x$Zmat_k,1,sum)}))
# })
# # Kt_mcmc_h[[h]] is a GxT matrix such that Kt_mcmc_h[[h]][,t] is the traceplot of the number
# # of feature in the h-th group at time t
# Kt_mcmc = Reduce(`+`, Kt_mcmc_h)
# par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
# for(t in 1:Ttot){
#   plot(Kt_mcmc[,t], type = "l", main = paste0("K_",t))
# }
# # Kt_mcmc is a GxT matrix such that Kt_mcmc[,t] is the traceplot of the 
# # total number of feature across all groups at time t
# 
# 
# # Final estimated filtered density in each group
# burnin = (G-1)/2
# filt_dens_h = lapply(G_th, function(x){
#   result <- Reduce(`+`, lapply(x$Particles_PG[(G-burnin+1):(G)], function(xx){xx$mean_k}) )
#   result = result / length((G-burnin+1):(G))
# })
# 
# for(h in 1:H){
#   par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
#   for(t in 1:Ttot){
#     res_t_vec = filt_dens_h[[h]][t,]
#     res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
#     plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
#               horizontal = F, main = paste0("Gr. ",h,"; Est. T = ",t))
#   }
# }
# 
# # Final estimated filtered density 
# filt_dens = Reduce(`+`, filt_dens_h )
# par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
# for(t in 1:Ttot){
#   res_t_vec = filt_dens[t,]
#   res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
#   plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
#             horizontal = F, main = paste0("Est. T = ",t))
# }


# PG with fix number of centers -------------------------------------------

seed = 42
set.seed(seed)

N  = 10
G  = 10 + 1
M1 = 0.25
use_VS = FALSE
H  = 4 # number of centers
sigma2_ker = 1e-4 # variance base measure

par_sig2X = set_par_invgamma(media = 0.1, var = 0.1)
par_sig2A = set_par_invgamma(media = 1,   var = 1) # I can not learn this
par_gamma = set_par_gamma(media = 0.1,    var = 0.5)
par_c     = set_par_gamma(media = 10,    var = 5)

# Gibbs sampler structures:
omega = 1
# Beta process (gamma, c, sigma)
BetaProcess_params = matrix(NA, nrow = G, ncol = 3)
colnames(BetaProcess_params) = c("gamma","c","sigma")

# Variance
sigma2_XA_mcmc = matrix(NA, nrow = G, ncol = 2)
colnames(sigma2_XA_mcmc) = c("sigma2_X","sigma2_A")

# Group-related quantities 
zetas_mcmc = vector("list", G)
H_mcmc <- Hstar_mcmc <- rep(NA,G)

# Auxiliaries
var_prop <- M1_mcmc <- Psurv_mcmc <- rep(0,G)

# Latent allocation process
Particles_PG = vector("list", G)
Ktot = rep(NA,G)


# Initialization
setPsurv = 0.2
setM1    = 0.1

var_prop[1] = 0.1 # initial adaptive variance
BetaProcess_params[1,1] = setM1/setPsurv  # gamma at first iteration
BetaProcess_params[1,2] = 1/setPsurv - 1  # c at first iteration
BetaProcess_params[,3]  = 0 # set sigma = 0 for all iterations
M1_mcmc[1]    = computeM1(BetaProcess_params[1,])
Psurv_mcmc[1] = computePsurv(BetaProcess_params[1,])

sigma2_XA_mcmc[1,1] = sigma2_X # initial value for sigma2_X
sigma2_XA_mcmc[1,2] = 0.1  # initial value for sigma2_A

H_mcmc[1] = H; Hstar_mcmc[1] = 0
zetas_mcmc[[1]] = lapply(1:H, function(h){zetas[h,]})
# zetas_mcmc[[1]] = lapply(1:H, function(h){rep(0,D)})

g = 1
Particles_PG[[1]] = CondSMC(X = X, N = N, D = D, Ttot = Ttot,
                            Bfix = rep(-1,Ttot), Particle_fix = NULL,
                            M1 = M1_mcmc[1], Psurv = Psurv_mcmc[1] ,
                            sigma2_X = sigma2_XA_mcmc[1,1],
                            sigma2_A = sigma2_XA_mcmc[1,2],
                            zeta = zetas_mcmc[[1]],
                            use_VS = use_VS)

# Amat          = A # set true values
# mean_vec_1T   = (Z%*%A) # set true values
# K             = 4 # set true values
# Nnew_1T       = c(1,1,0,1,0,1,0,0,0,0,0,0)
# Zmat          = Z # set true values
# Ktot[1]       = nrow(Amat)

g = 2
UpdateSig2X = FALSE; UpdateSig2A = FALSE
Updatec = TRUE; Updategamma = TRUE
updateGrAlloc = TRUE; updateCenters = TRUE; upNcenters = FALSE
pb = txtProgressBar(min = 2, max = G, initial = 2) 
for(g in 2:G){
  
  setTxtProgressBar(pb,g)
  
  # Extract sufficient statistics from Particles_PG
  mean_vec_1T = Particles_PG[[g-1]]$mean_k
  Zmat        = Particles_PG[[g-1]]$Zmat_k
  Amat        = Particles_PG[[g-1]]$Amat_k
  Ktot[g-1]   = ifelse(is.null(nrow(Amat)), 0, nrow(Amat))
  
  if(Ktot[g-1] == 0)
    cat("\n Ktot[",g-1,"] è 0, speriamo vada tutto bene \n")
  
  
  # Update sigma2_X
  if(UpdateSig2X){
    suff_stat_X = sum( apply(X-mean_vec_1T, 1, function(x){t(x)%*%x}) )
    sigma2_XA_mcmc[g,1] = r_fullcond_sigma2_X(par_sig2X[1], par_sig2X[2], D, Ttot, suff_stat_X)
  }else{
    sigma2_XA_mcmc[g,1] = sigma2_XA_mcmc[g-1,1]
  }
  
  # Update sigma2_A
  if(UpdateSig2A){
    suff_stat_A = 0
    for(h in 1:H){
      ind_h = which(Particles_PG[[g-1]]$gr_alloc_k == h)
      Amat_h = Amat[ind_h,]
      suff_stat_A = suff_stat_A + 
                    sum( apply(Amat_h, 1, 
                               function(x){t(x-zetas_mcmc[[g-1]][[h]])%*%(x-zetas_mcmc[[g-1]][[h]])}) )
    }
    # suff_stat_A = sum( apply(Amat, 1, function(x){t(x)%*%x}) )
    sigma2_XA_mcmc[g,2] = r_fullcond_sigma2_A(par_sig2A[1], 
                                              par_sig2A[2], 
                                              D, 
                                              Ktot[g-1], 
                                              suff_stat_A)
  }else{
    sigma2_XA_mcmc[g,2] = sigma2_XA_mcmc[g-1,2]
  }
  
  # Update_gamma
  if(Updategamma){
    BetaProcess_params[g,1] = r_fullcond_gamma(par_gamma[1], par_gamma[2], 
                                               Ktot[g-1], Ttot, 
                                               BetaProcess_params[g-1,2], 
                                               BetaProcess_params[g-1,3])
    # BetaProcess_params[g,1] = gamma # set true values
  }else{
    BetaProcess_params[g,1] = BetaProcess_params[g-1,1]
  }
  
  # Update c
  if(Updatec){
    suff_stat_Z = stat_suff_Z(Zmat = Zmat)
    BetaProcess_params[g,2] = r_fullcond_c(BetaProcess_params[g-1,2], 
                                           par_c[1], par_c[2],
                                           gamma = BetaProcess_params[g,1], 
                                           sigma = BetaProcess_params[g,3],
                                           NumTrailTot = suff_stat_Z[1], 
                                           NumFailTot  = suff_stat_Z[3], 
                                           NnewTot = Ktot[g-1], Ttot = Ttot,
                                           var_prop[g-1])
    # BetaProcess_params[g,2] = c # set true values
    # Update var_prop
    var_prop[g] = 0.01
  }else{
    BetaProcess_params[g,2] = BetaProcess_params[g-1,2]
    var_prop[g] = 0.01
  }
  
  # Update M1 and Psurv
  M1_mcmc[g]    = computeM1(BetaProcess_params[g,])
  Psurv_mcmc[g] = computePsurv(BetaProcess_params[g,])
  
  # cat("\n M1 = ",M1_mcmc[g],"; Psurv = ",Psurv_mcmc[g],"\n")
  
  
  # Update zetas
  temp = r_fullcond_zeta( data = Amat, 
                          gr_alloc_old = Particles_PG[[g-1]]$gr_alloc_k, 
                          centers_old = zetas_mcmc[[g-1]],  
                          sig2_A = sigma2_XA_mcmc[g,2],
                          sig2_ker = sigma2_ker,
                          updateGrAlloc = updateGrAlloc, 
                          updateCenters = updateCenters, 
                          upNcenters = upNcenters )
  zetas_mcmc[[g]] = temp$centers
  H_mcmc[g] = temp$H; Hstar_mcmc[g] = temp$Hstar
  gr_alloc_new = c(temp$gr_alloc)
    
  # Update cluster labels in conditioning particles
  Particle_fix = Particles_PG[[g-1]]$Path_k
  for(t in 1:Ttot){
    Particle_fix[[t]]$gr_alloc = gr_alloc_new[Particle_fix[[t]]$label_actives]
    Particle_fix[[t]]$gr_card  = tabulate(Particle_fix[[t]]$gr_alloc, 
                                          nbins = H_mcmc[g]+Hstar_mcmc[g])
  }
  
  # Run Conditional SMC
  Particles_PG[[g]] = CondSMC( X = X, N = N, D = D, Ttot = Ttot, 
                               Bfix = Particles_PG[[g-1]]$B_k,
                               Particle_fix = Particle_fix,
                               M1 = M1_mcmc[g], Psurv = Psurv_mcmc[g] ,
                               sigma2_X = sigma2_XA_mcmc[g,1],
                               sigma2_A = sigma2_XA_mcmc[g,2],
                               zeta = zetas_mcmc[[g]], 
                               use_VS = use_VS )
  
  H_mcmc[g] = length(unique(Particles_PG[[g-1]]$gr_alloc_k))
  
  
  if(g%%10 == 0){
    cat("\n ")
    cat("nrow(Amat) = ",nrow(Amat))
    
    # Plot mean_vec_1T
    par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
    for(t in 1:Ttot){
      res_t_vec = mean_vec_1T[t,]
      res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
      plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
                horizontal = F, main = paste0("mean T=",t) )
    }
    
    # Plot zetas
    par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
    for(i in 1:length(zetas_mcmc[[g]])){
      # plot figure
      res_t_vec = zetas_mcmc[[g]][[i]]
      res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
      plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
                horizontal = F, main = paste0("center[",i,"], g = ",g) )
    }
    
  }
  
}
close(pb)

# Inference sigma2_X and sigma2_A
quantiles   = apply(sigma2_XA_mcmc, 2, quantile, probs = c(0.0025,0.5,0.9975))
dens_thetas = apply(sigma2_XA_mcmc, 2, function(x){density(x)$y})

par(mfrow = c(1,1), mar = c(4,4,2,2), bty = "l")
plot(0,0,type = "n",  main = "Variances", xlab = "", ylab = "Dens.",
     xlim = range(quantiles), ylim = c(0,max(dens_thetas)) )
hist(sigma2_XA_mcmc[,1], col = "darkgreen", nclass = "FD", freq = F, add = T)
hist(sigma2_XA_mcmc[,2], col = "darkred", nclass = "FD", freq = F, add = T)
abline(v = apply(sigma2_XA_mcmc, 2, median), col = c("darkgreen"), lwd = 2, lty = 3 )
legend("topright", legend = colnames(sigma2_XA_mcmc) ,col = c("darkgreen","darkred"), pch = 16)



# Gamma
par(mfrow = c(1,1), mar = c(4,4,2,2), bty = "l")
summary(BetaProcess_params[,1])
plot(BetaProcess_params[,1], type = "l", main = "gamma")

# c
summary(BetaProcess_params[,2])
plot(BetaProcess_params[,2], type = "l", main = "c")

par(mfrow = c(1,3), bty = "l")
plot(M1_mcmc, type = "l", main = "M1")
plot(Psurv_mcmc, type = "l", main = "Prob. surv.");
plot(Ktot, type = "l", main = "Num. features")


# Traceplot for each time t
Kt_mcmc = t(sapply(Particles_PG, function(x){apply(x$Zmat_k,1,sum)}))
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  plot(Kt_mcmc[,t], type = "l", main = paste0("K_",t))
}

# Traceplot for each time t within each group h
# ..... TODO! .....


# Final estimated filtered density
burnin = (G-1)/2
result <- Reduce(`+`, lapply(Particles_PG[(G-burnin+1):(G)], function(x){x$mean_k}) )
result = result / length((G-burnin+1):(G))

par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  res_t_vec = result[t,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("Est. T = ",t))
}

# Final estimated filtered density in each group
# ..... TODO! .....


Afinal = Particles_PG[[G]]$Amat_k; dim(Afinal)
Zfinal = Particles_PG[[G]]$Zmat_k; dim(Zfinal)

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



save_img = FALSE
if(save_img)
  pdf("img/zeta_est_mcmc.pdf")
for(g in 1:G){
  Hg = length(zetas_mcmc[[g]])
  par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
  for(i in 1:Hg){
    res_t_vec = zetas_mcmc[[g]][[i]]
    res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
    plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
              horizontal = F, main = paste0("Zeta ",i,"; g = ",g))
  }
}
if(save_img)
  dev.off()