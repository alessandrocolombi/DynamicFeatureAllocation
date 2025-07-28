# Introduction ------------------------------------------------------------
rm(list = ls())
setwd("C:/Users/colom/DynamicFeatureAllocation/R/Scripts/postLA")
Rcpp::sourceCpp("../../../src/RcppFunctions.cpp")
source("../../Rfunctions.R" )
source("../../genera_img_ilaria.R")

# Custom functions --------------------------------------------------------


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


# Options -------------------------------------------

seed = 42
set.seed(seed)

N  = 500
G  = 10 + 1
M1 = 0.25
use_VS = TRUE
H  = 4 # number of centers
sigma2_ker = 1e-4 # variance base measure

par_sig2X = set_par_invgamma(media = 0.001, var = 0.1)
par_sig2A = set_par_invgamma(media = 0.01,   var = 2)
par_gamma = set_par_gamma(media = 0.5,    var = 0.001)
par_c     = set_par_gamma(media = 4,    var = 1)

# Gibbs sampler structures:

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
setPsurv = 0.5
setM1    = 0.1

var_prop[1] = 0.1 # initial adaptive variance
BetaProcess_params[1,1] = setM1/setPsurv  # gamma at first iteration
BetaProcess_params[1,2] = 1/setPsurv - 1  # c at first iteration
BetaProcess_params[,3]  = 0 # set sigma = 0 for all iterations
M1_mcmc[1]    = computeM1(BetaProcess_params[1,])
Psurv_mcmc[1] = computePsurv(BetaProcess_params[1,])

sigma2_XA_mcmc[1,1] = sig2_X # initial value for sigma2_X
sigma2_XA_mcmc[1,2] = 1e-5  # initial value for sigma2_A

H_mcmc[1] = H; Hstar_mcmc[1] = 0
zetas_mcmc[[1]] = lapply(1:H, function(h){zetas[h,]})
# zetas_mcmc[[1]] = lapply(1:H, function(h){rep(0,D)})



# First iteration (g=1) ---------------------------------------------------

sink("log.txt")
g = 1
Particles_PG[[1]] = CondSMC(X = X, N = N, D = D, Ttot = Ttot, 
                            Bfix = rep(-1,Ttot), Particle_fix = NULL,  
                            M1 = M1_mcmc[1], Psurv = Psurv_mcmc[1] ,
                            sigma2_X = sigma2_XA_mcmc[1,1],
                            sigma2_A = sigma2_XA_mcmc[1,2],
                            zeta = zetas_mcmc[[1]], 
                            proposal_Nnew_1T = NULL,
                            use_VS = use_VS)
sink()
temp = sapply(Particles_PG[[1]]$Path_k[1:Ttot],function(x) x$Nnew)
temp

dim(Particles_PG[[1]]$mean_k)
dim(Particles_PG[[1]]$Amat_k)
dim(Particles_PG[[1]]$Zmat_k)

media = Particles_PG[[1]]$mean_k
plot_multiple_imgs(media, plt_wnd = c(3,4))

Amat = Particles_PG[[1]]$Amat_k
plot_multiple_imgs(Amat, plt_wnd = c(3,4))

# Run (g >= 2) ---------------------------------------------------

g = 2
UpdateSig2X = FALSE; UpdateSig2A = FALSE
Updatec = FALSE; Updategamma = FALSE
updateGrAlloc = FALSE; updateCenters = FALSE; upNcenters = FALSE
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
  
  H_mcmc[g] = length(Particles_PG[[g]]$gr_alloc_k)
}
close(pb)


# Analysis ----------------------------------------------------------------


par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
plot(sigma2_XA_mcmc[,1], type = "l")

# Inference sigma2_X 
brn = (G-1)/2
quantiles_sig2X =  quantile(sigma2_XA_mcmc[brn:G,1], probs = c(0.0025,0.5,0.9975))
dens_sig2X = density(sigma2_XA_mcmc[brn:G,1])

par(mfrow = c(1,1), mar = c(4,4,2,2), bty = "l")
plot(0,0,type = "n",  
     main = "Variance Sig2X", xlab = "", ylab = "Dens.",
     xlim = c(0,max(dens_sig2X$x)),
     ylim = c(0,max(dens_sig2X$y)) )
lines(x = dens_sig2X$x, y = dens_sig2X$y, lwd = 3, col = "darkgreen")
abline(v = quantiles_sig2X[2], col = c("darkgreen"), lwd = 2, lty = 3 )
abline(v = sigma2_X, col = c("black"), lwd = 3, lty = 3 )

# Gamma
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

