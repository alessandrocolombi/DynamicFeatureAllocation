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




# Basic options -----------------------------------------------------------

seed = 42
set.seed(seed)

N  = 500
G  = 15 + 1
# M1 = M1init
use_VS = TRUE

par_sig2X = set_par_invgamma(media = 0.01, var = 0.0001)
par_sig2A = set_par_invgamma(media = 1,   var = 1) # I can not learn this
par_gamma = set_par_gamma(media = 0.4,    var = 2)
par_c     = set_par_gamma(media = 1.1,    var = 1)

# Plot prior sig2X
dens_prior = density( 1/rgamma(n = 10000, shape = par_sig2X[1], rate = par_sig2X[2]) )
par(mfrow = c(1,1), mar = c(2,3,2,1), bty = "l")
plot(x = dens_prior$x, y = dens_prior$y, type = "l", lwd = 3, xlim = c(0,0.1))

# Check M1 prior
Nprior = 10000
x = matrix(0,nrow = Nprior, ncol = 3)
x[,1] = rgamma(n = Nprior, shape = par_gamma[1], rate = par_gamma[2])
x[,2] = rgamma(n = Nprior, shape = par_c[1], rate = par_c[2])
x[,3] = rep(0,Nprior) 

M1Prior = apply(x,1,computeM1)
range(M1Prior); mean(M1Prior)
round(quantile(M1Prior, probs = c(0.025,0.5,0.995)),4)

PsurvPrior = apply(x,1,computePsurv)
range(PsurvPrior); mean(PsurvPrior)
round(quantile(PsurvPrior, probs = c(0.025,0.5,0.995)),4)

# Gibbs sampler structures
BetaProcess_params = matrix(NA, nrow = G, ncol = 3)
colnames(BetaProcess_params) = c("gamma","c","sigma")
sigma2_XA_mcmc = matrix(NA, nrow = G, ncol = 2)
colnames(sigma2_XA_mcmc) = c("sigma2_X","sigma2_A")
var_prop <- M1_mcmc <- Psurv_mcmc <- rep(0,G)

Particles_PG = vector("list", G)
Ktot = rep(NA,G)

# Initialization
var_prop[1] = 0.1 # initial adaptive variance
BetaProcess_params[1,1] = 1 # gamma = 1 at first iteration
BetaProcess_params[1,2] = 1 # c = 1 at first iteration
BetaProcess_params[,3]  = 0 # set sigma = 0 for all iterations
# BetaProcess_params[1,1:2] = c(gamma,c) # # set true values (if any)


# sigma2_XA_mcmc[1,] = c(0.05,0.5) # initial values for sigma2_X and sigma2_A
sigma2_XA_mcmc[1,] = c(sigma2_X,0.5) # set true values

M1_mcmc[1] = computeM1(BetaProcess_params[1,])
Psurv_mcmc[1] = computePsurv(BetaProcess_params[1,])

proposal_Nnew_1T = matrix(0,nrow = Ttot, ncol = G)
# PG init -----------------------------------------------------------------


Particles_PG[[1]] = Conditional_SeqMonteCarlo( X = X, N = N, D = D, Ttot = Ttot,
                                               Bfix = rep(-1,Ttot), Particle_fix = NULL,  
                                               M1 = M1_mcmc[1],Psurv = Psurv_mcmc[1],
                                               sigma2_A = sigma2_XA_mcmc[1,2], 
                                               sigma2_X = sigma2_XA_mcmc[1,1],
                                               proposal_Nnew_1T = NULL,
                                               use_VS = use_VS)

temp = sapply(Particles_PG[[1]]$Path_k[1:Ttot],function(x) x$Nnew)
proposal_Nnew_1T[,1] = temp
temp

# Survived
sapply(Particles_PG[[1]]$Path_k[1:Ttot], function(x) max(x$survived) )


dim(Particles_PG[[1]]$mean_k)
dim(Particles_PG[[1]]$Amat_k)
dim(Particles_PG[[1]]$Zmat_k)

media = Particles_PG[[1]]$mean_k
plot_multiple_imgs(media, plt_wnd = c(3,4))

Amat = Particles_PG[[1]]$Amat_k
plot_multiple_imgs(Amat, plt_wnd = c(3,4))


# Amat          = A # set true values
# mean_vec_1T   = (Z%*%A) # set true values
# K             = 4 # set true values
# Nnew_1T       = c(1,1,0,1,0,1,0,0,0,0,0,0)
# Zmat          = Z # set true values




# Run ---------------------------------------------------------------------
g = 2
UpdateSig2X = FALSE;UpdateSig2A = FALSE;
Updatec = TRUE; Updategamma = TRUE
UpdatePG = TRUE
pb = txtProgressBar(min = 2, max = G, initial = 2, style = 3) 
for(g in 2:G){
  
  # cat("\n ++++++ g = ",g," ++++++ \n")
  setTxtProgressBar(pb,g)
  if(UpdatePG){
    mean_vec_1T = Particles_PG[[g-1]]$mean_k
    Zmat        = Particles_PG[[g-1]]$Zmat_k
    Amat        = Particles_PG[[g-1]]$Amat_k
    Ktot[g-1]   = nrow(Particles_PG[[g-1]]$Amat_k)
  } else{
    Amat          = A # set true values
    mean_vec_1T   = (Z%*%A) # set true values
    K             = 4 # set true values
    Nnew_1T       = c(1,1,0,1,0,1,0,0,0,0,0,0)
    Zmat          = Z # set true values
    Ktot[g-1]     = nrow(Amat)
  }
  
  # Update sigma2_X
  if(UpdateSig2X){
    suff_stat_X = sum( apply(X-mean_vec_1T, 1, function(x){t(x)%*%x}) )
    sigma2_XA_mcmc[g,1] = r_fullcond_sigma2_X(par_sig2X[1], par_sig2X[2], D, Ttot, suff_stat_X)
    # sigma2_XA_mcmc[g,1] = sigma2_X # set true values
    # cat("\n sig2X = ",sigma2_XA_mcmc[g,1],"\n")
  }else{
    sigma2_XA_mcmc[g,1] = sigma2_XA_mcmc[g-1,1]
  }
  
  # Update sigma2_A
  if(UpdateSig2A){
    stop("Sei sicuro di voler aggiornare sig2A?")
    suff_stat_A = sum( apply(Amat, 1, function(x){t(x)%*%x}) )
    sigma2_XA_mcmc[g,2] = r_fullcond_sigma2_A(par_sig2A[1], par_sig2A[2], D, Ktot[g-1], suff_stat_A)
    # sigma2_XA_mcmc[g,2] = sigma2_A # set true values
  }else{
    sigma2_XA_mcmc[g,2] = sigma2_XA_mcmc[g-1,2]
  }
  
  # Update_gamma
  if(Updategamma){
    BetaProcess_params[g,1] = r_fullcond_gamma(par_gamma[1], par_gamma[2], 
                                               Ktot[g-1], Ttot, 
                                               BetaProcess_params[g-1,2], BetaProcess_params[g-1,3])
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
  
  # Run Conditional SMC
  if(UpdatePG){
    proposal_N = matrix(proposal_Nnew_1T[,1:g-1], nrow = Ttot, ncol = (g-1) )
    Particles_PG[[g]] = Conditional_SeqMonteCarlo( X = X, N = N, D = D, Ttot = Ttot,
                                                   Bfix = Particles_PG[[g-1]]$B_k,
                                                   Particle_fix = Particles_PG[[g-1]]$Path_k,  
                                                   M1 = M1_mcmc[g],Psurv = Psurv_mcmc[g],
                                                   sigma2_A = sigma2_XA_mcmc[g,2], 
                                                   sigma2_X = sigma2_XA_mcmc[g,1],
                                                   proposal_Nnew_1T = proposal_N,
                                                   use_VS = use_VS ) 
    # Update proposal_Nnew_1T:
    temp = sapply(Particles_PG[[g]]$Path_k[1:Ttot],function(x) x$Nnew)
    proposal_Nnew_1T[,g] = temp
  }
  
  
}
close(pb)
beepr::beep()
# stop("Fine sampler, fermo io")



# Analysis PG ----------------------------------------------------------------

K_mcmc = sapply(Particles_PG, function(x){ncol(x$Zmat_k)})
par(mfrow = c(1,3), bty = "l")
plot(M1_mcmc, type = "l", main = "M1")
plot(Psurv_mcmc, type = "l", main = "Prob. surv.");
plot(K_mcmc, type = "l", main = "Num. features")


Kt_mcmc = t(sapply(Particles_PG, function(x){apply(x$Zmat_k,1,sum)}))
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  plot(Kt_mcmc[,t], type = "l", main = paste0("K_=",t))
}

# # Survived
# par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
# for(t in 1:Ttot){
#   Nsurv_t_mcmc = sapply(Particles_PG, function(x) max(x$Path_k[[t]]$survived) )
#   plot(Nsurv_t_mcmc, type = "l", main = paste0("Surv. T=",t))
# }

# New 
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  Nnew_t_mcmc = sapply(Particles_PG, function(x) x$Path_k[[t]]$Nnew )
  plot(Nnew_t_mcmc, type = "l", main = paste0("T=",t))
}

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



# Analysis sig2X ----------------------------------------------------------

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





