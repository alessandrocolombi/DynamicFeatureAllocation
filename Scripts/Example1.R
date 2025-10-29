wd = "C:/Users/colom/DynamicFeatureAllocation/Scripts/"
# wd = "/home/lucia.paci/Lucia/Ale/DynamicFeatureAllocation/Scripts/"
setwd(wd)
source("./../R/Rfunctions.R")
Rcpp::sourceCpp("./../src/RcppFunctions.cpp")
mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)
library(fields)


# Generate data -----------------------------------------------------------

seed = 132332
set.seed(seed)

Ttot = 23
Ktrue = 3
V = 30

lambda_1 <- lambda_2 <- lambda_3 <- rep(0,V)
lambda_1[1:10]  =  1/length(1:10)
lambda_2[11:20] =  1/length(11:20)
lambda_3[21:30] =  1/length(21:30)

Lambda = cbind(lambda_1,lambda_2,lambda_3)

a = 0.5
Xi = matrix(0,nrow = Ktrue, ncol = Ttot)
## Zipfs decay
Xi[1,1:6]   = 400 * zipfs_decay(6,a)
Xi[2,7:14]  = 500 * zipfs_decay(8,a)
Xi[3,15:23] = 100 * zipfs_decay(9,-a)


Mean = Lambda %*% Xi

D <- matrix(rpois(V * Ttot, lambda = as.vector(Mean)), nrow = V, ncol = Ttot)


# -> Lambda is VxK
# -> Xi is KxT
# -> Mean is VxT
# -> D is VxT


# Visualize data ----------------------------------------------------------


ymax = 0.75
ymax2 = max(Xi)
par(mfrow = c(1,2), mar = c(3.5,3.5,2,1), mgp=c(2,0.5,0))
barplot( height = Lambda[,1], 
         names.arg = as.character(1:30),
         las = 1, col = "darkred", border = NA,
         xlab = "Word",
         main = "Lambda_1", ylab = "Prob.", ylim = c(0,ymax),
         cex.names = 0.5 )
barplot( height = Xi[1,], 
         names.arg = as.character(1:Ttot),
         las = 1, col = "darkred", border = NA,
         xlab = "Time",
         main = "Xi_1", ylab = "Intensity", ylim = c(0,ymax2),
         cex.names = 0.5 )
barplot( height = Lambda[,2], 
         names.arg = as.character(1:30),
         las = 1, col = "darkblue", border = NA,
         xlab = "Word",
         main = "Lambda_2", ylab = "Prob.", ylim = c(0,ymax),
         cex.names = 0.5)
barplot( height = Xi[2,], 
         names.arg = as.character(1:Ttot),
         las = 1, col = "darkblue", border = NA,
         xlab = "Word",
         main = "Xi_2", ylab = "Intensity", ylim = c(0,ymax2),
         cex.names = 0.5)
barplot( height = Lambda[,3], 
         names.arg = as.character(1:30),
         las = 1, col = "darkgreen", border = NA,
         xlab = "Word",
         main = "Lambda_3", ylab = "Prob.", ylim = c(0,ymax),
         cex.names = 0.5)
barplot( height = Xi[3,], 
         names.arg = as.character(1:Ttot),
         las = 1, col = "darkgreen", border = NA,
         xlab = "Word",
         main = "Xi_3", ylab = "Intensity", ylim = c(0,ymax2),
         cex.names = 0.5)

# Plot mean and D
par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:V, 
       t(Mean),   
       col = mycol,    
       xlab = "Time", 
       ylab = "Words",
       main = "Mean",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, V, length.out = min(V, 10)), 
     labels = round(seq(1, V, length.out = min(V, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:V, t(Mean),
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,                # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)

# D
par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:V, 
       t(D),   
       col = mycol,    
       xlab = "Time", 
       ylab = "Words",
       main = "D",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, V, length.out = min(V, 10)), 
     labels = round(seq(1, V, length.out = min(V, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:V, t(D),
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,                # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)



# Set options -------------------------------------------------------------
seed = 1234

H = 10 # number of atoms
delta = 0.01 # Dirichlet parameter

# Hyperparameters
a_phi=1;b_phi=1;
a_gamma=1;b_gamma=1;
a_sigma=1;b_sigma=1;
a_beta=1;b_beta=1; 
var_phi<-var_gamma<-var_sigma<-var_beta<-0.01

UpdateDitl=TRUE; UpdateS=TRUE; 
UpdateLambda=TRUE;UpdateXi=TRUE; 
UpdateU=TRUE;
UpdatePhi=TRUE; UpdateGamma=TRUE;
UpdateSigma=TRUE; UpdateBeta = TRUE; 
print = TRUE
param_DTM = set_param_DTM(H,delta, 
                          a_phi,b_phi,a_gamma,b_gamma,a_sigma,b_sigma,a_beta,b_beta, 
                          var_phi,var_gamma,var_sigma,var_beta,
                          UpdateDitl,UpdateS,UpdateLambda,UpdateXi, UpdateU,
                          UpdatePhi,UpdateGamma,UpdateSigma,UpdateBeta,seed,print)

# Initial values
Xi0 = matrix(sample(0:100,H*Ttot,TRUE), nrow = H, ncol = Ttot)
S0  = matrix(rgamma(n=H*Ttot,1,1), nrow = H, ncol = Ttot)
Lambda0 = vector("list",H)
Lambda0 = lapply(Lambda0, function(x){
  A = matrix(0, nrow = V, ncol = Ttot)
  A = apply(A,2,function(y){ a = rgamma(n=V, shape = delta, rate = 1); a/sum(a) })
})
phi0=1;gamma0=1;sigma0=0.5;beta0=1;


init_DTM = set_init_DTM(Xi0,Lambda0,S0,phi0,gamma0,sigma0,beta0)


# Run ---------------------------------------------------------------------
niter = 10000
nburn = 1
fit = GibbsSampler_DTM(niter,nburn,D,param_DTM,init_DTM)


# Diagnosis ---------------------------------------------------------------

phi_mcmc = fit$phi
plot(phi_mcmc, type = "l")
plot(fit$gamma, type = "l")

N_mean <- Reduce("+", fit$N)/length(fit$N)
dim(N_mean)
par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:H, 
       N_mean,   
       col = mycol,    
       xlab = "Time", 
       ylab = "Atoms",
       main = "<N_tl>",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, H, length.out = min(H, 10)), 
     labels = round(seq(1, H, length.out = min(H, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:H, N_mean,
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,                # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)


length(fit$Lambda)
length(fit$Lambda[[1]])


it = 10000
# Dl_final = Reduce("+",fit$Dl[[it]])
Dl_final_1 = fit$Dl[[it]][[1]]
par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:V, 
       t(Dl_final_1),   
       col = mycol,    
       xlab = "Time", 
       ylab = "Words",
       main = "Dl_final",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, V, length.out = min(H, 10)), 
     labels = round(seq(1, V, length.out = min(H, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:V, t(Dl_final),
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,                # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)


fit$Lambda[[it]][[1]]











