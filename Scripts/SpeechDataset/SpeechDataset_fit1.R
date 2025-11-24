wd = "C:/Users/colom/DynamicFeatureAllocation/Scripts/SpeechDataset"
# wd = "/home/lucia.paci/Lucia/Ale/DynamicFeatureAllocation/Scripts/SpeechDataset"
setwd(wd)
source("./../../R/Rfunctions.R")
Rcpp::sourceCpp("./../../src/RcppFunctions.cpp")
mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)
library(fields)


# Read data -----------------------------------------------------------

seed = 132332
set.seed(seed)
r = 20

# Read President names
Presidents_all <- read.csv("C:/Users/colom/DynamicFeatureAllocation/Scripts/data/Presidents_all.csv")
Presidents_all[,2]
data = read.table(paste0("../data/SpeechData_top",r,".txt"))
colnames(data) = Presidents_all[,2]

# Reduced sample
data = as.matrix(data)
# data = data[,41:60]

# Dimensions
V = nrow(data)
Ttot = ncol(data)


# -> Lambda is VxK
# -> Xi is KxT
# -> Mean is VxT
# -> D is VxT


# Visualize data ----------------------------------------------------------


# D
par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:V, 
       t(data),   
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
  1:Ttot, 1:V, t(data),
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,                # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)



# Set options -------------------------------------------------------------
seed = 22123

H = 30 # number of atoms
delta = 0.01 # Dirichlet parameter

# Hyperparameters
a_phi=1;b_phi=1;
a_gamma=1;b_gamma=1;
a_sigma=1;b_sigma=1;
a_beta=1;b_beta=1; 
var_phi<-var_gamma<-var_sigma<-var_beta<-0.01



# Initial values
Xi0 = matrix(sample(0:100,H*Ttot,TRUE), nrow = H, ncol = Ttot)
S0  = matrix(rgamma(n=H*Ttot,1,1), nrow = H, ncol = Ttot)
Lambda0 = vector("list",H)
Lambda0 = lapply(Lambda0, function(x){
  A = matrix(0, nrow = V, ncol = Ttot)
  A = apply(A,2,function(y){ a = rgamma(n=V, shape = delta, rate = 1); a/sum(a) })
  A
})



phi0=1;gamma0=1;sigma0=0.5;beta0=1;


init_DTM = set_init_DTM(Xi0,Lambda0,S0,phi0,gamma0,sigma0,beta0)

UpdateDitl=TRUE; UpdateS=TRUE; 
UpdateLambda=TRUE;UpdateXi=TRUE; 
UpdateU=TRUE;
UpdatePhi=FALSE; UpdateGamma=FALSE;
UpdateSigma=FALSE; UpdateBeta = FALSE; 
print = TRUE; JointAdp = FALSE
param_DTM = set_param_DTM(H,delta, 
                          a_phi,b_phi,a_gamma,b_gamma,a_sigma,b_sigma,a_beta,b_beta, 
                          var_phi,var_gamma,var_sigma,var_beta,
                          UpdateDitl,UpdateS,UpdateLambda,UpdateXi, UpdateU,
                          UpdatePhi,UpdateGamma,UpdateSigma,UpdateBeta,seed,
                          print,JointAdp)

# Run ---------------------------------------------------------------------
niter = 5000
nburn =    1
data  = as.matrix(data)
fit = GibbsSampler_DTM(niter,nburn,data,param_DTM,init_DTM)


# Diagnosis ---------------------------------------------------------------

it_start = 0 + nburn + 1     # da dove parto + burnin + val iniziale
it_end   = niter + nburn + 1 # ultima iterazione

## Xi ------------------------------------------------------------------

Xi_mean <- Reduce("+", fit$Xi[it_start:it_end])/length(fit$Xi[it_start:it_end])
par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:H, 
       t(Xi_mean),   
       col = mycol,    
       xlab = "Time", 
       ylab = "Atoms",
       main = "<Xi>",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, H, length.out = min(H, 10)), 
     labels = round(seq(1, H, length.out = min(H, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:H, Xi_mean,
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,              # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)


## S  ----------------------------------------------------------------------

S_mean <- Reduce("+", fit$S[it_start:it_end])/length(fit$S[it_start:it_end])
par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:H, 
       t(S_mean),   
       col = mycol,    
       xlab = "Time", 
       ylab = "Atoms",
       main = "<S>",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, H, length.out = min(H, 10)), 
     labels = round(seq(1, H, length.out = min(H, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:H, S_mean,
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,                # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)


## N_tl ------------------------------------------------------------------

N_mean <- Reduce("+", fit$N[it_start:it_end])/length(fit$N[it_start:it_end])
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



## D_l --------------------------------------------------------------------

Dl_means = lapply(1:H, function(l){
  temp = lapply(fit$Dl, function(x){x[[l]]})
  Reduce("+", temp[it_start:it_end])/length(temp[it_start:it_end])
})

for(l in 1:H){
  Dl_mean_l = Dl_means[[l]]
  par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
  image( 1:Ttot, 1:V, 
         t(Dl_mean_l),   
         col = mycol,    
         xlab = "Time", 
         ylab = "Words",
         main = paste0("<D_",l,">"),
         axes = FALSE )
  axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(2, at = seq(1, V, length.out = min(H, 10)), 
       labels = round(seq(1, V, length.out = min(H, 10))),
       cex.axis = 0.7)
  box()
  fields::image.plot(
    1:Ttot, 1:V, t(Dl_mean_l),
    col = mycol,
    legend.only = TRUE,
    horizontal = FALSE,
    legend.width = 1.2,            # controls legend thickness
    legend.shrink = 0.8,           # smaller legend
    legend.mar = 8.5,                # margin from image
    legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
  )
}




## Lambda --------------------------------------------------------------------

Lambda_means = lapply(1:H, function(l){
  temp = lapply(fit$Lambda, function(x){x[[l]]})
  Reduce("+", temp[it_start:it_end])/length(temp[it_start:it_end])
})


for(l in 1:H){
  Lambda_means_l = Lambda_means[[l]]
  par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
  image( 1:Ttot, 1:V, 
         t(Lambda_means_l),   
         col = mycol,    
         xlab = "Time", 
         ylab = "Words",
         main = paste0("<Lambda_",l,">"),
         axes = FALSE )
  axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(2, at = seq(1, V, length.out = min(V, 10)), 
       labels = round(seq(1, H, length.out = min(V, 10))),
       cex.axis = 0.7)
  box()
  fields::image.plot(
    1:Ttot, 1:V, Lambda_means_l,
    col = mycol,
    legend.only = TRUE,
    horizontal = FALSE,
    legend.width = 1.2,            # controls legend thickness
    legend.shrink = 0.8,           # smaller legend
    legend.mar = 8.5,                # margin from image
    legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
  )
  
}




## Hyperparameters ---------------------------------------------------------

par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( fit$phi, xlab = "Iter.", ylab = "phi", type = "l" )

par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( fit$gamma, xlab = "Iter.", ylab = "gamma", type = "l" )

par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( fit$sigma, xlab = "Iter.", ylab = "sigma", type = "l" )

par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( fit$beta, xlab = "Iter.", ylab = "beta", type = "l" )

par(mfrow = c(1,2), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( fit$t_sigma_gamma, xlab = "Iter.", ylab = "t", type = "l" )
plot( log(fit$t_sigma_gamma), xlab = "Iter.", ylab = "log(t)", type = "l" )


# Brutta ------------------------------------------------------------------

it_last = 1002
Xi_last = fit$Xi[[it_last]]
Lambda_last = fit$Lambda[[it_last]]

Xi_last[,1]
active = which(Xi_last[,1] > 0)
Lambda_active = Lambda_last[active]
Lambda_active[[1]]

length(fit)










