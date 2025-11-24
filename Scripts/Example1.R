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
seed = 22123

H = 10 # number of atoms
delta = 0.1 # Dirichlet parameter

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
# ---
# Set Lambda0 to the true one
for(t in 1:6){
  Lambda0[[1]][,t] = Lambda[,1]
}
for(t in 7:14){
  Lambda0[[2]][,t] = Lambda[,2]
}
for(t in 15:23){
  Lambda0[[3]][,t] = Lambda[,3]
}
# ---
# Set Xi0 to the true one
Xi0 = matrix(0, nrow = H, ncol = Ttot)
Xi0[1:3,] = Xi
# ---
# Set S0 to the true one
# load("Brutta_Strue.Rdat")
# S0 = S_mean
# ---


phi0=1;gamma0=1;sigma0=0.25;beta0=0.5;


init_DTM = set_init_DTM(Xi0,Lambda0,S0,phi0,gamma0,sigma0,beta0)

UpdateDitl=TRUE; UpdateS=TRUE; 
UpdateLambda=TRUE;UpdateXi=TRUE; 
UpdateU=TRUE;
UpdatePhi=FALSE; UpdateGamma=FALSE;
UpdateSigma=FALSE; UpdateBeta = FALSE; 
print = TRUE
JointAdp = FALSE
param_DTM = set_param_DTM(H,delta, 
                          a_phi,b_phi,a_gamma,b_gamma,a_sigma,b_sigma,a_beta,b_beta, 
                          var_phi,var_gamma,var_sigma,var_beta,
                          UpdateDitl,UpdateS,UpdateLambda,UpdateXi, UpdateU,
                          UpdatePhi,UpdateGamma,UpdateSigma,UpdateBeta,seed,
                          print,JointAdp)

# Run ---------------------------------------------------------------------
niter = 5000
nburn =    1

# sink("log.txt")
fit = GibbsSampler_DTM(niter,nburn,D,param_DTM,init_DTM)
# sink()

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

# Traceplot
l = 1; t = 1
x = sapply(1:length(fit$Xi), function(i) fit$Xi[[i]][l,t])
par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( x, xlab = "Iter.", ylab = paste0("Xi_",l,",",t), type = "l", ylim = c(0,max(x)) )


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

# Traceplot
l = 1; t = 1
x = sapply(1:length(fit$S), function(i) fit$S[[i]][l,t])
par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( x, xlab = "Iter.", ylab = paste0("S_",l,",",t), type = "l" )


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




## Mean values --------------------------------------------------------------------
meanRes <- matrix(0, nrow = Ttot, ncol = V)
for(it in 1:length(fit$Xi)) {
  
  Xi_it <- fit$Xi[[it]]          # H × Ttot
  Lambdas_it <- fit$Lambda[[it]] # list of H matrices (each V × Ttot)
  
  # Build Lambda_arr: H × V × Ttot (fully vectorized using simplify2array)
  Lambda_arr <- simplify2array(Lambdas_it)  # gives V × Ttot × H
  Lambda_arr <- aperm(Lambda_arr, c(3,1,2)) # reorder → H × V × Ttot
  
  # (2) Expand Xi: from H×Ttot into H×V×Ttot
  Xi_expanded <- array(Xi_it, dim = c(H, 1, Ttot))
  Xi_expanded <- Xi_expanded[, rep(1, V), , drop = FALSE]   # replicate across V
  # Now Xi_expanded is H × V × Ttot
  
  # Elementwise multiply and sum over H
  # result: V × Ttot after summing over axis=1, then transpose → Ttot × V
  res_it <- apply(Xi_expanded * Lambda_arr, c(2,3), sum)
  res_it <- t(res_it)   # now Ttot × V
  
  meanRes <- meanRes + res_it
}
meanRes <- meanRes / niter
dim(meanRes)

par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:V, 
       meanRes,   
       col = mycol,    
       xlab = "Time", 
       ylab = "Words",
       main = paste0("< Estimated Mean >"),
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, V, length.out = min(H, 10)), 
     labels = round(seq(1, V, length.out = min(H, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:V, meanRes,
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,                # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)
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

## U ------------------------------------------------------------------

U_mean <- Reduce("+", fit$U[it_start:it_end])/length(fit$U[it_start:it_end] )
par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:H, 
       t(U_mean),   
       col = mycol,    
       xlab = "Time", 
       ylab = "Atoms",
       main = "<U>",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, H, length.out = min(H, 10)), 
     labels = round(seq(1, H, length.out = min(H, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:H, U_mean,
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,              # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)

# Traceplot
l = 1; t = 1
x = sapply(1:length(fit$U), function(i) fit$U[[i]][l,t])
par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( x, xlab = "Iter.", ylab = paste0("U_",l,",",t), type = "l" )

## Traceplot Num. topics ------------------------------------------------------------------

# Num. topics at each time 
Kt_tr = sapply(fit$Xi, function(Xi_it) apply(Xi_it, 2, function(x) length(which(x > 0))) )
Kt_tr = t(Kt_tr)

par(mfrow = c(3,3), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
for(t in 1:Ttot){
  plot( Kt_tr[,t], xlab = "Iter.", ylab = paste0("K",t), type = "l" )
}

# Total number of distinct topics
K_tr = sapply(fit$Xi, function(Xi_it){
  idx_list = lapply(1:H, function(l) find_indices(Xi_it[l,]))
  total <- sum(vapply(idx_list, function(x) length(x[[1]]), integer(1)))
  total
} )

par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( K_tr, xlab = "Iter.", ylab = paste0("K"), type = "l" )


## Unique values for a given iteration ------------------------------------------------------------------
it = 4500
Xi_it = fit$Xi[[it]]
idx_list = lapply(1:H, function(l) find_indices(Xi_it[l,]))
for(l in 1:H){
  Num_act_l = length(idx_list[[l]]$idx_born) # how many distinct values
  
  if(Num_act_l > 0){
    for(jj in seq_along(idx_list[[l]]$idx_born) ){
      
      # Find unique value
      tj = idx_list[[l]]$idx_born[jj]
      lambda = fit$Lambda[[it]][[l]][,tj]
      
      # Find activity times
      tj_next = idx_list[[l]]$idx_born[jj+1]
      if(is.na(tj_next))
        tj_next = Inf
      
      time_act = c(tj,
                   idx_list[[l]]$idx_surv[which(idx_list[[l]]$idx_surv < tj_next)]
                   )
      active <- Xi_active <- rep(0,Ttot)
      active[ time_act ] = 1
      Xi_active[ time_act ] = fit$Xi[[it]][l,time_act]
      par(mfrow = c(1,2), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
      barplot( height = lambda, 
               names.arg = as.character(1:V),
               las = 1, col = "darkred", border = NA,
               xlab = "Word",
               main = paste0("Atom: ",l,",",jj), ylab = "Prob.", ylim = c(0,ymax),
               cex.names = 0.5 )
      barplot( height = Xi_active, 
               names.arg = as.character(1:Ttot),
               las = 1, col = "darkgreen", border = NA,
               xlab = "Time",
               main = "", ylab = "Intensity", ylim = c(0,max(Xi_active)),
               cex.names = 0.5 )
      # barplot( height = active, 
      #          names.arg = as.character(1:Ttot),
      #          las = 1, col = "black", border = NA,
      #          xlab = "Time",
      #          main = "", ylab = "", ylim = c(0,1.1),
      #          cex.names = 0.5 )
    }

  }

}

# Brutta ------------------------------------------------------------------
par(mfrow = c(1,3), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( fit$t_sigma_gamma[6000:6966], xlab = "Iter.", ylab = "t", type = "l" )
plot( fit$gamma[6000:6966], xlab = "Iter.", ylab = "gamma", type = "l" )
plot( fit$sigma[6000:6966], xlab = "Iter.", ylab = "sigma", type = "l" )

set_par_gamma(1,0.1)

it = 4280
fit$gamma[it]
fit$beta[it]
fit$sigma[it]
fit$t_sigma_gamma[it]

1/fit$sigma[it]

(fit$sigma[it]*10/fit$gamma[it])^20

msigma = mean(fit$sigma)
msigma
mgamma = mean(fit$gamma)
mgamma
mbeta = mean(fit$beta)
mbeta
mt = mean(fit$t_sigma_gamma)
mt










it = 388
t = 1
Xi_it = fit$Xi[[it]]
Lambdas_it = fit$Lambda[[it]]
Lambda_tl = vapply(Lambdas_it, function(x) x[,t], numeric(V) )
Lambda_tl = t(Lambda_tl)
sum <- as.numeric( t(Xi_it[, t]) %*% Lambda_tl )
round(sum,2)


res <- array(0, dim = c(length(fit$Xi), Ttot, V))

for(it in 1:niter) {
  Xi_it      <- fit$Xi[[it]]         # H × Ttot
  Lambdas_it <- fit$Lambda[[it]]     # list of H matrices (each V × Ttot)
  
  # Build Lambda_tl: H × V × Ttot
  # For each l=1..H: extract all columns at once
  Lambda_arr <- array(0, dim = c(H, V, Ttot))
  for(l in 1:H) {
    Lambda_arr[l, , ] <- Lambdas_it[[l]]  # V × Ttot
  }
  
  # For each t compute Xi_it[,t]ᵀ %*% Lambda_tl (vector length V)
  for(t in 1:Ttot) {
    Lambda_tl <- Lambda_arr[ , , t]          # H × V
    res[it, t, ] <- as.numeric( t(Xi_it[, t]) %*% Lambda_tl )
  }
}



