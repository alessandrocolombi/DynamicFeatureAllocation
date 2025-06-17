# Introduction ------------------------------------------------------------

setwd("C:/Users/colom/DynamicFeatureAllocation/R/Scripts")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../Rfunctions.R")
source("../genera_img_ilaria.R")

# Custom functions --------------------------------------------------------

computeM1 = function(x){
  x[1] * gamma(x[3]+1) * gamma(x[2]-x[3]+1) / gamma(x[2]+2)
}
computePsurv = function(x){
  (1-x[3])/(1+x[2])
}
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

# Sim. from model 2 ---------------------------------------------------------

Ttot = 24    # Time length
di   = 8
D    = di*di # Size of the problem

M = 4
zetas = rbind(A1,A2,A3,A4)
Psurv_true = 0.85
M1true = 1/M
c = 1/Psurv_true - 1
gamma = M1true * Psurv_true

sigma2_ker = 0.01;sigma2_A = 0.01;sigma2_X = 0.01



## Generate BetaBernoulli with M centers
base_list = list("active_values" = c(),"Nnew" = c(), "survived" = c())
Process   = lapply(1:M, function(x){lapply(1:Ttot, function(xx){base_list})})

# Time 1
t = 1
for(h in 1:M){
  Nnew   = rpois(n = 1, lambda = M1true)
  values = matrix(0,nrow = Nnew, ncol = D)
  for(j in 1:D){
    values[,j] = rnorm(n=Nnew, mean = zetas[h,j], sd = sqrt( sigma2_ker ) ) 
  }
  Process[[h]][[t]]$active_values   = values
  Process[[h]][[t]]$Nnew            = Nnew
  if(Nnew > 0)
    Process[[h]][[t]]$survived        = seq(1,Nnew)
}

# Time t (2...Ttot)
for(t in 2:Ttot){
  
  for(h in 1:M){
    # Update center h at time t
      
    # Thinning part:
    num_active_old = nrow(Process[[h]][[t-1]]$active_values)
    if( num_active_old > 0  ){
      survived = sample( c(0,1), size = num_active_old, 
                         replace = TRUE, prob = c(1-Psurv_true,Psurv_true) )
      n_surv = sum(survived)
    }else{
      survived = c()
      n_surv = 0
    }
      
    # Innovation:
    Nnew   = rpois(n = 1, lambda = M1true)
    values = matrix(0,nrow = Nnew, ncol = D)
    for(jj in 1:D){
      values[,jj] = rnorm(n=Nnew, mean = zetas[h,jj], sd = sqrt( sigma2_ker ) )
    }
      
    # Assemble thinning + innovation:
    if(n_surv){
      survived_values = Process[[h]][[t-1]]$active_values[which(survived > 0),]
      Process[[h]][[t]]$active_values = rbind( survived_values, values )
    }else{
      Process[[h]][[t]]$active_values =  values 
    }
      
    Process[[h]][[t]]$Nnew     = Nnew
    Process[[h]][[t]]$survived = which(survived > 0)
      
  }
    
}


# Find all features appeared at least once
A_m = lapply(1:M, function(x){c()})
cum_Ntot = lapply(1:M,function(x){c()})
Nnew_paths = lapply(1:M,function(x){c()})
for(h in 1:M){
  for(t in Ttot:1){
    Mat = Process[[ h ]][[t]]$active_values
    nuove = Process[[ h ]][[t]]$Nnew
    Nnew_paths[[h]] = c(nuove, Nnew_paths[[h]])
    cum_Ntot[[h]] = c(nuove, cum_Ntot[[h]])
    if(nuove > 0){
      if(nrow(Mat)-nuove+1 <= 0)
        stop("Unexprected behaviour")
      A_m[[h]] = rbind( Mat[(nrow(Mat)-nuove+1):nrow(Mat),] , A_m[[h]])
    }
  }
  cum_Ntot[[h]] = cumsum(cum_Ntot[[h]])
}

# Find survival path of each feature for each center
Z_m = vector("list",M)
for(h in 1:M){
  Z_m[[h]] = matrix(0, nrow = Ttot, ncol = nrow(A_m[[h]]))
  # Time 1:
  t = 1
  if(Nnew_paths[[h]][t] > 0)
    Z_m[[h]][t, 1:cum_Ntot[[h]][1] ] = 1
  for(t in 2:Ttot){
    # Time t (2...Ttot)
    if(Nnew_paths[[h]][t] > 0)
      Z_m[[h]][t, (cum_Ntot[[h]][t-1]+1):(cum_Ntot[[h]][t])] = 1
    Nsurvived_old = length(which(Z_m[[h]][t-1,] == 1))
    submat = matrix( Z_m[[h]][,which(Z_m[[h]][t-1,] == 1)],
                     nrow = Ttot, ncol = Nsurvived_old)
    submat[t, Process[[ h ]][[t]]$survived ] = 1
    Z_m[[h]][,which(Z_m[[h]][t-1,] == 1)] = submat
  }
}


mean_all = lapply(1:M, function(h){Z_m[[h]]%*%A_m[[h]]})
mean_all_sum = Reduce(`+`, mean_all)
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  mean_t_vec = mean_all_sum[t,]
  mean_t = matrix(data = mean_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img(mean_t, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = paste0("T = ",t))
}










# Sim. data from model 2  -------------------------------------------------

seed = 42
set.seed(seed)

Ttot = 24    # Time length
di   = 8
D    = di*di # Size of the problem

M = 4
zetas = rbind(A1,A2,A3,A4)
Psurv_true = 0.9
M1true = 1/M
c = 1/Psurv_true - 1
gamma = M1true * Psurv_true

sigma2_ker = 0.01;sigma2_A = 0.01;sigma2_X = 0.01


sim_data = sim_images_ale_mod2( zetas = zetas, Ti = Ttot, 
                                sig2X = sigma2_X, sig2A = sigma2_A,
                                Psurv = Psurv_true, M1 = M1true)

A = sim_data$A
Z = sim_data$Z
X = sim_data$X

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
  plot_img(X_t, center_value = 0, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = paste0("Obs. T = ",t))
}
