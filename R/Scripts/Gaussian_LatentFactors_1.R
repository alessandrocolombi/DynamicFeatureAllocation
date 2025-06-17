
# Introduction ------------------------------------------------------------

setwd("C:/Users/colom/DynamicFeatureAllocation/R/Scripts")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../Rfunctions.R")
source("../genera_img_ilaria.R")

# Custom functions --------------------------------------------------------
# Immagini Base -----------------------------------------------------------

di   = 8
D    = di*di # Size of the problem

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

# Data generation ---------------------------------------------------------
seed = 42
set.seed(seed)

Ttot = 24    # Time length
di   = 8
D    = di*di # Size of the problem
Psurv = 0.9

sigma2_X = 0.01
sigma2_A = 0.01

sim_data = sim_images_ale(img_base = img_base, Ti = Ttot, 
                          sig2X = sigma2_X, sig2A = sigma2_A,
                          Psurv = Psurv)

A = sim_data$A
Z = sim_data$Z
X = sim_data$X

Nfeat_tot = nrow(A)
par(mfrow = c(2,3), mar = c(2,2,2,1), bty = "l")
for(j in 1:Nfeat_tot){
  Mat_vec = A[j,]
  Mat = matrix(data = Mat_vec, nrow = di, ncol = di, byrow = F)
  plot_img(Mat, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = paste0("Feat. ",j,"th"))
}

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
  plot_img(X_t, 
           center_value = 0, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = paste0("Obs. T = ",t))
}


Psurv_true = 0.75
M1true = 0.75
c = 1/Psurv_true - 1
gamma = M1true * Psurv_true
sigma = 0
# c     = 1/0.1 - 1
# gamma = 1




# Set options -------------------------------------------------------------

seed = 42
set.seed(seed)

N   = 500 # Number of particles
M1  = computeM1(c(gamma,c,sigma))
mu0 = rep(0,D)
# Define structure to store particles during the filtering
base_list = list("active_values" = c(),"Nnew" = c(), "survived" = c())
Particles = lapply(1:N, function(x){lapply(1:Ttot, function(xx){base_list})})
# Particles[[k]][[t]]$active_values   : matrix with the active features at time t for particle k
# Particles[[k]][[t]]$Nnew     : number of new features at time t for particle k
# Particles[[k]][[t]]$survived : indeces of the features of parent particle survived from time t-1 to time t for particle k



# Define auxiliary structures
A <- B <- log_w <- w <- W <- matrix(NA, nrow = Ttot, ncol = N)
# The final row of A is empty
# A[1,]: vector with the parents of the particles at time 2, i.e., X[2,k] is generated conditionally to X[1,A[1,k]]


# Run filtering -----------------------------------------------------------

# Inputs :  N;D;X;c;M1;sigma2_A;sigma2_X;
# Outputs:  Particles


# Time 1
# a) Sample the particles
t  = 1
for(k in 1:N){
  Nnew   = rpois(n = 1, lambda = M1)
  # values = sample_A(Nnew, X[t,], mu0, sigma2_X, sigma2_A  )
  ## OLD
  scale_mean_post = sigma2_A / (sigma2_A + sigma2_X)
  scale_var_post  = sigma2_A*sigma2_X / (sigma2_A + sigma2_X)
  values = matrix(0,nrow = Nnew, ncol = D)
  for(j in 1:D){
    values[,j] = rnorm(n=Nnew, mean = scale_mean_post*X[t,j], sd = sqrt( scale_var_post ) )
  }
  ## END OLD
  Particles[[k]][[t]]$active_values   = values
  Particles[[k]][[t]]$Nnew            = Nnew
  if(Nnew > 0)
    Particles[[k]][[t]]$survived        = seq(1,Nnew)
}

# b) Compute the un-nomralized weights
for(k in 1:N){
  N_t_new = Particles[[k]][[t]]$Nnew
  # log_w[t,k] = log_dmarg_img( N_t_new, X[t,], mu0, sigma2_X, sigma2_A )
  ## OLD
  log_w[t,k] = log_marginal(X[t,],sigma2_A,sigma2_X,N_t_new)
  ## END OLD
}
w[t,] = exp( log_w[t,] - max(log_w[t,]) )
W[t,] = w[t,]/sum(w[t,])

# Time t (t = 2,...,Ttot)
for(t in 2:Ttot){
  cat("\n ++++++ t = ",t," ++++++ \n")
  # a) Sample the parent node
  A[t-1,] = sample(1:N, size = N, W[t-1,], replace = TRUE)
  
  # b) Sample the particles
  for(k in 1:N){
    j = A[t-1,k] # parent index
    
    # Thinning part:
    Psurv = (sigma+1)/(c+1)
    num_active = nrow(Particles[[j]][[t-1]]$active_values) 
    survived = sample(c(0,1), size = num_active, replace = TRUE, prob = c(1-Psurv,Psurv))
    
    # Innovation:
    Xstar  = X[t,] - colSums(Particles[[j]][[t-1]]$active_values)
    Nnew   = rpois(n = 1, lambda = M1)
    # values = sample_A(Nnew, Xstar, mu0, sigma2_X, sigma2_A  )
    
    ## OLD
    values = matrix(0,nrow = Nnew, ncol = D)
    for(jj in 1:D){
      values[,jj] = rnorm(n=Nnew, mean = scale_mean_post*Xstar[jj], sd = sqrt( scale_var_post ) )
    }
    ## END OLD
      
    # Assemble thinning + innovation:
    if(sum(survived > 0)){
      survived_values = Particles[[j]][[t-1]]$active_values[which(survived > 0),]
      Particles[[k]][[t]]$active_values = rbind( survived_values, values   )
    }else{
      Particles[[k]][[t]]$active_values =  values 
    }
    Particles[[k]][[t]]$Nnew     = Nnew
    Particles[[k]][[t]]$survived = which(survived > 0)
    
    # Compute the un-normalized weights
    # log_w[t,k] = log_dmarg_img( Nnew, Xstar, mu0, sigma2_X, sigma2_A ) 
    ## OLD
    log_w[t,k] = log_marginal(X[t,],sigma2_A,sigma2_X,N_t_new)
    ## END OLD
  }

  # c) Compute the normalized weights
  w[t,] = exp( log_w[t,] - max(log_w[t,]) )
  W[t,] = w[t,]/sum(w[t,])
}

# Find ancestral lineage and save final paths
B[Ttot,] = 1:N
for(k in 1:N){
  for(t in (Ttot-1):1){
    B[t,k] = A[ t, B[t+1,k] ] # B_t^k = A_t^{B_{t+1}^k}
  }
}


# Find all features appeared at least once
Features_particle = lapply(1:N, function(x){c()})
cum_Ntot = lapply(1:N,function(x){c()})
Nnew_paths = lapply(1:N,function(x){c()})
for(k in 1:N){
  for(t in Ttot:1){
    Mat = Particles[[ B[t,k] ]][[t]]$active_values
    nuove = Particles[[ B[t,k] ]][[t]]$Nnew
    Nnew_paths[[k]] = c(nuove, Nnew_paths[[k]])
    cum_Ntot[[k]] = c(nuove, cum_Ntot[[k]])
    if(nuove > 0){
      if(nrow(Mat)-nuove+1 <= 0)
        stop("Unexprected behaviour")
      Features_particle[[k]] = rbind( Mat[(nrow(Mat)-nuove+1):nrow(Mat),] , Features_particle[[k]])
    }
  }
  cum_Ntot[[k]] = cumsum(cum_Ntot[[k]])
}

cum_Ntot

K_finals = sapply(cum_Ntot, function(x){tail(x,1)})
plot(table(K_finals))


# Plot features of a specific particle
k = 5
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(j in 1:K_finals[k]){
  Mat_vec = Features_particle[[k]][j,]
  Mat = matrix(data = Mat_vec, nrow = di, ncol = di, byrow = F)
  ACutils::ACheatmap(Mat, center_value = NULL, 
                     use_x11_device = FALSE, 
                     horizontal = F, main = paste0("Feat. ",j,"th"))
}


# Find survival path of each feature for each particle
Survived_features = vector("list",N)
for(k in 1:N){
  Survived_features[[k]] = matrix(0, nrow = Ttot, ncol = nrow(Features_particle[[k]]))
  # Time 1:
  t = 1
  if(Nnew_paths[[k]][t] > 0)
    Survived_features[[k]][t, 1:cum_Ntot[[k]][1] ] = 1
  for(t in 2:Ttot){
    # Time t (2...Ttot)
    if(Nnew_paths[[k]][t] > 0)
      Survived_features[[k]][t, (cum_Ntot[[k]][t-1]+1):(cum_Ntot[[k]][t])] = 1
    # Survived_features[[k]]
    Nsurvived_old = length(which(Survived_features[[k]][t-1,] == 1))
    submat = matrix( Survived_features[[k]][,which(Survived_features[[k]][t-1,] == 1)],
                     nrow = Ttot, ncol = Nsurvived_old)
    submat
    submat[t, Particles[[ B[t,k] ]][[t]]$survived ] = 1
    Survived_features[[k]][,which(Survived_features[[k]][t-1,] == 1)] = submat
  }
}


# Plot the path of a single particle
k = 5
mean_k_vec = Survived_features[[k]]%*%Features_particle[[k]]
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  mean_kt_vec = mean_k_vec[t,]
  mean_t = matrix(data = mean_kt_vec, nrow = di, ncol = di, byrow = F)
  ACutils::ACheatmap(mean_t, center_value = 0.5, use_x11_device = FALSE, 
                     horizontal = F, main = paste0("Obs. T = ",t))
}

# Final estimated filtered density
BB = 10000
ind_finals = sample(1:N, size = BB, replace = TRUE, prob = W[Ttot,])
post_means_all = vector("list",N)
for(k in 1:N){
  post_means_all[[k]] = Survived_features[[k]]%*%Features_particle[[k]]
}

result <- Reduce(`+`, post_means_all[ind_finals])
result = result / BB

# pdf("img/Adam_Est.pdf")
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  res_t_vec = result[t,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  ACutils::ACheatmap(res_t, center_value = 0.5, use_x11_device = FALSE, 
                     horizontal = F, main = paste0("Est. T = ",t))
}
# dev.off()





# Example with SeqMonteCarlo functions ------------------------------------


seed = 42
set.seed(seed)

N     = 100 # Number of particles
Psurv_true = 0.5
M1true = 1
c = 1/Psurv_true - 1
gamma = M1true * Psurv_true
sigma = 0
fit_smc = SeqMonteCarlo(X,N,D,Ttot,M1,Psurv,sigma2_A,sigma2_X, use_VS = FALSE)

k = 1
dim(fit_smc$Particles_Zmat[[k]])
dim(fit_smc$Particles_Amat[[k]])
fit_smc$Final_Weights
dim(fit_smc$Final_Est)


result = fit_smc$Final_Est
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  res_t_vec = result[t,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  ACutils::ACheatmap(res_t, 
                     center_value = 0, col.lower = "grey95",col.upper = "grey10",
                     use_x11_device = FALSE, 
                     horizontal = F, main = paste0("Est. T = ",t))
}

Ktot_mcmc = sapply(fit_smc$Particles_Zmat, function(x){ncol(x)})
Kt_mcmc = t(sapply(fit_smc$Particles_Zmat, function(x){apply(x,1,sum)}))
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  plot(Kt_mcmc[,t], type = "l", main = paste0("K_",t))
}




