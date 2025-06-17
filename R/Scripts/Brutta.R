

t = 1;h=1
mean_t_vec = mu0
par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
mean_t = matrix(data = mean_t_vec, nrow = di, ncol = di, byrow = F)
plot_img(mean_t, center_value = 0.5, col.lower = "grey95",col.upper = "grey10",
         horizontal = F, main = paste0("Mean T = ",t))


log_dmarg_img( 0, X[t,], zetas[h,], sigma2_X, sigma2_ker )
log_dmarg_img( 1, X[t,], zetas[h,], sigma2_X, sigma2_ker )
log_dmarg_img( 2, X[t,], zetas[h,], sigma2_X, sigma2_ker )
log_dmarg_img( 5, X[t,], zetas[h,], sigma2_X, sigma2_ker )

Id = diag(D)
mvtnorm::dmvnorm(x = X[1,], mean = rep(0,D), sigma = sigma2_X*Id, log = TRUE)
mvtnorm::dmvnorm(x = X[1,], mean = zetas[1,], sigma = (sigma2_X+sigma2_ker)*Id, log = TRUE)


X = X; N = N; D = D; Ttot = Ttot;
mu0 = zetas_mcmc[[g]][];
Bfix = rep(-1,Ttot); Particle_fix = NULL;  
M1 = M1_mcmc[1];Psurv = Psurv_mcmc[1];
sigma2_A = sigma2_ker;
sigma2_X = sigma2_X_mcmc[g,1];
use_VS = use_VS

Particles_PG[[1]]$Nnew_paths_k
Particles_PG[[1]]$cum_Ntot_k

xx = Particles_PG[[1]]$mean_k
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  Xt_vec = xx[t,]
  X_t = matrix(data = Xt_vec, nrow = di, ncol = di, byrow = F)
  plot_img(X_t, 
           center_value = NULL, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = paste0("... prova ...",t))
}










# Final estimated filtered density in each group
filt_dens_h = lapply(G_th, function(x){
  result <- Reduce(`+`, lapply(x$Particles_PG[1], function(xx){xx$mean_k}) )
  result = result / 1
})

g = 4
filt_dens_h = lapply(G_th, function(x){x$Particles_PG[[g]]$mean_k})
for(h in 1:M){
  par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
  for(t in 1:Ttot){
    res_t_vec = filt_dens_h[[h]][t,]
    res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
    plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
              horizontal = F, main = paste0("Gr. ",h,"; Est. T = ",t))
  }
}

# Final estimated filtered density 
filt_dens = Reduce(`+`, filt_dens_h )
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  res_t_vec = filt_dens[t,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("Gr. ",h,"; Est. T = ",t))
}






h = 4; g = 1
N = 1000
temp = Conditional_SeqMonteCarlo( X = X, N = N, D = D, Ttot = Ttot,
                                   mu0 = zetas_mcmc[[g]][h,],
                                   Bfix = rep(-1,Ttot), Particle_fix = NULL,  
                                   M1 = M1_mcmc[1],Psurv = Psurv_mcmc[1],
                                   sigma2_A = sigma2_ker, 
                                   sigma2_X = sigma2_X_mcmc[g,1],
                                   use_VS = use_VS)

temp$Nnew_paths_k
temp$cum_Ntot_k

t = 1
mean_t_vec = X[t,]
par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
mean_t = matrix(data = mean_t_vec, nrow = di, ncol = di, byrow = F)
plot_img(mean_t, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
         horizontal = F, main = paste0("Mean T = ",t))




g
Ktot[1]
Kh = sapply(G_th, function(x){
  Kh = nrow(x$Particles_PG[[g-1]]$Amat_k)
  ifelse(is.null(Kh),0,Kh)
})
cgr_old = unlist(sapply(1:length(Kh), function(h){rep(h,Kh[h])}))
cgr_new = sample(1:M, size = Ktot[g-1], replace = TRUE)

Amat_new = rbind(
  Amat[which(cgr_new == 1),],
  Amat[which(cgr_new == 2),],
  Amat[which(cgr_new == 3),],
  Amat[which(cgr_new == 4),])


Zmat_new = cbind(
  Zmat[,which(cgr_new == 1)],
  Zmat[,which(cgr_new == 2)],
  Zmat[,which(cgr_new == 3)],
  Zmat[,which(cgr_new == 4)]
)


# aggiornare zeta_h | Amat_h è facile








Ktot_1T = nrow(Amat_k)
Zmat_k = matrix(0, nrow = Ttot, ncol = Ktot_1T)
t = 1
if(Nnew_paths_k[t] > 0){
  Zmat_k[t, 1:cum_Ntot_k[1] ] = 1
}
for(t in 2:Ttot){
  # Time t (2...Ttot)
  if(Nnew_paths_k[t] > 0)
    Zmat_k[t, (cum_Ntot_k[t-1]+1):(cum_Ntot_k[t])] = 1
  
  Nsurvived_old = length(which(Zmat_k[t-1,] == 1))
  submat = matrix( Zmat_k[,which(Zmat_k[t-1,] == 1)],
                   nrow = Ttot, ncol = Nsurvived_old)
  submat[t, Path_k[[t]]$survived ] = 1
  Zmat_k[,which(Zmat_k[t-1,] == 1)] = submat
}



Amat_k_h = Amat_k[which(gr_alloc_k == 1),]



















N = 50
X = X; N = N; D = D; Ttot = Ttot; 
Bfix = rep(-1,Ttot); Particle_fix = NULL;  
M1 = M1_mcmc[1]; Psurv = Psurv_mcmc[1] ;
sigma2_X = sigma2_XA_mcmc[1,1];
sigma2_A = sigma2_XA_mcmc[1,2];
zeta = zetas_mcmc[[1]]; 
use_VS = use_VS



vv <- c(1, 3, 5, 0, 4, 0, 2)

# Get the indices of non-zero and zero elements
nonzero_indices <- which(vv != 0)
zero_indices    <- which(vv == 0)

# Combine to get the new ordering
ordering <- c(nonzero_indices, zero_indices)

# Reorder vv (optional, just to check)
vv_sorted <- vv[ordering]


# Original list with 5 elements
my_list <- list("a", "b", "c", "d", "e")

# Check length before
length(my_list)
# [1] 5

# Remove the 3rd element
my_list[[3]] <- NULL

# Check length after
length(my_list)
# [1] 4

my_list <- list(a = 1, b = 2)
# Add a new element named "new"
my_list[["new"]] <- 3


# Example list
my_list <- as.list(letters[1:10])  # list of 10 elements: "a", "b", ..., "j"

# Index to move
h <- 3  # Move the 3rd element to the end

# Reorder the list
my_list <- my_list[c(setdiff(seq_along(my_list), h), h)]

# Check result
print(my_list)





vv = c(1,5,3,0,1,0,3,4)

data = Amat; 
Hold = H_mcmc[g-1]; Hstar_old = Hstar_mcmc[g-1];
gr_alloc_old = Particles_PG[[g-1]]$gr_alloc_k; 
centers_old = zetas_mcmc[[g-1]];  
sig2_A = sigma2_XA_mcmc[g,2];
sig2_ker = sigma2_ker;
updateGrAlloc = updateGrAlloc; 
updateCenters = updateCenters; 
upNcenters = upNcenters





k = 1
mean_t_vec = l3[9,]
par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
mean_t = matrix(data = mean_t_vec, nrow = di, ncol = di, byrow = F)
plot_img(mean_t, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
         horizontal = F, main = paste0("... prova ... "))

X = X; N = N; D = D; Ttot = Ttot; 
Bfix = Particles_PG[[g-1]]$B_k;
Particle_fix = Particle_fix;
M1 = M1_mcmc[g]; Psurv = Psurv_mcmc[g] 
sigma2_X = sigma2_XA_mcmc[g,1]
sigma2_A = sigma2_XA_mcmc[g,2]
zeta = zetas_mcmc[[g]]
use_VS = use_VS

Path_k = Particles_PG[[1]]$Path_k
Ktot = 0
t = 1
n_active = length(Path_k[[t]]$gr_alloc)
Path_k[[t]][["label_actives"]] <- c()  
if(Path_k[[t]]$Nnew > 0){
  Path_k[[t]]$label_actives =  c(Path_k[[t]]$label_actives,1:n_active)
  Ktot = Ktot + Path_k[[t]]$Nnew 
}
t = 2
for(t in 2:Ttot){
  n_active = length(Path_k[[t]]$gr_alloc)
  Path_k[[t]][["label_actives"]] <- c()
  # Survived
  n_surv = length(Path_k[[t]]$survived)
  if(n_surv > 0){
    Path_k[[t]]$label_actives = c( Path_k[[t]]$label_actives,
                                   Path_k[[t-1]]$label_actives[Path_k[[t]]$survived] )
  }
  # New
  if(Path_k[[t]]$Nnew > 0){
    
    Path_k[[t]]$label_actives = c( Path_k[[t]]$label_actives,
                                   (Ktot+1):(Ktot+Path_k[[t]]$Nnew) )
    Ktot = Ktot + Path_k[[t]]$Nnew 
  }
}


x = Particles_PG[[1]]$Path_k
Amat = Particles_PG[[1]]$Amat_k

Amat[x[[1]]$label_actives,]


x[[1]]$label_actives
x[[2]]$label_actives
x[[3]]$label_actives
x[[4]]$label_actives
x[[5]]$label_actives
x[[6]]$label_actives
x[[7]]$label_actives
x[[24]]$label_actives


g = 1
X = X; N = N; D = D; Ttot = Ttot; 
Bfix = rep(-1,Ttot);Particle_fix = NULL;
M1 = M1_mcmc[g]; Psurv = Psurv_mcmc[g] ;
sigma2_X = sigma2_XA_mcmc[g,1];
sigma2_A = sigma2_XA_mcmc[g,2];
zeta = zetas_mcmc[[g]]; 
use_VS = FALSE 



x = Particles_PG[[1]]
x$Nnew_paths_k
x$cum_Ntot_k
for(t in 1:Ttot){
  plot_v2m(x$mean_k[t,], di = 8, center_val_plot = NULL)
}
# plot_v2m(x$mean_k[1,], di = 8, center_val_plot = NULL)

plot_v2m(zetas_mcmc[[g]][[1]], di = 8, center_val_plot = NULL)
plot_v2m(zetas_mcmc[[g]][[2]], di = 8, center_val_plot = NULL)
plot_v2m(zetas_mcmc[[g]][[3]], di = 8, center_val_plot = NULL)
plot_v2m(zetas_mcmc[[g]][[4]], di = 8, center_val_plot = NULL)



g = 2
mean_vec_1T = Particles_PG[[g-1]]$mean_k
Amat        = Particles_PG[[g-1]]$Amat_k
plot_v2m(mean_vec_1T[20,], di = 8, center_val_plot = NULL)


g = 2
temp = r_fullcond_zeta( data = Amat, 
                        Hold = H_mcmc[g-1], Hstar_old = Hstar_mcmc[g-1],
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








sigma2_X = 0.001; sigma2_A = 0.001
log_dmarg_img(2,X[t,],zeta[[1]],sigma2_X,sigma2_A )

nnn_h = c(1,1,0,0)
www = log_dmarg_img(nnn_h[1],X[t,],zeta[[1]],sigma2_X,sigma2_A )+
      log_dmarg_img(nnn_h[2],X[t,],zeta[[2]],sigma2_X,sigma2_A )+
      log_dmarg_img(nnn_h[3],X[t,],zeta[[3]],sigma2_X,sigma2_A )+
      log_dmarg_img(nnn_h[4],X[t,],zeta[[4]],sigma2_X,sigma2_A )
www

plot_v2m(X[t,], di = 8, center_val_plot = NULL)
g = 51
plot_v2m(zetas_mcmc[[g]][[1]], di = 8, center_val_plot = NULL)
plot_v2m(zetas_mcmc[[g]][[2]], di = 8, center_val_plot = NULL)
plot_v2m(zetas_mcmc[[g]][[3]], di = 8, center_val_plot = NULL)
plot_v2m(zetas_mcmc[[g]][[4]], di = 8, center_val_plot = NULL)




x = c(25000,99999,0)
x[1] * exp( lgamma(x[3]+1) + lgamma(x[2]-x[3]+1) - lgamma(x[2]+2) )

for(t in 1:Ttot){
  plot_v2m(mean_vec_1T[t,], di = 8, center_val_plot = NULL)
}

log_dmarg_img(5,X[t,],zeta[[1]],sigma2_X,sigma2_A )




Particle_fix[[1]]$gr_card

lapply(Particle_fix, function(x){length(x$gr_card)})






r = 100
n1 <- n2 <- 1000
gamma_1 <- gamma_2 <- 1
-gamma_1*log(n1) - gamma_2*log(n2) + lgamma( gamma_1*(r+1)  ) + lgamma( gamma_2*(r+1)  ) - lgamma( gamma_1*(r)  ) - lgamma( gamma_2*(r)  )

gamma_1 <- gamma_2 <- 10
-gamma_1*log(n1) - gamma_2*log(n2) + lgamma( gamma_1*(r+1)  ) + lgamma( gamma_2*(r+1)  ) - lgamma( gamma_1*(r)  ) - lgamma( gamma_2*(r)  )

gamma_1 <- gamma_2 <- 100
-gamma_1*log(n1) - gamma_2*log(n2) + lgamma( gamma_1*(r+1)  ) + lgamma( gamma_2*(r+1)  ) - lgamma( gamma_1*(r)  ) - lgamma( gamma_2*(r)  )



gamma_1 <- gamma_2 <- 10000
-gamma_1*log(n1) - gamma_2*log(n2) + lgamma( gamma_1*(r+1)  ) + lgamma( gamma_2*(r+1)  ) - lgamma( gamma_1*(r)  ) - lgamma( gamma_2*(r)  )




r = 100
n1 <- n2 <- 1000
gamma_1 <- gamma_2 <- 1
-gamma_1*log(n1) - gamma_2*log(n2) + lgamma( gamma_1*(r+1)  ) + lgamma( gamma_2*(r+1)  ) - lgamma( gamma_1*(r)  ) - lgamma( gamma_2*(r)  )

gamma_1 <- gamma_2 <- 1e-1
-gamma_1*log(n1) - gamma_2*log(n2) + lgamma( gamma_1*(r+1)  ) + lgamma( gamma_2*(r+1)  ) - lgamma( gamma_1*(r)  ) - lgamma( gamma_2*(r)  )

gamma_1 <- gamma_2 <- 1e-2
-gamma_1*log(n1) - gamma_2*log(n2) + lgamma( gamma_1*(r+1)  ) + lgamma( gamma_2*(r+1)  ) - lgamma( gamma_1*(r)  ) - lgamma( gamma_2*(r)  )

gamma_1 <- gamma_2 <- 1e-6
-gamma_1*log(n1) - gamma_2*log(n2) + lgamma( gamma_1*(r+1)  ) + lgamma( gamma_2*(r+1)  ) - lgamma( gamma_1*(r)  ) - lgamma( gamma_2*(r)  )


gamma_1 <- gamma_2 <- 1e-10
-gamma_1*log(n1) - gamma_2*log(n2) + lgamma( gamma_1*(r+1)  ) + lgamma( gamma_2*(r+1)  ) - lgamma( gamma_1*(r)  ) - lgamma( gamma_2*(r)  )



r = 5
x = 1e-16
exp(lgamma(x*(r+1))-lgamma(x*r))
r/(r+1)

exp(lgamma(x*(r+1)))
exp(-lgamma(x*(r)))






ss = 0
for(j in 1:nrow(A)){
  ss = ss  + t(A[j,])%*%A[j,]
}
ss

ss = 0
for(j in 1:ncol(A)){
  ss = ss  + t(A[,j])%*%A[,j]
}
ss

par_sig2A = set_par_invgamma(media = 0.001,   var = 0.0001)
a = par_sig2A[1]; b = par_sig2A[2]
num = b + 0.5*ss
den = a + D*4*0.5 - 1

num/den


# I want to debug the variable selection procedure
X = X; N = N; D = D; Ttot = Ttot; Bfix = Particles_PG[[g-1]]$B_k;
Particle_fix = Particles_PG[[g-1]]$Path_k;  
M1 = M1_mcmc[g];Psurv = Psurv_mcmc[g];
sigma2_A = sigma2_XA_mcmc[g,2]
sigma2_X = sigma2_XA_mcmc[g,1]
use_VS = TRUE  
mu0 = rep(0,D)


t = 18
ActVal = Particles[[j]][[t-1]]$active_values
dim(ActVal)

par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
mean_t = matrix(data = values[3,], nrow = di, ncol = di, byrow = F)
plot_img(mean_t, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
         horizontal = F, main = paste0("... prova ... "))


Xobs = X[t,]
ActFeat = ActVal
sig2X = sigma2_X
Psurv = Psurv 

D = ncol(ActFeat)
p = nrow(ActFeat)
active = rep(1, p)
log_prop_move = 0
for(j in 1:p){
  active_prop_1 <- active_prop_0 <- active
  active_prop_1[j] = 1;active_prop_0[j] = 0
  mean_1 = active_prop_1 %*% ActFeat
  mean_0 = active_prop_0 %*% ActFeat
  par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
  mean_t = matrix(data = mean_1, nrow = di, ncol = di, byrow = F)
  plot_img(mean_t, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = paste0("... prova ... "))
  
  log_prob_1 <- log_dmvnorm(x = Xobs, mean = mean_1, Sigma = sig2X*diag(D)) - 2*D
  log_prob_0 <- log_dmvnorm(x = Xobs, mean = mean_0, Sigma = sig2X*diag(D))
  
  log_den = log_stable_sum( c(log_prob_1,log_prob_0), is_log = TRUE )
  log_prob_active_j = log_prob_1 - log_den  
  if(log(runif(1)) < log_prob_active_j){
    active_j = 1
    log_prop_move = log_prop_move + log_prob_active_j
  }else{
    active_j = 0
    log_prop_move = log_prop_move + log_prob_0 - log_den
  }
  
  active[j] = active_j
  
}
active


par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
mean_t = matrix(data = ActVal[1,]+ActVal[2,], nrow = di, ncol = di, byrow = F)
plot_img(mean_t, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
         horizontal = F, main = paste0("... prova ... "))





















