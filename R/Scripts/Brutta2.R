X = X; N = N; D = D; Ttot = Ttot;
Bfix = rep(-1,Ttot); Particle_fix = NULL;  
M1 = M1_mcmc[1];Psurv = Psurv_mcmc[1];
sigma2_A = sigma2_XA_mcmc[1,2]; 
sigma2_X = sigma2_XA_mcmc[1,1];
use_VS = TRUE; mu0 = rep(0,D)


survived_particles = A[t-1,]
selected_Nnew = sapply(Particles[survived_particles], function(x) x[[t-1]]$Nnew)



plot_multiple_imgs(ActFeat[9:12,], plt_wnd = c(2,2))

aaa = matrix(colSums(ActFeat[6:12,]), nrow = 1, ncol = D)
plot_multiple_imgs(aaa, plt_wnd = c(1,1))
aaa = matrix(colSums(ActFeat[1:5,]), nrow = 1, ncol = D)
plot_multiple_imgs(aaa, plt_wnd = c(1,1))



means = rbind(mean_0,mean_1)
plot_multiple_imgs(means, plt_wnd = c(1,2))

plot_multiple_imgs(values, plt_wnd = c(2,4))

aaa = matrix(colSums(values), nrow = 1, ncol = D)
plot_multiple_imgs(aaa, plt_wnd = c(1,1))

plot_multiple_imgs(aux, plt_wnd = c(1,1))

plot_multiple_imgs(Xstar, plt_wnd = c(1,1))

plot_multiple_imgs(Particles[[k]][[t]]$active_values, plt_wnd = c(3,3))
aaa = matrix(colSums(Particles[[k]][[t]]$active_values), nrow = 1, ncol = D)
plot_multiple_imgs(aaa, plt_wnd = c(1,1))


X_mat = matrix(X[3,],nrow = 1, ncol = D)
plot_multiple_imgs(X_mat, plt_wnd = c(1,1))


ActFeat = Particles[[k]][[2]]$active_values

ActFeat = Particles[[27]][[t-1]]$active_values




X = X; N = N; D = D; Ttot = Ttot;
Bfix = Particles_PG[[g-1]]$B_k;
Particle_fix = Particles_PG[[g-1]]$Path_k;  
M1 = M1_mcmc[g];Psurv = Psurv_mcmc[g];
sigma2_A = sigma2_XA_mcmc[g,2]; 
sigma2_X = sigma2_XA_mcmc[g,1];
proposal_Nnew_1T = proposal_N;
use_VS = use_VS
mu0 = rep(0,D)


sum(W[t,])
round(W[t,],3)
round(W[t,22],3)

plot_multiple_imgs(Xstar, plt_wnd = c(1,1))
aaa = matrix(matrix(Xstar,1,D), nrow = 1, ncol = D)
plot_multiple_imgs(aaa, plt_wnd = c(1,1))
ActFeat_h
aaa = matrix(matrix(ActFeat_h,1,D), nrow = 1, ncol = D)
plot_multiple_imgs(aaa, plt_wnd = c(1,1))
# check Bfix --------------------------------------------------------------

for(t in 1:Ttot){
  k = Bfix[t]
  ActFeat = Particles_PG[[g-1]]$Path_k[[t]]$active_values
  plot_multiple_imgs(ActFeat, plt_wnd = c(1,2))
}






