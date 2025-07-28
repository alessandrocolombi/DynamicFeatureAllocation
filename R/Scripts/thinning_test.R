# wd ----------------------------------------------------------------------
setwd("C:/Users/colom/DynamicFeatureAllocation/R/Scripts")


# funzioni ----------------------------------------------------------------

Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../Rfunctions.R" )
source("../genera_img_ilaria.R")


# setting -----------------------------------------------------------------

Nnew = 4
sig2X = 0.01
sig2A = 0.5 * 0.01

# Data generation ---------------------------------------------------------
seed = 2927
set.seed(seed)

Ttot = 24    # Time length
di   = 8
D    = di*di # Size of the problem
Psurv = 0.9

sigma2_X = sig2X
sigma2_A = sig2A

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

# pdf("img/Adam_Obs.pdf")
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  Xt_vec = X[t,]
  X_t = matrix(data = Xt_vec, nrow = di, ncol = di, byrow = F)
  plot_img(X_t, center_value = 0.5, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = paste0("Obs. T = ",t))
}
# dev.off()

# par(mfrow = c(1,1), mar = c(4,4,2,1), bty = "l")
# image(Z, col = c("black", "white"), 
#       xaxt = "n", yaxt = "n",
#       xlab="Features", ylab="Time") 


# Test --------------------------------------------------------------------
plot_multiple_imgs(X, plt_wnd = c(4,6))

t = 1
mu0 = rep(0,D)
#a) Sample A
values = sample_A(Nnew, X[t,], mu0, sig2X, sig2A  ) 
values = rbind(values , sample_A(Nnew, X[17,], mu0, sig2X, sig2A  ) )

dim(values)


plot_multiple_imgs(values, plt_wnd = c(2,4))

plot_multiple_imgs(matrix(colSums(values), nrow = 1, ncol = D), plt_wnd = c(1,1))


Psurv = 0.5
propose_thinning_Be(X[t,], values, sig2X, Psurv)

propose_thinning_VS(X[1,],values,sig2X, Psurv)
Xobs = X[1,]
ActFeat = values
perm = sample(1:nrow(values),size = nrow(values))
ActFeat = values[perm,]
plot_multiple_imgs(ActFeat, plt_wnd = c(2,4))

mat = rbind(mean_1,mean_0)
plot_multiple_imgs(mat, plt_wnd = c(1,2))

mean_final = active %*% ActFeat
mean_final_mat = matrix(mean_final,nrow = 1, ncol = D)
plot_multiple_imgs(mean_final_mat, plt_wnd = c(1,1))





