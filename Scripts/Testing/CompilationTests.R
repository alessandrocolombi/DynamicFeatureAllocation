wd = "C:/Users/colom/DynamicFeatureAllocation/Scripts/Testing/"
# wd = "/home/lucia.paci/Lucia/Ale/DynamicFeatureAllocation/Scripts/Testing/"
setwd(wd)

# sink("log.txt")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
# sink()

mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)

zipfs_decay = function(n,a){
  sapply(1:n,function(i){i^{-a}})
}

# test -----------------------------------------------------------
seed = 132332
set.seed(seed)

Ttot = 23
Ktrue = 3
V = 30

S <- matrix(rgamma(10*3, shape = 1, rate = 1), nrow = 10, ncol = 3)
lambda_1 <- lambda_2 <- lambda_3 <- rep(0,V)
lambda_1[1:10] =  S[,1]/sum(S[,1])
lambda_2[11:20] =  S[,2]/sum(S[,2])
lambda_3[21:30] =  S[,3]/sum(S[,3])

Lambda = cbind(lambda_1,lambda_2,lambda_3)

a = 1.001
Xi = matrix(0,nrow = Ktrue, ncol = Ttot)
Xi[1,1:6] =  400 * zipfs_decay(6,a)
Xi[2,7:14] = 500 * zipfs_decay(8,a)
Xi[3,15:23] = 40 * zipfs_decay(9,-a)

Mean = Lambda %*% Xi
D <- matrix(rpois(V * Ttot, lambda = as.vector(Mean)), nrow = V, ncol = Ttot)

classify_indices <- function(Z) {
  # ensure Z is numeric (0/1)
  Z <- as.numeric(Z)
  Ttot <- length(Z)
  
  idx_born = c()
  if(Z[1] > 0)
    idx_born = c(idx_born,1)
  idx_born <- c(idx_born,which(Z > 0 & (c(1, head(Z, -1)) == 0)))   # 1s preceded by 0 or start
  
  idx_noact <- which(Z == 0)                             # all zeros
  idx_surv <- which(Z > 0 & !(1:Ttot %in% idx_born))    # 1s that are not born
  
  list(idx_born = idx_born, idx_surv = idx_surv, idx_noact = idx_noact)
}
classify_indices(Xi)
aa = apply(Xi, 1, classify_indices)

delta = 0.5
Lambda_itl = vector("list",Ktrue)
Lambda_itl = lapply(Lambda_itl, function(x) {
  matrix( 0, nrow = V, ncol = Ttot )
})

for(l in 1:Ktrue){
  for(t in aa[[l]]$idx_born){
    Lambda_itl[[l]][,t] = Lambda[,l]
  }
  for(t in aa[[l]]$idx_surv){
    Lambda_itl[[l]][,t] = Lambda_itl[[l]][,t-1]
  }
  for(t in aa[[l]]$idx_noact){
    temp = rgamma(n = V, shape = delta, rate = 1)
    Lambda_itl[[l]][,t] = temp
  }
}



Z1 <- c(1,1,0,0,0,1,1,1)
Z2 <- c(1,1,0,0,0,0,0,0)
Z3 <- c(0,1,1,0,0,0,0,1)

res1 <- classify_indices(Z1)
res2 <- classify_indices(Z2)
res3 <- classify_indices(Z3)

res11 <- classify_indices_cpp(Z1)
res22 <- classify_indices_cpp(Z2)
res33 <- classify_indices_cpp(Z3)
# Print results
res1; res11

res2; res22
res3; res33


D_itl = vector("list",Ktrue)
D_itl = lapply(D_itl, function(x) {
  matrix( 1:(V*Ttot), nrow = V, ncol = Ttot )
})

zeta_list <- lapply(1:Ktrue, function(l) {
  Lambda[, l] %o% Xi[l,]   # outer product: V x T
})
zeta <- array(unlist(zeta_list), dim = c(V, Ttot, Ktrue))

D_itl = array(0,dim = c(V,Ttot,Ktrue))
for(i in 1:V){
  for(t in 1:Ttot){
    w_it = zeta[i,t,]
    if(sum(w_it)>0){
      w_it = w_it/sum(w_it)
      idx_it = sample(1:Ktrue, size = D[i,t], replace = TRUE, prob = w_it)
      D_itl[i,t,] = tabulate(idx_it,nbins = 3)
    }
  }
}

N_tl = apply(D_itl, c(2,3), sum)

D_itl_list = vector("list",3)
D_itl_list[[1]] = D_itl[,,1]
D_itl_list[[2]] = D_itl[,,2]
D_itl_list[[3]] = D_itl[,,3]

Test_sample_Ditl(Lambda_itl, Xi, D, seed)
Test_sample_Lambdaitl(D_itl_list, Xi, delta, seed)



D_itl[[1]]






