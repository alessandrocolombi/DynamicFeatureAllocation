wd = "C:/Users/colom/DynamicFeatureAllocation/Scripts/Testing/"
# wd = "/home/lucia.paci/Lucia/Ale/DynamicFeatureAllocation/Scripts/Testing/"
setwd(wd)

# sink("log.txt")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
# sink()



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



test_FC(Lambda,Xi,D,seed)



dim(Xi)



classify_indices <- function(Z) {
  # ensure Z is numeric (0/1)
  Z <- as.numeric(Z)
  Ttot <- length(Z)
  
  idx_born = c()
  if(Z[1] == 1)
    idx_born = c(idx_born,1)
  idx_born <- c(idx_born,which(Z == 1 & (c(1, head(Z, -1)) == 0)))   # 1s preceded by 0 or start
  
  idx_noact <- which(Z == 0)                             # all zeros
  idx_surv <- which(Z == 1 & !(1:Ttot %in% idx_born))    # 1s that are not born
  
  list(idx_born = idx_born, idx_surv = idx_surv, idx_noact = idx_noact)
}
Z1 <- c(1,1,0,0,0,1,1,1)
Z2 <- c(1,1,0,0,0,0,0,0)
Z3 <- c(0,1,1,0,0,0,0,1)

res1 <- classify_indices(Z1)
res2 <- classify_indices(Z2)
res3 <- classify_indices(Z3)

# Print results
res1
res2
res3


