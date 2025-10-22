wd = "C:/Users/colom/DynamicFeatureAllocation/Scripts/Testing/"
wd = "/home/lucia.paci/Lucia/Ale/DynamicFeatureAllocation/Scripts/Testing/"
setwd(wd)
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)

zipfs_decay = function(n,a){
  sapply(1:n,function(i){i^{-a}})
}

library(fields)
# Generate data -----------------------------------------------------------
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

# Plot and sizes:

# -> Lambda is VxK
# -> Xi is KxT
# -> Mean is VxT
# -> D is VxT

# Plot Lambda and Xi
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


# Sample D_itl ---------------------------------------------------------------
dim(Lambda)
dim(Xi)
zeta_list <- lapply(1:Ktrue, function(l) {
  Lambda[, l] %o% Xi[l,]   # outer product: V x T
})
# Combine into 3D array: V x T x K
zeta <- array(unlist(zeta_list), dim = c(V, Ttot, Ktrue))


dim(zeta)

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



# Plot Lambda, Xi, D_l
ymax = 0.75
ymax2 = max(Xi)

par(mfrow = c(1,3), mar = c(3.5,3.5,2,1), mgp=c(2,0.5,0))
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
image( 1:Ttot, 1:V, 
       t(D_itl[,,1]),   
       col = mycol,    
       xlab = "Time", 
       ylab = "Words",
       main = "D_1",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, V, length.out = min(V, 10)), 
     labels = round(seq(1, V, length.out = min(V, 10))),
     cex.axis = 0.7)
box()


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
image( 1:Ttot, 1:V, 
       t(D_itl[,,2]),   
       col = mycol,    
       xlab = "Time", 
       ylab = "Words",
       main = "D_1",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, V, length.out = min(V, 10)), 
     labels = round(seq(1, V, length.out = min(V, 10))),
     cex.axis = 0.7)
box()

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
image( 1:Ttot, 1:V, 
       t(D_itl[,,3]),   
       col = mycol,    
       xlab = "Time", 
       ylab = "Words",
       main = "D_3",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, V, length.out = min(V, 10)), 
     labels = round(seq(1, V, length.out = min(V, 10))),
     cex.axis = 0.7)
box()












