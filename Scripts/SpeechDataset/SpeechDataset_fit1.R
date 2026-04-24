wd = "C:/Users/colom/DynamicFeatureAllocation/Scripts/SpeechDataset"
# wd = "/home/lucia.paci/Lucia/Ale/DynamicFeatureAllocation/Scripts/SpeechDataset"
setwd(wd)
source("./../../R/Rfunctions.R")
Rcpp::sourceCpp("./../../src/RcppFunctions.cpp")
mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)
library(fields)
library(fangs)

# Read data -----------------------------------------------------------

seed = 132332
set.seed(seed)
r = 10

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

H = 20 # number of atoms
delta = 0.001 # Dirichlet parameter

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
niter = 2000
nburn = 5000
thin  =   10
data  = as.matrix(data)
fit = GibbsSampler_DTM(niter,nburn,thin,data,param_DTM,init_DTM)

# Atoms to path transformation ---------------------------------------------------------------

Lambda_fit = lapply(fit$Lambda_star_mcmc, function(Lam_it) t( do.call(rbind,Lam_it) ))
# --> list of length (V x K_it) matrices, one for every iteration
K_it = sapply(Lambda_fit, ncol)

topic_objs <- lapply(1:length(fit$Xi), function(it) { build_topic_matrices(fit$Xi[[it]], Ttot) })

# Diagnosis ---------------------------------------------------------------

it_start = 1     # da dove parto + burnin + val iniziale
it_end   = niter # dove finisco

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


## Active features matrix  ------------------------------------------------------------------

Zmat_list  = lapply(topic_objs[it_start:it_end], function(z) z$Activity)
Ximat_list = lapply(topic_objs[it_start:it_end], function(z) t(z$Xi_star) )

prova = fangs(Zmat_list)
View(prova)
est = prova$estimate
est

Kest = ncol(est)

par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Kest, 1:Ttot, 
       t(est),   
       col = mycol,    
       xlab = "Topics", 
       ylab = "Time",
       main = "<Z est>",
       axes = FALSE )
axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(1, at = seq(1, Kest, length.out = min(Kest, 10)), 
     labels = round(seq(1, Kest, length.out = min(Kest, 10))),
     cex.axis = 0.7)
box()
# fields::image.plot(
#   1:Ttot, 1:Kest, est,
#   col = mycol,
#   legend.only = TRUE,
#   horizontal = FALSE,
#   legend.width = 1.2,            # controls legend thickness
#   legend.shrink = 0.8,           # smaller legend
#   legend.mar = 8.5,                # margin from image
#   legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
# )

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


## Unique values for a given iteration ------------------------------------------------------------------
it = 1950
ymax_topic = max(Lambda_fit[[it]])
ymax_xi = max( topic_objs[[it]]$Xi_star )
for(l in 1:K_it[it]){
  par(mfrow = c(1,2), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
  barplot( height = Lambda_fit[[it]][,l], 
           names.arg = as.character(1:V),
           las = 1, col = "darkred", border = NA,
           xlab = "Word",
           main = paste0("Atom: ",l), ylab = "Prob.", ylim = c(0,ymax_topic),
           cex.names = 0.5 )
  barplot( height = topic_objs[[it]]$Xi_star[l,], 
           names.arg = as.character(1:Ttot),
           las = 1, col = "darkgreen", border = NA,
           xlab = "Time",
           main = "", ylab = "Intensity", ylim = c(0,ymax_xi),
           cex.names = 0.5 )
}

## Weighted pairwise similarity matrix ------------------------------------------------------------------

Ximat_list = lapply(topic_objs[it_start:it_end], function(z) t(z$Xi_star) )
Adj_mat_list = lapply(Ximat_list, function(x) x%*%t(x) )

Wpsm <- Reduce(`+`, Adj_mat_list)/length(Adj_mat_list)
Wpsm <- sqrt(Wpsm)

par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:Ttot, 
       Wpsm,   
       col = mycol,    
       xlab = "Time", 
       ylab = "Time",
       main = "Wpsm",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:Ttot, Wpsm,
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,                # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)

## Unweighted pairwise similarity matrix ------------------------------------------------------------------

Adj_bin_mat_list = lapply(Zmat_list, function(z) z%*%t(z) )

Upsm <- Reduce(`+`, Adj_bin_mat_list)/length(Adj_bin_mat_list)

par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:Ttot, 
       Upsm,   
       col = mycol,    
       xlab = "Time", 
       ylab = "Time",
       main = "Upsm",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:Ttot, Upsm,
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,                # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)



# Left ordering ---------------------------------------------------------


Z = Z_it
left_order <- function(Z) {
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  if (!all(Z %in% c(0, 1)))
    stop("Z must be a binary matrix (0/1).")
  
  # Create a binary "signature" for each column
  # Interpreted lexicographically: top row = most significant bit
  col_signature <- apply(Z, 2, function(col) {
    # Paste as "010110" etc.
    paste(col, collapse = "")
  })
  
  # # Convert to integers for sorting (safe up to ~50 rows)
  # col_value <- strtoi(col_signature, base = 2)
  # 
  # # Optional: use column sums as secondary key (common convention)
  # col_sum <- colSums(Z)
  
  # Order: decreasing lexicographic, then decreasing sum
  ord <- order(col_signature, decreasing = TRUE)
  
  return(Z[,ord])
}


par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
for(ii in 1995:2000){
  Z_it = Zmat_list[[ii]]
  Z_it = left_order(Z_it)
  image( 1:ncol(Z_it), 1:Ttot, 
         t(Z_it),   
         col = mycol,    
         xlab = "Topics", 
         ylab = "Time",
         main = " ",
         axes = FALSE )
  axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(1, at = seq(1, ncol(Z_it), length.out = min(Kest, 10)), 
       labels = round(seq(1, Kest, length.out = min(Kest, 10))),
       cex.axis = 0.7)
  box()
}

par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
# Finale
Z_it = est
Z_it = left_order(Z_it)
image( 1:ncol(Z_it), 1:Ttot, 
       t(Z_it),   
       col = mycol,    
       xlab = "Topics", 
       ylab = "Time",
       main = " ",
       axes = FALSE )
axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(1, at = seq(1, ncol(Z_it), length.out = min(Kest, 10)), 
     labels = round(seq(1, Kest, length.out = min(Kest, 10))),
     cex.axis = 0.7)
box()



sum(Z_it)


# Brutta iterazione -------------------------------------------------------


it = 2000
Z_it = Zmat_list[[it]]
Lambda_star_it = fit$Lambda_star_mcmc[[it]]
Xi_it = fit$Xi[[it]]

idx_list = lapply(1:H, function(l) find_indices(Xi_it[l,]))
for(l in 1:H){
  Num_act_l = length(idx_list[[l]]$idx_born) # how many distinct values
  if(Num_act_l > 0){
    for(jj in seq_along(idx_list[[l]]$idx_born) ){
      
      # Find unique value
      tj = idx_list[[l]]$idx_born[jj]
      # lambda = fit$Lambda[[it]][[l]][,tj]
      lambda = Lambda_star_it[[l]][jj,]
      
      # Find activity times
      tj_next = idx_list[[l]]$idx_born[jj+1]
      if(is.na(tj_next))
        tj_next = Ttot + 1
      
      time_act = c( tj,intersect(tj:tj_next, idx_list[[l]]$idx_surv) )
      active <- Xi_active <- rep(0,Ttot)
      active[ time_act ] = 1
      Xi_active[ time_act ] = fit$Xi[[it]][l,time_act]
      par(mfrow = c(1,2), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
      barplot( height = lambda, 
               names.arg = as.character(1:V),
               las = 1, col = "darkred", border = NA,
               xlab = "Word",
               main = paste0("Atom: ",l,",",jj), ylab = "Prob.", ylim = c(0,1),
               cex.names = 0.5 )
      barplot( height = Xi_active, 
               names.arg = as.character(1:Ttot),
               las = 1, col = "darkgreen", border = NA,
               xlab = "Time",
               main = "", ylab = "Intensity", ylim = c(0,max(Xi_active)),
               cex.names = 0.5 )
    }
    
  }
}




lambda = Lambda_star_it[[l]][jj,]

ind_lambda_sort = order(lambda, decreasing = TRUE)[1:20]
data[ind_lambda_sort, 1]

round(rbind(data[ind_lambda_sort, 1], lambda[ind_lambda_sort]),3)


