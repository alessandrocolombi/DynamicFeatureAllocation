wd = "C:/Users/colom/DynamicFeatureAllocation/Scripts/"
# wd = "/home/lucia.paci/Lucia/Ale/DynamicFeatureAllocation/Scripts/"
setwd(wd)

source("./../R/Rfunctions_centers.R")
Rcpp::sourceCpp("./../src/RcppFunctions_centers.cpp")

mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)
library(fields)
library(fangs)

# Generate data -----------------------------------------------------------

seed = 132332
set.seed(seed)

zipfs_decay = function(n,a){
  sapply(1:n,function(i){i^{-a}})
}

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
Xi[1,1:6]   = 400 * zipfs_decay(6,a)
Xi[2,7:14]  = 500 * zipfs_decay(8,a)
Xi[3,15:23] = 100 * zipfs_decay(9,-a)

Mean = Lambda %*% Xi
D <- matrix(rpois(V * Ttot, lambda = as.vector(Mean)), nrow = V, ncol = Ttot)
data = D
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


ymax_topic = max(Lambda)*1.2
ymax_xi = max(Xi)
for(k in 1:Ktrue){
  par(mfrow = c(1,2), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
  barplot( height = unname(Lambda[,k]),
           names.arg = as.character(1:V),
           las = 1, col = "darkred", border = NA,
           xlab = "Word", ylab = "Prob.", ylim = c(0,ymax_topic),
           cex.names = 0.5 )
  barplot( height = Xi[k,], 
           names.arg = as.character(1:Ttot),
           las = 1, col = "darkgreen", border = NA,
           xlab = "Time",
           main = "", ylab = "Intensity", ylim = c(0,ymax_xi),
           cex.names = 0.5 )
}

# Set options -------------------------------------------------------------

seed = 22123
set.seed(seed)

save_img = FALSE

H = 10              # number of atoms per center/process
M0 = Ktrue          # fixed number of centers

# Process hyperparameters. 
gamma = 1
sigma = 0.25
beta  = 0.5

# Center-clustering hyperparameters
a_phi_centers = 2
b_phi_centers = 1
delta0_centers = 1
omega = 1
a_omega = 1
b_omega = 1
var_phi_centers = 0.01
var_delta_centers = 0.01
mstar_max = 0       # bypassed in the fixed-M sampler


# Initial values ----------------------------------------------------------


Lambda_smooth = apply(Lambda, 2, smooth_simplex)

Xi0 = vector("list", M0)
S0 = vector("list", M0)
Lambda0 = vector("list", M0)
Zeta0 = vector("list", M0)

for(m in 1:M0){
  Xi0[[m]] = matrix(0, nrow = H, ncol = Ttot)
  Xi0[[m]][1,] = Xi[m,]

  S0[[m]] = matrix(rgamma(n=H*Ttot, shape=1, rate=1), nrow = H, ncol = Ttot)

  # Center parameters must be strictly positive.
  Zeta0[[m]] = 0.2 + 30 * Lambda_smooth[,m]

  Lambda0[[m]] = vector("list", H)
  for(l in 1:H){
    A = matrix(0, nrow = V, ncol = Ttot)
    for(t in 1:Ttot){
      A[,t] = rdirichlet_vec(Zeta0[[m]])
    }
    Lambda0[[m]][[l]] = A
  }

  # Give the first atom of each center a good starting value.
  for(t in 1:Ttot){
    Lambda0[[m]][[1]][,t] = Lambda_smooth[,m]
  }
}

init_DTM_centers = set_init_DTM_centers(Xi0,Lambda0,S0,Zeta0)


ymax_zeta0 = max(unlist(sapply(Zeta0,max)))*1.2
for(m in 1:M0){
  par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
  barplot( height = unname(Zeta0[[m]]),
           names.arg = as.character(1:V),
           las = 1, col = "darkred", border = NA,
           xlab = "Word", ylab = "Prob.", ylim = c(0,ymax_zeta0),
           cex.names = 0.5 )
}

# Sampler flags -----------------------------------------------------------

UpdateDitl = TRUE
UpdateS = TRUE
UpdateLambda = TRUE
UpdateXi = TRUE       # fixed here for this first minimal example
UpdateU = TRUE
UpdateCenters = TRUE
UpdateOmega = FALSE   
print = TRUE

param_DTM_centers = set_param_DTM_centers(H,gamma,sigma,beta,
                                          a_phi_centers,b_phi_centers,delta0_centers,
                                          omega,a_omega,b_omega,
                                          var_phi_centers,var_delta_centers,mstar_max,
                                          UpdateDitl,UpdateS,UpdateLambda,UpdateXi,UpdateU,
                                          UpdateCenters,UpdateOmega,
                                          seed,print)


# Run ---------------------------------------------------------------------

niter = 500
nburn = 100
thin  = 1

fit = GibbsSampler_DTM_centers(niter,nburn,thin,D,param_DTM_centers,init_DTM_centers)

cat("Run completed\n")
cat("Saved iterations:", length(fit$M), "\n")
cat("Unique M values:", paste(unique(fit$M), collapse = ", "), "\n")
cat("Unique Mstar values:", paste(unique(fit$Mstar), collapse = ", "), "\n")


# Atoms to path transformation ---------------------------------------------------------------

n_saved = length(fit$Xi)
M_fit = length(fit$Lambda_star_by_center[[1]])
stopifnot(M_fit == M0)
stopifnot(all(fit$M == M_fit))

center_objs = vector("list", M_fit)
names(center_objs) = paste0("center_", seq_len(M_fit))

for(m in seq_len(M_fit)){
  Lambda_fit_m = lapply(seq_len(n_saved), function(it) fit$Lambda_star_by_center[[it]][[m]])
  # --> list of length n_saved; each element is a V x K_{it,m} matrix
  K_it_m = vapply(Lambda_fit_m, ncol, integer(1))
  
  topic_objs_m = lapply(seq_len(n_saved), function(it) {
    build_topic_matrices(fit$Xi[[it]][[m]], Ttot)
  })
  
  center_objs[[m]] = list(
    Lambda_fit = Lambda_fit_m,
    K_it = K_it_m,
    topic_objs = topic_objs_m
  )
}

# For example, center_objs[[m]]$Lambda_fit[[it]] contains the atoms
# allocated to center m at iteration it, and center_objs[[m]]$topic_objs[[it]]
# contains the corresponding Activity and Xi_star matrices.

# Diagnosis ---------------------------------------------------------------

it_start = 1     # da dove parto + burnin + val iniziale
it_end   = niter # dove finisco


# Traceplot Num. topics ------------------------------------------------------------------

stopifnot(it_start >= 1, it_end <= n_saved, it_start <= it_end)
iter_an = seq(it_start,it_end)

# Num. topics at each time, separately for every center.
Kt_tr_by_center = lapply(seq_len(M_fit), function(m) {
  Kt_m = vapply(iter_an, function(it) {
    rowSums(center_objs[[m]]$topic_objs[[it]]$Activity)
  }, numeric(Ttot))
  t(Kt_m)
})
names(Kt_tr_by_center) = names(center_objs)

# Same object as in the non-center case, but aggregated over centers.
Kt_tr = Reduce("+", Kt_tr_by_center)

for(m in seq_len(M_fit)){
  par(mfrow = c(3,3), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
  for(t in 1:Ttot){
    plot( iter_an, Kt_tr_by_center[[m]][,t],
          xlab = "Iter.", ylab = paste0("K",t),
          main = paste0("Center ",m),
          type = "l" )
  }
}

par(mfrow = c(3,3), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
for(t in 1:Ttot){
  plot( iter_an, Kt_tr[,t],
        xlab = "Iter.", ylab = paste0("K",t),
        main = "Total",
        type = "l" )
}

# Total number of distinct topics, separately for every center.
K_tr_by_center = sapply(seq_len(M_fit), function(m) center_objs[[m]]$K_it[iter_an])
if(M_fit == 1){
  K_tr_by_center = matrix(K_tr_by_center, ncol = 1)
}
colnames(K_tr_by_center) = names(center_objs)

# Same object as in the non-center case, but aggregated over centers.
K_tr = rowSums(K_tr_by_center)

cols_center = hcl.colors(M_fit, palette = "Dark 3")

par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
matplot( iter_an, K_tr_by_center,
         xlab = "Iter.", ylab = "K",
         type = "l", lty = 1, col = cols_center )
lines( iter_an, K_tr, lwd = 2 )
legend( "topright",
        legend = c(paste0("Center ",seq_len(M_fit)), "Total"),
        col = c(cols_center, "black"),
        lty = 1,
        lwd = c(rep(1,M_fit), 2),
        bty = "n",
        cex = 0.7 )


# Active features matrix  ------------------------------------------------------------------

make_topic_objs_total = function(it){
  Activity_list = lapply(seq_len(M_fit), function(m) center_objs[[m]]$topic_objs[[it]]$Activity)
  Xi_star_list  = lapply(seq_len(M_fit), function(m) center_objs[[m]]$topic_objs[[it]]$Xi_star)
  
  list(
    Activity = do.call(cbind, Activity_list),
    Xi_star  = do.call(rbind, Xi_star_list)
  )
}

cor_from_ximat = function(x){
  if(ncol(x) == 0)
    return(matrix(0, nrow = Ttot, ncol = Ttot))
  
  x_bar = 1/ncol(x) * t(rep(1, ncol(x)) %*% t(x))
  Cov_it = 1/ncol(x) * (x %*% t(x)) - (x_bar %*% t(x_bar))
  
  sd_it = sqrt(diag(Cov_it))
  denom = sd_it %o% sd_it
  Cor_it = Cov_it
  
  idx = denom > 0
  Cor_it[idx] = Cov_it[idx] / denom[idx]
  Cor_it[!idx] = 0
  diag(Cor_it) = ifelse(sd_it > 0, 1, 0)
  
  Cor_it
}

make_feature_summary = function(topic_objs_list){
  Zmat_list  = lapply(topic_objs_list, function(z) z$Activity)
  Ximat_list = lapply(topic_objs_list, function(z) t(z$Xi_star) )
  
  nonempty_it = which(vapply(Zmat_list, ncol, integer(1)) > 0)
  
  if(length(nonempty_it) >= 2){
    n_fangs = length(nonempty_it)
    fangs_est = fangs(Zmat_list[nonempty_it],
                      nInit = min(16,n_fangs),
                      nSweet = max(1,min(4,n_fangs - 1)),
                      quiet = TRUE)
    est = fangs_est$estimate
    Zest_ordered = left_order(est)
    Kest = ncol(est)
  } else if(length(nonempty_it) == 1){
    fangs_est = NULL
    est = Zmat_list[[nonempty_it]]
    Zest_ordered = left_order(est)
    Kest = ncol(est)
  } else {
    fangs_est = NULL
    est = matrix(0L, nrow = Ttot, ncol = 0)
    Zest_ordered = est
    Kest = 0
  }
  
  Adj_mat_list = lapply(Ximat_list, function(x) x%*%t(x) )
  Wpsm = Reduce(`+`, Adj_mat_list)/length(Adj_mat_list)
  
  Cor_mat_list = lapply(Ximat_list, cor_from_ximat)
  Cor = Reduce(`+`, Cor_mat_list)/length(Cor_mat_list)
  
  Adj_bin_mat_list = lapply(Zmat_list, function(z) z%*%t(z) )
  Upsm = Reduce(`+`, Adj_bin_mat_list)/length(Adj_bin_mat_list)
  
  list(
    Zmat_list = Zmat_list,
    Ximat_list = Ximat_list,
    fangs_est = fangs_est,
    fangs_iter = iter_an[nonempty_it],
    est = est,
    Zest_ordered = Zest_ordered,
    Kest = Kest,
    Wpsm = Wpsm,
    Cor = Cor,
    Upsm = Upsm
  )
}

plot_Z_est = function(feature_obj, main_lab){
  Kest = feature_obj$Kest
  
  if(Kest == 0){
    par(mfrow = c(1,1), mar = c(3.5,3.5,2,1), mgp=c(2,0.5,0))
    plot.new()
    title(main = paste0("<Z est> - ",main_lab))
    text(0.5,0.5,"No active features")
    return(invisible(NULL))
  }
  
  par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
  image( 1:Kest, 1:Ttot,
         t(feature_obj$Zest_ordered),
         col = mycol,
         xlab = "Topics",
         ylab = "Time",
         main = paste0("<Z est> - ",main_lab),
         axes = FALSE )
  axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)),
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(1, at = seq(1, Kest, length.out = min(Kest, 10)),
       labels = round(seq(1, Kest, length.out = min(Kest, 10))),
       cex.axis = 0.7)
  box()
}

plot_pairwise_matrix = function(A, main_lab){
  par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
  image( 1:Ttot, 1:Ttot,
         A,
         col = mycol,
         xlab = "Time",
         ylab = "Time",
         main = main_lab,
         axes = FALSE )
  axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)),
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)),
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7)
  box()
  fields::image.plot(
    1:Ttot, 1:Ttot, A,
    col = mycol,
    legend.only = TRUE,
    horizontal = FALSE,
    legend.width = 1.2,
    legend.shrink = 0.8,
    legend.mar = 8.5,
    legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
  )
}

if(exists("Presidents_all") &&
   is.data.frame(Presidents_all) &&
   nrow(Presidents_all) >= Ttot &&
   ncol(Presidents_all) >= 2){
  pres_tick_idx = unique(round(seq(1,Ttot,length.out = min(Ttot,11))))
  pres_tick_lab = substr(as.character(Presidents_all[pres_tick_idx,2]), 1, 4)
} else {
  pres_tick_idx = unique(round(seq(1,Ttot,length.out = min(Ttot,10))))
  pres_tick_lab = as.character(pres_tick_idx)
}

save_img_current = exists("save_img") && isTRUE(save_img)
if(save_img_current && !dir.exists("img"))
  dir.create("img", recursive = TRUE)

plot_cor_matrix = function(Cor, main_lab, file_lab){
  if(save_img_current)
    pdf(paste0("img/SD_Corr_matrix_",file_lab,".pdf"), width=10, height=8)
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.08))
  par(mar = c(2.4,2.75,0.6,0), mgp=c(2,0.5,0), cex = 2, las = 1)
  image( 1:Ttot, 1:Ttot,
         Cor,
         col = mycol,
         xlab = "",
         ylab = "",
         main = main_lab,
         axes = FALSE,
         asp = 1 )
  axis(1, at = pres_tick_idx, labels = FALSE, tck = -0.015)
  axis(2, at = pres_tick_idx, labels = pres_tick_lab, tck = -0.015, las = 1, cex.axis = 0.9)
  text(
    x = pres_tick_idx,
    y = par("usr")[3] - 0.012 * diff(par("usr")[3:4]),
    labels = pres_tick_lab,
    srt = 45,
    adj = 1,
    xpd = NA,
    cex = 0.9
  )
  par(mar = c(2,0,0,0.5))
  fields::image.plot(
    1:Ttot, 1:Ttot, Cor,
    col = mycol,
    legend.only = TRUE,
    horizontal = FALSE,
    legend.width = 0.9,
    legend.shrink = 0.78,
    legend.mar = 2.2,
    legend.args = list(text = " ", side = 3, line = 0.2, cex = 0.8)
  )
  layout(1)
  
  if(save_img_current)
    dev.off()
}

topic_objs_total = lapply(iter_an, make_topic_objs_total)

features_by_center = lapply(seq_len(M_fit), function(m) {
  make_feature_summary(center_objs[[m]]$topic_objs[iter_an])
})
names(features_by_center) = names(center_objs)

features_total = make_feature_summary(topic_objs_total)

for(m in seq_len(M_fit)){
  plot_Z_est(features_by_center[[m]], paste0("Center ",m))
}
plot_Z_est(features_total, "Total")


# Weighted pairwise similarity matrix ------------------------------------------------------------------

for(m in seq_len(M_fit)){
  plot_pairwise_matrix(features_by_center[[m]]$Wpsm, paste0("Wpsm - Center ",m))
}
plot_pairwise_matrix(features_total$Wpsm, "Wpsm - Total")


# Wpsm as correlation ------------------------------------------------------------------

for(m in seq_len(M_fit)){
  plot_cor_matrix(features_by_center[[m]]$Cor,
                  paste0("Correlation - Center ",m),
                  paste0("center_",m))
}
plot_cor_matrix(features_total$Cor, "Correlation - Total", "total")


# Unweighted pairwise similarity matrix ------------------------------------------------------------------

for(m in seq_len(M_fit)){
  plot_pairwise_matrix(features_by_center[[m]]$Upsm, paste0("Upsm - Center ",m))
}
plot_pairwise_matrix(features_total$Upsm, "Upsm - Total")


# Monte Carlo mean of the empirical centers --------------------------------

for(m in seq_len(M_fit)){
  nonempty_it = iter_an[center_objs[[m]]$K_it[iter_an] > 0]
  
  if(length(nonempty_it) > 0){
    Lambda_center_it = lapply(nonempty_it, function(it) {
      rowMeans(center_objs[[m]]$Lambda_fit[[it]])
    })
    
    center_objs[[m]]$Lambda_mean = Reduce("+", Lambda_center_it)/length(Lambda_center_it)
    center_objs[[m]]$n_nonempty = length(nonempty_it)
  } else {
    center_objs[[m]]$Lambda_mean = rep(NA_real_, V)
    center_objs[[m]]$n_nonempty = 0
  }
}

ymax_center_mean = max(unlist(lapply(center_objs, function(x) x$Lambda_mean)), na.rm = TRUE)*1.2
if(is.finite(ymax_center_mean)){
  for(m in seq_len(M_fit)){
    if(center_objs[[m]]$n_nonempty > 0){
      par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
      barplot( height = unname(center_objs[[m]]$Lambda_mean),
               names.arg = as.character(1:V),
               las = 1, col = "darkred", border = NA,
               xlab = "Word", ylab = "Prob.",
               main = paste0("MC mean center ", m),
               ylim = c(0,ymax_center_mean),
               cex.names = 0.5 )
    }
  }
}


