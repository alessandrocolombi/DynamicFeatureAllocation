# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"DynamicFeatureAllocation/Scripts/SpeechDataset")
setwd(wd)

# Functions ---------------------------------------------------------------
source("./../../R/Rfunctions.R")
Rcpp::sourceCpp("./../../src/RcppFunctions.cpp")
library(fangs)

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
left_order_pair <- function(Zstar, Xi_star) {
  if (!is.matrix(Zstar)) Zstar <- as.matrix(Zstar)
  if (!is.matrix(Xi_star)) Xi_star <- as.matrix(Xi_star)
  
  if (!all(Zstar %in% c(0, 1)))
    stop("Zstar must be a binary matrix (0/1).")
  
  if (!all(dim(Zstar) == dim(Xi_star)))
    stop("Zstar and Xi_star must have the same dimensions.")
  
  # Binary signature of each column of Zstar
  col_signature <- apply(Zstar, 2, function(col) {
    paste(col, collapse = "")
  })
  
  # Lexicographic decreasing order
  ord <- order(col_signature, decreasing = TRUE)
  
  list(
    Zord = Zstar[, ord, drop = FALSE],
    Xiord = Xi_star[, ord, drop = FALSE],
    ord = ord
  )
}


mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)
tau_xi_plot = 0.5
# Read data -----------------------------------------------------------

seed = 132332
set.seed(seed)
r = 10

# Read President names
Presidents_all <- read.csv("../data/Presidents_all.csv")
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


# Set options list -------------------------------------------------------------


params_grid = matrix(0,nrow = 3, ncol = 4)
params_grid[1,] = c(1,0.1,0.1,0.1)
params_grid[2,] = c(1e-3,1,0.1,0.1)
params_grid[3,] = c(1e-2,1,1e-2,0.9)
params_grid = as.data.frame(params_grid, stringsAsFactors = FALSE)
colnames(params_grid) = c("delta","beta0","gamma0","sigma0")

cat("\n ---- Number of configurations to run: ",nrow(params_grid)," ---- \n")

save_all_chain = TRUE
seed = 22123
H = 20 # number of atoms
niter = 10000
nburn = 10000
thin  = 20

# Fixed hyperparameters / MCMC settings
a_phi = 1; b_phi = 1
a_sigma = 1; b_sigma = 1
prop_var_phi = 0.01

phi0 = 1

UpdateDitl = TRUE; UpdateS = TRUE
UpdateLambda = TRUE; UpdateXi = TRUE
UpdateU = TRUE
UpdatePhi = FALSE; UpdateGamma = FALSE
UpdateSigma = FALSE; UpdateBeta = FALSE
print = TRUE; JointAdp = FALSE

save_summary_dir = file.path(wd, "save_summary")
fit_files = list.files(save_summary_dir, pattern = "\\.rds$", full.names = TRUE)

if(length(fit_files) == 0)
  stop("No fit .rds files found in: ", save_dir)

cat("\nFound ", length(fit_files), " fit object(s) to summarize.\n", sep = "")


# Plot loop ---------------------------------------------------------------

i = 1
for(i in seq_along(fit_files)) {
  fit_file = fit_files[i]
  fit_name = basename(fit_file)
  summary_file = file.path(
    save_summary_dir,
    sub("\\.rds$", "_summary.rds", fit_name)
  )
  
  # cat("\n[", i, "/", length(fit_files), "] Reading ", fit_name, " ...\n", sep = "")
  fit_summary = readRDS(paste0("save_summary/cfg001_r10_delta_1e00_beta0_1em01_gamma0_1em01_sigma0_1em01_summary.rds"))
  fit_all = readRDS(paste0("save/cfg001_r10_delta_1e00_beta0_1em01_gamma0_1em01_sigma0_1em01.rds"))
  ## Plot Xi mean ------------------------------------------------------------
  Xi_mean = fit_summary$Xi_mean
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
  ## Plot S mean ------------------------------------------------------------
  S_mean = fit_summary$S_mean
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
  
  ## Plot N mean ------------------------------------------------------------
  N_mean <- fit_summary$N_mean
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
  ## Traceplot Num. topics ------------------------------------------------------------------
  
  # Num. topics at each time 
  Kt_tr = fit_summary$Kt_tr
  
  par(mfrow = c(3,3), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
  for(t in 1:Ttot){
    plot( Kt_tr[,t], xlab = "Iter.", ylab = paste0("K",t), type = "l" )
  }
  
  # Total number of distinct topics
  K_tr = fit_summary$K_it
  par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
  plot( K_tr, xlab = "Iter.", ylab = paste0("K"), type = "l" )
  
  ## Active features matrix  ------------------------------------------------------------------
  
  Zmat_list  = lapply(fit_summary$topic_objs, function(z) z$Activity)
  Ximat_list = lapply(fit_summary$topic_objs, function(z) t(z$Xi_star) )
  
  fangs_est = fangs(Zmat_list)
  est = fangs_est$estimate
  Zest_ordered = left_order(est)
  Kest = ncol(est)
  
  par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
  image( 1:Kest, 1:Ttot, 
         t(Zest_ordered),   
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
  ## Weighted pairwise similarity matrix ------------------------------------------------------------------
  
  Adj_mat_list = lapply(Ximat_list, function(x) x%*%t(x) )
  
  Wpsm <- Reduce(`+`, Adj_mat_list)/length(Adj_mat_list)
  
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
  
  ## Wpsm as correlation ------------------------------------------------------------------
  Cor_mat_list = lapply(Ximat_list, function(x){
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
  })
  Cor <- Reduce(`+`, Cor_mat_list)/length(Cor_mat_list)
  
  par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
  image( 1:Ttot, 1:Ttot, 
         Cor,   
         col = mycol,    
         xlab = "Time", 
         ylab = "Time",
         main = "Cor",
         axes = FALSE )
  axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7)
  box()
  fields::image.plot(
    1:Ttot, 1:Ttot, Cor,
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
  ## Mean values --------------------------------------------------------------------
  
  meanRes <- fit_summary$meanRes
  
  par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
  image( 1:Ttot, 1:V, 
         meanRes,   
         col = mycol,    
         xlab = "Time", 
         ylab = "Words",
         main = paste0("< Estimated Mean >"),
         axes = FALSE )
  axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(2, at = seq(1, V, length.out = min(H, 10)), 
       labels = round(seq(1, V, length.out = min(H, 10))),
       cex.axis = 0.7)
  box()
  fields::image.plot(
    1:Ttot, 1:V, meanRes,
    col = mycol,
    legend.only = TRUE,
    horizontal = FALSE,
    legend.width = 1.2,            # controls legend thickness
    legend.shrink = 0.8,           # smaller legend
    legend.mar = 8.5,                # margin from image
    legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
  )
  ## Lambda values --------------------------------------------------------------------
  it = 950
  ymax_topic = max(Lambda_fit[[it]])
  ymax_xi = max( fit_summary$topic_objs[[it]]$Xi_star )
  for(l in 1:nrow(fit_summary$topic_objs[[it]]$Xi_star)){
    par(mfrow = c(1,2), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
    barplot( height = Lambda_fit[[it]][,l], 
             names.arg = as.character(1:V),
             las = 1, col = "darkred", border = NA,
             xlab = "Word",
             main = paste0("Atom: ",l), ylab = "Prob.", ylim = c(0,ymax_topic),
             cex.names = 0.5 )
    barplot( height = fit_summary$topic_objs[[it]]$Xi_star[l,], 
             names.arg = as.character(1:Ttot),
             las = 1, col = "darkgreen", border = NA,
             xlab = "Time",
             main = "", ylab = "Intensity", ylim = c(0,ymax_xi),
             cex.names = 0.5 )
  }
  
  ## End: start reading new file
  
}









# Brutta ------------------------------------------------------------------

Ximat_list = lapply(fit_summary$topic_objs, function(z) t(z$Xi_star) )

t = 10
for(t in c(10000-10,10000)){
  plot_mat = Ximat_list[[t]]
  plot_mat[which(plot_mat < 10)] = 0
  par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
  image( 1:Ttot, 1:ncol(Ximat_list[[t]]),
         plot_mat,   
         col = c("blue",mycol),    
         ylab = "Topics", 
         xlab = "Time",
         main = paste0("Xi_star - iter. = ",t),
         axes = FALSE )
  axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(2, at = seq(1, ncol(Ximat_list[[t]]), length.out = min(ncol(Ximat_list[[t]]), 10)), 
       labels = round(seq(1, ncol(Ximat_list[[t]]), length.out = min(ncol(Ximat_list[[t]]), 10))),
       cex.axis = 0.7)
  box()
  fields::image.plot(
    1:Ttot, 1:ncol(Ximat_list[[t]]),
    plot_mat,
    col = mycol,
    legend.only = TRUE,
    horizontal = FALSE,
    legend.width = 1.2,            # controls legend thickness
    legend.shrink = 0.8,           # smaller legend
    legend.mar = 8.5,                # margin from image
    legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
  )
}


A = matrix(rpois(n = 20, 3), nrow = 4, ncol = 5)
A
A[which(A < 2)] = -1
A


# Brutta pair ordering ----------------------------------------------------


t = 10000
tau_xi_plot = 0.5
for(t in seq(niter-20,niter,by=1) ){
  Zstar = fit_summary$topic_objs[[t]]$Activity
  Xistar = t(fit_summary$topic_objs[[t]]$Xi_star)
  Z_Xi_ord = left_order_pair(Zstar,Xistar)
  Ktot_it = ncol(Z_Xi_ord$Zord)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.18))
  par(mar = c(3.5,3.5,2,2), mgp=c(2,0.5,0))
  image( 1:Ktot_it, 1:Ttot,
         t(Z_Xi_ord$Zord),   
         col = c(mycol),    
         xlab = "Topics", 
         ylab = "Time",
         main = paste0("Zord - iter. = ",t),
         axes = FALSE )
  axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(1, at = seq(1, Ktot_it, length.out = min(Ktot_it, 10)), 
       labels = round(seq(1, Ktot_it, length.out = min(Ktot_it, 10))),
       cex.axis = 0.7)
  box()
  plot_mat = Z_Xi_ord$Xiord
  max_plot_mat = max(plot_mat, na.rm = TRUE)
  par(mar = c(3.5,3.5,2,2), mgp=c(2,0.5,0))
  if(max_plot_mat <= tau_xi_plot)
    stop("Error, max_plot_mat can not be larger than tau_xi_plot")
  
  eps_break = .Machine$double.eps
  green_breaks = seq(tau_xi_plot, max_plot_mat + eps_break, length.out = length(mycol) + 1)
  breaks_xi = c(0, tau_xi_plot, green_breaks[-1])
  image( 1:Ktot_it, 1:Ttot,
           t(plot_mat),
           col = c("blue", mycol),
           breaks = breaks_xi,
           xlab = "Topics",
           ylab = "Time",
           main = paste0("Xiord - iter. = ",t),
           axes = FALSE )
  axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(1, at = seq(1, Ktot_it, length.out = min(Ktot_it, 10)), 
       labels = round(seq(1, Ktot_it, length.out = min(Ktot_it, 10))),
       cex.axis = 0.7)
  box()
}


# President over time topics --------------------------------------------------------

# t = 60 # Trump
t = 57
Kmax = max(K_tr)
Xi_Pres = lapply(fit_summary$topic_objs, function(x) x$Xi_star[, t])
Xi_Pres = lapply(Xi_Pres, function(x) {
  x = sort(x, decreasing = TRUE)
  if(length(x) < Kmax)
    x = c(x, rep(0, Kmax - length(x)))
  x
})
Xi_Pres_mat = do.call(rbind, Xi_Pres)

soglia = 0.1 * nrow(Xi_Pres_mat)
sel_col = which( colSums(Xi_Pres_mat) > soglia )
Xi_plot = Xi_Pres_mat[,sel_col]

par(mfrow = c(1,1), mar = c(3.5,3.5,2,2), mgp = c(2,0.5,0))
matplot( Xi_plot,type = "l", 
         lty = 1, lwd = 1,
         xlab = "Iteration", ylab = paste0("Xi[, ", t, "]"),
         main = "Ordered topic intensities" )

# Xi mass ordering --------------------------------------------------------

compromise_order_pair <- function(Zstar, Xi_star) {
  if (!is.matrix(Zstar)) Zstar <- as.matrix(Zstar)
  if (!is.matrix(Xi_star)) Xi_star <- as.matrix(Xi_star)
  
  if (!all(dim(Zstar) == dim(Xi_star)))
    stop("Zstar and Xi_star must have the same dimensions.")
  
  birth_time <- apply(Zstar, 2, function(col) {
    idx <- which(col > 0)
    if(length(idx) == 0) Inf else idx[1]
  })
  xi_mass <- colSums(Xi_star)
  xi_keys <- as.data.frame(-t(Xi_star))
  ord <- do.call(order, c(list(birth_time, -xi_mass), xi_keys))
  
  list(
    Zord = Zstar[, ord, drop = FALSE],
    Xiord = Xi_star[, ord, drop = FALSE],
    birth = birth_time[ord],
    mass = xi_mass[ord],
    ord = ord
  )
}
pad_ncol_right <- function(X, ncol_target) {
  if (!is.matrix(X)) X <- as.matrix(X)
  
  if (ncol(X) > ncol_target)
    stop("ncol_target must be at least ncol(X).")
  
  if (ncol(X) == ncol_target)
    return(X)
  
  cbind(X, matrix(0, nrow = nrow(X), ncol = ncol_target - ncol(X)))
}
tail_mass_matrix <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)
  
  apply(X, 2, function(col) rev(cumsum(rev(col))))
}

max_K_tr = max(K_tr)
Xiord_comp_list = lapply(fit_summary$topic_objs, function(obj) {
  Zstar = obj$Activity
  Xistar = t(obj$Xi_star)
  compromise_order_pair(Zstar, Xistar)$Xiord
})
Xiord_comp_pad_list = lapply(Xiord_comp_list, function(Xiord_it) {
  pad_ncol_right(Xiord_it, max_K_tr)
})
Xi_star_mean = Reduce(`+`, Xiord_comp_pad_list) / length(Xiord_comp_pad_list)
CumXi_comp_pad_list = lapply(Xiord_comp_pad_list, tail_mass_matrix)
CumXi_star_mean = Reduce(`+`, CumXi_comp_pad_list) / length(CumXi_comp_pad_list)

soglia = 0
sel_colums = which(colSums(Xi_star_mean) > soglia )
plot_mat = Xi_star_mean[,sel_colums]

max_plot_mat = max(plot_mat, na.rm = TRUE)
par(mfrow = c(1,1), mar = c(3.5,3.5,2,2), mgp = c(2,0.5,0))
green_breaks = seq(0, max_plot_mat, length.out = length(mycol)+1)
breaks_xi = green_breaks
image(1:ncol(plot_mat), 1:Ttot,
      t(plot_mat),
      col = mycol,
      breaks = breaks_xi,
      xlab = "Topics",
      ylab = "Time",
      main = "Xi_star_mean compromise-ord",
      axes = FALSE)
axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)),
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7)
axis(1, at = seq(1, ncol(plot_mat), length.out = min(ncol(plot_mat), 10)),
     labels = round(seq(1, ncol(plot_mat), length.out = min(ncol(plot_mat), 10))),
     cex.axis = 0.7)
box()

sel_colums = which(colSums(CumXi_star_mean) > soglia )
plot_mat = CumXi_star_mean[,sel_colums]

max_plot_mat = max(plot_mat, na.rm = TRUE)
par(mfrow = c(1,1), mar = c(3.5,3.5,2,2), mgp = c(2,0.5,0))
green_breaks = seq(0, max_plot_mat, length.out = length(mycol)+1)
breaks_xi = green_breaks
image(1:ncol(plot_mat), 1:Ttot,
      t(plot_mat),
      col = mycol,
      breaks = breaks_xi,
      xlab = "Topics",
      ylab = "Time",
      main = "CumXi_star_mean compromise-ord",
      axes = FALSE)
axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)),
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7)
axis(1, at = seq(1, ncol(plot_mat), length.out = min(ncol(plot_mat), 10)),
     labels = round(seq(1, ncol(plot_mat), length.out = min(ncol(plot_mat), 10))),
     cex.axis = 0.7)
box()





for(t in seq(niter-10,niter,by=1) ){
  Zstar = fit_summary$topic_objs[[t]]$Activity
  Xistar = t(fit_summary$topic_objs[[t]]$Xi_star)
  Z_Xi_comp = compromise_order_pair(Zstar, Xistar)
  Ktot_it = ncol(Z_Xi_comp$Zord)
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))
  par(mar = c(3.5,3.5,2,2), mgp = c(2,0.5,0))
  image(1:Ktot_it, 1:Ttot,
        t(Z_Xi_comp$Zord),
        col = c(mycol),
        xlab = "Topics",
        ylab = "Time",
        main = paste0("Z compromise-ord - iter. = ", t),
        axes = FALSE)
  axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)),
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7)
  axis(1, at = seq(1, Ktot_it, length.out = min(Ktot_it, 10)),
       labels = round(seq(1, Ktot_it, length.out = min(Ktot_it, 10))),
       cex.axis = 0.7)
  box()
  
  plot_mat = Z_Xi_comp$Xiord
  max_plot_mat = max(plot_mat, na.rm = TRUE)
  par(mar = c(3.5,3.5,2,2), mgp = c(2,0.5,0))
  if(max_plot_mat <= tau_xi_plot) {
    image(1:Ktot_it, 1:Ttot,
          t(plot_mat),
          col = "blue",
          zlim = c(0, 1),
          xlab = "Topics",
          ylab = "Time",
          main = paste0("Xi compromise-ord - iter. = ", t),
          axes = FALSE)
  } else {
    eps_break = .Machine$double.eps
    green_breaks = seq(tau_xi_plot, max_plot_mat + eps_break, length.out = length(mycol) + 1)
    breaks_xi = c(0, tau_xi_plot, green_breaks[-1])
    image(1:Ktot_it, 1:Ttot,
          t(plot_mat),
          col = c("blue", mycol),
          breaks = breaks_xi,
          xlab = "Topics",
          ylab = "Time",
          main = paste0("Xi compromise-ord - iter. = ", t),
          axes = FALSE)
  }
  axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)),
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7)
  axis(1, at = seq(1, Ktot_it, length.out = min(Ktot_it, 10)),
       labels = round(seq(1, Ktot_it, length.out = min(Ktot_it, 10))),
       cex.axis = 0.7)
  box()
}



# Read full object --------------------------------------------------------
fit = readRDS("save/cfg001_r10_delta_1e00_beta0_1em01_gamma0_1em01_sigma0_1em01.rds")
View(fit)




