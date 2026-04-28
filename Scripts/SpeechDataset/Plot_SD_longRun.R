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

mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)
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

i = 1
for(i in seq_along(fit_files)) {
  fit_file = fit_files[i]
  fit_name = basename(fit_file)
  summary_file = file.path(
    save_summary_dir,
    sub("\\.rds$", "_summary.rds", fit_name)
  )
  
  cat("\n[", i, "/", length(fit_files), "] Reading ", fit_name, " ...\n", sep = "")
  fit_summary = readRDS(fit_file)

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
  ## Mean values --------------------------------------------------------------------
  
  meanRes <- matrix(0, nrow = Ttot, ncol = V)
  for(it in it_start:it_end ) {
    
    meanRes = meanRes + t(Lambda_fit[[it]]%*%topic_objs[[it]]$Xi_star)
    
  }
  meanRes <- meanRes / niter
  dim(meanRes)
  
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






