# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[4] # <--- modify here
wd = paste0(choose_wd,"DynamicFeatureAllocation/Scripts/SpeechDataset")
setwd(wd)

# Functions ---------------------------------------------------------------
source("./../../R/Rfunctions.R")
Rcpp::sourceCpp("./../../src/RcppFunctions.cpp")

library(parallel)

avail_cores = parallel::detectCores(logical = TRUE)
if(is.na(avail_cores))
  avail_cores = 1L
n_cores = 4 # <---

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

delta_all   = c(1e-4,1)  
beta0_all   = c(1e-2,1e-1,1)
gamma0_all  = c(1e-2,1e-1)
sigma0_all  = c(0.1,0.5,0.9)

params_grid = expand.grid(
  delta = delta_all,
  beta0 = beta0_all,
  gamma0 = gamma0_all,
  sigma0 = sigma0_all,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

params_grid = matrix(0,nrow = 3, ncol = 4)
params_grid[1,] = c(1,0.1,0.1,0.1)
params_grid[2,] = c(1e-3,1,0.1,0.1)
params_grid[3,] = c(1e-2,1,1e-2,0.9)
params_grid = as.data.frame(params_grid, stringsAsFactors = FALSE)
colnames(params_grid) = c("delta","beta0","gamma0","sigma0")

cat("\n ---- Number of configurations to run: ",nrow(params_grid)," ---- \n")
# Parallel MCMC runner ----------------------------------------------------

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

save_dir = file.path(wd, "save")
if(!dir.exists(save_dir))
  stop("save folder does not exist: ", save_dir)
save_summary_dir = file.path(wd, "save_summary")
if(!dir.exists(save_summary_dir))
  dir.create(save_summary_dir, recursive = TRUE, showWarnings = FALSE)

fit_files = list.files(save_dir, pattern = "\\.rds$", full.names = TRUE)

if(length(fit_files) == 0)
  stop("No fit .rds files found in: ", save_dir)

cat("\nFound ", length(fit_files), " fit object(s) to summarize.\n", sep = "")

summarize_fit = function(fit, H, Ttot, niter_expected = NULL) {
  n_saved = length(fit$Xi)
  if(n_saved == 0)
    stop("fit$Xi is empty")
  
  it_start = 1
  it_end = n_saved
  
  Lambda_fit = lapply(fit$Lambda_star_mcmc, function(Lam_it) t(do.call(rbind, Lam_it)))
  K_it = sapply(Lambda_fit, ncol)
  topic_objs = lapply(seq_len(n_saved), function(it) build_topic_matrices(fit$Xi[[it]], Ttot))
  
  V = if(length(Lambda_fit) > 0) nrow(Lambda_fit[[1]]) else NA_integer_
  meanRes = matrix(0, nrow = Ttot, ncol = V)
  for(it in it_start:it_end) {
    meanRes = meanRes + t(Lambda_fit[[it]] %*% topic_objs[[it]]$Xi_star)
  }
  meanRes = meanRes / length(it_start:it_end)
  
  Xi_mean = Reduce("+", fit$Xi[it_start:it_end]) / length(fit$Xi[it_start:it_end])
  S_mean = Reduce("+", fit$S[it_start:it_end]) / length(fit$S[it_start:it_end])
  N_mean = Reduce("+", fit$N[it_start:it_end]) / length(fit$N[it_start:it_end])
  
  Kt_tr = sapply(fit$Xi, function(Xi_it) {
    apply(Xi_it, 2, function(x) length(which(x > 0)))
  })
  Kt_tr = t(Kt_tr)
  
  list(
    meta = list(
      n_saved = n_saved,
      H = H,
      Ttot = Ttot,
      V = V,
      niter_expected = niter_expected
    ),
    Lambda_fit = Lambda_fit,
    topic_objs = topic_objs,
    K_it = K_it,
    meanRes = meanRes,
    Xi_mean = Xi_mean,
    S_mean = S_mean,
    N_mean = N_mean,
    Kt_tr = Kt_tr
  )
}

summary_info = vector("list", length(fit_files))

for(i in seq_along(fit_files)) {
  fit_file = fit_files[i]
  fit_name = basename(fit_file)
  summary_file = file.path(
    save_summary_dir,
    sub("\\.rds$", "_summary.rds", fit_name)
  )
  
  cat("\n[", i, "/", length(fit_files), "] Loading ", fit_name, " ...\n", sep = "")
  fit = readRDS(fit_file)
  
  fit_summary = summarize_fit(fit, H = H, Ttot = Ttot, niter_expected = niter)
  saveRDS(fit_summary, file = summary_file)
  
  summary_size_mb = file.info(summary_file)$size / 1024^2
  under_limit = summary_size_mb < 100
  
  cat("Saved summary: ", basename(summary_file), "\n", sep = "")
  cat("Size: ", round(summary_size_mb, 2), " MB\n", sep = "")
  cat("Is size < 100 MB? ", under_limit, "\n", sep = "")
  
  summary_info[[i]] = data.frame(
    fit_file = fit_name,
    summary_file = basename(summary_file),
    size_mb = round(summary_size_mb, 2),
    less_than_100mb = under_limit,
    stringsAsFactors = FALSE
  )
  
  rm(fit, fit_summary)
  gc(verbose = FALSE)
}

summary_info = do.call(rbind, summary_info)
write.csv(
  summary_info,
  file = file.path(save_summary_dir, paste0("summary_sizes_r", r, ".csv")),
  row.names = FALSE
)

cat("\nSummary size check:\n")
print(summary_info)
