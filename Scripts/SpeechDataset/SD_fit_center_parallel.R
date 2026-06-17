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
source("./../../R/Rfunctions_centers.R")
Rcpp::sourceCpp("./../../src/RcppFunctions_centers.cpp")

library(parallel)

avail_cores = parallel::detectCores(logical = TRUE)
if(is.na(avail_cores))
  avail_cores = 1L
n_cores = min(4, avail_cores) # <---

# Read data ---------------------------------------------------------------

seed_data = 132332
set.seed(seed_data)
r = 10

Presidents_all <- read.csv("../data/Presidents_all.csv")
Presidents_all[,2]
data = read.table(paste0("../data/SpeechData_top",r,".txt"))
colnames(data) = Presidents_all[,2]
data = as.matrix(data)
# data = data[,41:60]

V = nrow(data)
Ttot = ncol(data)

# Model and tuning options ------------------------------------------------

seed = 22123

# Fixed dimensions of the dynamic model with centers.
H = 10
M0 = 3

# Static Poisson-NMF initialization.
static_nstart = 20
static_niter = 500
eps_init = 1e-8
zeta_floor = 0.2
zeta_strength = 30

# Grid: quantities that vary across parallel runs.
gamma_all = c(0.001,0.01,0.1,1,10)
delta0_centers_all = c(0.001,0.1,1,10)

params_grid = expand.grid(
  gamma = gamma_all,
  delta0_centers = delta0_centers_all,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Fixed process hyperparameters.
# The C++ sampler uses gamma/M internally because M is fixed here.
sigma = 0.1
beta  = 0.1

# Fixed center-clustering hyperparameters.
a_phi_centers = 1
b_phi_centers = 1
omega = 1
a_omega = 1
b_omega = 1
var_phi_centers = 0.01
var_delta_centers = 0.01
mstar_max = 0

# Update flags.
UpdateDitl = TRUE
UpdateS = TRUE
UpdateLambda = TRUE
UpdateXi = TRUE
UpdateU = TRUE
UpdateCenters = TRUE
UpdateOmega = FALSE
print = TRUE

# MCMC control.
niter = 100#5000
nburn = 10#5000
thin  = 1#10

# Output control.
save_all_chain = TRUE
save_dir = file.path(wd, "centers_save")
log_dir = file.path(save_dir, "logs")
output_dir = save_dir

if(identical(Sys.getenv("SD_CENTER_PARALLEL_SMOKE"), "1")){
  n_cores = 1
  params_grid = params_grid[1,,drop = FALSE]
  static_nstart = 2
  static_niter = 5
  niter = 2
  nburn = 0
  thin = 1
  save_all_chain = FALSE
  save_dir = file.path(tempdir(), "centers_parallel_smoke")
  log_dir = file.path(save_dir, "logs")
  output_dir = save_dir
}

if(!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if(!dir.exists(log_dir))
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n ---- Number of configurations to run: ",nrow(params_grid)," ---- \n")

# Static initialization helpers ------------------------------------------

fit_static_poisson_nmf = function(D, K, nstart = 10, niter = 500,
                                  eps = 1e-8, seed = NULL){
  if(!is.null(seed))
    set.seed(seed)
  
  V = nrow(D)
  Ttot = ncol(D)
  doc_totals = colSums(D)
  
  best_loss = Inf
  best_Lambda = NULL
  best_Xi = NULL
  
  for(s in seq_len(nstart)){
    Lambda = matrix(rexp(V*K, rate = 1), nrow = V, ncol = K)
    Lambda = sweep(Lambda, 2, colSums(Lambda), "/")
    
    Xi = matrix(runif(K*Ttot, min = 0.5, max = 1.5), nrow = K, ncol = Ttot)
    Xi = sweep(Xi, 2, colSums(Xi), "/")
    Xi = sweep(Xi, 2, pmax(doc_totals, eps), "*")
    
    for(iter in seq_len(niter)){
      Mean = Lambda %*% Xi + eps
      Xi = Xi * (t(Lambda) %*% (D/Mean)) /
        matrix(colSums(Lambda), nrow = K, ncol = Ttot)
      Xi = pmax(Xi, eps)
      
      Mean = Lambda %*% Xi + eps
      Lambda = Lambda * ((D/Mean) %*% t(Xi)) /
        matrix(rowSums(Xi), nrow = V, ncol = K, byrow = TRUE)
      Lambda = pmax(Lambda, eps)
      
      lambda_scale = colSums(Lambda)
      Lambda = sweep(Lambda, 2, lambda_scale, "/")
      Xi = sweep(Xi, 1, lambda_scale, "*")
    }
    
    Mean = Lambda %*% Xi + eps
    loss_mat = Mean - D
    idx_pos = D > 0
    loss_mat[idx_pos] = loss_mat[idx_pos] + D[idx_pos]*log(D[idx_pos]/Mean[idx_pos])
    loss = sum(loss_mat)
    
    if(loss < best_loss){
      best_loss = loss
      best_Lambda = Lambda
      best_Xi = Xi
    }
  }
  
  center_time = apply(best_Xi, 1, function(x) weighted.mean(seq_len(Ttot), x + eps))
  ord = order(center_time)
  
  list(
    Lambda = best_Lambda[,ord,drop = FALSE],
    Xi = best_Xi[ord,,drop = FALSE],
    loss = best_loss,
    center_time = center_time[ord]
  )
}

build_center_init = function(static_init, H, Ttot, zeta_floor, zeta_strength){
  M0 = nrow(static_init$Xi)
  center_profiles = static_init$Lambda
  Xi_static = static_init$Xi
  
  Xi0 = vector("list", M0)
  S0 = vector("list", M0)
  Lambda0 = vector("list", M0)
  Zeta0 = vector("list", M0)
  
  for(m in seq_len(M0)){
    Xi0[[m]] = matrix(0L, nrow = H, ncol = Ttot)
    Xi0[[m]][1,] = as.integer(round(Xi_static[m,]))
    
    S0[[m]] = matrix(rgamma(n = H*Ttot, shape = 1, rate = 1),
                     nrow = H, ncol = Ttot)
    
    Zeta0[[m]] = zeta_floor + zeta_strength * center_profiles[,m]
    
    Lambda0[[m]] = vector("list", H)
    for(l in seq_len(H)){
      Lambda0[[m]][[l]] = replicate(Ttot, rdirichlet_vec(Zeta0[[m]]))
    }
    
    for(t in seq_len(Ttot)){
      Lambda0[[m]][[1]][,t] = center_profiles[,m]
    }
  }
  
  set_init_DTM_centers(Xi0,Lambda0,S0,Zeta0)
}

format_tag_value = function(x){
  out = format(signif(x, 6), scientific = TRUE, trim = TRUE)
  out = gsub("\\+","",out)
  out = gsub("-","m",out)
  out = gsub("\\.","p",out)
  out
}

make_config_tag = function(cfg, cfg_id, r, M0, H){
  paste0(
    "cfg",sprintf("%03d",cfg_id),
    "_center_r",r,
    "_M",M0,
    "_H",H,
    "_gamma_",format_tag_value(cfg$gamma),
    "_delta0_",format_tag_value(cfg$delta0_centers)
  )
}

# One shared static initialization for all configurations.
static_init = fit_static_poisson_nmf(data, M0,
                                     nstart = static_nstart,
                                     niter = static_niter,
                                     eps = eps_init,
                                     seed = seed)

cat("Static Poisson-NMF center initialization completed\n")
cat("KL objective:", signif(static_init$loss, 5), "\n")
cat("Center times:", paste(round(static_init$center_time, 2), collapse = ", "), "\n")
cat("Center masses:", paste(round(rowSums(static_init$Xi), 2), collapse = ", "), "\n")

# Parallel MCMC runner ----------------------------------------------------

project_dir = normalizePath(file.path(wd, "..", ".."), winslash = "/", mustWork = TRUE)
rfunctions_path = file.path(project_dir, "R", "Rfunctions_centers.R")
rcpp_path = file.path(project_dir, "src", "RcppFunctions_centers.cpp")

run_single_config = function(cfg, cfg_id, data, H, Ttot, r, M0,
                              static_init, zeta_floor, zeta_strength,
                              sigma, beta,
                              a_phi_centers, b_phi_centers,
                              omega, a_omega, b_omega,
                              var_phi_centers, var_delta_centers, mstar_max,
                              UpdateDitl, UpdateS, UpdateLambda, UpdateXi, UpdateU,
                              UpdateCenters, UpdateOmega,
                              print, seed,
                              niter, nburn, thin,
                              output_dir, log_dir, save_dir, save_all_chain,
                              static_nstart, static_niter, eps_init){
  cfg = as.list(as.data.frame(cfg, stringsAsFactors = FALSE))
  cfg$gamma = as.numeric(cfg$gamma)
  cfg$delta0_centers = as.numeric(cfg$delta0_centers)
  
  if(any(!is.finite(c(cfg$gamma, cfg$delta0_centers))) ||
     cfg$gamma <= 0 || cfg$delta0_centers <= 0)
    stop("Invalid configuration: gamma and delta0_centers must be positive finite scalars")
  
  tag = make_config_tag(cfg, cfg_id, r, M0, H)
  pdf_file = file.path(output_dir, paste0(tag, ".pdf"))
  log_file = file.path(log_dir, paste0(tag, ".log"))
  fit_file = file.path(save_dir, paste0(tag, ".rds"))
  setup_file = file.path(save_dir, paste0(tag, "_setup.rds"))
  chain_seed = seed + cfg_id
  
  log_open = FALSE
  pdf_open = FALSE
  on.exit({
    if(pdf_open)
      try(grDevices::dev.off(), silent = TRUE)
    if(log_open)
      try(sink(), silent = TRUE)
  }, add = TRUE)
  
  result = tryCatch({
    set.seed(chain_seed)
    
    init_DTM_centers = build_center_init(static_init, H, Ttot,
                                         zeta_floor, zeta_strength)
    
    param_DTM_centers = set_param_DTM_centers(H,cfg$gamma,sigma,beta,
                                              a_phi_centers,b_phi_centers,cfg$delta0_centers,
                                              omega,a_omega,b_omega,
                                              var_phi_centers,var_delta_centers,mstar_max,
                                              UpdateDitl,UpdateS,UpdateLambda,UpdateXi,UpdateU,
                                              UpdateCenters,UpdateOmega,
                                              chain_seed,print)
    
    tuning_options = list(
      r = r,
      seed = chain_seed,
      H = H,
      M0 = M0,
      static_nstart = static_nstart,
      static_niter = static_niter,
      eps_init = eps_init,
      zeta_floor = zeta_floor,
      zeta_strength = zeta_strength,
      gamma = cfg$gamma,
      sigma = sigma,
      beta = beta,
      a_phi_centers = a_phi_centers,
      b_phi_centers = b_phi_centers,
      delta0_centers = cfg$delta0_centers,
      omega = omega,
      a_omega = a_omega,
      b_omega = b_omega,
      var_phi_centers = var_phi_centers,
      var_delta_centers = var_delta_centers,
      mstar_max = mstar_max,
      UpdateDitl = UpdateDitl,
      UpdateS = UpdateS,
      UpdateLambda = UpdateLambda,
      UpdateXi = UpdateXi,
      UpdateU = UpdateU,
      UpdateCenters = UpdateCenters,
      UpdateOmega = UpdateOmega,
      niter = niter,
      nburn = nburn,
      thin = thin
    )
    
    sink(log_file, split = FALSE)
    log_open = TRUE
    
    cat("Start dynamic centers MCMC\n")
    cat("tag:", tag, "\n")
    cat("gamma:", cfg$gamma, "\n")
    cat("delta0_centers:", cfg$delta0_centers, "\n")
    cat("niter:", niter, "nburn:", nburn, "thin:", thin, "\n")
    
    fit = GibbsSampler_DTM_centers(niter,nburn,thin,data,param_DTM_centers,init_DTM_centers)
    
    cat("Run completed\n")
    cat("Saved iterations:", length(fit$M), "\n")
    cat("Unique M values:", paste(unique(fit$M), collapse = ", "), "\n")
    cat("Unique Mstar values:", paste(unique(fit$Mstar), collapse = ", "), "\n")
    
    if(save_all_chain){
      saveRDS(fit, file = fit_file)
      saveRDS(list(
        tuning_options = tuning_options,
        param_DTM_centers = param_DTM_centers,
        init_DTM_centers = init_DTM_centers,
        static_init = static_init
      ), file = setup_file)
    } else {
      grDevices::pdf(pdf_file, width = 11, height = 8.5)
      pdf_open = TRUE
      par(mar = c(1,1,1,1))
      plot.new()
      text(0.02, 0.98,
           labels = paste("Run completed",
                          paste0("tag: ",tag),
                          paste0("gamma: ",cfg$gamma),
                          paste0("delta0_centers: ",cfg$delta0_centers),
                          paste0("saved iterations: ",length(fit$M)),
                          sep = "\n"),
           adj = c(0,1))
    }
    
    rm(fit, init_DTM_centers, param_DTM_centers)
    gc(verbose = FALSE)
    
    if(pdf_open){
      try(grDevices::dev.off(), silent = TRUE)
      pdf_open = FALSE
    }
    if(log_open){
      try(sink(), silent = TRUE)
      log_open = FALSE
    }
    
    list(
      config_id = cfg_id,
      tag = tag,
      gamma = cfg$gamma,
      delta0_centers = cfg$delta0_centers,
      fit_file = fit_file,
      setup_file = setup_file,
      pdf_file = pdf_file,
      log_file = log_file,
      status = "success",
      error_message = NA_character_
    )
  }, error = function(e){
    if(log_open){
      try(sink(), silent = TRUE)
      log_open = FALSE
    }
    
    writeLines(
      c(
        paste0("MCMC failed for ",tag),
        paste0("r = ",r),
        paste0("seed = ",chain_seed),
        paste0("gamma = ",cfg$gamma),
        paste0("delta0_centers = ",cfg$delta0_centers),
        conditionMessage(e)
      ),
      con = log_file
    )
    
    list(
      config_id = cfg_id,
      tag = tag,
      gamma = cfg$gamma,
      delta0_centers = cfg$delta0_centers,
      fit_file = fit_file,
      setup_file = setup_file,
      pdf_file = pdf_file,
      log_file = log_file,
      status = "failed",
      error_message = conditionMessage(e)
    )
  })
  
  result
}

cl = parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)

cat(sprintf("[%s] START whole script: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "SD_fit_center_parallel.R"))
flush.console()

parallel::clusterExport(
  cl,
  varlist = c(
    "data", "H", "Ttot", "r", "M0",
    "static_init", "zeta_floor", "zeta_strength",
    "sigma", "beta",
    "a_phi_centers", "b_phi_centers",
    "omega", "a_omega", "b_omega",
    "var_phi_centers", "var_delta_centers", "mstar_max",
    "UpdateDitl", "UpdateS", "UpdateLambda", "UpdateXi", "UpdateU",
    "UpdateCenters", "UpdateOmega",
    "print", "seed",
    "niter", "nburn", "thin",
    "output_dir", "log_dir", "save_dir", "save_all_chain",
    "static_nstart", "static_niter", "eps_init",
    "rfunctions_path", "rcpp_path",
    "params_grid",
    "format_tag_value", "make_config_tag",
    "build_center_init", "run_single_config"
  ),
  envir = environment()
)

for(worker_id in seq_along(cl)){
  parallel::clusterCall(cl[worker_id], function(rfunctions_path, rcpp_path){
    source(rfunctions_path)
    Rcpp::sourceCpp(rcpp_path)
    NULL
  }, rfunctions_path = rfunctions_path, rcpp_path = rcpp_path)
}

results = parallel::parLapplyLB(cl, seq_len(nrow(params_grid)), function(cfg_id){
  run_single_config(
    cfg = params_grid[cfg_id,,drop = FALSE],
    cfg_id = cfg_id,
    data = data,
    H = H,
    Ttot = Ttot,
    r = r,
    M0 = M0,
    static_init = static_init,
    zeta_floor = zeta_floor,
    zeta_strength = zeta_strength,
    sigma = sigma,
    beta = beta,
    a_phi_centers = a_phi_centers,
    b_phi_centers = b_phi_centers,
    omega = omega,
    a_omega = a_omega,
    b_omega = b_omega,
    var_phi_centers = var_phi_centers,
    var_delta_centers = var_delta_centers,
    mstar_max = mstar_max,
    UpdateDitl = UpdateDitl,
    UpdateS = UpdateS,
    UpdateLambda = UpdateLambda,
    UpdateXi = UpdateXi,
    UpdateU = UpdateU,
    UpdateCenters = UpdateCenters,
    UpdateOmega = UpdateOmega,
    print = print,
    seed = seed,
    niter = niter,
    nburn = nburn,
    thin = thin,
    output_dir = output_dir,
    log_dir = log_dir,
    save_dir = save_dir,
    save_all_chain = save_all_chain,
    static_nstart = static_nstart,
    static_niter = static_niter,
    eps_init = eps_init
  )
})

results_df = do.call(rbind, lapply(results, as.data.frame))

summary_file = file.path(output_dir, paste0("parallel_centers_summary_r",r,".csv"))

write.csv(
  results_df,
  file = summary_file,
  row.names = FALSE
)

fmt_path = function(x){
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

cat("\nSaved objects / output manifest\n")
cat("summary_csv: ", fmt_path(summary_file), "\n", sep = "")
cat("output_dir:  ", fmt_path(output_dir), "\n", sep = "")
cat("log_dir:     ", fmt_path(log_dir), "\n", sep = "")
cat("save_dir:    ", fmt_path(save_dir), "\n", sep = "")

for(i in seq_len(nrow(results_df))){
  cat("\nConfiguration ", i, "/", nrow(results_df), ": ",
      results_df$tag[i], "\n", sep = "")
  cat("status:    ", results_df$status[i], "\n", sep = "")
  cat("fit_rds:   ", fmt_path(results_df$fit_file[i]), "\n", sep = "")
  cat("setup_rds: ", fmt_path(results_df$setup_file[i]), "\n", sep = "")
  cat("log_file:  ", fmt_path(results_df$log_file[i]), "\n", sep = "")
  if(!save_all_chain || identical(results_df$status[i], "failed"))
    cat("pdf_file:  ", fmt_path(results_df$pdf_file[i]), "\n", sep = "")
}

cat(sprintf("[%s] END whole script: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "SD_fit_center_parallel.R"))
flush.console()

results_df
