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

delta_all     = c(1e-3,1e-2)  # Dirichlet parameter
mu_gamma_all  = c(0.1,10) 
var_gamma_all = c(1,10)
mu_beta_all   = c(0.1,10)
var_beta_all  = c(1,10)
sigma0_all    = c(0.25,0.5,0.75)

params_grid = expand.grid(
  delta = delta_all,
  mu_gamma = mu_gamma_all,
  var_gamma = var_gamma_all,
  mu_beta = mu_beta_all,
  var_beta = var_beta_all,
  sigma0 = sigma0_all,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
names(params_grid) = c("delta","mu_gamma","var_gamma","mu_beta","var_beta","sigma0")


params_grid = params_grid[c(10,13,19),]
# Parallel MCMC runner ----------------------------------------------------

library(parallel)

seed = 22123
H = 20 # number of atoms
niter = 200
nburn = 5
thin  = 2

# Fixed hyperparameters / MCMC settings
a_phi = 1; b_phi = 1
a_sigma = 1; b_sigma = 1
prop_var_phi = 0.01
prop_var_gamma = 0.01
prop_var_sigma = 0.01
prop_var_beta = 0.01

phi0 = 1; gamma0 = 1; beta0 = 1

UpdateDitl = TRUE; UpdateS = TRUE
UpdateLambda = TRUE; UpdateXi = TRUE
UpdateU = TRUE
UpdatePhi = FALSE; UpdateGamma = TRUE
UpdateSigma = FALSE; UpdateBeta = TRUE
print = TRUE; JointAdp = FALSE

output_dir = file.path(wd, "img", "parallel")
if(!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_dir = file.path(output_dir, "log")
if(!dir.exists(log_dir))
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

project_dir = normalizePath(file.path(wd, "..", ".."), winslash = "/", mustWork = TRUE)
rfunctions_path = file.path(project_dir, "R", "Rfunctions.R")
rcpp_path = file.path(project_dir, "src", "RcppFunctions.cpp")

format_tag_value = function(x) {
  out = format(signif(x, 6), scientific = TRUE, trim = TRUE)
  out = gsub("\\+", "", out)
  out = gsub("-", "m", out)
  out = gsub("\\.", "p", out)
  out
}

make_config_tag = function(cfg, cfg_id, r) {
  paste0(
    "cfg", sprintf("%03d", cfg_id),
    "_r", r,
    "_delta_", format_tag_value(cfg$delta),
    "_muGamma_", format_tag_value(cfg$mu_gamma),
    "_varGamma_", format_tag_value(cfg$var_gamma),
    "_muBeta_", format_tag_value(cfg$mu_beta),
    "_varBeta_", format_tag_value(cfg$var_beta),
    "_sigma0_", format_tag_value(cfg$sigma0)
  )
}

plot_traceplot_num_topics = function(fit, H, Ttot, main_prefix = "") {
  Kt_tr = sapply(fit$Xi, function(Xi_it) {
    apply(Xi_it, 2, function(x) length(which(x > 0)))
  })
  Kt_tr = t(Kt_tr)
  
  topics_per_page = 9
  topic_chunks = split(seq_len(Ttot), ceiling(seq_len(Ttot) / topics_per_page))
  
  for(chunk in topic_chunks) {
    par(mfrow = c(3,3), mar = c(3,3,2,1), mgp = c(2,0.5,0), bty = "l")
    for(t in chunk) {
      plot(Kt_tr[,t], xlab = "Iter.", ylab = paste0("K", t), type = "l",
           main = paste0(main_prefix, "Time ", t))
    }
    if(length(chunk) < topics_per_page) {
      for(dummy in seq_len(topics_per_page - length(chunk))) {
        plot.new()
      }
    }
  }
  
  K_tr = sapply(fit$Xi, function(Xi_it) {
    idx_list = lapply(seq_len(H), function(l) find_indices(Xi_it[l,]))
    sum(vapply(idx_list, function(x) length(x$idx_born), integer(1)))
  })
  
  par(mfrow = c(1,1), mar = c(3,3,2,1), mgp = c(2,0.5,0), bty = "l")
  plot(K_tr, xlab = "Iter.", ylab = "K", type = "l",
       main = paste0(main_prefix, "Total number of distinct topics"))
}

run_single_config = function(cfg, cfg_id, data, V, Ttot, H, r, output_dir, log_dir, seed,
                              niter, nburn, thin,
                              a_phi, b_phi, a_sigma, b_sigma,
                              prop_var_phi, prop_var_gamma, prop_var_sigma, prop_var_beta,
                              phi0, gamma0, beta0,
                              UpdateDitl, UpdateS, UpdateLambda, UpdateXi, UpdateU,
                              UpdatePhi, UpdateGamma, UpdateSigma, UpdateBeta,
                              print, JointAdp) {
  cfg = as.list(cfg)
  cfg$delta = as.numeric(cfg$delta)
  cfg$mu_gamma = as.numeric(cfg$mu_gamma)
  cfg$var_gamma = as.numeric(cfg$var_gamma)
  cfg$mu_beta = as.numeric(cfg$mu_beta)
  cfg$var_beta = as.numeric(cfg$var_beta)
  cfg$sigma0 = as.numeric(cfg$sigma0)
  
  tag = make_config_tag(cfg, cfg_id, r)
  pdf_file = file.path(output_dir, paste0(tag, ".pdf"))
  log_file = file.path(log_dir, paste0(tag, ".log"))
  chain_seed = seed + cfg_id
  log_open = FALSE
  pdf_open = FALSE
  on.exit({
    if(pdf_open) {
      try(grDevices::dev.off(), silent = TRUE)
    }
    if(log_open) {
      try(sink(), silent = TRUE)
    }
  }, add = TRUE)
  
  result = tryCatch({
    set.seed(chain_seed)
    
    delta = cfg$delta
    mu_gamma = cfg$mu_gamma
    var_gamma = cfg$var_gamma
    mu_beta = cfg$mu_beta
    var_beta = cfg$var_beta
    sigma0 = cfg$sigma0
    
    ab_gamma = set_par_gamma(mu_gamma, var_gamma)
    a_gamma = ab_gamma[1]; b_gamma = ab_gamma[2]
    ab_beta = set_par_gamma(mu_beta, var_beta)
    a_beta = ab_beta[1]; b_beta = ab_beta[2]
    
    Xi0 = matrix(sample(0:100, H * Ttot, TRUE), nrow = H, ncol = Ttot)
    S0 = matrix(rgamma(n = H * Ttot, 1, 1), nrow = H, ncol = Ttot)
    Lambda0 = vector("list", H)
    Lambda0 = lapply(Lambda0, function(x) {
      A = matrix(0, nrow = V, ncol = Ttot)
      apply(A, 2, function(y) {
        a = rgamma(n = V, shape = delta, rate = 1)
        a / sum(a)
      })
    })
    
    init_DTM = set_init_DTM(Xi0, Lambda0, S0, phi0, gamma0, sigma0, beta0)
    param_DTM = set_param_DTM(
      H, delta,
      a_phi, b_phi, a_gamma, b_gamma, a_sigma, b_sigma, a_beta, b_beta,
      prop_var_phi, prop_var_gamma, prop_var_sigma, prop_var_beta,
      UpdateDitl, UpdateS, UpdateLambda, UpdateXi, UpdateU,
      UpdatePhi, UpdateGamma, UpdateSigma, UpdateBeta,
      chain_seed, print, JointAdp
    )
    
    sink(log_file, split = FALSE)
    log_open = TRUE
    
    fit = GibbsSampler_DTM(niter, nburn, thin, data, param_DTM, init_DTM)
    
    grDevices::pdf(pdf_file, width = 12, height = 8)
    pdf_open = TRUE
    
    plot_traceplot_num_topics(
      fit = fit,
      H = H,
      Ttot = Ttot,
      main_prefix = paste0(tag, " | ")
    )
    
    rm(fit, Xi0, S0, Lambda0, init_DTM, param_DTM)
    gc(verbose = FALSE)
    try(grDevices::dev.off(), silent = TRUE)
    pdf_open = FALSE
    try(sink(), silent = TRUE)
    log_open = FALSE
    
    list(
      config_id = cfg_id,
      tag = tag,
      pdf_file = pdf_file,
      log_file = log_file,
      status = "success",
      error_message = NA_character_
    )
  }, error = function(e) {
    try({
      grDevices::pdf(pdf_file, width = 11, height = 8.5)
      par(mar = c(1,1,1,1))
      plot.new()
      text(
        0.02, 0.98,
        labels = paste(
          "MCMC failed",
          paste0("tag: ", tag),
          paste0("r: ", r),
          paste0("delta: ", cfg$delta),
          paste0("mu_gamma: ", cfg$mu_gamma),
          paste0("var_gamma: ", cfg$var_gamma),
          paste0("mu_beta: ", cfg$mu_beta),
          paste0("var_beta: ", cfg$var_beta),
          paste0("sigma0: ", cfg$sigma0),
          paste0("seed: ", chain_seed),
          "",
          paste(strwrap(conditionMessage(e), width = 100), collapse = "\n"),
          sep = "\n"
        ),
        adj = c(0,1)
      )
      dev.off()
    }, silent = TRUE)
    
    writeLines(
      c(
        paste0("MCMC failed for ", tag),
        paste0("r = ", r),
        paste0("seed = ", chain_seed),
        conditionMessage(e)
      ),
      con = log_file
    )
    
    list(
      config_id = cfg_id,
      tag = tag,
      pdf_file = pdf_file,
      log_file = log_file,
      status = "failed",
      error_message = conditionMessage(e)
    )
  })
  
  result
}

avail_cores = parallel::detectCores(logical = TRUE)
if(is.na(avail_cores))
  avail_cores = 1L
n_cores = 3 # <---

cl = parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)

parallel::clusterExport(
  cl,
  varlist = c(
    "data", "V", "Ttot", "H", "r", "output_dir", "log_dir", "seed",
    "niter", "nburn", "thin",
    "a_phi", "b_phi", "a_sigma", "b_sigma",
    "prop_var_phi", "prop_var_gamma", "prop_var_sigma", "prop_var_beta",
    "phi0", "gamma0", "beta0",
    "UpdateDitl", "UpdateS", "UpdateLambda", "UpdateXi", "UpdateU",
    "UpdatePhi", "UpdateGamma", "UpdateSigma", "UpdateBeta",
    "print", "JointAdp",
    "rfunctions_path", "rcpp_path",
    "params_grid",
    "format_tag_value", "make_config_tag",
    "plot_traceplot_num_topics", "run_single_config"
  ),
  envir = environment()
)

parallel::clusterEvalQ(cl, {
  source(rfunctions_path)
  Rcpp::sourceCpp(rcpp_path)
  NULL
})

results = parallel::parLapplyLB(cl, seq_len(nrow(params_grid)), function(cfg_id) {
  run_single_config(
    cfg = params_grid[cfg_id, , drop = FALSE],
    cfg_id = cfg_id,
    data = data,
    V = V,
    Ttot = Ttot,
    H = H,
    r = r,
    output_dir = output_dir,
    log_dir = log_dir,
    seed = seed,
    niter = niter,
    nburn = nburn,
    thin = thin,
    a_phi = a_phi,
    b_phi = b_phi,
    a_sigma = a_sigma,
    b_sigma = b_sigma,
    prop_var_phi = prop_var_phi,
    prop_var_gamma = prop_var_gamma,
    prop_var_sigma = prop_var_sigma,
    prop_var_beta = prop_var_beta,
    phi0 = phi0,
    gamma0 = gamma0,
    beta0 = beta0,
    UpdateDitl = UpdateDitl,
    UpdateS = UpdateS,
    UpdateLambda = UpdateLambda,
    UpdateXi = UpdateXi,
    UpdateU = UpdateU,
    UpdatePhi = UpdatePhi,
    UpdateGamma = UpdateGamma,
    UpdateSigma = UpdateSigma,
    UpdateBeta = UpdateBeta,
    print = print,
    JointAdp = JointAdp
  )
})

results_df = do.call(rbind, lapply(results, as.data.frame))
results_df = cbind(params_grid[results_df$config_id, , drop = FALSE], results_df)

write.csv(
  results_df,
  file = file.path(output_dir, paste0("parallel_summary_r", r, ".csv")),
  row.names = FALSE
)

results_df

