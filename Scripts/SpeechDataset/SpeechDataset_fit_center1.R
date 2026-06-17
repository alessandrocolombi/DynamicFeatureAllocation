wd = "C:/Users/colom/DynamicFeatureAllocation/Scripts/SpeechDataset"
# wd = "/home/lucia.paci/Lucia/Ale/DynamicFeatureAllocation/Scripts/SpeechDataset"
setwd(wd)
source("./../../R/Rfunctions_centers.R")
Rcpp::sourceCpp("./../../src/RcppFunctions_centers.cpp")
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


# Model and tuning options ------------------------------------------------

seed = 22123
set.seed(seed)

# Fixed dimensions of the dynamic model with centers.
H = 10
M0 = 3

# Static Poisson-NMF initialization.
static_nstart = 20
static_niter = 500
eps_init = 1e-8
zeta_floor = 0.2
zeta_strength = 30

# Process hyperparameters.
gamma = 0.1
sigma = 0.1
beta  = 0.1

# Center-clustering hyperparameters.
a_phi_centers = 2
b_phi_centers = 1
delta0_centers = 1
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
niter = 2000
nburn = 1000
thin  = 5

# Output control.
run_mcmc = TRUE
save_fit = TRUE
save_dir = "save/centers"


# Static initialization for centers ---------------------------------------

# This block ignores the dynamic structure and fits a static Poisson NMF,
# data ~= Lambda_static %*% Xi_static. The columns of Lambda_static are used
# as initial center profiles, while Xi_static initializes the first atom in
# each center-specific latent process.
set.seed(seed)

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

static_init = fit_static_poisson_nmf(data, M0,
                                     nstart = static_nstart,
                                     niter = static_niter,
                                     eps = eps_init, seed = seed)

center_profiles = static_init$Lambda
Xi_static = static_init$Xi
center_time = static_init$center_time

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
  
  # Give the first atom of each center the static lexical profile.
  for(t in seq_len(Ttot)){
    Lambda0[[m]][[1]][,t] = center_profiles[,m]
  }
}

init_DTM_centers = set_init_DTM_centers(Xi0,Lambda0,S0,Zeta0)

cat("Static Poisson-NMF center initialization completed\n")
cat("KL objective:", signif(static_init$loss, 5), "\n")
cat("Center times:", paste(round(center_time, 2), collapse = ", "), "\n")
cat("Center masses:", paste(round(rowSums(Xi_static), 2), collapse = ", "), "\n")

par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
matplot(seq_len(Ttot), t(Xi_static),
        xlab = "Time", ylab = "Static intensity",
        type = "l", lty = 1, lwd = 2)
legend("topright",
       legend = paste0("Center ",seq_len(M0)),
       col = seq_len(M0), lty = 1, lwd = 2,
       bty = "n", cex = 0.8)

ymax_center_profile = max(center_profiles)*1.2
for(m in seq_len(M0)){
  par(mfrow = c(1,1), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
  barplot( height = unname(center_profiles[,m]),
           names.arg = rownames(data),
           las = 2, col = "darkred", border = NA,
           xlab = "Word", ylab = "Prob.",
           main = paste0("Initial center profile ",m),
           ylim = c(0,ymax_center_profile),
           cex.names = 0.65 )
}


# Dynamic model run -------------------------------------------------------

param_DTM_centers = set_param_DTM_centers(H,gamma,sigma,beta,
                                          a_phi_centers,b_phi_centers,delta0_centers,
                                          omega,a_omega,b_omega,
                                          var_phi_centers,var_delta_centers,mstar_max,
                                          UpdateDitl,UpdateS,UpdateLambda,UpdateXi,UpdateU,
                                          UpdateCenters,UpdateOmega,
                                          seed,print)

tuning_options = list(
  r = r,
  seed = seed,
  H = H,
  M0 = M0,
  static_nstart = static_nstart,
  static_niter = static_niter,
  eps_init = eps_init,
  zeta_floor = zeta_floor,
  zeta_strength = zeta_strength,
  gamma = gamma,
  sigma = sigma,
  beta = beta,
  a_phi_centers = a_phi_centers,
  b_phi_centers = b_phi_centers,
  delta0_centers = delta0_centers,
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

format_tag_value = function(x){
  out = format(x, scientific = TRUE, trim = TRUE)
  out = gsub("\\+","",out)
  out = gsub("-","m",out)
  out = gsub("\\.","p",out)
  out
}

fit_tag = paste0("center1_r",r,
                 "_M",M0,
                 "_H",H,
                 "_gamma",format_tag_value(gamma),
                 "_sigma",format_tag_value(sigma),
                 "_beta",format_tag_value(beta))

if(run_mcmc){
  cat("Start dynamic centers MCMC\n")
  cat("niter:", niter, "nburn:", nburn, "thin:", thin, "\n")
  
  fit = GibbsSampler_DTM_centers(niter,nburn,thin,data,param_DTM_centers,init_DTM_centers)
  
  cat("Run completed\n")
  cat("Saved iterations:", length(fit$M), "\n")
  cat("Unique M values:", paste(unique(fit$M), collapse = ", "), "\n")
  cat("Unique Mstar values:", paste(unique(fit$Mstar), collapse = ", "), "\n")
  
  if(save_fit){
    if(!dir.exists(save_dir))
      dir.create(save_dir, recursive = TRUE)
    
    fit_file = file.path(save_dir, paste0(fit_tag,".rds"))
    setup_file = file.path(save_dir, paste0(fit_tag,"_setup.rds"))
    
    saveRDS(fit, file = fit_file)
    saveRDS(list(
      tuning_options = tuning_options,
      param_DTM_centers = param_DTM_centers,
      init_DTM_centers = init_DTM_centers,
      static_init = static_init,
      center_profiles = center_profiles,
      Xi_static = Xi_static,
      Presidents_all = Presidents_all
    ), file = setup_file)
    
    cat("Saved fit:", fit_file, "\n")
    cat("Saved setup:", setup_file, "\n")
  }
} else {
  cat("run_mcmc is FALSE: dynamic centers MCMC was not launched.\n")
}

