# Single-chain analysis for center-model runs -----------------------------

# Select one fitted chain and inspect centers, feature traces, and
# center-update diagnostics. This script is intentionally separate from the
# batch summary script.

# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here on the VM if needed
wd = paste0(choose_wd,"DynamicFeatureAllocation/Scripts/SpeechDataset")
setwd(wd)

# Data --------------------------------------------------------------------

seed_data = 132332
set.seed(seed_data)
r = 10

Presidents_all <- read.csv("../data/Presidents_all.csv")
data = read.table(paste0("../data/SpeechData_top",r,".txt"))
colnames(data) = Presidents_all[,2]
data = as.matrix(data)

V = nrow(data)
Ttot = ncol(data)
vocab = rownames(data)
if(is.null(vocab))
  vocab = as.character(seq_len(V))

# Select chain ------------------------------------------------------------

fit_basename = paste0("cfg003_center_r10_M3_H10_gamma_1em03_delta0_1em02_phi_1e00",".rds")
downloads_dir = "C:/Users/colom/Downloads"
local_save_dir = file.path(wd, "centers_save")

fit_file = file.path(downloads_dir, fit_basename)
if(!file.exists(fit_file))
  fit_file = file.path(local_save_dir, fit_basename)

setup_file = file.path(dirname(fit_file), sub("\\.rds$","_setup.rds", basename(fit_file)))
if(!file.exists(setup_file))
  setup_file = file.path(local_save_dir, sub("\\.rds$","_setup.rds", fit_basename))

if(!file.exists(fit_file))
  stop("fit_file does not exist: ", fit_file)

# Fallback initialization settings if setup is missing.
seed = 22123
M0 = 7
static_nstart = 20
static_niter = 500
eps_init = 1e-8
zeta_floor = 0.2
zeta_strength = 30

# Helpers ----------------------------------------------------------------

normalize_positive = function(x){
  x = as.numeric(x)
  s = sum(x)
  if(!is.finite(s) || s <= 0)
    stop("Cannot normalize a vector with non-positive or invalid sum")
  x/s
}

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

get_initial_zeta = function(setup_file){
  if(file.exists(setup_file)){
    setup = readRDS(setup_file)
    if(!is.null(setup$init_DTM_centers$Zeta0))
      return(setup$init_DTM_centers$Zeta0)
    if(!is.null(setup$static_init$Lambda) &&
       !is.null(setup$tuning_options$zeta_floor) &&
       !is.null(setup$tuning_options$zeta_strength)){
      return(lapply(seq_len(ncol(setup$static_init$Lambda)), function(m) {
        setup$tuning_options$zeta_floor +
          setup$tuning_options$zeta_strength*setup$static_init$Lambda[,m]
      }))
    }
    stop("setup file found, but it does not contain Zeta0 or static_init")
  }
  
  message("setup file not found. Reconstructing Zeta0 from the static Poisson-NMF initialization.")
  static_init = fit_static_poisson_nmf(data, M0,
                                       nstart = static_nstart,
                                       niter = static_niter,
                                       eps = eps_init,
                                       seed = seed)
  lapply(seq_len(M0), function(m) zeta_floor + zeta_strength*static_init$Lambda[,m])
}

compute_center_check = function(fit, Zeta0){
  if(is.null(fit$Zeta) || length(fit$Zeta) == 0)
    stop("fit$Zeta is missing or empty")
  
  n_saved = length(fit$Zeta)
  M_fit = length(fit$Zeta[[1]])
  if(length(Zeta0) != M_fit)
    stop("length(Zeta0) does not match the number of fitted centers")
  
  Zeta_mean = lapply(seq_len(M_fit), function(m){
    vals_m = lapply(seq_len(n_saved), function(it) as.numeric(fit$Zeta[[it]][[m]]))
    Reduce("+", vals_m)/n_saved
  })
  
  list(
    meta = list(
      fit_file = fit_file,
      setup_file = if(file.exists(setup_file)) setup_file else NA_character_,
      n_saved = n_saved,
      M_fit = M_fit,
      V = length(Zeta_mean[[1]])
    ),
    Zeta0 = Zeta0,
    Zeta_mean = Zeta_mean,
    Zeta0_norm = lapply(Zeta0, normalize_positive),
    Zeta_mean_norm = lapply(Zeta_mean, normalize_positive)
  )
}

compute_center_update_diagnostics = function(fit){
  if(is.null(fit$centers_aux) || length(fit$centers_aux) == 0)
    return(NULL)
  
  nonempty_aux = which(vapply(fit$centers_aux, function(x) length(x) > 0, logical(1)))
  if(length(nonempty_aux) == 0)
    return(list(update_recorded = FALSE))
  
  M_fit = length(fit$centers_aux[[nonempty_aux[1]]]$accept_zeta)
  
  accept_mat = t(vapply(nonempty_aux, function(it) {
    as.numeric(fit$centers_aux[[it]]$accept_zeta)
  }, numeric(M_fit)))
  
  cluster_size_mat = t(vapply(nonempty_aux, function(it) {
    as.numeric(fit$centers_aux[[it]]$cluster_size)
  }, numeric(M_fit)))
  
  log_acc_mat = t(vapply(nonempty_aux, function(it) {
    as.numeric(fit$centers_aux[[it]]$log_acc_zeta)
  }, numeric(M_fit)))
  
  list(
    update_recorded = TRUE,
    iterations_with_aux = nonempty_aux,
    accept_mat = accept_mat,
    cluster_size_mat = cluster_size_mat,
    log_acc_mat = log_acc_mat,
    accept_rate = colMeans(accept_mat),
    ever_accepted = colSums(accept_mat) > 0,
    mean_cluster_size = colMeans(cluster_size_mat),
    min_log_acc = apply(log_acc_mat, 2, min),
    median_log_acc = apply(log_acc_mat, 2, median),
    max_log_acc = apply(log_acc_mat, 2, max)
  )
}

compute_feature_trace = function(fit){
  if(is.null(fit$Lambda_star) || length(fit$Lambda_star) == 0)
    stop("fit$Lambda_star is missing or empty")
  if(is.null(fit$Lambda_star_by_center) || length(fit$Lambda_star_by_center) == 0)
    stop("fit$Lambda_star_by_center is missing or empty")
  
  n_saved = length(fit$Lambda_star)
  M_fit = length(fit$Lambda_star_by_center[[1]])
  
  K_total = vapply(seq_len(n_saved), function(it) {
    ncol(fit$Lambda_star[[it]])
  }, integer(1))
  
  K_by_center = t(vapply(seq_len(n_saved), function(it) {
    vapply(seq_len(M_fit), function(m) {
      ncol(fit$Lambda_star_by_center[[it]][[m]])
    }, integer(1))
  }, integer(M_fit)))
  
  colnames(K_by_center) = paste0("center_",seq_len(M_fit))
  
  list(
    K_total = K_total,
    K_by_center = K_by_center
  )
}

plot_feature_trace = function(feature_trace){
  iter = seq_along(feature_trace$K_total)
  M_fit = ncol(feature_trace$K_by_center)
  
  par(mfrow = c(1,1), mar = c(3,3,2,1), mgp=c(2,0.5,0), bty = "l")
  plot(iter, feature_trace$K_total,
       xlab = "Saved iteration",
       ylab = "K",
       main = "Total number of features",
       type = "l")
  
  par(mfrow = c(M_fit,1), mar = c(2.5,3,1.5,1), mgp=c(2,0.5,0), bty = "l")
  for(m in seq_len(M_fit)){
    plot(iter, feature_trace$K_by_center[,m],
         xlab = "Saved iteration",
         ylab = paste0("K_",m),
         main = paste0("Center ",m),
         type = "l")
  }
}

plot_center_check = function(center_check, vocab, top_n = V, order_vocab = FALSE){
  M_fit = center_check$meta$M_fit
  V = center_check$meta$V
  top_n = min(top_n, V)
  
  for(m in seq_len(M_fit)){
    init_m = center_check$Zeta0_norm[[m]]
    est_m = center_check$Zeta_mean_norm[[m]]
    
    if(order_vocab){
      idx = order(pmax(init_m, est_m), decreasing = TRUE)[seq_len(top_n)]
    } else {
      idx = seq_len(top_n)
    }
    ymax = max(init_m[idx], est_m[idx])*1.2
    
    par(mfrow = c(1,1), mar = c(5,3,3,1), mgp=c(2,0.5,0), bty = "l")
    barplot(
      height = rbind(initial = init_m[idx], estimated = est_m[idx]),
      beside = TRUE,
      names.arg = vocab[idx],
      las = 2,
      col = c("grey70","darkred"),
      border = NA,
      ylab = "Normalized center",
      main = paste0("Center ",m,": initial vs estimated"),
      ylim = c(0,ymax),
      cex.names = 0.8
    )
    legend("topright",
           legend = c("initial", "estimated"),
           fill = c("grey70","darkred"),
           border = NA,
           bty = "n")
  }
}

plot_wordcloud_panel = function(words, freq, col = "darkred",
                                scale = c(4,0.7), seed = 1){
  ord = order(freq, decreasing = TRUE)
  words = words[ord]
  freq = freq[ord]
  
  if(requireNamespace("wordcloud", quietly = TRUE)){
    set.seed(seed)
    wordcloud::wordcloud(
      words = words,
      freq = freq,
      random.order = FALSE,
      rot.per = 0,
      colors = col,
      scale = scale
    )
  } else {
    plot.new()
    set.seed(seed)
    x = runif(length(words), 0.08, 0.92)
    y = runif(length(words), 0.10, 0.90)
    cex = 0.65 + 2.35*sqrt(freq/max(freq))
    text(x, y, labels = words, cex = cex, col = col)
    mtext("install package 'wordcloud' for a true word cloud",
          side = 1, line = -1, cex = 0.7)
  }
}

plot_center_wordcloud_hist = function(center_check, vocab,
                                      centers = seq_len(center_check$meta$M_fit),
                                      top_n_cloud = min(80, length(vocab)),
                                      hist_top_n = length(vocab),
                                      col = "darkred",
                                      seed = 1){
  V = center_check$meta$V
  top_n_cloud = min(top_n_cloud, V)
  hist_top_n = min(hist_top_n, V)
  
  old_par = par(no.readonly = TRUE)
  on.exit({
    layout(1)
    par(old_par)
  })
  
  for(m in centers){
    est_m = center_check$Zeta_mean_norm[[m]]
    idx_cloud = order(est_m, decreasing = TRUE)[seq_len(top_n_cloud)]
    idx_hist = order(est_m, decreasing = TRUE)[seq_len(hist_top_n)]
    
    layout(matrix(c(1,2), nrow = 1), widths = c(1.05,1))
    
    par(mar = c(0.5,0.2,2,0.2))
    plot_wordcloud_panel(
      words = vocab[idx_cloud],
      freq = est_m[idx_cloud],
      col = col,
      seed = seed + m
    )
    title(main = paste0("Center ",m), line = 0.2)
    
    par(mar = c(3.2,2.8,2,0.6), mgp = c(1.8,0.45,0), bty = "l")
    barplot(
      height = est_m[idx_hist],
      names.arg = rep("", hist_top_n),
      col = col,
      border = NA,
      ylim = c(0, max(est_m[idx_hist])*1.08),
      ylab = "Mass",
      main = "Estimated center",
      axes = FALSE
    )
    axis(2, las = 1, cex.axis = 0.8)
  }
}

# Load selected chain -----------------------------------------------------

cat("Reading fit file:\n", normalizePath(fit_file, winslash = "/", mustWork = TRUE), "\n", sep = "")
fit = readRDS(fit_file)



View(fit)

m = 1
Kt_tr = sapply(fit$Xi, function(Xi_it_list) apply(Xi_it_list[[m]], 2, function(x) length(which(x > 0))) )
Kt_tr = t(Kt_tr)

it = 5000
Kt_tr[it,]
fit$Xi[[it]][[m]]
fit$S[[it]][[m]]
fit$U[[it]][[m]]


par(mfrow = c(1,2), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
plot( fit$t_sigma_gamma, xlab = "Iter.", ylab = "t", type = "l" )
plot( log(fit$t_sigma_gamma), xlab = "Iter.", ylab = "log(t)", type = "l" )


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


 