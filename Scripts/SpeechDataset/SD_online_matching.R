# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc, wd_unicatt, wd_g100, wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here if needed
wd = paste0(choose_wd, "DynamicFeatureAllocation/Scripts/SpeechDataset")
setwd(wd)

# Libraries ---------------------------------------------------------------
#
# The online matching itself does not require any special package.
# Everything below is written with base R.

# Plot options -------------------------------------------------------------
eps_break = .Machine$double.eps
mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)
save_img = FALSE
width = 12; height = 6
# Input / output ----------------------------------------------------------
#
# This script reads one compact summary object from `save_summary/`,
# builds an online consensus template across iterations, and saves the
# result back into `save_summary/`.
save_summary_dir = file.path(wd, "save_summary")
if(!dir.exists(save_summary_dir))
  stop("save_summary folder does not exist: ", save_summary_dir)

summary_file_name = NULL

summary_files = list.files(
  save_summary_dir,
  pattern = "_summary\\.rds$",
  full.names = TRUE
)

if(length(summary_files) == 0)
  stop("No *_summary.rds file found in: ", save_summary_dir)

if(is.null(summary_file_name)) {
  summary_file = summary_files[1]
} else {
  summary_file = file.path(save_summary_dir, summary_file_name)
  if(!file.exists(summary_file))
    stop("Summary file not found: ", summary_file)
}

cat("\nLoading summary object:\n", basename(summary_file), "\n", sep = "")
fit_summary = readRDS(summary_file)

# Load the full-chain fit as well, because Lambda_star is only available
# in the raw MCMC output.
save_dir = file.path(wd, "save")
if(!dir.exists(save_dir))
  stop("save folder does not exist: ", save_dir)

fit_file = file.path(
  save_dir,
  sub("_summary\\.rds$", ".rds", basename(summary_file))
)
if(!file.exists(fit_file))
  stop("Matching raw fit file not found: ", fit_file)

cat("Loading full-chain fit object:\n", basename(fit_file), "\n", sep = "")
fit_all = readRDS(fit_file)
Lambda_fit = lapply(fit_all$Lambda_star_mcmc, function(Lam_it) t(do.call(rbind, Lam_it)))

# Basic checks ------------------------------------------------------------
if(is.null(fit_summary$topic_objs) || length(fit_summary$topic_objs) == 0)
  stop("fit_summary$topic_objs is missing or empty.")
if(length(Lambda_fit) != length(fit_summary$topic_objs))
  stop("Lambda_fit and topic_objs have different lengths.")

Ttot = nrow(fit_summary$topic_objs[[1]]$Activity)
Niter = length(fit_summary$topic_objs)
V = nrow(Lambda_fit[[1]])

# Matching weights / threshold -------------------------------------------
#
# The distance combines:
#   - birth/death discrepancy
#   - Euclidean distance on normalized Xi
#
# A new column is created if the best distance is above `match_threshold`.
# Since support is fully determined by birth/death under the model
# assumptions, we do not include a separate Hamming term on Z.
w_birth = 1
w_death = 1
w_xi = 1
match_threshold = 1000

# Helpers -----------------------------------------------------------------

# Extract the paired objects for one iteration:
#   - Zstar      : Ttot x K
#   - Xistar     : Ttot x K
#   - Lambda_it  : V x K
#
# Column k of Lambda_it is paired with column k of Xistar and Zstar.
get_topic_triplet = function(topic_obj, Lambda_it) {
  Zstar = topic_obj$Activity
  Xistar = t(topic_obj$Xi_star)
  
  if(!all(dim(Zstar) == dim(Xistar)))
    stop("Zstar and Xistar do not have matching dimensions.")
  if(ncol(Lambda_it) != ncol(Zstar))
    stop("Lambda_it does not have the same number of topic columns as Zstar/Xistar.")
  
  list(Zstar = Zstar, Xistar = Xistar, Lambda = Lambda_it)
}

# Normalize Xi so the Xi distance reflects SHAPE more than scale.
normalize_xi = function(x) {
  sx = sum(x)
  if(sx <= 0) {
    rep(0, length(x))
  } else {
    x / sx
  }
}

# Extract birth/death times from one binary topic support.
extract_birth_death = function(z_col) {
  idx = which(z_col > 0)
  
  if(length(idx) == 0) {
    list(birth = Inf, death = -Inf)
  } else {
    list(birth = idx[1], death = idx[length(idx)])
  }
}

# Distance between one incoming topic and one template topic.
#
# Template Z columns are probabilistic (running means), so we threshold
# them at 0.5 to recover a binary activity pattern for the Hamming part.
topic_distance = function(template_z_prob, template_xi_mean, z_new, xi_new,
                          w_birth, w_death, w_xi) {
  template_z_bin = as.integer(template_z_prob >= 0.5)
  
  bd_template = extract_birth_death(template_z_bin)
  bd_new = extract_birth_death(z_new)
  
  db = abs(bd_template$birth - bd_new$birth)
  dd = abs(bd_template$death - bd_new$death)
  
  x_template_norm = normalize_xi(template_xi_mean)
  x_new_norm = normalize_xi(xi_new)
  dxi = sqrt(sum((x_template_norm - x_new_norm)^2))
  
  w_birth * db +
    w_death * dd +
    w_xi * dxi
}

# Add a new template column. We keep:
#   - running mean of Z
#   - running mean of Xi
#   - running mean of tail-mass Xi
#   - running mean of Lambda (conditional on being matched/present)
#   - match count
#   - first / last iteration where the column was seen
append_template_column = function(template, z_new, xi_new, lambda_new, it) {
  cum_xi_new = rev(cumsum(rev(xi_new)))
  
  template$Z_mean = cbind(template$Z_mean, z_new)
  template$Xi_mean = cbind(template$Xi_mean, xi_new)
  template$CumXi_mean = cbind(template$CumXi_mean, cum_xi_new)
  template$Lambda_mean = cbind(template$Lambda_mean, lambda_new)
  template$count = c(template$count, 1L)
  template$first_seen = c(template$first_seen, it)
  template$last_seen = c(template$last_seen, it)
  
  template
}

# Update an existing template column via running averages.
update_template_column = function(template, k, z_new, xi_new, lambda_new, it) {
  m = template$count[k]
  cum_xi_new = rev(cumsum(rev(xi_new)))
  
  template$Z_mean[, k] =
    (m * template$Z_mean[, k] + z_new) / (m + 1)
  template$Xi_mean[, k] =
    (m * template$Xi_mean[, k] + xi_new) / (m + 1)
  template$CumXi_mean[, k] =
    (m * template$CumXi_mean[, k] + cum_xi_new) / (m + 1)
  template$Lambda_mean[, k] =
    (m * template$Lambda_mean[, k] + lambda_new) / (m + 1)
  
  template$count[k] = m + 1L
  template$last_seen[k] = it
  
  template
}

# Initialize the online template from iteration 1.
initialize_template = function(topic_triplet) {
  Z0 = topic_triplet$Zstar
  Xi0 = topic_triplet$Xistar
  Lambda0 = topic_triplet$Lambda
  K0 = ncol(Z0)
  
  CumXi0 = apply(Xi0, 2, function(col) rev(cumsum(rev(col))))
  if(!is.matrix(CumXi0))
    CumXi0 = matrix(CumXi0, ncol = K0)
  
  list(
    Z_mean = Z0,
    Xi_mean = Xi0,
    CumXi_mean = CumXi0,
    Lambda_mean = Lambda0,
    count = rep(1L, K0),
    first_seen = rep(1L, K0),
    last_seen = rep(1L, K0)
  )
}

# One online update step:
# process one iteration column-by-column, matching each new column to the
# current template if the best distance is small enough; otherwise append.
online_update = function(template, topic_triplet, it,
                         w_birth, w_death, w_xi,
                         match_threshold) {
  Zstar = topic_triplet$Zstar
  Xistar = topic_triplet$Xistar
  Lambda_it = topic_triplet$Lambda
  K_new = ncol(Zstar)
  
  matched_template = rep(FALSE, length(template$count))
  assignment = rep(NA_integer_, K_new)
  assignment_dist = rep(NA_real_, K_new)
  
  for(j in seq_len(K_new)) {
    K_template = length(template$count)
    dist_vec = rep(Inf, K_template)
    
    # We do not allow two new columns from the same iteration to update
    # the same template slot in the same pass.
    for(k in seq_len(K_template)) {
      if(!matched_template[k]) {
        dist_vec[k] = topic_distance(
          template_z_prob = template$Z_mean[, k],
          template_xi_mean = template$Xi_mean[, k],
          z_new = Zstar[, j],
          xi_new = Xistar[, j],
          w_birth = w_birth,
          w_death = w_death,
          w_xi = w_xi
        )
      }
    }
    
    best_k = which.min(dist_vec)
    best_d = dist_vec[best_k]
    
    if(length(best_k) == 0 || !is.finite(best_d) || best_d > match_threshold) {
      template = append_template_column(
        template = template,
        z_new = Zstar[, j],
        xi_new = Xistar[, j],
        lambda_new = Lambda_it[, j],
        it = it
      )
      assignment[j] = length(template$count)
      assignment_dist[j] = NA_real_
      matched_template = c(matched_template, TRUE)
    } else {
      template = update_template_column(
        template = template,
        k = best_k,
        z_new = Zstar[, j],
        xi_new = Xistar[, j],
        lambda_new = Lambda_it[, j],
        it = it
      )
      assignment[j] = best_k
      assignment_dist[j] = best_d
      matched_template[best_k] = TRUE
    }
  }
  
  list(
    template = template,
    assignment = assignment,
    assignment_dist = assignment_dist
  )
}

# Online matching ---------------------------------------------------------
 
it_start = 1 #Niter/2
Lsaved_iter = Niter - it_start

# We start from iteration 1, use it as the initial template, then process
# all remaining iterations one by one.
topic_triplet_1 = get_topic_triplet(fit_summary$topic_objs[[it_start]], Lambda_fit[[it_start]])
template = initialize_template(topic_triplet_1)

n_saved = 1
history = vector("list", Lsaved_iter)
history[[n_saved]] = list(
  template_size = length(template$count),
  assignment = seq_len(length(template$count)),
  assignment_dist = rep(0, length(template$count))
)
n_saved = n_saved + 1

Ncols = length(template$count)
cat("\n","Ncols = ",Ncols,"\n")

pb = txtProgressBar(min = it_start + 1, max = Niter, style = 3)
for(it in (it_start+1):Niter) {
  if(Ncols != length(template$count)){
    Ncols = length(template$count)
    cat("\n","Ncols = ",Ncols,"\n")
  }
  
  topic_triplet_it = get_topic_triplet(fit_summary$topic_objs[[it]], Lambda_fit[[it]])
  
  step_out = online_update(
    template = template,
    topic_triplet = topic_triplet_it,
    it = it,
    w_birth = w_birth,
    w_death = w_death,
    w_xi = w_xi,
    match_threshold = match_threshold
  )
  
  template = step_out$template
  history[[n_saved]] = list(
    template_size = length(template$count),
    assignment = step_out$assignment,
    assignment_dist = step_out$assignment_dist
  )
  n_saved = n_saved + 1
  setTxtProgressBar(pb, it)
}
close(pb)

# Reordering final consensus columns -------------------------------------
#
# Once the online pass is finished, we order the final template by:
#   1. decreasing match count
#   2. increasing mean birth time
#   3. decreasing total Xi mass
#
# This gives a readable final consensus ordering.
template_birth = apply(template$Z_mean, 2, function(col) {
  idx = which(col >= 0.5)
  if(length(idx) == 0) Inf else idx[1]
})
template_mass = colSums(template$Xi_mean)

ord_final = order(-template$count, template_birth, -template_mass)

template$Z_mean = template$Z_mean[, ord_final, drop = FALSE]
template$Xi_mean = template$Xi_mean[, ord_final, drop = FALSE]
template$CumXi_mean = template$CumXi_mean[, ord_final, drop = FALSE]
template$Lambda_mean = template$Lambda_mean[, ord_final, drop = FALSE]
template$count = template$count[ord_final]
template$first_seen = template$first_seen[ord_final]
template$last_seen = template$last_seen[ord_final]


# Final rescaling over ALL iterations ------------------------------------
#
# Up to this point, the running means are conditional on the topic being
# matched. For example, if count[l] = 1, then Xi_mean[, l] is just the one
# observed Xi vector for that topic.
#
# For posterior summaries, we also want the unconditional mean over all
# Niter iterations, where an unmatched topic contributes a zero column.
# If a column is seen count[l] times, then:
#   unconditional mean = conditional mean * count[l] / Niter
#
# We keep both versions:
#   - *_conditional : average only over matched occurrences
#   - the rescaled template fields below: average over ALL iterations
template$Z_mean_conditional = template$Z_mean
template$Xi_mean_conditional = template$Xi_mean
template$CumXi_mean_conditional = template$CumXi_mean
template$Lambda_mean_conditional = template$Lambda_mean

count_scale = template$count / Lsaved_iter
template$Z_mean = sweep(template$Z_mean, 2, count_scale, `*`)
template$Xi_mean = sweep(template$Xi_mean, 2, count_scale, `*`)
template$CumXi_mean = sweep(template$CumXi_mean, 2, count_scale, `*`)

# Save result -------------------------------------------------------------
out = list(
  meta = list(
    summary_file = basename(summary_file),
    Ttot = Ttot,
    Niter = Niter,
    weights = c(
      w_birth = w_birth,
      w_death = w_death,
      w_xi = w_xi
    ),
    match_threshold = match_threshold
  ),
  template = template,
  history = history,
  Z_online_mean = template$Z_mean,
  Xi_online_mean = template$Xi_mean,
  CumXi_online_mean = template$CumXi_mean,
  Lambda_online_mean = template$Lambda_mean_conditional,
  Z_online_mean_conditional = template$Z_mean_conditional,
  Xi_online_mean_conditional = template$Xi_mean_conditional,
  CumXi_online_mean_conditional = template$CumXi_mean_conditional,
  Lambda_online_mean_conditional = template$Lambda_mean_conditional
)

out_file = file.path(
  save_summary_dir,
  sub("_summary\\.rds$", "_online_matching.rds", basename(summary_file))
)

saveRDS(out, file = out_file)

cat("\nSaved online matching object:\n", basename(out_file), "\n", sep = "")
cat("Final number of consensus columns: ", length(template$count), "\n", sep = "")
cat("Match counts per column:\n")
print(template$count)

# Plots and analysis -------------------------------------------------------------
if(FALSE){
  res = readRDS("save_summary/cfg001_r10_delta_1e00_beta0_1em01_gamma0_1em01_sigma0_1em01_online_matching.rds")
  
  ## Xi mean plot -------------------------------------------------------------
  soglia = 25
  plot_mat = res$Xi_online_mean
  sel_colums = which(colSums(plot_mat) > soglia )
  plot_mat = plot_mat[,sel_colums]
  
  max_plot_mat = max(plot_mat, na.rm = TRUE)
  K_plot = ncol(plot_mat)
  
  par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
  green_breaks = seq(soglia, max_plot_mat + eps_break, length.out = length(mycol) + 1)
  breaks_xi = c(0, soglia, green_breaks[-1])
  image( 1:K_plot, 1:Ttot,
         t(plot_mat),
         col = c("blue", mycol),
         breaks = breaks_xi,
         xlab = "Topics",
         ylab = "Time",
         main = "< Xi >",
         axes = FALSE )
  axis(2, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(1, at = seq(1, K_plot, length.out = min(K_plot, 10)), 
       labels = round(seq(1, K_plot, length.out = min(K_plot, 10))),
       cex.axis = 0.7)
  box()
  fields::image.plot(
    1:K_plot, 1:Ttot,
    t(plot_mat),
    col = c("blue", mycol),
    breaks = breaks_xi,
    legend.only = TRUE,
    horizontal = FALSE,
    legend.width = 1.2,            # controls legend thickness
    legend.shrink = 0.8,           # smaller legend
    legend.mar = 8.5,                # margin from image
    legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
  )
  ## Xi and Lambda paired plot (all) -------------------------------------------------------------
  ymax_topic = max( res$Lambda_online_mean )
  ymax_xi = max( res$Xi_online_mean )
  Kest = ncol(res$Xi_online_mean)
  
  if(save_img)
    pdf("img/Lambda_Xi_paired_all.pdf",width = width, height = height)
  for(l in 1:Kest){
    par(mfrow = c(1,2), mar = c(3,3,1,1), mgp=c(2,0.5,0), bty = "l")
    barplot( height = res$Lambda_online_mean[,l], 
             names.arg = as.character(1:V),
             las = 1, col = "darkred", border = NA,
             xlab = "Word",
             main = paste0("Topic: ",l), ylab = "Prob.", ylim = c(0,ymax_topic),
             cex.names = 0.5 )
    barplot( height = res$Xi_online_mean[,l], 
             names.arg = as.character(1:Ttot),
             las = 1, col = "darkgreen", border = NA,
             xlab = "Time",
             main = paste0("Xi: ",l), ylab = "Intensity", ylim = c(0,ymax_xi),
             cex.names = 0.5 )
  }
  if(save_img)
    dev.off()
  ## Xi and Lambda paired plot (selected) -------------------------------------------------------------
  sel_cols = which(apply(res$Xi_online_mean, 2, max) > 25)
  sel_Lambda = res$Lambda_online_mean[,sel_cols]
  sel_Xi = res$Xi_online_mean[,sel_cols]
  
  ymax_topic = max( sel_Lambda )
  ymax_xi = max( sel_Xi )
  K_sel = ncol( sel_Xi )
  
  par(mfrow = c(2,5), mar = c(1,4,1,1), mgp=c(2,0.5,0), bty = "l")
  # L'ho salvato a mano come SD_TopicXi_all_1 e SD_TopicXi_all_2
  # Dimensioni 20x8
  for(l in 1:K_sel){
    barplot( height = sel_Lambda[,l], 
             names.arg = "",#as.character(1:V),
             las = 1, col = "darkred", border = NA,
             # xlab = "Word",ylab = "Prob.",
             ylim = c(0,ymax_topic),
             cex.axis = 2,
             cex.names = 0.5 )
    title(
      main = paste0("k = ", l),
      cex.main = 2.5,
      line = -2
    )
    barplot( height = sel_Xi[,l], 
             names.arg = "",  #as.character(1:Ttot),
             las = 1, col = "darkgreen", border = NA,
             # xlab = "Time", ylab = "Intensity",
             ylim = c(0,ymax_xi),
             cex.axis = 2,
             cex.names = 1.5 )
    title(
      main = paste0("k = ", l),
      cex.main = 2.5,
      line = -2
    )
  }
  
  ## Cloud + Lambda + Xi (selected) ----------------------------------------------------------------
  r = 10
  Presidents_all <- read.csv("../data/Presidents_all.csv")
  data = read.table(paste0("../data/SpeechData_top",r,".txt"))
  colnames(data) = Presidents_all[,2]
  data = as.matrix(data)
  vocab = rownames(data)
  if(is.null(vocab))
    vocab = as.character(seq_len(nrow(data)))
  
  if(!requireNamespace("wordcloud", quietly = TRUE))
    stop("Package 'wordcloud' is required for the analysis plots.")
  
  if(save_img)
    pdf("img/Lambda_Xi_wordcloud_sel.pdf", width = width, height = height)
  
  for(l in seq_len(K_sel)) {
    
    top_idx = which( sel_Lambda[, l] > quantile(sel_Lambda[, l],0.9) )
    par(mfrow = c(1,3), mar = c(3,3,2,1), mgp = c(2,0.5,0), bty = "l")
    # Word cloud: word size proportional to Lambda weight.
    wordcloud::wordcloud(
      words = vocab[top_idx],
      freq = sel_Lambda[top_idx, l],
      scale = c(4, 0.8),
      min.freq = min(sel_Lambda[top_idx, l]),
      max.words = 50,
      random.order = FALSE,
      ordered.colors = TRUE,
      rot.per = 0,
      colors = rep("darkred", length(top_idx))
    )
    title(main = paste0("Word cloud: ", l), line = -1)
    
    barplot(
      height = sel_Lambda[, l],
      names.arg = vocab,
      las = 1, col = "darkred", border = NA,
      xlab = "Word",
      main = paste0("Topic: ", l), ylab = "Prob.",
      ylim = c(0, ymax_topic),
      cex.names = 0.5
    )
    barplot(
      height = sel_Xi[, l],
      names.arg = as.character(1:Ttot),
      las = 1, col = "darkgreen", border = NA,
      xlab = "Time",
      main = paste0("Xi: ", l), ylab = "Intensity",
      ylim = c(0, ymax_xi),
      cex.names = 0.5
    )
  }
  
  if(save_img)
    dev.off()
  
  ## Cloud only (all words) ----------------------------------------------------------------
  top_idx = 1:V
  
  if(save_img)
    pdf("img/_wordcloud_full_sel.pdf", width = width, height = height)
  
  for(l in seq_len(K_sel)) {
    
    par(mfrow = c(1,1), mar = c(3,3,2,1), mgp = c(2,0.5,0), bty = "l")
    wordcloud::wordcloud(
      words = vocab[top_idx],
      freq = sel_Lambda[top_idx, l],
      scale = c(4, 0.8),
      min.freq = min(sel_Lambda[top_idx, l]),
      max.words = 50,
      random.order = FALSE,
      ordered.colors = TRUE,
      rot.per = 0,
      colors = rep("darkred", length(top_idx))
    )
  }
  
  if(save_img)
    dev.off()

  ## Normalized xi_t ---------------------------------------------------------
  sel_cols = which(apply(res$Xi_online_mean, 2, max) > soglia)
  par(mfrow = c(5,6), mar = c(1,3,2,1), mgp=c(2,0.5,0), bty = "l")
  for(t in 1:Ttot){
    xi_t = res$Xi_online_mean[t,sel_cols]
    xi_t_norm = unname(xi_t/sum(xi_t))
    barplot( height = xi_t_norm, 
             names.arg = as.character(1:length(sel_cols)),
             las = 1, col = "darkgreen", border = NA,
             xlab = "", ylab = "Prop.", 
             main = Presidents_all[t,2], #paste0("t ",t), 
             ylim = c(0,0.6),
             cex.names = 0.5 )
    
  }
  
  
}

