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

# Basic checks ------------------------------------------------------------
if(is.null(fit_summary$topic_objs) || length(fit_summary$topic_objs) == 0)
  stop("fit_summary$topic_objs is missing or empty.")

Ttot = nrow(fit_summary$topic_objs[[1]]$Activity)
Niter = length(fit_summary$topic_objs)

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

# Extract Ttot x K matrices from one saved topic object.
get_topic_pair = function(topic_obj) {
  Zstar = topic_obj$Activity
  Xistar = t(topic_obj$Xi_star)
  
  if(!all(dim(Zstar) == dim(Xistar)))
    stop("Zstar and Xistar do not have matching dimensions.")
  
  list(Zstar = Zstar, Xistar = Xistar)
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
#   - match count
#   - first / last iteration where the column was seen
append_template_column = function(template, z_new, xi_new, it) {
  cum_xi_new = rev(cumsum(rev(xi_new)))
  
  template$Z_mean = cbind(template$Z_mean, z_new)
  template$Xi_mean = cbind(template$Xi_mean, xi_new)
  template$CumXi_mean = cbind(template$CumXi_mean, cum_xi_new)
  template$count = c(template$count, 1L)
  template$first_seen = c(template$first_seen, it)
  template$last_seen = c(template$last_seen, it)
  
  template
}

# Update an existing template column via running averages.
update_template_column = function(template, k, z_new, xi_new, it) {
  m = template$count[k]
  cum_xi_new = rev(cumsum(rev(xi_new)))
  
  template$Z_mean[, k] =
    (m * template$Z_mean[, k] + z_new) / (m + 1)
  template$Xi_mean[, k] =
    (m * template$Xi_mean[, k] + xi_new) / (m + 1)
  template$CumXi_mean[, k] =
    (m * template$CumXi_mean[, k] + cum_xi_new) / (m + 1)
  
  template$count[k] = m + 1L
  template$last_seen[k] = it
  
  template
}

# Initialize the online template from iteration 1.
initialize_template = function(topic_pair) {
  Z0 = topic_pair$Zstar
  Xi0 = topic_pair$Xistar
  K0 = ncol(Z0)
  
  CumXi0 = apply(Xi0, 2, function(col) rev(cumsum(rev(col))))
  if(!is.matrix(CumXi0))
    CumXi0 = matrix(CumXi0, ncol = K0)
  
  list(
    Z_mean = Z0,
    Xi_mean = Xi0,
    CumXi_mean = CumXi0,
    count = rep(1L, K0),
    first_seen = rep(1L, K0),
    last_seen = rep(1L, K0)
  )
}

# One online update step:
# process one iteration column-by-column, matching each new column to the
# current template if the best distance is small enough; otherwise append.
online_update = function(template, topic_pair, it,
                         w_birth, w_death, w_xi,
                         match_threshold) {
  Zstar = topic_pair$Zstar
  Xistar = topic_pair$Xistar
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
 
it_start = Niter/2
Lsaved_iter = Niter - it_start

# We start from iteration 1, use it as the initial template, then process
# all remaining iterations one by one.
topic_pair_1 = get_topic_pair(fit_summary$topic_objs[[it_start]])
template = initialize_template(topic_pair_1)

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
  
  topic_pair_it = get_topic_pair(fit_summary$topic_objs[[it]])
  
  step_out = online_update(
    template = template,
    topic_pair = topic_pair_it,
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

count_scale = template$count / Lsaved_iter
template$Z_mean = sweep(template$Z_mean, 2, count_scale, `*`)
template$Xi_mean = sweep(template$Xi_mean, 2, count_scale, `*`)
template$CumXi_mean = sweep(template$CumXi_mean, 2, count_scale, `*`)

## Plot -------------------------------------------------------------
soglia = 0.5
plot_mat = template$Xi_mean #template$Xi_mean
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
       main = "Xi_aligned_mean",
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
  Z_online_mean_conditional = template$Z_mean_conditional,
  Xi_online_mean_conditional = template$Xi_mean_conditional,
  CumXi_online_mean_conditional = template$CumXi_mean_conditional
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
