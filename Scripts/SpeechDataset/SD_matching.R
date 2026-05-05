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
# We use `clue::solve_LSAP()` to solve the linear assignment problem
# (Hungarian matching). The script stops with a clear error if `clue`
# is not installed.
if(!requireNamespace("clue", quietly = TRUE))
  stop("Package 'clue' is required. Please install it with install.packages('clue').")


eps_break = .Machine$double.eps
# Input / output ----------------------------------------------------------
#
# This script works on the compact summary objects saved in `save_summary/`.
# By default, it uses the first *_summary.rds file found in that folder.
# If you want a specific file, set `summary_file_name` manually.
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

# Matching strategy -------------------------------------------------------
#
# We align every iteration to the FINAL iteration.
# Each column is a topic-specific trajectory represented by the pair:
#   - Zstar[, k]  : binary activity pattern
#   - Xistar[, k] : intensity path
#
# Matching is done by minimizing a cost that combines:
#   1. distance between birth times
#   2. distance between death times
#   3. Euclidean distance between normalized Xi trajectories
#
# Since the support is fully determined by birth/death under the model
# constraints, we do not include an extra Hamming term on Zstar.
w_birth = 1
w_death = 1
w_xi = 1
unmatch_cost = 25

# Helpers -----------------------------------------------------------------

# Convert one topic object into the two Ttot x K matrices used below.
get_topic_pair = function(topic_obj) {
  Zstar = topic_obj$Activity
  Xistar = t(topic_obj$Xi_star)
  
  if(!all(dim(Zstar) == dim(Xistar)))
    stop("Zstar and Xistar do not have matching dimensions.")
  
  list(Zstar = Zstar, Xistar = Xistar)
}

# For Xi comparison we use normalized trajectories, so that the distance
# is about SHAPE rather than raw total mass.
normalize_xi = function(x) {
  sx = sum(x)
  if(sx <= 0) {
    rep(0, length(x))
  } else {
    x / sx
  }
}

# Extract the "life interval" of one topic column.
# Because topics should not reappear after death, these two times are
# meaningful structural summaries.
extract_birth_death = function(z_col) {
  idx = which(z_col > 0)
  
  if(length(idx) == 0) {
    list(birth = Inf, death = -Inf)
  } else {
    list(birth = idx[1], death = idx[length(idx)])
  }
}

# Build the K_ref x K_cur cost matrix.
# Lower cost means "better match".
build_cost_matrix = function(refZ, refXi, curZ, curXi,
                             w_birth, w_death, w_xi) {
  K_ref = ncol(refZ)
  K_cur = ncol(curZ)
  
  cost = matrix(0, nrow = K_ref, ncol = K_cur)
  
  ref_bd = lapply(seq_len(K_ref), function(k) extract_birth_death(refZ[, k]))
  cur_bd = lapply(seq_len(K_cur), function(k) extract_birth_death(curZ[, k]))
  
  ref_xi_norm = lapply(seq_len(K_ref), function(k) normalize_xi(refXi[, k]))
  cur_xi_norm = lapply(seq_len(K_cur), function(k) normalize_xi(curXi[, k]))
  
  for(i in seq_len(K_ref)) {
    for(j in seq_len(K_cur)) {
      db = abs(ref_bd[[i]]$birth - cur_bd[[j]]$birth)
      dd = abs(ref_bd[[i]]$death - cur_bd[[j]]$death)
      dxi = sqrt(sum((ref_xi_norm[[i]] - cur_xi_norm[[j]])^2))
      
      cost[i, j] = w_birth * db +
                   w_death * dd +
                   w_xi * dxi
    }
  }
  
  cost
}

# Align one iteration to the reference iteration.
#
# We keep the reference number of columns, K_ref.
# If a reference topic is unmatched, we leave that aligned column at zero.
align_to_reference = function(refZ, refXi, curZ, curXi,
                              w_birth, w_death, w_xi,
                              unmatch_cost) {
  K_ref = ncol(refZ)
  K_cur = ncol(curZ)
  
  cost_main = build_cost_matrix(
    refZ = refZ, refXi = refXi,
    curZ = curZ, curXi = curXi,
    w_birth = w_birth,
    w_death = w_death,
    w_xi = w_xi
  )
  
  # We append one dummy column per reference topic. Assigning row i
  # to its own dummy column means "leave reference slot i unmatched".
  big_penalty = 1e8
  cost_dummy = matrix(big_penalty, nrow = K_ref, ncol = K_ref)
  diag(cost_dummy) = unmatch_cost
  
  cost_aug = cbind(cost_main, cost_dummy)
  assignment = clue::solve_LSAP(cost_aug)
  
  Z_aligned = matrix(0L, nrow = nrow(curZ), ncol = K_ref)
  Xi_aligned = matrix(0, nrow = nrow(curXi), ncol = K_ref)
  
  matched_current = rep(NA_integer_, K_ref)
  matched_cost = rep(NA_real_, K_ref)
  
  for(i in seq_len(K_ref)) {
    j = assignment[i]
    matched_cost[i] = cost_aug[i, j]
    
    # Real topic match
    if(j <= K_cur) {
      Z_aligned[, i] = curZ[, j]
      Xi_aligned[, i] = curXi[, j]
      matched_current[i] = j
    }
  }
  
  list(
    Z_aligned = Z_aligned,
    Xi_aligned = Xi_aligned,
    matched_current = matched_current,
    matched_cost = matched_cost,
    assignment = assignment,
    cost_aug = cost_aug
  )
}

# Reference iteration -----------------------------------------------------
#
# The user requested the FINAL iteration as reference.
ref_it = Niter
ref_pair = get_topic_pair(fit_summary$topic_objs[[ref_it]])
refZ = ref_pair$Zstar
refXi = ref_pair$Xistar
K_ref = ncol(refZ)

cat("\nReference iteration: ", ref_it, "\n", sep = "")
cat("Reference number of topics: ", K_ref, "\n", sep = "")

# Align all iterations ----------------------------------------------------
Nstart = Niter - 100
Lsummary = Niter - Nstart
aligned_list = vector("list", Niter-Nstart)

for(it in 1:Lsummary) {
  cur_pair = get_topic_pair(fit_summary$topic_objs[[it]])
  
  aligned_list[[it]] = align_to_reference(
    refZ = refZ,
    refXi = refXi,
    curZ = cur_pair$Zstar,
    curXi = cur_pair$Xistar,
    w_birth = w_birth,
    w_death = w_death,
    w_xi = w_xi,
    unmatch_cost = unmatch_cost
  )
}

# Posterior summaries after alignment ------------------------------------
#
# After alignment, every iteration has the same number of columns K_ref.
# This makes posterior averaging meaningful.
Z_aligned_list = lapply(aligned_list, function(x) x$Z_aligned)
Xi_aligned_list = lapply(aligned_list, function(x) x$Xi_aligned)

Z_aligned_mean = Reduce(`+`, Z_aligned_list) / length(Z_aligned_list)
Xi_aligned_mean = Reduce(`+`, Xi_aligned_list) / length(Xi_aligned_list)




suppressWarnings({
  mycol = hcl.colors(n = 100, palette = "Greens", rev = TRUE)
})

soglia = 15.5
sel_colums = which(colSums(Xi_aligned_mean) > soglia )
plot_mat = Xi_aligned_mean[,sel_colums]
if(!is.matrix(plot_mat))
  plot_mat = matrix(plot_mat, ncol = length(sel_colums))

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


# A compact summary of how often each aligned slot is occupied.
slot_occupancy = sapply(Z_aligned_list, function(Z_it) colSums(Z_it) > 0)
if(is.vector(slot_occupancy))
  slot_occupancy = matrix(slot_occupancy, nrow = 1)
slot_match_prob = colMeans(t(slot_occupancy))

# Save result -------------------------------------------------------------
out = list(
  meta = list(
    summary_file = basename(summary_file),
    ref_it = ref_it,
    K_ref = K_ref,
    Ttot = Ttot,
    Niter = Niter,
    weights = c(
      w_birth = w_birth,
      w_death = w_death,
      w_xi = w_xi,
      unmatch_cost = unmatch_cost
    )
  ),
  refZ = refZ,
  refXi = refXi,
  aligned_list = aligned_list,
  Z_aligned_list = Z_aligned_list,
  Xi_aligned_list = Xi_aligned_list,
  Z_aligned_mean = Z_aligned_mean,
  Xi_aligned_mean = Xi_aligned_mean,
  slot_match_prob = slot_match_prob
)

out_file = file.path(
  save_summary_dir,
  sub("_summary\\.rds$", "_aligned_ref_last.rds", basename(summary_file))
)

saveRDS(out, file = out_file)

cat("\nSaved aligned object:\n", basename(out_file), "\n", sep = "")
cat("Number of aligned slots: ", K_ref, "\n", sep = "")
cat("Average slot occupancy:\n")
print(round(slot_match_prob, 3))
