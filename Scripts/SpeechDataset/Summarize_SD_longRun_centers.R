# Batch summary for center-model runs ------------------------------------

# This script reads many center-chain .rds files and produces:
# 1) a compact CSV summary;
# 2) LaTeX code for the same table, printed to console.

# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here on the VM if needed
wd = paste0(choose_wd,"DynamicFeatureAllocation/Scripts/SpeechDataset")
setwd(wd)

# Paths -------------------------------------------------------------------

downloads_dir = "C:/Users/colom/Downloads"
local_save_dir = file.path(wd, "centers_save")
summary_dir = file.path(wd, "centers_summary")

if(!dir.exists(summary_dir))
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_file = file.path(summary_dir, "centers_config_feature_summary.csv")

# Helpers ----------------------------------------------------------------

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
  
  list(
    update_recorded = TRUE,
    accept_mat = accept_mat,
    accept_rate = colMeans(accept_mat)
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
  
  list(
    K_total = K_total,
    K_by_center = K_by_center
  )
}

parse_tag_value = function(x){
  as.numeric(gsub("m","-",x,fixed = TRUE))
}

parse_center_fit_name = function(fit_name){
  cfg = sub("^cfg([0-9]+)_.*$","\\1",fit_name)
  gamma_tag = sub("^.*_gamma_([^_]+)_delta0_.*\\.rds$","\\1",fit_name)
  delta0_tag = sub("^.*_delta0_([^_]+)\\.rds$","\\1",fit_name)
  
  data.frame(
    config_id = as.integer(cfg),
    gamma = parse_tag_value(gamma_tag),
    delta0_centers = parse_tag_value(delta0_tag),
    stringsAsFactors = FALSE
  )
}

summarize_center_fit_file = function(fit_file){
  fit_name = basename(fit_file)
  meta = parse_center_fit_name(fit_name)
  
  tryCatch({
    fit = readRDS(fit_file)
    
    feature_trace = compute_feature_trace(fit)
    center_update_diag = compute_center_update_diagnostics(fit)
    
    K_by_center_mean = colMeans(feature_trace$K_by_center)
    M_fit = length(K_by_center_mean)
    
    accept_by_center = rep(NA_real_, M_fit)
    if(!is.null(center_update_diag) &&
       !identical(center_update_diag$update_recorded, FALSE)){
      accept_by_center = 100*center_update_diag$accept_rate
    }
    
    out = cbind(
      meta,
      data.frame(
        Ktot = mean(feature_trace$K_total),
        K1 = K_by_center_mean[1],
        K2 = K_by_center_mean[2],
        K3 = K_by_center_mean[3],
        acc_center1_percent = accept_by_center[1],
        acc_center2_percent = accept_by_center[2],
        acc_center3_percent = accept_by_center[3],
        stringsAsFactors = FALSE
      )
    )
    
    rm(fit, feature_trace, center_update_diag)
    gc(verbose = FALSE)
    out
  }, error = function(e){
    cbind(
      meta,
      data.frame(
        Ktot = NA_real_,
        K1 = NA_real_,
        K2 = NA_real_,
        K3 = NA_real_,
        acc_center1_percent = NA_real_,
        acc_center2_percent = NA_real_,
        acc_center3_percent = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  })
}

find_center_fit_files = function(){
  fit_dirs = unique(c(downloads_dir, local_save_dir))
  fit_files = unlist(lapply(fit_dirs, function(d) {
    if(!dir.exists(d))
      return(character(0))
    list.files(d,
               pattern = "^cfg[0-9]+_center_r10_M3_H10_gamma_.*_delta0_.*\\.rds$",
               full.names = TRUE)
  }))
  
  fit_files = fit_files[!grepl("_setup\\.rds$", fit_files)]
  fit_files = fit_files[!duplicated(basename(fit_files))]
  sort(fit_files)
}

format_latex_number = function(x){
  ifelse(is.na(x), "--", sprintf("%.3f", x))
}

make_latex_table = function(tab){
  lines = c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\scriptsize",
    "\\begin{tabular}{rrrrrrrrrr}",
    "\\hline",
    "config & $\\gamma$ & $\\delta_0$ & $K_{\\mathrm{tot}}$ & $K_1$ & $K_2$ & $K_3$ & acc$_1$ (\\%) & acc$_2$ (\\%) & acc$_3$ (\\%) \\\\",
    "\\hline"
  )
  
  body = apply(tab, 1, function(row){
    paste0(
      as.integer(row[["config_id"]]), " & ",
      format_latex_number(as.numeric(row[["gamma"]])), " & ",
      format_latex_number(as.numeric(row[["delta0_centers"]])), " & ",
      format_latex_number(as.numeric(row[["Ktot"]])), " & ",
      format_latex_number(as.numeric(row[["K1"]])), " & ",
      format_latex_number(as.numeric(row[["K2"]])), " & ",
      format_latex_number(as.numeric(row[["K3"]])), " & ",
      format_latex_number(as.numeric(row[["acc_center1_percent"]])), " & ",
      format_latex_number(as.numeric(row[["acc_center2_percent"]])), " & ",
      format_latex_number(as.numeric(row[["acc_center3_percent"]])), " \\\\"
    )
  })
  
  c(
    lines,
    body,
    "\\hline",
    "\\end{tabular}",
    "\\caption{Summary of feature counts and center-update acceptance rates across center-model configurations.}",
    "\\label{tab:center-config-summary}",
    "\\end{table}"
  )
}

# Run batch summary -------------------------------------------------------

fit_files = find_center_fit_files()
if(length(fit_files) == 0)
  stop("No center fit .rds files found in Downloads or centers_save")

config_feature_summary = do.call(rbind, lapply(fit_files, summarize_center_fit_file))
config_feature_summary = config_feature_summary[order(config_feature_summary$config_id),]

output_cols = c("config_id","gamma","delta0_centers",
                "Ktot","K1","K2","K3",
                "acc_center1_percent","acc_center2_percent","acc_center3_percent")
config_feature_summary = config_feature_summary[,output_cols]

numeric_cols = c("gamma","delta0_centers",
                 "Ktot","K1","K2","K3",
                 "acc_center1_percent","acc_center2_percent","acc_center3_percent")
config_feature_summary[numeric_cols] =
  lapply(config_feature_summary[numeric_cols], function(x) round(x, 3))

write.csv(config_feature_summary, file = summary_file, row.names = FALSE)

latex_table_code = make_latex_table(config_feature_summary)
cat(paste(latex_table_code, collapse = "\n"), "\n")
