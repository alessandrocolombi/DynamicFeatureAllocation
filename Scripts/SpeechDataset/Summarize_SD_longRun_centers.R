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

M_select = 7 # <--- choose 3 or 7
if(!M_select %in% c(3,7))
  stop("M_select must be either 3 or 7")

downloads_dir = "C:/Users/colom/Downloads"
local_save_dir = file.path(wd, "centers_save")
summary_dir = file.path(wd, "centers_summary")

if(!dir.exists(summary_dir))
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_file = file.path(summary_dir,
                         paste0("centers_config_feature_summary_M",M_select,".csv"))

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

make_center_summary_cols = function(prefix, M){
  paste0(prefix, seq_len(M))
}

make_accept_summary_cols = function(M){
  paste0("acc_center", seq_len(M), "_percent")
}

summarize_center_fit_file = function(fit_file, M_expected){
  fit_name = basename(fit_file)
  meta = parse_center_fit_name(fit_name)
  
  tryCatch({
    fit = readRDS(fit_file)
    
    feature_trace = compute_feature_trace(fit)
    center_update_diag = compute_center_update_diagnostics(fit)
    
    K_by_center_mean = colMeans(feature_trace$K_by_center)
    M_fit = length(K_by_center_mean)
    if(M_fit != M_expected)
      stop("Expected M = ",M_expected," centers, but fit has M = ",M_fit)
    
    accept_by_center = rep(NA_real_, M_fit)
    if(!is.null(center_update_diag) &&
       !identical(center_update_diag$update_recorded, FALSE)){
      accept_by_center = 100*center_update_diag$accept_rate
    }
    
    K_cols = as.data.frame(as.list(K_by_center_mean[seq_len(M_expected)]))
    names(K_cols) = make_center_summary_cols("K", M_expected)
    
    acc_cols = as.data.frame(as.list(accept_by_center[seq_len(M_expected)]))
    names(acc_cols) = make_accept_summary_cols(M_expected)
    
    out = cbind(
      meta,
      data.frame(
        Ktot = mean(feature_trace$K_total),
        stringsAsFactors = FALSE
      ),
      K_cols,
      acc_cols
    )
    
    rm(fit, feature_trace, center_update_diag)
    gc(verbose = FALSE)
    out
  }, error = function(e){
    K_cols = as.data.frame(as.list(rep(NA_real_, M_expected)))
    names(K_cols) = make_center_summary_cols("K", M_expected)
    
    acc_cols = as.data.frame(as.list(rep(NA_real_, M_expected)))
    names(acc_cols) = make_accept_summary_cols(M_expected)
    
    cbind(
      meta,
      data.frame(
        Ktot = NA_real_,
        stringsAsFactors = FALSE
      ),
      K_cols,
      acc_cols
    )
  })
}

find_center_fit_files = function(M){
  fit_dirs = unique(c(downloads_dir, local_save_dir))
  fit_pattern = paste0("^cfg[0-9]+_center_r10_M",M,
                       "_H10_gamma_.*_delta0_.*\\.rds$")
  
  fit_files = unlist(lapply(fit_dirs, function(d) {
    if(!dir.exists(d))
      return(character(0))
    list.files(d,
               pattern = fit_pattern,
               full.names = TRUE)
  }))
  
  fit_files = fit_files[!grepl("_setup\\.rds$", fit_files)]
  fit_files = fit_files[!duplicated(basename(fit_files))]
  sort(fit_files)
}

format_latex_number = function(x){
  ifelse(is.na(x), "--", sprintf("%.3f", x))
}

make_latex_table = function(tab, M){
  K_cols = make_center_summary_cols("K", M)
  acc_cols = make_accept_summary_cols(M)
  
  header = c(
    "config",
    "$\\gamma$",
    "$\\delta_0$",
    "$K_{\\mathrm{tot}}$",
    paste0("$K_",seq_len(M),"$"),
    paste0("acc$_",seq_len(M),"$ (\\%)")
  )
  
  lines = c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\scriptsize",
    paste0("\\begin{tabular}{",paste(rep("r", length(header)), collapse = ""),"}"),
    "\\hline",
    paste0(paste(header, collapse = " & "), " \\\\"),
    "\\hline"
  )
  
  body = apply(tab, 1, function(row){
    fields = c(
      as.integer(row[["config_id"]]),
      format_latex_number(as.numeric(row[["gamma"]])),
      format_latex_number(as.numeric(row[["delta0_centers"]])),
      format_latex_number(as.numeric(row[["Ktot"]])),
      vapply(K_cols, function(x) format_latex_number(as.numeric(row[[x]])), character(1)),
      vapply(acc_cols, function(x) format_latex_number(as.numeric(row[[x]])), character(1))
    )
    paste0(paste(fields, collapse = " & "), " \\\\")
  })
  
  c(
    lines,
    body,
    "\\hline",
    "\\end{tabular}",
    paste0("\\caption{Summary of feature counts and center-update acceptance rates across center-model configurations with $M = ",M,"$.}"),
    paste0("\\label{tab:center-config-summary-M",M,"}"),
    "\\end{table}"
  )
}

# Run batch summary -------------------------------------------------------

fit_files = find_center_fit_files(M_select)
if(length(fit_files) == 0)
  stop("No center fit .rds files found for M = ",M_select,
       " in Downloads or centers_save")

config_feature_summary = do.call(rbind, lapply(fit_files, summarize_center_fit_file,
                                               M_expected = M_select))
config_feature_summary = config_feature_summary[order(config_feature_summary$config_id),]

output_cols = c("config_id","gamma","delta0_centers",
                "Ktot",
                make_center_summary_cols("K", M_select),
                make_accept_summary_cols(M_select))
config_feature_summary = config_feature_summary[,output_cols]

numeric_cols = c("gamma","delta0_centers",
                 "Ktot",
                 make_center_summary_cols("K", M_select),
                 make_accept_summary_cols(M_select))
config_feature_summary[numeric_cols] =
  lapply(config_feature_summary[numeric_cols], function(x) round(x, 3))

write.csv(config_feature_summary, file = summary_file, row.names = FALSE)

latex_table_code = make_latex_table(config_feature_summary, M_select)
cat(paste(latex_table_code, collapse = "\n"), "\n")
