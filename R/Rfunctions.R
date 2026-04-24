# Set parameters ----------------------------------------------------------
set_param_DTM = function(H,delta,a_phi,b_phi,a_gamma,b_gamma,a_sigma,b_sigma,a_beta,b_beta, 
                         var_phi,var_gamma,var_sigma,var_beta,
                         UpdateDitl,UpdateS,UpdateLambda,UpdateXi, UpdateU,
                         UpdatePhi,UpdateGamma,UpdateSigma,UpdateBeta,
                         seed, print, JointAdp){
  return(list(
    "H" = H,
    "delta" = delta,
    "a_phi" = a_phi,
    "b_phi" = b_phi,
    "a_gamma" = a_gamma,
    "b_gamma" = b_gamma,
    "a_sigma" = a_sigma,
    "b_sigma" = b_sigma,
    "a_beta" = a_beta,
    "b_beta" = b_beta,
    "var_phi" = var_phi,
    "var_gamma" = var_gamma,
    "var_sigma" = var_sigma,
    "var_beta" = var_beta,
    "UpdateDitl" = UpdateDitl,
    "UpdateS" = UpdateS,
    "UpdateLambda" = UpdateLambda,
    "UpdateXi" = UpdateXi,
    "UpdateU" = UpdateU,
    "UpdatePhi" = UpdatePhi,
    "UpdateGamma" = UpdateGamma,
    "UpdateSigma" = UpdateSigma,
    "UpdateBeta" = UpdateBeta,
    "seed" = seed,
    "JointAdp" = JointAdp,
    "print" = print
  ))
}

set_init_DTM = function(Xi0,Lambda0,S0,phi0,gamma0,sigma0,beta0){
  return(list(
    "Xi0" = Xi0,
    "Lambda0" = Lambda0,
    "S0" = S0,
    "phi0" = phi0,
    "gamma0" = gamma0,
    "sigma0" = sigma0,
    "beta0" = beta0
  ))
}

# Gibbs sampler -----------------------------------------------------------

GibbsSampler_DTM = function(niter,nburn,thin,data,param_DTM,init_DTM){
  
  if(!is.matrix(data))
    stop("data must be a matrix of size V x T")
  V = nrow(data)
  Ttot = ncol(data)
  
  # Checks
  if(param_DTM$H <= 1 || param_DTM$a_phi <=0 || param_DTM$b_phi <=0 ||param_DTM$a_gamma <=0 || param_DTM$b_gamma <=0 || param_DTM$a_sigma <=0 || param_DTM$b_sigma <=0 || param_DTM$a_beta <=0 || param_DTM$b_beta <=0 || 
     param_DTM$var_phi <=0 || param_DTM$var_gamma <=0 || param_DTM$var_sigma <=0 || param_DTM$var_beta <=0 )
    stop("Invalid param_DTM")
  
  H = param_DTM$H
  if(nrow(init_DTM$Xi0) != H || ncol(init_DTM$Xi0) != Ttot)
    stop("Invalid Xi0")
  if(length(init_DTM$Lambda0)!= H || nrow(init_DTM$Lambda0[[1]]) != V || ncol(init_DTM$Lambda0[[1]]) != Ttot  )
    stop("Invalid Lambda0")
  if(nrow(init_DTM$S0) != H || ncol(init_DTM$S0) != Ttot)
    stop("Invalid S0")
  if( init_DTM$phi0 <= 0 || init_DTM$gamma0 <= 0 || init_DTM$sigma0 <= 0 || init_DTM$beta0 <= 0 )
    stop("Invalid values of hyperparameters")
  
  res = GibbsSampler_DTM_c(niter,nburn,thin,data,H,V,Ttot,param_DTM,init_DTM)
  res
}

# Utilities ---------------------------------------------------------------

find_indices <- function(Z_l) {
  Ttot <- length(Z_l)
  
  idx_born  <- integer(0)
  idx_surv  <- integer(0)
  idx_noact <- integer(0)
  
  for (t in seq_len(Ttot)) {
    if (Z_l[t] > 0) {
      # active
      if (t == 1 || Z_l[t - 1] == 0) {
        idx_born <- c(idx_born, t)
      } else {
        idx_surv <- c(idx_surv, t)
      }
    } else {
      # no activity
      idx_noact <- c(idx_noact, t)
    }
  }
  
  list(
    idx_born  = idx_born,
    idx_surv  = idx_surv,
    idx_noact = idx_noact
  )
}

zipfs_decay = function(n,a){
  sapply(1:n,function(i){i^{-a}})
}

set_par_invgamma = function(media, var)
{
  a = 2 + media^2/var
  b = media*(a-1)
  out = c(a,b)
  return(out)
}
set_par_gamma = function(media, var)
{
  b = media/var
  a = media*b
  out = c(a,b)
  return(out)
}


build_topic_matrices <- function(Xi_it, Ttot) {
  # Xi_it: H x Ttot
  paths <- apply(Xi_it, 1, find_indices)
  
  K_now <- sum(vapply(paths, function(p) length(p$idx_born), integer(1)))
  
  Activity <- matrix(0L, nrow = Ttot, ncol = K_now)
  Xi_star  <- matrix(0,  nrow = K_now, ncol = Ttot)
  
  topic_id <- 0
  
  for (l in seq_along(paths)) {
    p <- paths[[l]]
    
    if (length(p$idx_born) > 0) {
      for (j in seq_along(p$idx_born)) {
        topic_id <- topic_id + 1
        
        tj <- p$idx_born[j]
        tj_next <- p$idx_born[j + 1]
        
        if (is.na(tj_next)) {
          tj_next <- Ttot + 1
        }
        
        time_act <- c(tj, intersect(tj:tj_next, p$idx_surv))
        
        Activity[time_act, topic_id] <- 1L
        Xi_star[topic_id, time_act] <- Xi_it[l, time_act]
      }
    }
  }
  
  list(
    Activity = Activity,   # Ttot x K_now
    Xi_star  = Xi_star    # K_now x Ttot
    # paths    = paths,
    # K        = K_now
  )
}
## Examples: (a) single iteration: build_topic_matrices(fit$Xi[[it]], Ttot)
## Examples: (a) all iteration: lapply(1:length(fit$Xi), function(it) { build_topic_matrices(fit$Xi[[it]], Ttot) })

