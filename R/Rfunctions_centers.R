# Parameters for DTM with random centers ---------------------------------

set_param_DTM_centers = function(H,gamma,sigma,beta,
                                 a_phi_centers,b_phi_centers,delta0_centers,
                                 omega,a_omega,b_omega,
                                 var_phi_centers,var_delta_centers,mstar_max,
                                 UpdateDitl,UpdateS,UpdateLambda,UpdateXi,UpdateU,
                                 UpdateCenters,UpdateOmega,
                                 seed,print,
                                 phi_centers = 1){
  return(list(
    "H" = H,
    "gamma" = gamma,
    "sigma" = sigma,
    "beta" = beta,
    "phi_centers" = phi_centers,
    "a_phi_centers" = a_phi_centers,
    "b_phi_centers" = b_phi_centers,
    "delta0_centers" = delta0_centers,
    "omega" = omega,
    "a_omega" = a_omega,
    "b_omega" = b_omega,
    "var_phi_centers" = var_phi_centers,
    "var_delta_centers" = var_delta_centers,
    "mstar_max" = mstar_max,
    "UpdateDitl" = UpdateDitl,
    "UpdateS" = UpdateS,
    "UpdateLambda" = UpdateLambda,
    "UpdateXi" = UpdateXi,
    "UpdateU" = UpdateU,
    "UpdateCenters" = UpdateCenters,
    "UpdateOmega" = UpdateOmega,
    "seed" = seed,
    "print" = print
  ))
}

set_init_DTM_centers = function(Xi0,Lambda0,S0,Zeta0,U0=NULL,omega0=NULL,Mstar0=NULL){
  out = list(
    "Xi0" = Xi0,
    "Lambda0" = Lambda0,
    "S0" = S0,
    "Zeta0" = Zeta0
  )
  if(!is.null(U0))
    out$U0 = U0
  if(!is.null(omega0))
    out$omega0 = omega0
  return(out)
}

# Gibbs sampler with random centers --------------------------------------

GibbsSampler_DTM_centers = function(niter,nburn,thin,data,param_DTM_centers,init_DTM_centers){

  if(!is.matrix(data))
    stop("data must be a matrix of size V x T")
  V = nrow(data)
  Ttot = ncol(data)

  H = param_DTM_centers$H
  if(is.null(param_DTM_centers$phi_centers))
    param_DTM_centers$phi_centers = 1
  
  if(H <= 1 || param_DTM_centers$gamma <= 0 || param_DTM_centers$beta <= 0 ||
     param_DTM_centers$sigma <= 0 || param_DTM_centers$sigma >= 1 ||
     param_DTM_centers$phi_centers <= 0 ||
     param_DTM_centers$delta0_centers <= 0 || param_DTM_centers$omega <= 0 ||
     param_DTM_centers$a_omega <= 0 || param_DTM_centers$b_omega <= 0 ||
     param_DTM_centers$var_delta_centers <= 0 ||
     param_DTM_centers$mstar_max < 0)
    stop("Invalid param_DTM_centers")

  if(!is.list(init_DTM_centers$Xi0) || !is.list(init_DTM_centers$Lambda0) ||
     !is.list(init_DTM_centers$S0) || !is.list(init_DTM_centers$Zeta0))
    stop("Xi0, Lambda0, S0 and Zeta0 must be lists")

  M0 = length(init_DTM_centers$Xi0)
  if(M0 <= 0)
    stop("Invalid number of initial centers")
  if(length(init_DTM_centers$Lambda0) != M0 || length(init_DTM_centers$S0) != M0)
    stop("Xi0, Lambda0 and S0 must have the same length")
  if(length(init_DTM_centers$Zeta0) != M0)
    stop("fixed-M sampler requires length(Zeta0) == length(Xi0)")
  if(!is.null(init_DTM_centers$Mstar0) && init_DTM_centers$Mstar0 != 0)
    stop("Mstar0 is bypassed in the fixed-M sampler and must be NULL or 0")

  for(m in seq_len(M0)){
    if(nrow(init_DTM_centers$Xi0[[m]]) != H || ncol(init_DTM_centers$Xi0[[m]]) != Ttot)
      stop("Invalid Xi0")
    if(nrow(init_DTM_centers$S0[[m]]) != H || ncol(init_DTM_centers$S0[[m]]) != Ttot)
      stop("Invalid S0")
    if(length(init_DTM_centers$Lambda0[[m]]) != H)
      stop("Invalid Lambda0")
    for(l in seq_len(H)){
      if(nrow(init_DTM_centers$Lambda0[[m]][[l]]) != V ||
         ncol(init_DTM_centers$Lambda0[[m]][[l]]) != Ttot)
        stop("Invalid Lambda0")
    }
  }

  if(!is.null(init_DTM_centers$U0)){
    if(length(init_DTM_centers$U0) != M0)
      stop("Invalid U0")
    for(m in seq_len(M0)){
      if(nrow(init_DTM_centers$U0[[m]]) != H || ncol(init_DTM_centers$U0[[m]]) != Ttot)
        stop("Invalid U0")
    }
  }

  for(m in seq_along(init_DTM_centers$Zeta0)){
    if(length(init_DTM_centers$Zeta0[[m]]) != V || any(init_DTM_centers$Zeta0[[m]] <= 0))
      stop("Invalid Zeta0")
  }

  GibbsSampler_DTM_centers_c(niter,nburn,thin,data,H,V,Ttot,param_DTM_centers,init_DTM_centers)
}



# Helpers -----------------------------------------------------------------

smooth_simplex = function(x, eps = 1e-4){
  x = x + eps
  x/sum(x)
}

rdirichlet_vec = function(alpha){
  x = rgamma(length(alpha), shape = alpha, rate = 1)
  x/sum(x)
}

find_indices <- function(Z_l) {
  Ttot <- length(Z_l)
  
  idx_born  <- integer(0)
  idx_surv  <- integer(0)
  idx_noact <- integer(0)
  
  for (t in seq_len(Ttot)) {
    if (Z_l[t] > 0) {
      if (t == 1 || Z_l[t - 1] == 0) {
        idx_born <- c(idx_born, t)
      } else {
        idx_surv <- c(idx_surv, t)
      }
    } else {
      idx_noact <- c(idx_noact, t)
    }
  }
  
  list(
    idx_born  = idx_born,
    idx_surv  = idx_surv,
    idx_noact = idx_noact
  )
}

build_topic_matrices <- function(Xi_it, Ttot) {
  # Xi_it: H x Ttot matrix for one center at one MCMC iteration.
  # Every positive run of a row is converted into one topic path.
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
    Activity = Activity,
    Xi_star  = Xi_star
  )
}

left_order <- function(Z) {
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  if (!all(Z %in% c(0, 1)))
    stop("Z must be a binary matrix (0/1).")
  
  col_signature <- apply(Z, 2, function(col) paste(col, collapse = ""))
  ord <- order(col_signature, decreasing = TRUE)
  
  Z[,ord, drop = FALSE]
}

