# Set parameters ----------------------------------------------------------
set_param_DTM = function(H,delta,a_phi,b_phi,a_gamma,b_gamma,a_sigma,b_sigma,a_beta,b_beta, 
                         var_phi,var_gamma,var_sigma,var_beta,
                         UpdateDitl,UpdateS,UpdateLambda,UpdateXi, UpdateU,
                         UpdatePhi,UpdateGamma,UpdateSigma,UpdateBeta,
                         seed, print){
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

GibbsSampler_DTM = function(niter,nburn,data,param_DTM,init_DTM){
  
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
  
  res = GibbsSampler_DTM_c(niter,nburn,data,H,V,Ttot,param_DTM,init_DTM)
  res
}

# Utilities ---------------------------------------------------------------

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