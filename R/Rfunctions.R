HelloR = function(){
  cat("\n Hello from R function \n")
}



# Utilities ---------------------------------------------------------------


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


plot_img = function(Mat, 
                    center_value = NULL, 
                    col.upper = "grey5", col.center = "grey45", col.lower = "grey95", 
                    col.n_breaks = 59, 
                    remove_diag = FALSE, 
                    main = " ", x_label = " ", y_label = " ", 
                    horizontal = TRUE) 
{
  if (any(is.na(Mat))) {
    cat("\n NA values have been removed from Matrix  \n")
  }
  if (col.n_breaks%%2 == 0) {
    warning("col.n_breaks is even but it has to be odd. Adding 1")
    col.n_breaks = col.n_breaks + 1
  }
  colorTable = fields::designer.colors(col.n_breaks, c(col.lower, 
                                                       col.center, col.upper))
  col_length = (col.n_breaks + 1)/2
  if (remove_diag) {
    diag(Mat) = NA
  }
  min_val = min(Mat, na.rm = T)
  max_val = max(Mat, na.rm = T)
  p_row = dim(Mat)[1]
  p_col = dim(Mat)[2]
  if (!is.null(center_value)) {
    brks = c(seq(min_val, center_value - 1e-04, l = col_length), 
             seq(center_value + 1e-04, max_val, l = col_length))
  } else {
    brks = seq(min_val, max_val, l = 2 * col_length)
  }
  if(brks[1] == tail(brks,1) ){
    brks = seq(min(brks)- 1e-10,max(brks)+1e-10,length.out = length(brks))
  }
  colnames(Mat) = 1:p_col
  rownames(Mat) = 1:p_row

  # par(mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))
  par(mar = c(0,0,2,4), mgp = c(1,0, 0))
  if (horizontal) {
    fields::image.plot(Mat, axes = F, horizontal = T, main = main, 
                       col = colorTable, breaks = brks, xlab = x_label, 
                       legend.width = 0.5,    # width of the legend bar
                       bigplot = c(0.05, 0.80, 0.1, 0.9),  # Adjusts the main plot to leave less room
                       smallplot = c(0.82, 0.84, 0.1, 0.9),
                       ylab = y_label)
  }
  else {
    fields::image.plot(Mat, axes = F, horizontal = FALSE, 
                       main = main, col = colorTable, breaks = brks, xlab = x_label, 
                       legend.width = 0.5,    # width of the legend bar
                       bigplot = c(0.05, 0.80, 0.1, 0.9),  # Adjusts the main plot to leave less room
                       smallplot = c(0.82, 0.84, 0.1, 0.9),
                       ylab = y_label)
  }
  box()
}



plot_v2m = function(x, di = 8, center_val_plot = NULL){
  par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
  mean_t = matrix(data = x, nrow = di, ncol = di, byrow = F)
  plot_img(mean_t, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
           horizontal = F, main = " ")
}

plot_multiple_imgs = function(mat, plt_wnd = c(3,4),center_value = NULL){
  rows = nrow(mat)
  cols = ncol(mat)
  di = sqrt(cols)
  par(mfrow = plt_wnd, mar = c(2,2,2,1), bty = "l")
  for(i in 1:rows){
    mean_t = matrix(data = mat[i,], nrow = di, ncol = di, byrow = F)
    plot_img(mean_t, center_value = NULL, col.lower = "grey95",col.upper = "grey10",
             horizontal = F, main = paste0(" i = ",i))  
  }
}


# This function maps a positive integer number in a relative number (integer with sign)
map_to_Z = function(n){
  if(n %% 2 == 0)
    n/2
  else
   -(n+1)/2
}
# This function maps a relative number (integer with sign) in a positive integer number 
map_to_N = function(z){
  if( z < 0 )
    (2*abs(z) - 1)
  else
    2*abs(z)
}

# BetaBernoulli - Moving images --------------------------------------------------
computeM1 = function(x){
  x[1] * exp( lgamma(x[3]+1) + lgamma(x[2]-x[3]+1) - lgamma(x[2]+2) )
}
computePsurv = function(x){
  (1-x[3])/(1+x[2])
}

stat_suff_Z = function(Zmat){
  Ttot = nrow(Zmat)
  K    = ncol(Zmat)
  if(K == 0)
    return(c(0,0,0))
  fun_temp  = function(x){
    x = x[1:(Ttot-1)]
    res = c(0,0)
    time_birth = which(x == 1)[1]
    if( is.na(time_birth) ){
      # j-th features appears at time T-1
      return(res)
    }
    if(time_birth >= 1){
      x = x[time_birth:length(x)]
      time_death = which(x == 0)[1]
      if( is.na(time_death) )
        time_death = length(x)
      x = x[1:time_death]
      res[1] = length(x) - 1 # first appearance is innovation, not thinning
      res[2] = sum(x) - 1    # first appearance is innovation, not thinning
      return(res)
    }else{
      stop("Error in stat_suff_Z")
    }
  }
  mat_temp = apply(Zmat, 2, fun_temp)
  row.names(mat_temp) = c("trials","successes")
  res = colSums(t(mat_temp))
  res = c(res,res[1]-res[2])
  names(res)[3] = "failures"
  return(res)
}




log_marginal = function(x,sigma2_A,sigma2_X,Nnew){
  D = length(x)
  eval = 0
  sigma_tot = sigma2_X + sigma2_A*Nnew
  for(j in 1:D){
    eval = eval + log_dnorm(x[j],mean = 0, sd = sqrt(sigma_tot))
  }
  eval
}


## Static quantities

# gamma
r_fullcond_gamma = function(a_gam, b_gam, NnewTot, Ttot, c, sigma){
  return(rgamma(1, shape = a_gam + NnewTot, 
                   rate  = b_gam + Ttot * (gamma(c+sigma+1)*gamma(1+sigma))/(gamma(c+2))
                )
        )
}

# c
log_fullcond_c = function(x, a_c, b_c, 
                          gamma, sigma, 
                          NumTrailTot, NumFailTot, NnewTot, Ttot) 
{
  if( (x+sigma+1) < 0 )
    stop("Error in log_fullcond_c")
  res = (a_c + NumFailTot - 1)*log(x)
  res = res + NumTrailTot * log(1+x)
  res = res + NnewTot * ( lgamma(x+sigma+1) - lgamma(x+2) )
  res = res - b_c*x - gamma*gamma(1+sigma)*Ttot*gamma(x+sigma+1)/lgamma(x+2)
  return(res)
}
r_fullcond_c = function(c_old, a_c, b_c, 
                        gamma, sigma,
                        NumTrailTot, NumFailTot, NnewTot, Ttot,
                        var_prop){

  # Metropolis Hastings
  log_c_old = log(c_old)
  log_c_new = log_c_old + rnorm( n = 1, mean = 0, sd = sqrt(var_prop) )
  c_new     = exp(log_c_new)
  
  if( (c_new+sigma) <= 0 )
    return(c_old)
  
  full_cond_old = log_fullcond_c(c_old, a_c, b_c, 
                                 gamma, sigma, 
                                 NumTrailTot, NumFailTot, NnewTot, Ttot)
  full_cond_new = log_fullcond_c(c_new, a_c, b_c, 
                                 gamma, sigma, 
                                 NumTrailTot, NumFailTot, NnewTot, Ttot)
  
  log_acceptance = full_cond_new - full_cond_old + log_c_new - log_c_old
  log_acceptance = min(log_acceptance, 0)
  log_u = log( runif(1) )
  if(log_u < log_acceptance)
    return(c_new)
  else
    return(c_old)
}


# sigma2_x
r_fullcond_sigma2_X = function(a_x, b_x, D, Ttot, suff_stat)
{
  shape_post = a_x + 0.5*(D*Ttot)
  rate_post  = b_x + 0.5*suff_stat
  res = rgamma(1,shape = shape_post, rate = rate_post)
  1/res
}

# sigma2_A
r_fullcond_sigma2_A = function(a_A, b_A, D, K, suff_stat)
{
  shape_post = a_A + 0.5*(D*K)
  rate_post  = b_A + 0.5*suff_stat
  res = rgamma(1,shape = shape_post, rate = rate_post)
  1/res
}




# Proposal in thinning process
propose_thinning_Be = function(Xobs, ActFeat, sig2X, Psurv){
  D = ncol(ActFeat)
  p = nrow(ActFeat)
  survived = sample(c(0,1), size = p, replace = TRUE, prob = c(1-Psurv,Psurv))
  n_surv = sum(survived)
  log_prop_move = n_surv * log(Psurv) + (p-n_surv) * log(1-Psurv)
  return(list("survived" = survived, "log_prop_move" = log_prop_move))
}
propose_thinning_VS = function(Xobs, ActFeat, sig2X, Psurv){
  D = ncol(ActFeat)
  p = nrow(ActFeat)
  active = rep(1, p)
  log_prop_move = 0
  for(j in 1:p){
    active_prop_1 <- active_prop_0 <- active
    active_prop_1[j] = 1;active_prop_0[j] = 0
    mean_1 = active_prop_1 %*% ActFeat
    mean_0 = active_prop_0 %*% ActFeat
    log_prob_1 <- log_dmvnorm(x = Xobs, mean = mean_1, Sigma = sig2X*diag(D)) - 2*D # score if I keep it (adjusted)
    log_prob_0 <- log_dmvnorm(x = Xobs, mean = mean_0, Sigma = sig2X*diag(D)) # score if I drop it
    log_den = log_stable_sum( c(log_prob_1,log_prob_0), is_log = TRUE )
    log_prob_active_j = log_prob_1 - log_den  
    if(log(runif(1)) < log_prob_active_j){
      active_j = 1
      log_prop_move = log_prop_move + log_prob_active_j
    }else{
      active_j = 0
      log_prop_move = log_prop_move + log_prob_0 - log_den
    }
    
    active[j] = active_j
    
  }
  return(list("survived" = active, "log_prop_move" = log_prop_move))
}
choose_propose_thinning <- function(use_VS = TRUE) {
  if(use_VS){
    return(propose_thinning_VS)
  }else{
    return(propose_thinning_Be)
  }
  
}




SeqMonteCarlo = function(X,N,D,Ttot,
                         M1,Psurv,sigma2_A,sigma2_X,
                         mu0 = rep(0,D),
                         use_VS = TRUE, 
                         seed0 = 42){
  
  thinning_func = choose_propose_thinning(use_VS)
  # Define structure to store particles during the filtering
  base_list = list("active_values" = c(),"Nnew" = c(), "survived" = c())
  Particles = lapply(1:N, function(x){lapply(1:Ttot, function(xx){base_list})})
  # Particles[[k]][[t]]$active_values   : matrix with the active features at time t for particle k
  # Particles[[k]][[t]]$Nnew     : number of new features at time t for particle k
  # Particles[[k]][[t]]$survived : indeces of the features of parent particle survived from time t-1 to time t for particle k
  
  
  # Define auxiliary structures
  A <- B <- log_w <- w <- W <- matrix(NA, nrow = Ttot, ncol = N)
  # The final row of A is empty
  # A[1,]: vector with the parents of the particles at time 2, i.e., X[2,k] is generated conditionally to X[1,A[1,k]]
  

  # Time 1
  # a) Sample the particles
  t  = 1
  scale_mean_post = sigma2_A / (sigma2_A + sigma2_X)
  scale_var_post  = sigma2_A*sigma2_X / (sigma2_A + sigma2_X)
  for(k in 1:N){
    seed_kt = ceiling(seed0 * k * t  * runif(n=1,min = 1, max = 100))
    Nnew   = rpois(n = 1, lambda = M1)
    values = sample_A(Nnew, X[t,], mu0, sigma2_X, sigma2_A, seed_kt  )    
    ## OLD
    # values = matrix(0,nrow = Nnew, ncol = D)
    # for(j in 1:D){
    #   values[,j] = rnorm(n=Nnew, mean = scale_mean_post*X[t,j], sd = sqrt( scale_var_post ) ) 
    # }
    ## END OLD
    Particles[[k]][[t]]$active_values   = values
    Particles[[k]][[t]]$Nnew            = Nnew
    if(Nnew > 0)
      Particles[[k]][[t]]$survived        = seq(1,Nnew)
    
    # Compute the unnormalized weights
    ## OLD
    # log_w[t,k] = log_marginal(X[t,],sigma2_A,sigma2_X,Nnew)
    ## END OLD
    log_w[t,k] = log_dmarg_img( Nnew, X[t,], mu0, sigma2_X, sigma2_A )
  }
  
  # b) Compute the unomralized weights
  w[t,] = exp( log_w[t,] - max(log_w[t,]) )
  W[t,] = w[t,]/sum(w[t,])
  
  # Time t (t = 2,...,Ttot)
  pb = txtProgressBar(min = 1, max = Ttot, initial = 1)
  for(t in 2:Ttot){
    setTxtProgressBar(pb,t)
    # cat("\n ++++++ t = ",t," ++++++ \n")
    # a) Sample the parent node
    A[t-1,] = sample(1:N, size = N, W[t-1,], replace = TRUE)
    
    # b) Sample the particles
    for(k in 1:N){
      seed_kt = ceiling(seed0 * k * t * runif(n=1,min=1,max = 100))
      j = A[t-1,k] # parent index
      
      # Thinning part:
      num_active_old = nrow(Particles[[j]][[t-1]]$active_values)
      if( num_active_old > 0  ){
        thinning = thinning_func( Xobs = X[t,], 
                                  ActFeat = Particles[[j]][[t-1]]$active_values,
                                  sig2X = sigma2_X, Psurv = Psurv )
        
        survived = thinning$survived
        n_surv   = sum(survived)
        log_prop_thin = thinning$log_prop_move
        log_prod_bern = n_surv * log(Psurv) + (num_active_old-n_surv) * log(1-Psurv)
      }else{
        survived = c()
        n_surv = 0
        log_prop_thin <- log_prod_bern <- 0
      }
      
      # Innovation:
      Xstar  = X[t,] - colSums(Particles[[j]][[t-1]]$active_values)
      Nnew   = rpois(n = 1, lambda = M1)
      values = sample_A(Nnew, Xstar, mu0, sigma2_X, sigma2_A, seed_kt  )
      ## OLD
      # values = matrix(0,nrow = Nnew, ncol = D)
      # for(jj in 1:D){
      #   values[,jj] = rnorm(n=Nnew, mean = scale_mean_post*Xstar[jj], sd = sqrt( scale_var_post ) )
      # }
      ## END OLD
      
      
      # Assemble thinning + innovation:
      if(sum(survived > 0)){
        survived_values = Particles[[j]][[t-1]]$active_values[which(survived > 0),]
        Particles[[k]][[t]]$active_values = rbind( survived_values, values   )
      }else{
        Particles[[k]][[t]]$active_values =  values 
      }
      Particles[[k]][[t]]$Nnew     = Nnew
      Particles[[k]][[t]]$survived = which(survived > 0)
      
      # Compute the un-normalized weights
      # OLD: log_w[t,k] = log_marginal(Xstar,sigma2_A,sigma2_X,Nnew) + log_prod_bern - log_prop_thin
      log_w[t,k] = log_dmarg_img(Nnew,Xstar,mu0,sigma2_X,sigma2_A)  + log_prod_bern - log_prop_thin
    }
    
    # c) Compute the normalized weights
    w[t,] = exp( log_w[t,] - max(log_w[t,]) )
    W[t,] = w[t,]/sum(w[t,])
  }
  close(pb)
  
  
  # Find ancestral lineage and save final paths
  B[Ttot,] = 1:N
  for(k in 1:N){
    for(t in (Ttot-1):1){
      B[t,k] = A[ t, B[t+1,k] ] # B_t^k = A_t^{B_{t+1}^k}
    }
  }
  
  
  # Find all features appeared at least once
  Features_particle = lapply(1:N, function(x){c()})
  cum_Ntot = lapply(1:N,function(x){c()})
  Nnew_paths = lapply(1:N,function(x){c()})
  for(k in 1:N){
    for(t in Ttot:1){
      Mat = Particles[[ B[t,k] ]][[t]]$active_values
      nuove = Particles[[ B[t,k] ]][[t]]$Nnew
      Nnew_paths[[k]] = c(nuove, Nnew_paths[[k]])
      cum_Ntot[[k]] = c(nuove, cum_Ntot[[k]])
      if(nuove > 0){
        if(nrow(Mat)-nuove+1 <= 0)
          stop("Unexprected behaviour")
        Features_particle[[k]] = rbind( Mat[(nrow(Mat)-nuove+1):nrow(Mat),] , Features_particle[[k]])
      }
    }
    cum_Ntot[[k]] = cumsum(cum_Ntot[[k]])
  }
  
  # Find survival path of each feature for each particle
  Survived_features = vector("list",N)
  for(k in 1:N){
    Survived_features[[k]] = matrix(0, nrow = Ttot, ncol = nrow(Features_particle[[k]]))
    # Time 1:
    t = 1
    if(Nnew_paths[[k]][t] > 0)
      Survived_features[[k]][t, 1:cum_Ntot[[k]][1] ] = 1
    for(t in 2:Ttot){
      # Time t (2...Ttot)
      if(Nnew_paths[[k]][t] > 0)
        Survived_features[[k]][t, (cum_Ntot[[k]][t-1]+1):(cum_Ntot[[k]][t])] = 1
      Nsurvived_old = length(which(Survived_features[[k]][t-1,] == 1))
      submat = matrix( Survived_features[[k]][,which(Survived_features[[k]][t-1,] == 1)],
                       nrow = Ttot, ncol = Nsurvived_old)
      submat[t, Particles[[ B[t,k] ]][[t]]$survived ] = 1
      Survived_features[[k]][,which(Survived_features[[k]][t-1,] == 1)] = submat
    }
  }
  
  # Final estimated filtered density
  BB = 10000
  ind_finals = sample(1:N, size = BB, replace = TRUE, prob = W[Ttot,])
  post_means_all = vector("list",N)
  for(k in 1:N){
    post_means_all[[k]] = Survived_features[[k]]%*%Features_particle[[k]]
  }
  
  result <- Reduce(`+`, post_means_all[ind_finals])
  result = result / BB
  
  
  res = list("Particles_Zmat" = Survived_features,
             "Particles_Amat" = Features_particle,
             "Final_Weights" = W[Ttot,],
             "Final_Est" = result)
  res
}



# Set Bfix = rep(-1,Ttot);Particle_fix = NULL to run without conditioning
Conditional_SeqMonteCarlo = function( X,N,D,Ttot,
                                      Bfix, Particle_fix,  
                                      M1,Psurv,sigma2_A,sigma2_X,
                                      mu0 = rep(0,D),
                                      proposal_Nnew_1T = NULL,
                                      use_VS = TRUE,
                                      seed0 = 42)
{
  
  thinning_func = choose_propose_thinning(use_VS)
  # Define structure to store particles during the filtering
  base_list = list("active_values" = c(),"Nnew" = c(), "survived" = c())
  Particles = lapply(1:N, function(x){lapply(1:Ttot, function(xx){base_list})})
  # Particles[[k]][[t]]$active_values   : matrix with the active features at time t for particle k
  # Particles[[k]][[t]]$Nnew     : number of new features at time t for particle k
  # Particles[[k]][[t]]$survived : indices of the features of parent particle survived from time t-1 to time t for particle k
  
  
  # Define auxiliary structures
  A <- B <- log_w <- w <- W <- matrix(NA, nrow = Ttot, ncol = N)
  # The final row of A is empty
  # A[1,]: vector with the parents of the particles at time 2, i.e., X[2,k] is generated conditionally to X[1,A[1,k]]
  
  ### ---> Start filtering
  
  # Time 1
  # a) Sample the particles
  t  = 1
  # cat("\n ++++++ t = ",t," ++++++ \n")
  scale_mean_post = sigma2_A / (sigma2_A + sigma2_X)
  scale_var_post  = sigma2_A*sigma2_X / (sigma2_A + sigma2_X)
  
  # cat("\n Part. fissa: ",Bfix[t],"\n")
  for(k in 1:N){
    seed_kt = ceiling(seed0 * k * t * runif(n=1, min = 1, max = 100) )
    # cat("\n Part. #",k,"; ")
    if( (!is.null(Bfix)) && (k == Bfix[t]) ){
      # cat(" Fissa!!! ")
      # Set the k-th particle to the value we are conditioning on
      Particles[[k]][[t]] = Particle_fix[[t]]
      Nnew = Particles[[k]][[t]]$Nnew
    }else{
      
      if(!is.null(proposal_Nnew_1T)){
        cat("\n","Dentro alla nuova proposal","\n")
        disc_adapt_prop = proposal_Nnew_1T[t,]
        disc_adapt_prop = table(disc_adapt_prop)
        Kproposal_ind = sample(1:length(disc_adapt_prop), size = 1, prob = disc_adapt_prop)
        Kprop1 = as.numeric(names(disc_adapt_prop)[Kproposal_ind])
        Kprop2 = rpois(n = 1, lambda = M1)
      }else{
        Kprop1 <- Kprop2 <- rpois(n = 1, lambda = M1)
      }

      # Draw a new particle
      Nnew   = ifelse(runif(1)<0.75, Kprop1, Kprop2)
      values = sample_A(Nnew, X[t,], mu0, sigma2_X, sigma2_A, seed_kt  ) 
      Particles[[k]][[t]]$active_values   = values
      Particles[[k]][[t]]$Nnew            = Nnew
      if(Nnew > 0)
        Particles[[k]][[t]]$survived      = seq(1,Nnew)
    }
      # Compute the unnormalized weights
      log_w[t,k] = log_dmarg_img(Nnew,X[t,],mu0,sigma2_X,sigma2_A)
  }
  
  # b) Compute the nomralized weights
  w[t,] = exp( log_w[t,] - max(log_w[t,]) )
  W[t,] = w[t,]/sum(w[t,])
  
  t = 2
  # Time t (t = 2,...,Ttot)
  for(t in 2:Ttot){
    # cat("\n ++++++ t = ",t," ++++++ \n")
    # a) Sample the parent node
    A[t-1,] = sample(1:N, size = N, W[t-1,], replace = TRUE)
    if( all(Bfix > 0) )
      A[t-1, Bfix[t] ] = Bfix[t-1] # Cond. SMC
    
    # b) Sample the particles
    # cat("\n Part. fissa: ",Bfix[t],"\n")
    for(k in 1:N){
      seed_kt = ceiling(seed0 * k * t * runif(n=1, min = 1, max = 100) )
      # cat("\n Part. #",k,"; ")
      j = A[t-1,k] # parent index
      # cat(" genitore",j)
      if((!is.null(Bfix)) && (k == Bfix[t])){
        # Set the k-th particle to the value we are conditioning on
        # cat(" Fissa!!! ")
        Particles[[k]][[t]] = Particle_fix[[t]]
        Nnew = Particles[[k]][[t]]$Nnew
        log_w[t,k] = Particles[[k]][[t]]$log_w
      }else{
        # Update particle k at time t
        
        # Thinning part:
        num_active_old = nrow(Particles[[j]][[t-1]]$active_values)
        if( num_active_old > 0  ){
          thinning = thinning_func( Xobs = X[t,], 
                                    ActFeat = Particles[[j]][[t-1]]$active_values,
                                    sig2X = sigma2_X, Psurv = Psurv )
          
          survived = thinning$survived
          n_surv   = sum(survived)
          # cat("; n_surv = ",n_surv)
          log_prop_thin = thinning$log_prop_move
          log_prod_bern = n_surv * log(Psurv) + (num_active_old-n_surv) * log(1-Psurv)
          
        }else{
          survived = c()
          n_surv = 0
          log_prop_thin <- log_prod_bern <- 0
        }
        
        # Innovation: 
        Xstar  = X[t,] 
        if(n_surv > 0){
          survived_values = Particles[[j]][[t-1]]$active_values[which(survived > 0),]
          if(!is.matrix(survived_values))
            survived_values = matrix(survived_values, nrow = n_surv, ncol = D, byrow = T)
          
          Xstar  = Xstar - rep(1,n_surv)%*%survived_values
        }
        # if(!is.null( nrow(Particles[[j]][[t-1]]$active_values) ) )
          # Xstar  = Xstar - colSums(Particles[[j]][[t-1]]$active_values)

        if(!is.null(proposal_Nnew_1T)){
          disc_adapt_prop = proposal_Nnew_1T[t,]
          disc_adapt_prop = table(disc_adapt_prop)
          Kproposal_ind = sample(1:length(disc_adapt_prop), size = 1, prob = disc_adapt_prop)
          Kprop1 = as.numeric(names(disc_adapt_prop)[Kproposal_ind])
          Kprop2 = rpois(n = 1, lambda = M1)
        }else{
          Kprop1 <- Kprop2 <- rpois(n = 1, lambda = M1)
        }
        
        # Draw a new particle
        Nnew   = ifelse(runif(1)<0.75, Kprop1, Kprop2)
        
        # cat("; Nnew = ",Nnew)
        values = sample_A(Nnew, c(Xstar), mu0, sigma2_X, sigma2_A, seed_kt)
        
        # Assemble thinning + innovation:
        if(n_surv > 0){
          survived_values = Particles[[j]][[t-1]]$active_values[which(survived > 0),]
          Particles[[k]][[t]]$active_values = rbind( survived_values, values   )
        }else{
          Particles[[k]][[t]]$active_values =  values 
        }
        
        Particles[[k]][[t]]$Nnew     = Nnew
        Particles[[k]][[t]]$survived = which(survived > 0)
        
        # Compute the un-normalized weights
        log_w[t,k] = log_dmarg_img(Nnew,Xstar,mu0,sigma2_X,sigma2_A) + log_prod_bern - log_prop_thin
      }
    }
    
    # c) Compute the normalized weights
    w[t,] = exp( log_w[t,] - max(log_w[t,]) )
    W[t,] = w[t,]/sum(w[t,])
  }
  
  ### End of filtering <--- 
  
  # Select the chosen particle
  kfinal = sample(1:N, size = 1, replace = TRUE, prob = W[Ttot,])
  
  # Find ancestral lineage and save final paths
  B_k = rep(0,Ttot)
  B_k[Ttot] = kfinal
  for(t in (Ttot-1):1){
    B_k[t] = A[ t, B_k[t+1] ] # B_t^k = A_t^{B_{t+1}^k}
  }
  
  # Save selected path
  Path_k = vector("list", Ttot)
  for(t in 1:Ttot){
    Path_k[[t]] = Particles[[ B_k[t] ]][[t]]
    Path_k[[t]]$log_w = log_w[ t, B_k[t] ]
  }

  # Find all features appeared at least once
  Amat_k <- cum_Ntot_k <- Nnew_paths_k <- c()
  for(t in Ttot:1){
    Mat = Path_k[[t]]$active_values
    nuove = Path_k[[t]]$Nnew
    Nnew_paths_k = c(nuove, Nnew_paths_k)
    cum_Ntot_k   = c(nuove, cum_Ntot_k)
    if(nuove > 0){
      if(nrow(Mat)-nuove+1 <= 0)
        stop("Unexprected behaviour")
      Amat_k = rbind( Mat[(nrow(Mat)-nuove+1):nrow(Mat),] , Amat_k)
    }
  }
  cum_Ntot_k = cumsum(cum_Ntot_k)
  
  
  # Find survival path of each feature for each particle
  if(is.null(nrow(Amat_k))){
    res = list("Zmat_k" = matrix(0,nrow = Ttot, ncol = 0),
               "Amat_k" = matrix(0,nrow = 0, ncol = D),
               "mean_k" = rep(0,D),
               "Path_k" = Path_k,
               "kfinal" = kfinal,
               "B_k"    = B_k,
               "Nnew_paths_k" = Nnew_paths_k,
               "cum_Ntot_k" = cum_Ntot_k)
    return(res)
  }else{
    
    Zmat_k = matrix(0, nrow = Ttot, ncol = nrow(Amat_k))
    t = 1
    if(Nnew_paths_k[t] > 0){
      Zmat_k[t, 1:cum_Ntot_k[1] ] = 1
    }
    for(t in 2:Ttot){
      # Time t (2...Ttot)
      if(Nnew_paths_k[t] > 0)
        Zmat_k[t, (cum_Ntot_k[t-1]+1):(cum_Ntot_k[t])] = 1
      
      Nsurvived_old = length(which(Zmat_k[t-1,] == 1))
      submat = matrix( Zmat_k[,which(Zmat_k[t-1,] == 1)],
                       nrow = Ttot, ncol = Nsurvived_old)
      submat[t, Path_k[[t]]$survived ] = 1
      Zmat_k[,which(Zmat_k[t-1,] == 1)] = submat
    }
    
    res = list("Zmat_k" = Zmat_k,
               "Amat_k" = Amat_k,
               "mean_k" = Zmat_k%*%Amat_k,
               "Path_k" = Path_k,
               "kfinal" = kfinal,
               "B_k"    = B_k,
               "Nnew_paths_k" = Nnew_paths_k,
               "cum_Ntot_k" = cum_Ntot_k)
    res
    
  }
}









# Model 2 -----------------------------------------------------------------

# Set Bfix = rep(-1,Ttot);Particle_fix = NULL to run without conditioning
CondSMC = function( X,N,D,Ttot, 
                    Bfix, Particle_fix,  
                    M1,Psurv,
                    sigma2_A,sigma2_X,zeta, # these form theta
                    proposal_Nnew_1T = NULL,
                    use_VS = FALSE,
                    seed0 = 42){
  
  thinning_func = choose_propose_thinning(use_VS)
  # Define structure to store particles during the filtering
  base_list = list("active_values" = c(),"Nnew" = c(), 
                   "survived" = c(), 
                   "gr_alloc" = c(), "gr_card" = c())
  Particles = lapply(1:N, function(x){lapply(1:Ttot, function(xx){base_list})})
  # Particles[[k]][[t]]$active_values   : matrix with the active features at time t for particle k
  # Particles[[k]][[t]]$Nnew     : number of new features at time t for particle k
  # Particles[[k]][[t]]$survived : indices of the features of parent particle survived from time t-1 to time t for particle k
  
  H = length(zeta) # number of centers
  
  # Define auxiliary structures
  A <- B <- log_w_mat <- w <- W <- matrix(NA, nrow = Ttot, ncol = N)
  # The final row of A is empty
  # A[1,]: vector with the parents of the particles at time 2, i.e., X[2,k] is generated conditionally to X[1,A[1,k]]
  
  ### ---> Start filtering
  # Time 1
  # a) Sample the particles
  t  = 1
  # cat("\n --------- \n t = ",t,"\n")
  for(k in 1:N){
    # cat("\n --------- \n")
    if( (!is.null(Bfix)) && (k == Bfix[t]) ){
      # Set the k-th particle to the value we are conditioning on
      Particles[[k]][[t]] = Particle_fix[[t]]
      log_w = Reduce(`+`, lapply(1:H, function(h){
        log_dmarg_img( Particles[[k]][[t]]$gr_card[h],
                       X[t,],zeta[[h]],sigma2_X,sigma2_A ) 
      }))
    }else{
      # Free particle, draw new values
      Nnew = 0; log_w = 0; gr_alloc = c(); gr_card = rep(0,H)
      values = matrix(0,nrow = 0, ncol = D)
      
      # For each group h=1,...,H
      for(h in 1:H){
        seed_kth = ceiling(seed0 * k * t * h * runif(n=1,min = 1,max = 100))
        # Discrete adaptive proposal for Nnew_h
        if(!is.null(proposal_Nnew_1T)){
          cat("\n Dentro alla nuova proposal! \n")
          disc_adapt_prop = proposal_Nnew_1T[t,]
          disc_adapt_prop = table(disc_adapt_prop)
          Kproposal_ind = sample(1:length(disc_adapt_prop), size = 1, prob = disc_adapt_prop)
          Kprop1 = as.numeric(names(disc_adapt_prop)[Kproposal_ind])
          Kprop2 = rpois(n = 1, lambda = M1)
        }else{
          Kprop1 <- Kprop2 <- rpois(n = 1, lambda = M1)
        }
        Nnew_h = ifelse(runif(1)<0.75, Kprop1, Kprop2)
        values_h = sample_A(   Nnew_h, X[t,], zeta[[h]], sigma2_X, sigma2_A, seed_kth ) # new featues in group h
        log_w_h = log_dmarg_img(Nnew_h,X[t,], zeta[[h]], sigma2_X, sigma2_A ) # log un-norm weight in group h
        # cat("\n h = ",h,"; Nnew_h = ",Nnew_h,"; log_w_h = ",log_w_h)
        # Assemble the different groups
        Nnew = Nnew + Nnew_h
        gr_card[h] = Nnew_h
        if(Nnew_h > 0){
          gr_alloc = c(gr_alloc,rep(h,Nnew_h)) # group allocation vector for each feature 
          values = rbind(values,values_h)
        }
        log_w = log_w + log_w_h
      }

      Particles[[k]][[t]]$active_values = values
      Particles[[k]][[t]]$Nnew          = Nnew
      Particles[[k]][[t]]$gr_card       = gr_card
      Particles[[k]][[t]]$gr_alloc      = gr_alloc
      # cat("\n t = ",t,"; k = ",k,"; #new = ",Nnew,"; alloc: ",Particles[[k]][[t]]$gr_alloc,"\n")
      if(Nnew > 0){
        Particles[[k]][[t]]$survived  = seq(1,Nnew)
      }
    }
    # Compute the unnormalized weights
    log_w_mat[t,k] = log_w
    
    # °°° Finish updating particle k
  }
  # b) Compute the nomralized weights
  w[t,] = exp( log_w_mat[t,] - max(log_w_mat[t,]) )
  W[t,] = w[t,]/sum(w[t,])
  
  # Time t (t = 2,...,Ttot)
  t = 2
  for(t in 2:Ttot){
    # cat("\n --------- \n t = ",t,"\n")
    # a) Sample the parent node
    A[t-1,] = sample(1:N, size = N, W[t-1,], replace = TRUE)
    if( all(Bfix > 0) )
      A[t-1, Bfix[t] ] = Bfix[t-1] # Cond. SMC
    
    # b) Sample the particles
    for(k in 1:N){
      
      # cat("\n --------- \n")
      
      j = A[t-1,k] # parent index
      Xstar  = X[t,] 
      log_prod_bern <- log_prop_thin <- 0
      if((!is.null(Bfix)) && (k == Bfix[t])){
        # Compute Xstar
        if(  !is.null( nrow(Particles[[j]][[t-1]]$active_values))  )
          Xstar  = Xstar - colSums(Particles[[j]][[t-1]]$active_values)
        
        # Set the k-th particle to the value we are conditioning on
        Particles[[k]][[t]] = Particle_fix[[t]]
        log_w = Reduce(`+`, lapply(1:H, function(h){
          log_dmarg_img( Particles[[k]][[t]]$gr_card[h],
                         Xstar,zeta[[h]],sigma2_X,sigma2_A ) 
        }))
      }else{
        
        # Update particle k at time t
        
        # Thinning, for each group h = 1, ..., H
        ind_survived = c(); log_prop_thin = 0; log_prod_bern = 0
        gr_alloc_surv = c(); gr_card_new = rep(0,H)
        for(h in 1:H){
          # Thinning part (group h)
          num_active_old_h = Particles[[j]][[t-1]]$gr_card[h]
          if( num_active_old_h > 0  ){
            
            ind_h = which(Particles[[j]][[t-1]]$gr_alloc == h)
            ActFeat_h = matrix(Particles[[j]][[t-1]]$active_values[ind_h,],
                               ncol = D)# active features in group h at time t
            thinning = thinning_func( Xobs = X[t,], 
                                      ActFeat = ActFeat_h,
                                      sig2X = sigma2_X, Psurv = Psurv )
            
            survived_h = thinning$survived # this is a 0/1 vector  
            survived_ind_h = ind_h[which(survived_h == 1)] # this is a vector with the survived indices
            n_surv_h   = sum(survived_h)
            # cat("\n h = ",h,"; #old = ",num_active_old_h,"; #new = ",n_surv_h)
            ind_survived = c(ind_survived, survived_ind_h) # collect ind_survived across all different groups
            gr_alloc_surv = c(gr_alloc_surv, rep(h,n_surv_h) ) # collect group allocation of survived atoms
            gr_card_new[h] = gr_card_new[h] + n_surv_h # update number of features in group h
            
            log_prop_thin = log_prop_thin + 
                            thinning$log_prop_move
            log_prod_bern = log_prod_bern + 
                            n_surv_h * log(Psurv) + 
                            (num_active_old_h-n_surv_h) * log(1-Psurv)
          }
        }
        
        # Compute Xstar (the remaining part of the image)
        n_surv = length(ind_survived) # total number of survived features
        if(n_surv > 0){
          # At least one feature survived
          if( length(gr_alloc_surv) != length(ind_survived) )
            stop("Error in CondSMC. The number of survived atoms mismatches the group allocation survived vector")
          survived_values = Particles[[j]][[t-1]]$active_values[ind_survived,]
          if(!is.matrix(survived_values))
            survived_values = matrix(survived_values, nrow = n_surv, ncol = D, byrow = T)
          
          # Subtract the sum of survived features
          Xstar  = Xstar - rep(1,n_surv)%*%survived_values
        }
        
        # Innovation: ... For each group h=1,...,H
        Nnew = 0; log_w = 0; gr_alloc = c(); 
        values = matrix(0,nrow = 0, ncol = D)
        for(h in 1:H){
          seed_kth = ceiling(seed0 * k * t * h * runif(n=1,min = 1, max = 100))
          
          # Discrete adaptive proposal for Nnew_h
          if(!is.null(proposal_Nnew_1T)){
            disc_adapt_prop = proposal_Nnew_1T[t,]
            disc_adapt_prop = table(disc_adapt_prop)
            Kproposal_ind = sample(1:length(disc_adapt_prop), size = 1, prob = disc_adapt_prop)
            Kprop1 = as.numeric(names(disc_adapt_prop)[Kproposal_ind])
            Kprop2 = rpois(n = 1, lambda = M1)
          }else{
            Kprop1 <- Kprop2 <- rpois(n = 1, lambda = M1)
          }
          Nnew_h = ifelse(runif(1)<0.75, Kprop1, Kprop2)
          
          # values_h = sample_A(Nnew_h, X[t,], zeta[[h]], sigma2_X, sigma2_A ) # new featues in group h
          # log_w_h = log_dmarg_img(Nnew_h,X[t,],zeta[[h]],sigma2_X,sigma2_A ) # log un-norm weight in group h
          
          values_h = sample_A(    Nnew_h, Xstar, zeta[[h]], sigma2_X, sigma2_A, seed_kth ) # new featues in group h
          log_w_h = log_dmarg_img(Nnew_h, Xstar, zeta[[h]], sigma2_X, sigma2_A ) # log un-norm weight in group h
          
          # cat("\n h = ",h,"; Nnew_h = ",Nnew_h,"; log_w_h = ",log_w_h)
          # Assemble the different groups
          Nnew = Nnew + Nnew_h
          gr_card_new[h] = gr_card_new[h] + Nnew_h # update number of features in group h
          if(Nnew_h > 0){
            gr_alloc = c(gr_alloc,rep(h,Nnew_h)) # group allocation vector for each feature 
            values = rbind(values,values_h)
          }
          log_w = log_w + log_w_h
        }
        
        # Assemble thinning + innovation:
        # n_surv = length(ind_survived) # già calcolato
        if(n_surv > 0){
          Particles[[k]][[t]]$active_values = rbind( survived_values, values   )
          Particles[[k]][[t]]$gr_alloc  = c(gr_alloc_surv, gr_alloc)
        }else{
          # None of the features survived
          Particles[[k]][[t]]$active_values = values 
          Particles[[k]][[t]]$gr_alloc      = gr_alloc 
        }
        
        Particles[[k]][[t]]$Nnew     = Nnew
        Particles[[k]][[t]]$survived = ind_survived
        Particles[[k]][[t]]$gr_card  = gr_card_new
        
      }
      
      # cat("\n t = ",t,";k = ",k,";#new = ",Nnew,";n_surv = ",n_surv,";alloc: ",Particles[[k]][[t]]$gr_alloc,"\n")
      # Compute the un-normalized weights
      log_w_mat[t,k] = log_w + log_prod_bern - log_prop_thin
    }
    
    # c) Compute the normalized weights
    w[t,] = exp( log_w_mat[t,] - max(log_w_mat[t,]) )
    W[t,] = w[t,]/sum(w[t,])
  }
  ### End of filtering <--- 
  
  # Select the chosen particle
  kfinal = sample(1:N, size = 1, replace = TRUE, prob = W[Ttot,])
  
  # Find ancestral lineage and save final paths
  B_k = rep(0,Ttot)
  B_k[Ttot] = kfinal
  for(t in (Ttot-1):1){
    B_k[t] = A[ t, B_k[t+1] ] # B_t^k = A_t^{B_{t+1}^k}
  }
  
  # Save selected path
  Path_k = vector("list", Ttot)
  for(t in 1:Ttot){
    Path_k[[t]] = Particles[[ B_k[t] ]][[t]]
  }
  
  # Find all features appeared at least once
  Amat_k <- cum_Ntot_k <- Nnew_paths_k <- gr_alloc_k <- c()
  for(t in Ttot:1){
    
    Mat = Path_k[[t]]$active_values
    nuove = Path_k[[t]]$Nnew
    Nnew_paths_k = c(nuove, Nnew_paths_k)
    cum_Ntot_k   = c(nuove, cum_Ntot_k)
    if(nuove > 0){
      if(nrow(Mat)-nuove+1 <= 0)
        stop("Unexprected behaviour")
      Amat_k = rbind( Mat[(nrow(Mat)-nuove+1):nrow(Mat),] , Amat_k)
      gr_alloc_k = c(tail(Path_k[[t]]$gr_alloc, nuove), gr_alloc_k)
    }
  }
  cum_Ntot_k = cumsum(cum_Ntot_k)
  # length(gr_alloc_k) = nrow(Amat_k); 
  # gr_alloc_k[j] = h iff Amat_k[j,] belong to h-th group 
  
  
  
  # Find survival path of each feature for each particle
  # and label the active particles
  if(is.null(nrow(Amat_k))){
    res = list("Zmat_k" = matrix(0,nrow = Ttot, ncol = 0),
               "Amat_k" = matrix(0,nrow = 0, ncol = D),
               "mean_k" = rep(0,D),
               "Path_k" = Path_k,
               "kfinal" = kfinal,
               "B_k"    = B_k,
               "Nnew_paths_k" = Nnew_paths_k,
               "cum_Ntot_k" = cum_Ntot_k,
               "gr_alloc_k" = gr_alloc_k)
    return(res)
  }else{
    
    t = 1
    Ktot = 0
    Zmat_k_old = matrix(0, nrow = Ttot, ncol = nrow(Amat_k))
    Zmat_k = matrix(0,nrow = Ttot, ncol = nrow(Amat_k))
    Path_k[[t]][["label_actives"]] <- c()  
    
    n_active = length(Path_k[[t]]$gr_alloc)
    if(Path_k[[t]]$Nnew > 0){
      Path_k[[t]]$label_actives =  c(Path_k[[t]]$label_actives,1:n_active)
      Ktot = Ktot + Path_k[[t]]$Nnew 
      Zmat_k_old[t, 1:cum_Ntot_k[1] ] = 1
      Zmat_k[t,Path_k[[t]]$label_actives] = 1
    }
    for(t in 2:Ttot){
      # Time t (2...Ttot)
      #a) labels
      n_active = length(Path_k[[t]]$gr_alloc)
      Path_k[[t]][["label_actives"]] <- c()
      ## Survived
      n_surv = length(Path_k[[t]]$survived)
      if(n_surv > 0){
        Path_k[[t]]$label_actives = c( Path_k[[t]]$label_actives,
                                       Path_k[[t-1]]$label_actives[Path_k[[t]]$survived] )
      }
      ## New
      if(Path_k[[t]]$Nnew > 0){
        Path_k[[t]]$label_actives = c( Path_k[[t]]$label_actives,
                                       (Ktot+1):(Ktot+Path_k[[t]]$Nnew) )
        Ktot = Ktot + Path_k[[t]]$Nnew 
      }
      #b) Z matrix
      Zmat_k[t,Path_k[[t]]$label_actives] = 1
      
      if(Nnew_paths_k[t] > 0)
        Zmat_k_old[t, (cum_Ntot_k[t-1]+1):(cum_Ntot_k[t])] = 1
      
      Nsurvived_old = length(which(Zmat_k_old[t-1,] == 1))
      submat = matrix( Zmat_k_old[,which(Zmat_k_old[t-1,] == 1)],
                       nrow = Ttot, ncol = Nsurvived_old)
      submat[t, Path_k[[t]]$survived ] = 1
      Zmat_k_old[,which(Zmat_k_old[t-1,] == 1)] = submat
    }
    
    if(Ktot != nrow(Amat_k))
      stop("Something wrong with the total number of groups")
    
    res = list("Zmat_k" = Zmat_k,
               "Zmat_k_old" = Zmat_k_old,
               "Amat_k" = Amat_k,
               "mean_k" = Zmat_k%*%Amat_k,
               "Path_k" = Path_k,
               "kfinal" = kfinal,
               "B_k"    = B_k,
               "Nnew_paths_k" = Nnew_paths_k,
               "cum_Ntot_k" = cum_Ntot_k,
               "gr_alloc_k" = gr_alloc_k)
    res
    
  }
  
}




r_fullcond_zeta = function( data, 
                            gr_alloc_old, # size: nrow(data); 
                            centers_old,  # size; length(table(gr_alloc_data))
                            sig2_A, sig2_ker, Lambda,
                            updateGrAlloc = FALSE, 
                            updateCenters = FALSE, 
                            upNcenters = FALSE )
{
  
  
  D = ncol(data)
  n = nrow(data)
  if( n != length(gr_alloc_old) )
    stop("The number of data does not match with the length of the group labels")
  if(sig2_A < 0 || sig2_ker < 0)
    stop("Negative variance")
  
  # re-order input labels
  Hold = length(table(gr_alloc_old)) # num. allocated groups
  Hstar_old = length(centers_old) - Hold # num. non allocared groups
  if(Hstar_old < 0)
    stop("Hstar_old is negative, this is impossible")
  
  n_h_old = tabulate(gr_alloc_old, nbins = (Hold+Hstar_old))
  counter = 1; h = 1
  while(counter <= length(n_h_old)){
    if(n_h_old[h] == 0){
      # h-th group is non-allocated:
      # move h-th cluster in final position and decrease the cluster label
      # of all subsequent clusters
      gr_alloc_old[gr_alloc_old > h] <- gr_alloc_old[gr_alloc_old > h] - 1
      n_h_old = c(n_h_old[-h],n_h_old[h])
      centers_old <- centers_old[c(setdiff(seq_along(centers_old), h), h)]
      h = h - 1
    }
    counter = counter + 1
    h = h + 1
  }
  
  # Update the number of centers
  if(upNcenters){
    z = map_to_Z(Hstar_old); # map the current value of Hstar in Z
    support_proposal = c((-5):(-1),1:5)
    Uproposal = sample( support_proposal, size = 1 )  #draw uniform proposal from -proposal to proposal expect 0
    temp = ifelse(runif(1) < 0.5, z,-z) # proposal step
    zproposal = Uproposal + temp # complete proposal step
    Hproposal = map_to_N(zproposal)
    
    # Compute log probability of accepting the move
    log_alpha = ( Hproposal - Hstar_old )*log(Lambda) + 
                lfactorial(Hstar_old) - lfactorial(Hproposal) + 
                (n-1) * ( log(Hstar_old+Hold) - log(Hproposal+Hold) )
    # Acceptance-rejection move
    Hstar_temp = ifelse(runif(1) < exp(log_alpha), Hproposal, Hstar_old )
    # cat("\n current: ",Hstar_old,"; proposed: ",Hproposal,"; prob = ",log_alpha,"\n")
    # Check
    if(Hstar_temp < 0)
      stop("Hstar_temp can not be negative")
  }else{
    Hstar_temp = Hstar_old
  }
  
  # Update the center values
  if(updateCenters){
    centers_new = lapply(1:(Hold+Hstar_temp), function(h){rep(NA,D)})
    # Allocated components
    for(h in 1:Hold){
      ind_h = which(gr_alloc_old == h) # ind. of elements in group h
      n_h = length(ind_h) # group cardinality
      if(is.null(n_h) || n_h == 0)
        stop("Empty allocated component")
      
      data_h = matrix(data[ind_h,],ncol = D)
      stat_suff_h = apply(data_h, 2, sum)
      scale_var_post = (sig2_A*sig2_ker)/(sig2_A+ n_h*sig2_ker)
      scale_mean_post = (sig2_ker)/(sig2_A+ n_h*sig2_ker)
      mean_post = scale_mean_post * stat_suff_h
      for(j in 1:D){
        centers_new[[h]][j] = rnorm(n = 1, mean = mean_post[j], sd = sqrt(scale_var_post))
      }
    }
    # Non Allocated components
    if(Hstar_temp > 0){
      for(h in (Hold+1):(Hold+Hstar_temp)){
        centers_new[[h]] = rnorm(n = D, mean = 0, sd = sig2_ker )
      }
    }
  }else{
    centers_new = centers_old
  }
  
  #check
  if(length(centers_new) != (Hold+Hstar_temp) )
    stop("The number of new centers does not match with the length of the vector of new centers")
  
  # Update the group allocations
  if(updateGrAlloc){
    gr_alloc_new = t(apply(data, 1, function(x){
                      logw = sapply(centers_new, function(z){
                                log_dmvnorm(x,z,diag(sig2_A,D))
                             })
                      w = exp( logw - max(logw) ) # unnormalized
                      sample(1:(Hold+Hstar_temp), size = 1, prob = w)
                    }))
    
   # Update H and re-order the groups
    n_h_new = tabulate(gr_alloc_new, nbins = (Hold+Hstar_temp))
    Hnew = length(which(n_h_new > 0)) # find num. allocated components
    Hstar_new = length(which(n_h_new == 0)) # find num. non-allocated components
    
    counter = 1; h = 1
    while(counter <= length(n_h_new)){
      if(n_h_new[h] == 0){
        # h-th group is non-allocated:
        # move h-th cluster in final position and decrease the cluster label
        # of all subsequent clusters
        gr_alloc_new[gr_alloc_new > h] <- gr_alloc_new[gr_alloc_new > h] - 1
        n_h_new = c(n_h_new[-h],n_h_new[h])
        centers_new <- centers_new[c(setdiff(seq_along(centers_new), h), h)]
        h = h - 1
      }
      counter = counter + 1
      h = h + 1
    }
  }else{
    Hnew = Hold
    Hstar_new = Hstar_old
    gr_alloc_new = gr_alloc_old
  }
  
  return(list( "H" = Hnew, "Hstar" = Hstar_new, 
               "centers" = centers_new, 
               "gr_alloc" = gr_alloc_new) )

}

