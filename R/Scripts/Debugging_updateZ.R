# Introduction ------------------------------------------------------------
rm(list = ls())
setwd("C:/Users/colom/DynamicFeatureAllocation/R/Scripts")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../Rfunctions.R")
source("../genera_img_ilaria.R")

# Custom functions --------------------------------------------------------


# Immagini Base -----------------------------------------------------------

di   = 8
D    = di*di # Size of the problem

par(mfrow = c(3,2), mar = c(2,2,2,1), bty = "l")
A1 = rep(0, D)
idx = c(2, 9:11, 18)
A1[idx] = rep(1, length(idx))
image((matrix(data = A1[1:D], nrow = di, ncol = di, byrow = F)), 
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 

A2 = rep(0, D)
idx = c(5:8, 13,16, 21,24, 29:32)
A2[idx] = rep(1, length(idx))
image((matrix(data = A2[1:D], nrow = di, ncol = di, byrow = F)),
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 


A3 = rep(0, D)
idx = c(33, 41:42, 49:51, 57:60)
A3[idx] = rep(1, length(idx))
image((matrix(data = A3[1:D], nrow = di, ncol = di, byrow = F)), 
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 

A4 = rep(0, D)
idx = c(38:40, 47, 55, 63)
A4[idx] = rep(1, length(idx))
image((matrix(data = A4[1:D], nrow = di, ncol = di, byrow = F)), 
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 

A5 = A2[D:1]
image((matrix(data = A5[1:D], nrow = di, ncol = di, byrow = F)), 
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 

A6 = c(rep(c(rep(0,3), 1,1,rep(0,3)), di))
image((matrix(data = A6[1:D], nrow = di, ncol = di, byrow = F)), 
      col = c("black", "white"), 
      xaxt="n", yaxt="n") 



img_base = matrix(NA,nrow = 6, ncol = D)
img_base[1,] = A1; img_base[2,] = A2
img_base[3,] = A3; img_base[4,] = A4
img_base[5,] = A5; img_base[6,] = A6

# Sim. data from zetas  -------------------------------------------------

seed = 42
set.seed(seed)

Ttot = 24    # Time length
di   = 8
D    = di*di # Size of the problem

H = 4
zetas = rbind(A1,A2,A3,A4)

sig2_A = 0.0001; sigma2_X = 0.01
sig2_ker = 0.001


par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
for(h in 1:H){
  res_t_vec = zetas[h,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("Est. T = ",t))
}


Nh = 10
data = matrix(0,nrow = 0, ncol = D)
for(h in 1:H){
  x = zetas[h,]
  data = rbind(data, mvtnorm::rmvnorm(n = Nh, mean = x, sigma = sig2_A*diag(D) ))
}
dim(data)

par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
for(i in 1:nrow(data)){
  # plot figure
  res_t_vec = data[i,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
}


gr_alloc_true = c()
for(h in 1:H)
  gr_alloc_true = c(gr_alloc_true, rep(h,Nh))

# Fit - 1 iter --------------------------------------------------------------------
Amat = data
gr_alloc_old = gr_alloc_true #sample(1:H, size = nrow(Amat), replace = TRUE)
centers_old  = lapply(1:H, function(h){rep(0,D)})
# centers_old  = lapply(1:H, function(h){zetas[h,]})
# Update zetas
temp = r_fullcond_zeta( data = Amat, 
                        gr_alloc_old = gr_alloc_old,
                        centers_old = centers_old,
                        sig2_A = sig2_A,
                        sig2_ker = sig2_ker,
                        updateGrAlloc = TRUE, 
                        updateCenters = TRUE, 
                        upNcenters = FALSE )
centers_new = temp$centers
H_new = temp$H; Hstar_new = temp$Hstar
gr_alloc_new = c(temp$gr_alloc)


par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
for(i in 1:length(centers_new)){
  # plot figure
  res_t_vec = centers_new[[i]]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
}

# Fit - multiple iterations --------------------------------------------------------------------
# set.seed(1234) # giusto
set.seed(112233) # sbagliato

Amat = data
gr_alloc_old =  sample(1:H, size = nrow(Amat), replace = TRUE)
centers_old  = lapply(1:H, function(h){rep(0,D)})

# gr_alloc_old =  gr_alloc_true
# centers_old  = lapply(1:H, function(h){zetas[h,]})

G = 10

pdf("img/debug_zetas.pdf")
pb = txtProgressBar(min = 1, max = G, initial = 1) 
for(g in 1:G){
  setTxtProgressBar(pb,g)
  # Update zetas
  temp = r_fullcond_zeta( data = Amat, 
                          gr_alloc_old = gr_alloc_old,
                          centers_old = centers_old,
                          sig2_A = sig2_A,
                          sig2_ker = sig2_ker,
                          updateGrAlloc = TRUE, 
                          updateCenters = TRUE, 
                          upNcenters = FALSE )
  centers_old = temp$centers
  gr_alloc_old = c(temp$gr_alloc)
  
  par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
  for(i in 1:length(centers_old)){
    # plot figure
    res_t_vec = centers_old[[i]]
    res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
    plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
              horizontal = F, main = paste0("center[",i,"], g = ",g) )
  }
}
close(pb)
dev.off()


centers_new = temp$centers
H_new = temp$H; Hstar_new = temp$Hstar
gr_alloc_new = c(temp$gr_alloc)


par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
for(i in 1:length(centers_new)){
  # plot figure
  res_t_vec = centers_new[[i]]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
}
# Debugging 1 ---------------------------------------------------------------

set.seed(112233) # sbagliato
Amat = data
gr_alloc_old =  sample(1:H, size = nrow(Amat), replace = TRUE)
centers_old  = lapply(1:H, function(h){rep(0,D)})


data = Amat
gr_alloc_old = gr_alloc_old
centers_old = centers_old
sig2_A = sig2_A
sig2_ker = sig2_ker
updateGrAlloc = TRUE 
updateCenters = TRUE 
upNcenters = FALSE 


  
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
    Hstar_temp = rpois(n=1,lambda = 1) # A caso, da modificare
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
      
      par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
      res_t_vec = mean_post 
      res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
      plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
                horizontal = F, main = paste0("mean_post ",h," ") )
      
      
      
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
    
    # Allocate items to new cluster
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
    
    par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
    for(i in 1:length(centers_new)){
      # plot figure
      res_t_vec = centers_new[[i]]
      res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
      plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
                horizontal = F, main = paste0("center[",i,"], g = ",g) )
    }
    
    
    counter = 1; h = 1
    while(counter <= length(n_h_new)){
      if(n_h_new[h] == 0){
        # h-th group is non-allocated:
        # move h-th cluster in final position and decrease the cluster label
        # of all subsequent clusters
        gr_alloc_new[gr_alloc_new > h] <- gr_alloc_new[gr_alloc_new > h] - 1
        n_h_new = c(n_h_new[-h],n_h_new[h])
        centers_new <- centers_new[ c(setdiff(seq_along(centers_new), h), h) ]
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
  


# repeat
centers_old = centers_new
gr_alloc_old = c(gr_alloc_new)





par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
# plot figure
res_t_vec = mean_post 
res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("mean_post ",h," ") )



par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
for(i in 1:nrow(data_h)){
  # plot figure
  res_t_vec = data_h[i,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
}


par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
for(i in 1:length(centers_old)){
  # plot figure
  res_t_vec = centers_old[[i]]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("center[",i,"], g = ",g) )
}
# Debugging 2 ---------------------------------------------------------------

data = Amat
gr_alloc_old = c(1,2,3,4)#Particles_PG[[g-1]]$gr_alloc_k
centers_old = zetas_mcmc[[g-1]]
sig2_A = sigma2_XA_mcmc[g,2]
sig2_ker = sigma2_ker
updateGrAlloc = TRUE 
updateCenters = TRUE 
upNcenters = FALSE 



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
  Hstar_temp = rpois(n=1,lambda = 1) # A caso, da modificare
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
    
    par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
    res_t_vec = mean_post 
    res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
    plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
              horizontal = F, main = paste0("mean_post ",h," ") )
    
    
    
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
  
  # Allocate items to new cluster
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
  
  # par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
  # for(i in 1:length(centers_new)){
  #   # plot figure
  #   res_t_vec = centers_new[[i]]
  #   res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  #   plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
  #             horizontal = F, main = paste0("center[",i,"], g = ",g) )
  # }
  # 
  
  counter = 1; h = 1
  while(counter <= length(n_h_new)){
    if(n_h_new[h] == 0){
      # h-th group is non-allocated:
      # move h-th cluster in final position and decrease the cluster label
      # of all subsequent clusters
      gr_alloc_new[gr_alloc_new > h] <- gr_alloc_new[gr_alloc_new > h] - 1
      n_h_new = c(n_h_new[-h],n_h_new[h])
      centers_new <- centers_new[ c(setdiff(seq_along(centers_new), h), h) ]
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

# repeat
centers_old = centers_new
gr_alloc_old = c(gr_alloc_new)





par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
# plot figure
res_t_vec = mean_post 
res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
          horizontal = F, main = paste0("mean_post ",h," ") )



par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
for(i in 1:nrow(data_h)){
  # plot figure
  res_t_vec = data_h[i,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
}


par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
for(i in 1:length(centers_new)){
  # plot figure
  res_t_vec = centers_new[[i]]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("center[",i,"], g = ",g) )
}


# Debugging 3 ---------------------------------------------------------------

data = Afinal
gr_alloc_old =  sample(1:H, size = nrow(Afinal), replace = TRUE)
centers_old  = lapply(1:H, function(h){rep(0,D)})

sig2_A = sigma2_XA_mcmc[g,2]
sig2_ker = sigma2_ker
updateGrAlloc = TRUE 
updateCenters = TRUE 
upNcenters = FALSE 

G = 1000

pdf("img/debug_zetas_Afinal.pdf")
pb = txtProgressBar(min = 1, max = G, initial = 1) 
for(g in 1:G){
  setTxtProgressBar(pb,g)
  # Update zetas
  temp = r_fullcond_zeta( data = Amat, 
                          gr_alloc_old = gr_alloc_old,
                          centers_old = centers_old,
                          sig2_A = sig2_A,
                          sig2_ker = sig2_ker,
                          updateGrAlloc = TRUE, 
                          updateCenters = TRUE, 
                          upNcenters = FALSE )
  centers_old = temp$centers
  gr_alloc_old = c(temp$gr_alloc)
  
  par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
  for(i in 1:length(centers_old)){
    # plot figure
    res_t_vec = centers_old[[i]]
    res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
    plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
              horizontal = F, main = paste0("center[",i,"], g = ",g) )
  }
}
close(pb)
dev.off()


centers_new = temp$centers
H_new = temp$H; Hstar_new = temp$Hstar
gr_alloc_new = c(temp$gr_alloc)


par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
for(i in 1:length(centers_new)){
  # plot figure
  res_t_vec = centers_new[[i]]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
}