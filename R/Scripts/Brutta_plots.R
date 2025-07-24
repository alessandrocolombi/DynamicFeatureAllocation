
# Show mean_vec_1T
par(mfrow = c(3,4), mar = c(2,2,2,1), bty = "l")
for(t in 1:Ttot){
  res_t_vec = mean_vec_1T[t,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("mean T=",t) )
}


# Plot Amat & Zmat
for(i in 1:nrow(Amat)){
  # plot figure
  par(mfrow = c(1,2), mar = c(2,2,2,1), bty = "l")
  res_t_vec = Amat[i,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
  
  # plot survival path
  plot(x = 1:Ttot, y = Zmat[,i], pch = 16, xlab = "", ylab = "")
  segments(x0 = 1:Ttot, x1 = 1:Ttot, y0 = rep(0,Ttot), y1 = Zmat[,i])
}


# Plot centers old
par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
for(i in 1:length(zetas_mcmc[[g-1]])){
  # plot figure
  res_t_vec = zetas_mcmc[[g-1]][[i]]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("center[",i,"], g = ",g) )
}

# Plot centers new
par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
for(i in 1:length(zetas_mcmc[[g]])){
  # plot figure
  res_t_vec = zetas_mcmc[[g]][[i]]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("center[",i,"], g = ",g) )
}


# Plot Amat 
par(mfrow = c(2,2), mar = c(2,2,2,1), bty = "l")
for(i in 1:nrow(Amat)){
  # plot figure
  res_t_vec = Amat[i,]
  res_t = matrix(data = res_t_vec, nrow = di, ncol = di, byrow = F)
  plot_img( res_t,center_value = NULL, col.lower = "grey95",col.upper = "grey10",
            horizontal = F, main = paste0("A[",i,"]") )
}
