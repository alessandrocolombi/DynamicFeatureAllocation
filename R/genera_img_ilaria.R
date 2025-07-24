#***** LINEAR GAUSSIAN LATENT FEATURE MODEL: simulation of data *****#
sim_images = function(epsX, Ti){
  set.seed(1989)
  # rm(list=ls())
  di = 8
  D = di*di
  # 0 nero, 1 bianco
  
  # NUOVI FEATURES:
  A1 = rep(0, D )
  idx = c(2, 9:11, 18)
  A1[idx] = rep(1, length(idx))
  postscript("A1.eps", hor=F)
  image((matrix(data = A1[1:D], nrow = di, ncol = di, byrow = F)), col = c("black", "white"), 
        xaxt="n", yaxt="n") 
  dev.off()
  
  A2 = rep(0, D)
  idx = c(5:8, 13,16, 21,24, 29:32)
  A2[idx] = rep(1, length(idx))
  postscript("A2.eps", hor=F)
  image((matrix(data = A2[1:D], nrow = di, ncol = di, byrow = F)), col = c("black", "white"), 
        xaxt="n", yaxt="n") 
  dev.off()
  
  
  A3 = rep(0, D)
  idx = c(33, 41:42, 49:51, 57:60)
  A3[idx] = rep(1, length(idx))
  postscript("A3.eps", hor=F)
  image((matrix(data = A3[1:D], nrow = di, ncol = di, byrow = F)), col = c("black", "white"), 
        xaxt="n", yaxt="n") 
  dev.off()
  
  A4 = rep(0, D)
  idx = c(38:40, 47, 55, 63)
  A4[idx] = rep(1, length(idx))
  postscript("A4.eps", hor=F)
  image((matrix(data = A4[1:D], nrow = di, ncol = di, byrow = F)), col = c("black", "white"), 
        xaxt="n", yaxt="n") 
  dev.off()
  
  A5 = A2[D:1]
  postscript("A5.eps", hor=F)
  image((matrix(data = A5[1:D], nrow = di, ncol = di, byrow = F)), col = c("black", "white"), 
        xaxt="n", yaxt="n") 
  dev.off()
  
  A6 = c(rep(c(rep(0,3), 1,1,rep(0,3)), di))
  postscript("A6.eps", hor=F)
  image((matrix(data = A6[1:D], nrow = di, ncol = di, byrow = F)), col = c("black", "white"), 
        xaxt="n", yaxt="n") 
  dev.off()
  
  if(Ti > 50)
  {
    K = 6
    A = rbind(A1, A2, A3, A4, A5, A6)
  } 
  else
  {
    K = 4
    A = rbind(A1, A2, A3, A4)    
  }
  
  dim(A)
  # Genero Z(1), ..., Z(T)
  
  ### SIMULATE Z
  Z = matrix(nrow = Ti, ncol = K)
  for(t in 1:Ti){
    if(t == 1){
      Zt = c(1, rep(0, K-1))
    }else{
      Zt = numeric(K)
      for(i in 1:K){
        if(Z[t-1,i]==1){
          Zt[i] = rbinom(1, 1, 0.9)
        }else{
          Zt[i] = 0
        }
      }
      # mettere la dipendenza temporale!
      if(t == (floor(Ti/K)-1)) Zt[2] = 1
      if(t == 2*(floor(Ti/K)-1)) Zt[3] = 1
      if(t == 3*(floor(Ti/K)-1)) Zt[4] = 1
      if(Ti>50){
        if(t == 4*(floor(Ti/K)-1)) Zt[5] = 1
        if(t == 5*(floor(Ti/K)-1)) Zt[6] = 1        
      }
    }
    Z[t, ] = Zt
  }
  
  
  # Simulate images:
  require("mvtnorm")
  X = matrix(nrow = Ti, ncol = D)
  # x11()
  # par(mfrow=c(4,5))
  for(t in 1:Ti){
    Zt = Z[t, ]
    aux = Zt%*%A # è la media
    X[t,] = rmvnorm(1, mean = aux, sigma = epsX*diag(D))
    postscript(paste("image",t,".eps", sep = ""), hor=F)
    #  x11()
    image((matrix(data = X[t,][1:D], nrow = di, ncol = di, byrow = F)), 
          col = gray.colors(20), xaxt="n", yaxt="n")
    dev.off() 
  }
  return(list(X = X, Z = Z, A = A) )
}


sim_images_ale = function(img_base, Ti, sig2X, sig2A, Psurv){
  # img_base is KxD
  A = img_base + rnorm( n = nrow(img_base)*ncol(img_base), sd = sqrt(sig2A) )
  K = ifelse(Ttot > 50, 6, 4)
  A = A[1:K,]
  
  # Genero Z(1), ..., Z(T)
  starts = seq(1,Ttot,by = floor(Ttot/K))
  Z = matrix(0,nrow = Ti, ncol = K)
  for(j in 1:K){
    life_span = rgeom(n=1,prob = 1-Psurv)
    end = min(starts[j]+life_span, Ttot)
    Z[starts[j]:end,j] = 1
  }
  
  # Simulate images:
  require("mvtnorm")
  X = matrix(nrow = Ti, ncol = D)
  for(t in 1:Ti){
    Zt = Z[t, ]
    aux = Zt%*%A # è la media
    X[t,] = rmvnorm(1, mean = aux, sigma = sig2X*diag(D))
  }
  return(list(X = X, Z = Z, A = A) )
}



sim_images_ale_mod2 = function(zetas, Ti, sig2X, sig2A, Psurv, M1){
  
  # img_base is KxD
  A = zetas + rnorm( n = nrow(zetas)*ncol(zetas), sd = sqrt(sig2A) )
  K = ifelse(Ttot > 50, 6, 4)
  A = A[1:K,]
  
  # Genero Z(1), ..., Z(T)
  starts = seq(1,Ttot,by = floor(Ttot/K))
  Z = matrix(0,nrow = Ti, ncol = K)
  for(j in 1:K){
    Z[starts[j],j] = 1
    for(t in (starts[j]+1):Ttot){
      if(t > Ttot)
        break
      if(Z[t-1,j] == 1)
        Z[t,j] = sample(c(0,1), size = 1, prob = c(1-Psurv,Psurv))
      if(Z[t-1,j] == 0){
        Pinn = ppois(0,M1,lower.tail = FALSE)
        Z[t,j] = sample(c(0,1), size = 1, prob = c(1-Pinn,Pinn))
      }
    }
  }
  
  # Simulate images:
  require("mvtnorm")
  X = matrix(nrow = Ti, ncol = D)
  for(t in 1:Ti){
    Zt = Z[t, ]
    aux = Zt%*%A # è la media
    X[t,] = rmvnorm(1, mean = aux, sigma = sig2X*diag(D))
  }
  
  return(list(X = X, Z = Z, A = A) )
}