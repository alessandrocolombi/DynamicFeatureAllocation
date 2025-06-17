seed = 42
set.seed(seed)

N    = 5 # Number of particles
Ttot = 3 # Time length
D    = 2 # Size of the problem

base_list = list("active_values" = c(),"Nnew" = c(), "survived" = c())

Particles = lapply(1:N, function(x){lapply(1:Ttot, function(xx){base_list})})
# Particles[[k]][[t]]$active_values   : matrix with the active features at time t for particle k
# Particles[[k]][[t]]$Nnew     : number of new features at time t for particle k
# Particles[[k]][[t]]$survived : indeces of the features of parent particle survived from time t-1 to time t for particle k



# Define basic structures
A <- B <- log_w <- w <- W <- matrix(NA, nrow = Ttot, ncol = N)

# The final row of A is empty
# A[1,]: vector with the parents of the particles at time 2, i.e., X[2,k] is generated conditionally to X[1,A[1,k]]


# I want to replicate the path in the Figure of Andrieu et at. (2014)
A[1,] = c(3,1,4,3,4)
A[2,] = c(2,4,2,1,3)


# Time 1

# a) Sample the particles
t = 1
for(k in 1:N){
  Nnew   = rpois(n = 1, lambda = 10)
  values = mvtnorm::rmvnorm(n = Nnew, mean = rep(0,D)) # (Nnew x D) matrix
  Particles[[k]][[t]]$active_values   = values
  Particles[[k]][[t]]$Nnew     = Nnew
  Particles[[k]][[t]]$survived = seq(1,Nnew)
}

# b) Compute the un-nomralized weights
log_w[1,] = 0
w[1,] = exp( log_w[1,] - max(log_w[1,]) )
W[1,] = w[1,]/sum(w[1,])

# Time t (t = 1,...,Ttot)
for(t in 2:Ttot){
  # a) Sample the parent node
  # A[n-1,] = sample(1:N, size = N, W[n-1,], replace = TRUE)
  
  # b) Sample the particles
  for(k in 1:N){
    j = A[t-1,k] # parent index
    
    # Thinning part:
    num_active = nrow(Particles[[j]][[t-1]]$active_values) 
    survived = sample(c(0,1), size = num_active, replace = TRUE, prob = c(0.3,0.7))
    
    # Innovation:
    Nnew   = rpois(n = 1, lambda = 3)
    values = mvtnorm::rmvnorm(n = Nnew, mean = rep(0,D)) # (Nnew x D) matrix
    
    if(sum(survived > 0)){
      survived_values = Particles[[j]][[t-1]]$active_values[which(survived > 0),]
      Particles[[k]][[t]]$active_values = rbind( survived_values, values   )
    }else{
      Particles[[k]][[t]]$active_values =  values 
    }
    Particles[[k]][[t]]$Nnew     = Nnew
    Particles[[k]][[t]]$survived = which(survived > 0)
  }
  
  # c) Compute the un-normalized weights
  log_w[t,] = 0
  w[t,] = exp( log_w[t,] - max(log_w[t,]) )
  W[t,] = w[t,]/sum(w[t,])
  
}

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
for(k in 1:N){
  for(t in Ttot:1){
    Mat = Particles[[ B[t,k] ]][[t]]$active_values
    nuove = Particles[[ B[t,k] ]][[t]]$Nnew
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
  Survived_features[[k]][t, Particles[[ B[t,k] ]][[t]]$survived ] = 1
  for(t in 2:Ttot){
    # Time t (2...Ttot)
    Survived_features[[k]][t, (cum_Ntot[[k]][t-1]+1):(cum_Ntot[[k]][t])] = 1
    Survived_features[[k]]
    submat = Survived_features[[k]][,which(Survived_features[[k]][t-1,] == 1)]
    submat
    submat[t, Particles[[ B[t,k] ]][[t]]$survived ] = 1
    Survived_features[[k]][,which(Survived_features[[k]][t-1,] == 1)] = submat
  }
}





k = 4
B[,k]

# Time 1:
t = 1
Particles[[ B[t,k] ]][[t]]

# Time 2:
t = 2
Particles[[ B[t,k] ]][[t]]

# Time 3:
t = 3
Particles[[ B[t,k] ]][[t]]

Features_particle[[k]]
cum_Ntot[[k]]



Survived_features[[k]]



