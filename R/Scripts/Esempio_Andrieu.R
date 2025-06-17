
# Define custom functions -------------------------------------------------

mean_g = function(xn){
  xn^2 / 20
}
mean_f = function(xold,n){
  xold/2 + 25*(xold)/(1+xold^2) + 8*cos(1.2*n)
}
eval_log_g = function(y,x,theta){
  dnorm(x = y, mean = mean_g(x), sd = sqrt(theta[1]) ,log = T)
}
eval_log_f = function(xnew,xold,n,theta){
  dnorm(x = new, mean = mean_f(xold,n), sd = sqrt(theta[2]) ,log = T)
}
d_invgamma = function(val,a,b){
  -2*log(val) + dgamma(x = 1/val, shape = a, rate = b)
}



# Data generation ---------------------------------------------------------


set.seed(12245)
theta_true = c(10,10)
T_tot = 101

y_true_1T <- x_true_1T <- rep(0,T_tot)

y_true_1T[1] = NA
x_true_1T[1] = rnorm(n = 1, mean = 0, sd = sqrt(5) )

for(n in 2:T_tot){
  y_true_1T[n] = rnorm(n = 1, mean = mean_g(x_true_1T[n-1]),   sd = sqrt(theta_true[1]))
  x_true_1T[n] = rnorm(n = 1, mean = mean_f(x_true_1T[n-1],n), sd = sqrt(theta_true[2]))
}


y_true_1T = y_true_1T[-1]
x_true_1T = x_true_1T[1:(T_tot-1)]
T_tot = length(y_true_1T)


par(mar = c(4,4,2,2), bty = "l")
plot( x = 1:T_tot, y = y_true_1T, type = "b", pch = 16, lty = 1, lwd = 2, xlab = "Time", ylab = "Y_n")

par(mar = c(4,4,2,2), bty = "l")
plot( x = 1:T_tot, y = x_true_1T, type = "b", pch = 16, lty = 1, lwd = 2, xlab = "Time", ylab = "X_n")


par(mar = c(4,4,2,2), bty = "l")
plot( x = x_true_1T, y = y_true_1T, type = "p", pch = 16, lty = 1, lwd = 2, xlab = "X_n", ylab = "Y_n")



# 1) SMC with known parameters --------------------------------------------

data = y_true_1T
N = 2000
X <- A <- B <- log_w <- w <- W <- matrix(NA, nrow = T_tot, ncol = N)

# TxN matrix: 
## X[t,] -> values of the N particles at time t
## X[,i] -> values across time

# Time 1

# a) Sample the particles
X[1,] = rnorm( n = N, mean = 0, sd = sqrt(5) )

# b) Compute the un-nomralized weights
log_w[1,] = eval_log_g(data[1], X[1,], theta_true[1] )

w[1,] = exp( log_w[1,] - max(log_w[1,]) )
W[1,] = w[1,]/sum(w[1,])

# Time n (n = 2...T)

for(n in 2:T_tot){
  # a) Sample the parent node
  A[n-1,] = sample(1:N, size = N, W[n-1,], replace = TRUE)
  
  # b) Sample the particles
  X_parents = X[n-1,A[n-1,]]
  X[n,] = rnorm( n = N, mean = mean_f(X_parents,n = n), sd = sqrt(theta_true[2]) )
  
  # c) Compute the un-normalized weights
  log_w[n,] = eval_log_g(data[n], X[n,], theta_true[1] )
  
  w[n,] = exp( log_w[n,] - max(log_w[n,]) )
  W[n,] = w[n,]/sum(w[n,])
  
}

# Find ancestral lineage and save final paths
X_1T_final = matrix(0,nrow = T_tot, ncol = N)
B[T_tot,] = 1:N
X_1T_final[T_tot,] = X[T_tot, B[T_tot,] ]


for(k in 1:N){
  for(n in (T_tot-1):1){
    B[n,k] = A[ n, B[n+1,k] ] # B_n^k = A_n^{B_{n+1}^k}
  }
  
  temp = matrix(0,nrow = (T_tot-1):1, ncol = 2)
  temp[,1] = (T_tot-1):1
  temp[,2] = B[(T_tot-1):1,k]
  X_1T_final[ (T_tot-1):1 , k ] = X[ temp ] # X_{1:T}^k = (X_1^{B_1^k},X_2^{B_2^k},...,X_T^k)
}
W_final = W[T_tot,]


# X_1T_final is a TxN matrix
## X_1T_final[,k] is a possible sequence of values
## W_final[k] is its associated weight
## P(X_{1:t} \in dX_{1:T}) = \sum_{k=1}^N W_final[k]*\delta_{X_1T_final[,k]}



par(mar = c(4,4,2,2), bty = "l")
matplot(X_1T_final, type = "l" )




# 2) Particle Gibbs ----------------------------------------------------------

data = y_true_1T
a <- b <- 0.01 # invgamma parameters

N = 200 # Particles number
G = 1000+1 # Gibbs iterations

theta_PG = matrix(NA,nrow = G, ncol = 2) 
#--> Gx2 matrix, theta_PG[g,i] is the value of the i-th unknown variable during the g-th iteration of the Gibbs-sampling

X_1T_PG  = matrix(NA,nrow = G, ncol = T_tot)
#--> GxT matrix, X_1T_PG[g,] is the sequence across all times (1:T) during the g-th iteration of the Gibbs-sampling

B_1T_PG  = matrix(NA,nrow = G, ncol = T_tot)
#--> GxT matrix, B_1T_PG[g,] is the sequence across all times (1:T) during the g-th iteration of the Gibbs-sampling


## Inizializzo:
theta_PG[1,] = c(1,1)
X_1T_PG[1, ] = rep(1,T_tot)
B_1T_PG[1, ] = sample(1:N, size = T_tot, replace = TRUE)



# Set true values (for debugging)
# X_1T_PG = t(apply(X_1T_PG, 1, function(x){x_true_1T}))


for(g in 2:G){
  cat("\n ++++++ it = ",g," ++++++ \n")
  # i) sample \theta(g) | y_{1:T}, X_{1:T}(g-1)
  post_shape   = c(a + T_tot/2, a + (T_tot-1)/2)
  temp_rate_1  = sum( (data - mean_g(X_1T_PG[g-1, ]))^2 )
  
  means_f      = mean_f( X_1T_PG[g-1, (1:T_tot-1) ], 2:T_tot )
  temp_rate_2  = sum( ( X_1T_PG[g-1, (2:T_tot) ]  - means_f )^2 )
  
  post_rate    = b + 0.5*c(temp_rate_1, temp_rate_2)
  inv_theta    = rgamma(n=2, shape = post_shape, rate = post_rate)
  
  theta_PG[g,] = 1/inv_theta
  
  # # ii) Conditional SMC
  X <- A <- B <- log_w <- w <- W <- matrix(NA, nrow = T_tot, ncol = N)

  # TxN matrix:
  ## X[t,] -> values of the N particles at time t
  ## X[,i] -> values across time

  # Time 1

  # a) Sample the particles
  X[1,] = rnorm( n = N, mean = 0, sd = sqrt(5) )
  X[1, B_1T_PG[g-1,1] ] = X_1T_PG[ g-1, 1 ]

  #--> this is a conditional Gibbs sampler. For ease of implementation, I'm first drawing all the values
  #    and then I fix the first one to the chosen one

  # b) Compute the un-nomralized weights
  log_w[1,] = eval_log_g(data[1], X[1,], theta_PG[g,1] )

  w[1,] = exp( log_w[1,] - max(log_w[1,]) )
  W[1,] = w[1,]/sum(w[1,])

  # Time n (n = 2...T)

  for(n in 2:T_tot){

    # a) Sample the parent node
    A[n-1,] = sample(1:N, size = N, W[n-1,], replace = TRUE)
    A[n-1, B_1T_PG[g-1, n] ] = B_1T_PG[g-1, n-1] # Cond. SMC

    # b) Sample the particles
    X_parents = X[n-1,A[n-1,]]
    X[n,] = rnorm( n = N, mean = mean_f(X_parents,n = n), sd = sqrt(theta_PG[g,2]) )
    X[n, B_1T_PG[g-1,n] ] = X_1T_PG[ g-1, n ]
    
    # c) Compute the un-normalized weights
    log_w[n,] = eval_log_g(data[n], X[n,], theta_PG[g,1] )

    w[n,] = exp( log_w[n,] - max(log_w[n,]) )
    W[n,] = w[n,]/sum(w[n,])

  }

  # Find ancestral lineage and save final paths
  X_1T_final = matrix(0,nrow = T_tot, ncol = N)
  B[T_tot,] = 1:N
  X_1T_final[T_tot,] = X[T_tot, B[T_tot,] ]


  for(k in 1:N){
    for(n in (T_tot-1):1){
      B[n,k] = A[ n, B[n+1,k] ] # B_n^k = A_n^{B_{n+1}^k}
    }

    temp = matrix(0,nrow = (T_tot-1):1, ncol = 2)
    temp[,1] = (T_tot-1):1
    temp[,2] = B[(T_tot-1):1,k]
    X_1T_final[ (T_tot-1):1 , k ] = X[ temp ] # X_{1:T}^k = (X_1^{B_1^k},X_2^{B_2^k},...,X_T^k)
  }
  W_final = W[T_tot,]


  # X_1T_final is a TxN matrix; W_final is a vector of length N
  ## X_1T_final[,k] is a possible sequence of values
  ## W_final[k] is its associated weight
  ## P(X_{1:t} \in dX_{1:T}) = \sum_{k=1}^N W_final[k]*\delta_{X_1T_final[,k]}

  #iii) Update X_1T_PG and B_1T_PG
  knew = sample(1:N, size = 1, prob = W_final)
  X_1T_PG[g,] = X_1T_final[,knew]
  B_1T_PG[g,] = B[,knew]
}


quantiles   = apply(theta_PG, 2, quantile, probs = c(0.0025,0.5,0.9975))
dens_thetas = apply(theta_PG, 2, function(x){density(x)$y})

par(mfrow = c(1,1), mar = c(4,4,2,2), bty = "l")
plot(0,0,type = "n",  main = "Thetas", xlab = "", ylab = "Prob.", xlim = range(quantiles), ylim = c(0,max(dens_thetas)) )
hist(theta_PG[,1], col = "darkgreen", nclass = "FD", freq = F, add = T)
hist(theta_PG[,2], col = "darkred", nclass = "FD", freq = F, add = T)
abline(v = apply(theta_PG, 2, median), col = c("darkgreen","darkred"), lwd = 2, lty = 3 )







par(mar = c(4,4,2,2), bty = "l")
matplot(t(X_1T_PG[(G-100):G,]), type = "l", xlab = "Time", ylab = "X_n" )



