#set.seed(310198)
#Functions with respect to the Black-Scholes model,
#with relevant parameters as defaults

#price process simulation (at maturity)
#we use the GBM-formula and the standard sampling routines for the normal dist.
S <- function(N, S_0 = 100, r = 0.05, sigma = 0.1, T = 0.2){
  return(S_0 * exp((r - sigma^2/2)*T + sigma * sqrt(T) * rnorm(N)))
}

#discounted European put option payoff
f <- function(S = 100, K = 120, r = 0.05, T = 0.2){
  return(exp(-r*T) * pmax(K-S,0))
}

#price process simulation (path) with discrete (uniform) monitoring:
#the function computes a matrix with N rows, each row has m+1 entries
#corresponding to a path of the price process. Essentially, this is Alg. 3.4
#from the lecture notes. For efficiency, the normal increments are generated
#simultaneously and all calculations have been vectorized.
S_path <- function(N, m = 50, S_0 = 100, r = 0.05, sigma = 0.1,
                   T = 0.2){
  dt <- T / m
  t <- seq(0,T,dt)
  X <- matrix(rnorm((m+1)*N, 0, sqrt(dt)), nrow = N)
  X[,1] <- numeric(length = N)
  X <- t(apply(X,1,cumsum))
  S <- S_0 * exp(t( (r-sigma^2/2) * t + t(sigma * X) ))
  return(S)
}

#discounted put option payoff (Asian), can handle multiple paths as matrix
#generated from previous function
g <- function(S_path, K = 120, r = 0.05, T = 0.2){
  S <- S_path[,-1]
  if(is.null(dim(S))) return(exp(-r*T) * pmax(K - mean(S),0))
  return(as.vector(exp(-r*T) * pmax(K - rowMeans(S),0)))
}

#auxiliary function w for BS-formula wrt. volatility
w <- function(sigma = 0.1, S = 100, K = 120, r = 0.05, T = 0.2){
  return((log(K/S) - (r-sigma^2/2) * T) / sigma / sqrt(T))
}

#pricing formula (Black-Scholes formula)
Put <- function(sigma = 0.1, S = 100, K = 120, r = 0.05, T = 0.2){
  return(exp(-r*T) * K * pnorm(w(sigma, S, K, r, T)) - 
           S * pnorm(w(sigma, S, K, r, T) - sigma * sqrt(T)))
}

#pricing formula wrt. sigma as objective for root-finder
I <- function(sigma) Put(sigma)

#drift-shift likelihood ratio (single path vector), computes likelihood ratio
#as derived in Lemma 3.2 to a single vector of length m
lratio <- function(x, S_0 = 100, r = 0.05, r_IS = 0.1, sigma = 0.1, T = 0.2){
  m <- length(x)
  dt <- T/m
  x_end <- x[length(x)]
  out <- exp(-(r_IS-r) * (2 * log(x_end/S_0) - m * dt * (r + r_IS - sigma^2)) 
             / (2 * sigma^2))
  return(out)
}

#drift-shift likelihood ratio (vectorized), employs the previous function to
#compute likelihood-ratios for multiple vectors of length m simultaneously
lratio_vectorized <- function(X, S_0 = 100, r = 0.05, r_IS = 0.1, sigma = 0.1,
                              T = 0.2){
  if(is.null(dim(X))) return(as.numeric(lapply(X, lratio, S_0 = S_0, r = r,
                                               r_IS = r_IS, sigma = sigma, T = T)))
  return(apply(X, 1, lratio, S_0 = S_0, r = r, r_IS = r_IS, sigma = sigma, T = T))
}