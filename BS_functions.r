#Functions with respect to the Black-Scholes model,
#with relevant parameters as defaults...
library(mvtnorm)

#price process simulation (maturity)
S <- function(N, S_0 = 100, r = 0.05, sigma = 0.1, T = 0.2){
  return(S_0 * exp((r - sigma^2/2)*T + sigma * sqrt(T) * rnorm(N)))
}

#discounted put option payoff
f <- function(S = 100, K = 120, r = 0.05, T = 0.2){
  return(exp(-r*T) * pmax(K-S,0))
}

#price process simulation (path)
S_path <- function(N, m = 50, S_0 = 100, r = 0.05, sigma = 0.1,
                   T = 0.2){
  dt <- T / m
  t <- seq(0,T,dt)
  X <- matrix(rnorm((m+1)*N, 0, sqrt(dt)), nrow = N)
  X[,1] <- numeric(length = N)
  X <- t(apply(X,1,cumsum))
  S <- S_0 * exp((r-sigma/2) * t + sigma * X)
  return(S)
}

#discounted put option payoff (asian)
g <- function(S_path, K = 120, r = 0.05, T = 0.2){
  S <- S_path[,-1]
  if(is.null(dim(S))) return(exp(-r*T) * pmax(K - mean(S),0))
  return(as.vector(exp(-r*T) * pmax(K - rowMeans(S),0)))
}

#auxiliary function w for BS-formula wrt. volatility
w <- function(sigma = 0.1, S = 100, K = 120, r = 0.05, T = 0.2){
  return((log(K/S) - (r-sigma^2/2) * T) / sigma / sqrt(T))
}

#pricing formula
Put <- function(sigma = 0.1, S = 100, K = 120, r = 0.05, T = 0.2){
  return(exp(-r*T) * K * pnorm(w(sigma, S, K, r, T)) - 
           S * pnorm(w(sigma, S, K, r, T) - sigma * sqrt(T)))
}

#pricing formula wrt. sigma
I <- function(sigma) Put(sigma)

#price process pdf (single path vector)
p <- function(x, S_0 = 100, r = 0.05, sigma = 0.1, T = 0.2){
  m <- length(x)
  dt <- T/m
  S <- c(S_0, x)
  factors <- c()
  for(i in 1:length(x)){
    d <- (log(S[i+1]/S[i]) - (r - sigma^2/2) * dt) / (sigma * sqrt(dt))
    factor <- dnorm(d) / (S[i+1] * sigma * sqrt(dt))
    factors <- c(factors, factor)
  }
  return(prod(factors))
}

#price process pdf (vectorized)
p_vectorized <- function(X, S_0 = 100, r = 0.05, sigma = 0.1, T = 0.2){
  if(is.null(dim(X))) return(as.numeric(lapply(X, p, S_0 = S_0, r = r,
                                              sigma = sigma, T = T)))
  return(apply(X, 1, p, S_0 = S_0, r = r, sigma = sigma, T = T))
}

#S_test <- S_path(100000)[,-1]
#p_vectorized(S_test)
