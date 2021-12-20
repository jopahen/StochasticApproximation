#Functions with respect to the Black-Scholes model,
#with relevant parameters as defaults...

#discounted put option payoff
f <- function(S = 100, K = 120, r = 0.05, T = 0.2){
  return(exp(-r*T) * max(K-S,0))
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
