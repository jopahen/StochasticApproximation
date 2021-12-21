#We implement the RM-algorithm to find the IV with regard to the put option price
source("BS_functions.r")

RM_IV <- function(n = 10000, N = 100, I = 22, sigma_0 = 0.2, alpha_0 = 2/(120+100),
                  rho = 1){
  sigma <- sigma_0
  sigma_new <- sigma - alpha_0 * (mean(f(S(N, sigma = sigma))) - I)
  iter <- 0
  err <- abs(sigma_new - sigma)
  while(iter < n){
    sigma <- sigma_new
    alpha <- alpha_0 / (iter+1)^rho
    sigma_new <- sigma - alpha_0 * (mean(f(S(N, sigma = sigma))) - I)
    err <- abs(sigma_new - sigma)
    iter <- iter + 1
  }
  return(sigma_new)
}

#sigma_IV_RM <- RM_IV(alpha_0 = 2/(K+S_0))
