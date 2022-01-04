#We implement the RM-algorithm to find the IV with regard to the put option price
source("BS_functions.r") #load auxiliary functions

#note that we provide default input data according to the hint for alpha_0
RM_IV <- function(n = 1000, N = 10000, I = 22, sigma_0 = 0.2, alpha_0 = 2/(120+100),
                  rho = 1){
  sigma <- sigma_0
  #RM-update:
  sigma_new <- sigma - alpha_0 * (mean(f(S(N, sigma = sigma))) - I) #crude MC-est.
  iter <- 0
  #iterate...
  while(iter < n){
    sigma <- sigma_new
    alpha <- alpha_0 / (iter+1)^rho
    sigma_new <- sigma - alpha_0 * (mean(f(S(N, sigma = sigma))) - I) #crude MC-est.
    iter <- iter + 1
  }
  return(sigma_new)
}

sigma_IV_RM <- RM_IV()
