#We implement the RM-algorithm to find the IV with regard to the Asian
#put option price
library(tictoc)
source("BS_functions.r")

RM_IV_Asian <- function(n = 1000, N = 1000, I = 22, sigma_0 = 0.1, alpha_0 = 2/(120+100),
                  rho = 1, batch_sd = 50, sd_monitor = FALSE){
  sigma <- sigma_0
  sigma_new <- sigma - alpha_0 * (mean(g(S_path(N, sigma = sigma))) - I)
  sigmas <- sigma_new
  batch_sds <- c()
  iter <- 0
  err <- abs(sigma_new - sigma)
  while(iter < n){
    sigma <- sigma_new
    alpha <- alpha_0 / (iter+1)^rho
    sigma_new <- sigma - alpha_0 * (mean(g(S_path(N, sigma = sigma))) - I)
    sigmas <- c(sigmas, sigma_new)
    err <- abs(sigma_new - sigma)
    iter <- iter + 1
    if(iter > batch_sd){
      batch_sds <- c(batch_sds, sd(sigmas[(iter - batch_sd):iter]))
    }
    if(length(batch_sds) > 1 & sd_monitor){
      if(batch_sds[length(batch_sds)] > batch_sds[length(batch_sds)-1]) break
    }
    #print(iter)
  }
  return(list(sigma = sigma_new, sigmas = sigmas, batch_sds = batch_sds))
}

#we implement a high-iteration MC-pricer to validate accuracy
Put_Asian_pricer <- function(N = 10^5, S_0 = 100, r = 0.05, sigma = 0.713,
                             K = 120, T = 0.2){
  return(mean(g(S_path(N, S_0 = S_0, r = r, sigma = sigma, T = T))))
}

tic()
sigma_IV_RM_Asian <- RM_IV_Asian(sd_monitor = TRUE)
x <- toc()
#plot(sigma_IV_RM_Asian$sigmas)
#plot(sigma_IV_RM_Asian$batch_sds, type = "l")
#Put_Asian_pricer(sigma = sigma_IV_RM_Asian$sigma)


