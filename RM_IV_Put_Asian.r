#We implement the RM-algorithm to find the IV with regard to the Asian
#put option price
source("BS_functions.r")

RM_IV_Asian <- function(n = 200, N = 10000, I = 49.3, sigma_0 = 1, alpha_0 = 2/(150+100),
                  rho = 0.8, K = 150, batch_sd = 50, sd_monitor = FALSE){
  sigma <- sigma_0
  #RM-update (w/ crude MC-estimator):
  sigma_new <- sigma - alpha_0 * (mean(g(S_path(N, sigma = sigma), K = K)) - I)
  sigmas <- sigma_new
  batch_sds <- c()
  iter <- 0
  #iterate...
  while(iter < n){
    sigma <- sigma_new
    alpha <- alpha_0 / (iter+1)^rho
    sigma_new <- sigma - alpha_0 * (mean(g(S_path(N, sigma = sigma), K = K)) - I)
    sigmas <- c(sigmas, sigma_new)
    iter <- iter + 1
    #batch-sd monitoring routine: when batch-sd starts to go rise again,
    #terminate the algorithm (only if sd_monitor flag is TRUE)
    if(iter > batch_sd){
      batch_sds <- c(batch_sds, sd(sigmas[(iter - batch_sd):iter]))
    }
    if(length(batch_sds) > 1 & sd_monitor){
      if(batch_sds[length(batch_sds)] > batch_sds[length(batch_sds)-1]) break
    }
    print(iter)
  }
  return(list(sigma = sigma_new, sigmas = sigmas, batch_sds = batch_sds))
}

#we implement a high-iteration MC-pricer to validate accuracy
Put_Asian_pricer <- function(N = 10^5, S_0 = 100, r = 0.05, sigma = 0.8,
                             K = 120, T = 0.2){
  price <- mean(g(S_path(N, S_0 = S_0, r = r, sigma = sigma, T = T), K = K, r = r, T = T))
  se <- sd(g(S_path(N, S_0 = S_0, r = r, sigma = sigma, T = T), K = K, r = r, T = T))
  price_CI_lower <- price - qnorm(0.975) * se / sqrt(N)
  price_CI_upper <- price + qnorm(0.975) * se / sqrt(N)
  return(list(price = price, CI_lower = price_CI_lower,
              CI_upper = price_CI_upper, se = se))
}

#plot results
sigma_IV_RM_Asian <- RM_IV_Asian(sd_monitor = FALSE)
plot(sigma_IV_RM_Asian$sigmas, type = "l", main = "Evolution of RM-iterations",
     ylab = "sigma", xlab = "Iterations", col = "blue")
plot(sigma_IV_RM_Asian$batch_sds, type = "l",
     main = "Evolution of batch-standard errors (batch size = 50)",
     xlab = "Iterations - batch size",
     ylab = "batch std. err.", col = "red")
Put_Asian_pricer_IS(sigma = sigma_IV_RM_Asian$sigma)


