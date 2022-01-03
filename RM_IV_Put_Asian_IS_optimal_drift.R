#We implement the RM-algorithm to find the IV with regard to the Asian
#put option price, using importance sampling with optimal drift
source("BS_functions.r")
source("optimal_drift.r")

Put_Asian_pricer_IS <- function(N = 10^5, S_0 = 100, r = 0.05, r_IS = -0.5,
                                sigma = 0.8, K = 150, T = 0.2){
  ISample <- S_path(N, S_0 = S_0, r = r_IS, sigma = sigma, T = T)
  frac <- g(ISample, K = K, r = r, T = T) * lratio_vectorized(ISample[,-1], S_0 = S_0, r = r, r_IS = r_IS, sigma = sigma, T = T)
  price <- mean(frac)
  se <- sd(frac)
  price_CI_lower <- price - qnorm(0.975) * se / sqrt(N)
  price_CI_upper <- price + qnorm(0.975) * se / sqrt(N)
  return(list(price = price, CI_lower = price_CI_lower,
              CI_upper = price_CI_upper, se = se))
}

RM_IV_Asian_IS_OD <- function(n = 1000, N = 1000, I = 49.3, sigma_0 = 1, alpha_0 = 2/(150+100),
                           rho = 1, K = 150, batch_sd = 100, sd_monitor = FALSE){
  sigma <- sigma_0
  drifts <- optimal_r(sigma)
  sigma_new <- sigma - alpha_0 * (Put_Asian_pricer_IS(N, K = K, sigma = sigma, r_IS = drifts)$price - I)
  sigmas <- sigma_new
  batch_sds <- c()
  iter <- 0
  err <- abs(sigma_new - sigma)
  while(iter < n){
    sigma <- sigma_new
    r_IS <- optimal_r(sigma)
    alpha <- alpha_0 / (iter+1)^rho
    sigma_new <- sigma - alpha_0 * (Put_Asian_pricer_IS(N, K = K, sigma = sigma, r_IS = r_IS)$price - I)
    sigmas <- c(sigmas, sigma_new)
    drifts <- c(drifts, r_IS)
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
  return(list(sigma = sigma_new, sigmas = sigmas, batch_sds = batch_sds, drifts = drifts))
}

RM <- RM_IV_Asian_IS_OD(sd_monitor = FALSE)
plot(RM$sigmas, type = "l")
plot(RM$batch_sds, type = "l")
plot(RM$drifts, type = "l")
Put_Asian_pricer_IS(sigma = RM$sigma)
