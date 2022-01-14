set.seed(310198)
#We implement the RM-algorithm to find the IV with regard to the Asian
#put option price, using importance sampling
source("BS_functions.r")

#MC pricer with IS using modified drift
Put_Asian_pricer_IS <- function(N = 10^5, S_0 = 100, r = 0.05, r_IS = -0.5,
                                sigma = 0.8, K = 150, T = 0.2){
  ISample <- S_path(N, S_0 = S_0, r = r_IS, sigma = sigma, T = T)
  frac <- g(ISample, K = K, r = r, T = T) * 
    lratio_vectorized(ISample[,-1], S_0 = S_0, r = r, r_IS = r_IS,
                      sigma = sigma, T = T)
  price <- mean(frac)
  se <- sd(frac)
  price_CI_lower <- price - qnorm(0.975) * se / sqrt(N)
  price_CI_upper <- price + qnorm(0.975) * se / sqrt(N)
  return(list(price = price, CI_lower = price_CI_lower,
              CI_upper = price_CI_upper, se = se))
}

RM_IV_Asian_IS <- function(n = 500, N = 10000, I = 49.3, sigma_0 = 1,
                           alpha_0 = 2/(150+100), rho = 0.8, K = 150,
                           batch_sd = 100, sd_monitor = FALSE, tol = 10^-4){
  sigma <- sigma_0
  #this time use MC-pricer with IS
  sigma_new <- sigma - alpha_0 * 
    (Put_Asian_pricer_IS(N, K = K, sigma = sigma)$price - I)
  sigmas <- sigma_new
  batch_sds <- c()
  iter <- 0
  learn_flag <- FALSE
  k <- 1
  while(iter < n){
    sigma <- sigma_new
    alpha <- alpha_0
    #adjust learning rate only when batch-sd starts to rise again
    if(learn_flag){
      alpha <- alpha_0 / k^rho
      k <- k+1
      if(alpha < tol) break
    }
    sigma_new <- sigma - alpha * 
      (Put_Asian_pricer_IS(N, K = K, sigma = sigma)$price - I)
    sigmas <- c(sigmas, sigma_new)
    iter <- iter + 1
    if(iter > batch_sd){
      batch_sds <- c(batch_sds, sd(sigmas[(iter - batch_sd):iter]))
    }
    if(length(batch_sds) > 1 & sd_monitor){
      if(batch_sds[length(batch_sds)] > batch_sds[length(batch_sds)-1]){
        learn_flag <- TRUE
      }
    }
    #print(iter)
  }
  return(list(sigma = sigma_new, sigmas = sigmas, batch_sds = batch_sds))
}

##diagnostic plots/validation:
#sigma_IV_RM_Asian_IS <- RM_IV_Asian_IS(sd_monitor = TRUE)
#plot(sigma_IV_RM_Asian_IS$sigmas, type = "l", main = "Evolution of RM-iterations",
#     ylab = "sigma", xlab = "Iterations", col = "blue")
#plot(sigma_IV_RM_Asian_IS$batch_sds, type = "l",
#     main = "Evolution of batch standard errors (batch size = 100)",
#     xlab = "Iterations - batch size",
#     ylab = "batch std. err.", col = "red")
#Put_Asian_pricer_IS(sigma = sigma_IV_RM_Asian_IS$sigma)