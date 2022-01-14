set.seed(31011998)
#We implement the RM-algorithm to find the IV with regard to the Asian
#put option price, having introduced a stopping procedure for the algorithm
source("BS_functions.r") #get auxiliary functions

#n = maximum iterations
#sd_monitor = flag to determine, if the stopping criterion should be used,
#otherwise the algorithm runs without adjusting learning rate!
#all other parameters according to notation in report
RM_IV_Asian <- function(n = 500, N = 10000, I = 49.3, sigma_0 = 1, alpha_0 = 2/(150+100),
                  rho = 0.8, K = 150, batch_sd = 100, sd_monitor = FALSE, tol = 10^-4){
  sigma <- sigma_0
  #RM-update (w/ crude MC-estimator):
  sigma_new <- sigma - alpha_0 * (mean(g(S_path(N, sigma = sigma), K = K)) - I)
  sigmas <- sigma_new
  batch_sds <- c()
  iter <- 0
  learn_flag <- FALSE
  k <- 1
  #iterate...
  while(iter < n){
    sigma <- sigma_new
    alpha <- alpha_0
    #adjust learning rate only when batch-sd starts to rise again
    if(learn_flag){
      alpha <- alpha_0 / k^rho
      k <- k+1
      if(alpha < tol) break
    } 
    sigma_new <- sigma - alpha * (mean(g(S_path(N, sigma = sigma), K = K)) - I)
    sigmas <- c(sigmas, sigma_new)
    iter <- iter + 1
    #batch-sd monitoring routine: when batch-sd starts to go rise again,
    #let the learning rate kick in...
    if(iter > batch_sd){
      batch_sds <- c(batch_sds, sd(sigmas[(iter - batch_sd):iter]))
    }
    if(length(batch_sds) > 1 & sd_monitor){
      if(batch_sds[length(batch_sds)] > batch_sds[length(batch_sds)-1]){
        learn_flag <- TRUE
      } 
    }
    #print(iter) #print to see that everything works (comment out if desired)
  }
  return(list(sigma = sigma_new, sigmas = sigmas, batch_sds = batch_sds))
}

#we implement a MC-pricer to validate accuracy
Put_Asian_pricer <- function(N = 10^5, S_0 = 100, r = 0.05, sigma = 0.8,
                             K = 120, T = 0.2){
  price <- mean(g(S_path(N, S_0 = S_0, r = r, sigma = sigma, T = T), K = K,
                  r = r, T = T))
  se <- sd(g(S_path(N, S_0 = S_0, r = r, sigma = sigma, T = T), K = K, r = r,
             T = T))
  price_CI_lower <- price - qnorm(0.975) * se / sqrt(N)
  price_CI_upper <- price + qnorm(0.975) * se / sqrt(N)
  return(list(price = price, CI_lower = price_CI_lower,
              CI_upper = price_CI_upper, se = se))
}

##plot results
#sigma_IV_RM_Asian <- RM_IV_Asian(sd_monitor = TRUE)
#plot(sigma_IV_RM_Asian$sigmas, type = "l", main = "Evolution of RM-iterations",
#     ylab = "sigma", xlab = "Iterations", col = "blue")
#plot(sigma_IV_RM_Asian$batch_sds, type = "l",
#     main = "Evolution of batch standard errors
#     (batch size = 100)",
#     xlab = "Iterations - batch size",
#     ylab = "batch std. err.", col = "red")

##give out results
#sigma_IV_RM_Asian$sigma
##check accuracy with high-iteration MC-pricer
#Put_Asian_pricer(sigma = sigma_IV_RM_Asian$sigma, K = 150)