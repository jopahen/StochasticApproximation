#We implement the RM-algorithm to find the IV with regard to the put option price
source("BS_functions.r") #load auxiliary functions

#note that we provide default input data according to the hint for alpha_0
RM_IV <- function(n = 1000, N = 10000, I = 22, sigma_0 = 0.2, alpha_0 = 2/(120+100),
                  rho = 1, lrate_delay = 200){
  sigma <- sigma_0
  sigmas <- sigma_0
  #RM-update:
  sigma_new <- sigma - alpha_0 * (mean(f(S(N, sigma = sigma))) - I) #crude MC-est.
  iter <- 0
  #iterate...
  while(iter < n){
    sigma <- sigma_new
    alpha <- alpha_0
    #learning rate is kept constant until lrate_delay iterations
    #to give the algorithm time to evolve to the area
    #of the solution, then it is modulated as proposed in the project description
    if(iter > lrate_delay) alpha <- alpha_0 / (iter+1-lrate_delay)^rho
    sigma_new <- sigma - alpha * (mean(f(S(N, sigma = sigma))) - I) #crude MC-est.
    sigmas <- c(sigmas, sigma_new)
    iter <- iter + 1
  }
  return(list(sigma = sigma_new, sigmas = sigmas))
}

##plot results
#sigma_IV_RM <- RM_IV()
#plot(sigma_IV_RM$sigmas, type = "l", col = "blue", ylab = "sigma",
#     main = "RM-iterates")
