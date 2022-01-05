source("BS_IV_rootfind.r") #compute reference IV via root-finder (Brent)
source("RM_IV_Put.r") #load RM-algorithm
#We run the analysis of the RM-algorithm for calculating the BS implied vol.
#select n = 1000 iterations, N = 1, 3, 10, 31, 100, 316, 1000, 3162, 10000
#samples per iteration
#rho = 0.8 or 1
n <- 1000 #number of iterations of RM-algorithm
Ns <- floor(10^seq(0,4,0.5)) #sample sizes for MC-estimator to be tested
rho <- 1 #can be adjusted to 0.8 or 1 as desired
Nsims <- 20 #number of simulations to estimate MSE

MSE <- numeric(length = length(Ns))

for(i in 1:length(Ns)){
  N <- Ns[i]
  sims <- numeric(length = Nsims)
  for(j in 1:Nsims){
    sims[j] <- RM_IV(n,N, rho = rho)$sigma
  }
  MSE[i] <- mean((sims - sigma_IV)^2)
}

#plot simulation results
par(pty="s")
plot(log10(n*Ns), log10(MSE), pch = 20, main = "MSE decay for RM-estimates of
     implied volatility, rho = 1", xlab = "lg(N * n)", ylab = "lg(MSE)",
     asp = 1) #plot MSEs on log-scales
abline(a = 2, b = -1, col = "red", lty = "dashed") #O(1/N) for comparison

