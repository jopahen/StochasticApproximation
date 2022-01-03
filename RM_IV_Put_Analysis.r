source("BS_IV_rootfind.r")
source("RM_IV_Put.r")
#We run the analysis of the RM-algorithm for calculating the BS implied vol.
#select n = 1000 iterations, N = 1, 10, 100, 1000 samples per iteration
#rho = 0.8 or 1
n <- 1000
Ns <- floor(10^seq(0,3,0.25))
rho <- 0.8
Nsims <- 20

MSE <- numeric(length = length(Ns))

for(i in 1:length(Ns)){
  N <- Ns[i]
  sims <- numeric(length = Nsims)
  for(j in 1:Nsims){
    sims[j] <- RM_IV(n,N, rho = rho)
  }
  MSE[i] <- mean((sims - sigma_IV)^2)
  print(MSE[i])
}

#plot simulation results
plot(log10(n*Ns), log10(MSE))
abline(a = 2, b = -1, col = "red", lty = "dashed")

