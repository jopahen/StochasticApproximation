set.seed(310198)
#Performance analysis of RM for BS-IV of Asian put option
library(tictoc) #package for time measurement
source("RM_IV_Put_Asian.r") #load RM-alg. with crude MC-pricer
source("RM_IV_Put_Asian_IS.r") #load RM-alg. with constant IS-drift MC-pricer
source("RM_IV_Put_Asian_IS_optimal_drift.r") #...with optimal IS-drift MC-pricer
#IMPORTANT: comment out diagnostics in the respective scripts before running
#this analysis!

Nsims <- 10
sigmas_crude <- numeric(length = length(Nsims))
times_crude <- numeric(length = length(Nsims))
sigmas_IS <- numeric(length = length(Nsims))
times_IS <- numeric(length = length(Nsims))
sigmas_opt_IS <- numeric(length = length(Nsims))
times_opt_IS <- numeric(length = length(Nsims))

#run RM-algorithms Nsims times and record estimates + runtimes
print("Starting performance analysis...")
for(i in 1:Nsims){
  tic()
  sigmas_crude[i] <- RM_IV_Asian(sd_monitor = TRUE)$sigma
  t <- toc()
  times_crude[i] <- t$toc - t$tic
  
  tic()
  sigmas_IS[i] <- RM_IV_Asian_IS(sd_monitor = TRUE)$sigma
  t <- toc()
  times_IS[i] <- t$toc - t$tic
  
  tic()
  sigmas_opt_IS[i] <- RM_IV_Asian_IS_OD(sd_monitor = TRUE)$sigma
  t <- toc()
  times_opt_IS[i] <- t$toc - t$tic
  print(paste("Iteration", i, "done.", sep = " "))
}

#average runtimes for each algorithm
mean(times_crude)
mean(times_IS)
mean(times_opt_IS)

#standard errors of the estimates
sd(sigmas_crude)
sd(sigmas_IS)
sd(sigmas_opt_IS)

#plot RM-estimates for all 3 algorithms
sds <- rbind(sigmas_crude,sigmas_IS,sigmas_opt_IS)
matplot(sds, col = "blue", pch = 20,
        axes = FALSE, ylab = "implied volatility", main = "RM-estimates")
axis(2)
axis(side=1, at=1:nrow(sds),
     labels=c("crude MC", "constant IS-drift", "optimal IS-drift"))
points(rowMeans(sds), pch = 19, col = "red")
