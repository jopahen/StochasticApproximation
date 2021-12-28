#Analysis of sd-monitored early stopping RM for BS-IV of Asian put option
library(tictoc)
source("RM_IV_Put_Asian.r")

Nsims <- 20
batch_sizes <- c(20, 30, 50, 70, 100, 250, 500, 1000)
MSEs <- numeric(length = length(batch_sizes))
av_runtimes <- numeric(length = length(batch_sizes))
I <- 22

for(i in 1:length(batch_sizes)){
  price_err <- c()
  times <- c()
  for(nsim in 1:Nsims){
    tic()
    RM <- RM_IV_Asian(batch_sd = batch_sizes[i], sd_monitor = TRUE)
    t <- toc()
    price_err <- c(price_err, Put_Asian_pricer(sigma = RM$sigma)$price)
    times <- c(times, t$toc - t$tic)
  }
  MSEs[i] <- mean((I-price_err)^2)
  av_runtimes[i] <- mean(times)
  print(paste("batch size", batch_sizes[i],"done.", sep = " "))
}

plot(batch_sizes, log(MSEs), type = "l")
points(batch_sizes, log(MSEs), pch = 20)

plot(batch_sizes, av_runtimes, type = "l")
points(batch_sizes, av_runtimes, pch = 20)

