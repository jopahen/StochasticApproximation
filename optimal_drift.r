#finding optimal drift
library(splines)
source("BS_functions.r")

#objective (expression (18) in project description)
#to be minimized with respect to drift r_IS (log-transformed)
objective <- function(r_IS, S = S_sim, sigma = 0.1, r = 0.05,
                      K = 150, T = 0.2){
  objective <- numeric(length = length(r_IS))
  for(i in 1:length(r_IS)){
    objective[i] <- log(mean(g(S, K, r, T)^2 * 
                               lratio_vectorized(S[,-1], S[1,1], r,
                                                 r_IS = r_IS[i], sigma, T)))
  }
  return(objective)
}

#simulate paths and optimize r_IS for given range of volatilities (sigmas)
N <- 10000
sigmas <- seq(0.1,5,0.1)
optimal_rs <- numeric(length = length(sigmas))

for(i in 1:length(sigmas)){
  S_sim <- S_path(N, S_0 = 100, r = 0.05, sigma = sigmas[i], T = 0.2)
  optimal_rs[i] <- optimize(objective, interval = c(-6,1),
                            sigma = sigmas[i])$minimum
}

#plot optimal r_IS wrt. sigmas
plot(sigmas, optimal_rs, pch = 20)

#spline interpolation (degree 5) of optimal drifts based on sigma
data <- data.frame(sigmas, optimal_rs)
spline_fit <- lm(optimal_rs ~ bs(sigmas, degree = 5), data = data)

optimal_r <- function(sigma){
  return(as.numeric(predict(spline_fit, newdata = list(sigmas = sigma))))
}

#plot spline interpolation
x <- seq(0,5,0.01)
lines(x, optimal_r(x), type = "l", col = "red")
