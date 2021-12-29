#experiments regarding importance sampling for the MC-pricer
source("BS_functions.r")

Put_Asian_pricer <- function(N = 10^5, S_0 = 100, r = 0.05, sigma = 0.8,
                             K = 120, T = 0.2){
  price <- mean(g(S_path(N, S_0 = S_0, r = r, sigma = sigma, T = T), K = K, r = r, T = T))
  se <- sd(g(S_path(N, S_0 = S_0, r = r, sigma = sigma, T = T), K = K, r = r, T = T))
  price_CI_lower <- price - qnorm(0.975) * se / sqrt(N)
  price_CI_upper <- price + qnorm(0.975) * se / sqrt(N)
  return(list(price = price, CI_lower = price_CI_lower,
              CI_upper = price_CI_upper, se = se))
}

Put_Asian_pricer_IS <- function(N = 10^5, S_0 = 100, r = 0.05, r_IS = 0.1,
                                sigma = 0.8, K = 120, T = 0.2){
  ISample <- S_path(N, S_0 = S_0, r = r_IS, sigma = sigma, T = T)
  frac <- g(ISample, K = K, r = r, T = T) * lratio_vectorized(ISample[,-1], S_0 = S_0, r = r, r_IS = r_IS, sigma = sigma, T = T)
  #frac <- g(ISample, K = K, r = r, T = T) * p_vectorized(ISample[,-1], S_0 = S_0, r = r, sigma = sigma, T = T) / p_vectorized(ISample[,-1], S_0 = S_0, r = r_IS, sigma = sigma, T = T)
  price <- mean(frac)
  se <- sd(frac)
  price_CI_lower <- price - qnorm(0.975) * se / sqrt(N)
  price_CI_upper <- price + qnorm(0.975) * se / sqrt(N)
  return(list(price = price, CI_lower = price_CI_lower,
              CI_upper = price_CI_upper, se = se))
}

test1 <- Put_Asian_pricer(K = 150, sigma = 0.8)
test2 <- Put_Asian_pricer_IS(K = 150, sigma = 0.8, r_IS = -0.5)

rs <- seq(-2,1,0.1)
ses <- numeric(length = length(rs))
prices <- numeric(length = length(rs))
for(i in 1:length(rs)){
  Pricer <- Put_Asian_pricer_IS(K = 150, sigma = 0.8, r_IS = rs[i])
  ses[i] <- Pricer$se
  prices[i] <- Pricer$price
}
plot(rs,ses, type = "l")
plot(rs,prices, type = "l")
