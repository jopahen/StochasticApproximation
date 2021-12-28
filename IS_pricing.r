#experiments regarding importance sampling for the MC-pricer
source("BS_functions.r")

Put_Asian_pricer <- function(N = 10^4, S_0 = 100, r = 0.05, sigma = 0.713,
                             K = 120, T = 0.2){
  price <- mean(g(S_path(N, S_0 = S_0, r = r, sigma = sigma, T = T), K = K, r = r, T = T))
  se <- sd(g(S_path(N, S_0 = S_0, r = r, sigma = sigma, T = T), K = K, r = r, T = T))
  price_CI_lower <- price - qnorm(0.975) * se / sqrt(N)
  price_CI_upper <- price + qnorm(0.975) * se / sqrt(N)
  return(list(price = price, CI_lower = price_CI_lower,
              CI_upper = price_CI_upper, se = se/sqrt(N)))
}

Put_Asian_pricer_IS <- function(N = 10^4, S_0 = 100, r = 0.05, r_IS = 0.1,
                                sigma = 0.713, K = 120, T = 0.2){
  ISample <- S_path(N, S_0 = S_0, r = r_IS, sigma = sigma, T = T)
  frac <- g(ISample, K = K, r = r, T = T) * p_vectorized(ISample, S_0 = S_0, r = r, sigma = sigma, T = T) / p_vectorized(ISample, S_0 = S_0, r = r_IS, sigma = sigma, T = T)
  price <- mean(frac)
  se <- sd(frac)
  price_CI_lower <- price - qnorm(0.975) * se / sqrt(N)
  price_CI_upper <- price + qnorm(0.975) * se / sqrt(N)
  return(list(price = price, CI_lower = price_CI_lower,
              CI_upper = price_CI_upper, se = se/sqrt(N)))
}

test1 <- Put_Asian_pricer()
test2 <- Put_Asian_pricer_IS(r_IS = -0.1)

rs <- seq(-1.5,1,0.1)
ses <- numeric(length = length(rs))
for(i in 1:length(rs)){
  ses[i] <- Put_Asian_pricer_IS(r_IS = rs[i])$se
}
plot(rs,ses, type = "l")
