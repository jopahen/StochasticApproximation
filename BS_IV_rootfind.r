#We use the Brent root finding algorithm to calculate the implied volatility
#in the given setting:
source("BS_functions.r")

library(pracma) #package implementing the Brent root-finding algorithm

objective <- function(sigma) I(sigma) - 22 #define objective to find root

sigma_IV <- brent(objective, 0.1, 2)$root #run Brent w/ search area [0.1,2]