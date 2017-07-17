## NashSutcliffe
# Function is intended to calculate Nash-Sutcliffe Efficiency given 
# vectors of simulated and observed data.
#
# Example data:
# x <- seq(1,100,1)
# observed <- 0.75*(x^2)
# simulated <- 0.75*(x^2)*rnorm(100, mean=1, sd=0.1)
# out.NSE <- NashSutcliffe(simulated, observed)

NashSutcliffe <- function(sim, obs) {
  if (length(sim) != length(obs)) stop("vectors not the same size")
  return(1-sum((obs-sim)^2)/sum((obs-mean(obs))^2)  )
}

## RMSE
# Function is intended to calculate RMSE given vectors of simulated and observed data.
#
# Example data (use above):
# out.RMSE <- RMSE(simulated, observed)

RMSE <- function(sim, obs) {
  if (length(sim) != length(obs)) stop("vectors not the same size")
  return(mean((obs-sim)^2)^0.5)
}

NRMSE <- function(sim, obs){
  if (length(sim) != length(obs)) stop("vectors not the same size")
  return((mean((obs-sim)^2)^0.5)/abs(max(obs)-min(obs)))
}

## R2
# Function is intended to calculate R^2, the coefficient of determination.
#
# out.R2 <- R2(simulated, observed)

R2 <- function(sim, obs) {
  if (length(sim) != length(obs)) stop("vectors not the same size")
  return((sum((obs-mean(obs))*(sim-mean(sim)))/
    ((sum((obs-mean(obs))^2)^0.5)*(sum((sim-mean(sim))^2)^0.5)))^2)
}