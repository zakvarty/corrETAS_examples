ETAS_RMSE <- function(mcmc, true_etas_values, burnin){
  squared_errors <- apply(mcmc$etas[-(1:burnin),], 1, function(x){(x - true_etas_values)^2})
  rmse <- apply(squared_errors, 1, function(x){sqrt(mean(x))})
  return(rmse)
}
