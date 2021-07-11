plot_ETAS_chains <- function(mcmc, add_truth = TRUE, true_values){
  for(i in 1:5){
    plot(
      y = mcmc$etas[[i]],
      x = seq_along(mcmc$etas[[i]]),
      xlab = "Index",
      ylab = paste0(etas_names[i]),
      type = "l", col = 1 + i)
    if(add_truth){
      abline(h = true_values[i], lwd = 2)
    }
  }
}

plot_mag_chains <- function(mcmc, add_truth = TRUE, true_values, is_corr = FALSE){
  for(i in 1:(4 + is_corr)){
    plot(
      y = mcmc$mag[[i]],
      x = seq_along(mcmc$mag[[i]]),
      xlab = "Index",
      ylab = paste0(mag_names_corr[i]),
      type = "l", col = 1 + i)
    if(add_truth){
      abline(h = true_values[i], lwd = 2)
    }
  }
}
