get_return_levels <- function(MCMC, return_periods, m_0, event_type = 0){

  rl_df <- data.frame(temp = rep(NA, length(MCMC$etas$mu)))

  for(i in seq_along(return_periods)){
    return_levels <- corrETAS::qgpd(
      p = 1 - 1/return_periods[i],
      nu = MCMC$mag[[1 + 2 * event_type]],
      shape = MCMC$mag[[2 + 2 * event_type]],
      mu = m_0)

    rl_df[,i] <- return_levels
  }

  colnames(rl_df) <- paste0("rl_",return_periods)

  return(rl_df)
}

plot_sampled_rls <- function(sampled_rls, alpha = 0.05, catalogue, event_type = 0, line_colours = c("lightblue", "blue", "lightblue"), add = FALSE, ...){
  return_periods <- as.integer(substr(x = colnames(sampled_rls), start = 4, stop = 100))

  point_estimates <- apply(sampled_rls, MARGIN = 2,FUN =  mean)
  CIs <- apply(sampled_rls, MARGIN = 2, FUN =  quantile, probs = c(1 - alpha/2, alpha / 2))

  if(!add){
    plot(
      x = rep(return_periods, each = 2),
      y = CIs,
      log = "xy",
      xlab = "return period",
      ylab = "return level",
      type = "n",
      ...
    )
  }

  lines(x = return_periods, y = CIs[1,], col = line_colours[1], lwd = 2)
  lines(x = return_periods, y = CIs[2,], col = line_colours[3], lwd = 2)
  lines(x = return_periods, y = point_estimates, col = line_colours[2], lwd = 2)

  if(event_type == 0){
    mag_subset <- cat$ms[cat$b == 0]
  } else if (event_type == 1){
    mag_subset <- cat$ms[cat$b > 0]
  } else if (event_type == 2){
    mag_subset <- cat$ms
  }

  points(
    x = (return_periods[return_periods < length(mag_subset)]),
    y = quantile(x = mag_subset, probs= 1 - 1/(return_periods[return_periods < length(mag_subset)]))
  )

}




get_return_levels_single <- function(MAG_MCMC, return_periods, m_0){

  rl_df <- data.frame(temp = rep(NA, length(MAG_MCMC$nu)))

  for(i in seq_along(return_periods)){
    return_levels <- corrETAS::qgpd(
      p = 1 - 1/return_periods[i],
      nu = MAG_MCMC$nu,
      shape = MAG_MCMC$xi,
      mu = m_0)

    rl_df[,i] <- return_levels
  }

  colnames(rl_df) <- paste0("rl_",return_periods)

  return(rl_df)
}

plot_sampled_rls_single <- function(sampled_rls, alpha = 0.05, catalogue, line_colours = c("lightblue","blue","lightblue"), line_type = 1, add = FALSE, ...){
  return_periods <- as.integer(substr(x = colnames(sampled_rls), start = 4, stop = 100))

  point_estimates <- apply(sampled_rls, MARGIN = 2,FUN =  mean, na.rm = TRUE)
  CIs <- apply(sampled_rls, MARGIN = 2, FUN =  quantile, probs = c(1 - alpha/2, alpha / 2), na.rm = TRUE)

  if(add == FALSE){
    plot(
      x = rep(return_periods, each = 2),
      y = CIs,
      log = "xy",
      xlab = "return period",
      ylab = "return level",
      type = "n",
      ...
    )
  }
  lines(x = return_periods, y = CIs[1,], col = line_colours[1], lwd = 2, lty = line_type)
  lines(x = return_periods, y = CIs[2,], col = line_colours[3], lwd = 2, lty = line_type)
  lines(x = return_periods, y = point_estimates, col = line_colours[2], lwd = 2, lty = line_type)

  mag_subset <- catalogue$ms

  points(
    x = (return_periods[return_periods < length(mag_subset)]),
    y = quantile(x = mag_subset, probs= 1 - 1/(return_periods[return_periods < length(mag_subset)]))
  )

}

qgpd_mix<- function(p, nu_0, xi_0, nu_1, xi_1,m_0, prop_trig, lower = 0 , upper = 20,...){
  objective <- function(x_p, p, nu_0, xi_0, nu_1, xi_1,m_0, prop_trig){
    p_0 <- corrETAS::pgpd(q = x_p, nu = nu_0, shape = xi_0, mu = m_0)
    p_1 <- corrETAS::pgpd(q = x_p, nu = nu_1, shape = xi_1, mu = m_0)
    (p_0 * (1 - prop_trig) + p_1 * prop_trig - p)^2
  }
  quantile <- optimise(f = objective,
                       p = p,
                       nu_0 = nu_0,
                       xi_0 = xi_0,
                       nu_1 = nu_1,
                       xi_1 = xi_1,
                       m_0 = m_0,
                       prop_trig = prop_trig,
                       upper = upper,
                       lower = lower,
                       ...)$minimum
  return(quantile)
}
# qgpd_mix(p = 0.5, nu_0 = 0.5, xi_0 = -0.1, nu_1 = 0.1, xi_1 = 1, m_0 = 1.5, prop_trig = 0.5)

mcmc_return_levels <- function(MCMC, return_period, m_0, ...){
  purrr::pmap_dbl(
    .f = qgpd_mix,
    .l = list(
      p = 1 - 1/return_period,
      m_0 = m_0,
      #nu_0 = MCMC$mag$nu_m0,
      #xi_0 = MCMC$mag$xi_m0,
      #nu_1 = MCMC$mag$nu_m1,
      #xi_1 = MCMC$mag$xi_m1,
      nu_0 = MCMC$mag$nu_m0,
      xi_0 = MCMC$mag$xi_m0,
      nu_1 = MCMC$mag$nu_m1,
      xi_1 = MCMC$mag$xi_m1,
      prop_trig = rowMeans(MCMC$b > 0),
      ...
    )
  )
}
