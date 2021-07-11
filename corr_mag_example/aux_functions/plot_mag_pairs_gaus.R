plot_mag_pairs_gaus <- function(m, b,mag_pars, m_0, text_loc = c(-2,2), text_cex = 1,...){
  nu_mb <- mag_pars[1]
  xi_mb <- mag_pars[2]
  nu_ma <- mag_pars[3]
  xi_ma <- mag_pars[4]
  sig_mb <- nu_mb / (1 + xi_mb)
  sig_ma <- nu_ma / (1 + xi_ma)

  triggered_indices <- which(b>0)
  parent_indices <- b[triggered_indices]
  parent_triggered <- (b[parent_indices]>0)

  m_trig <- m[triggered_indices]
  z_trig <- seq_along(m_trig)
  for(i in seq_along(z_trig)){
   z_trig[i] <-  corrETAS::gpd_to_gaus(x = m_trig[i], sig = sig_ma, xi = xi_ma, mu = m_0)
  }

  m_par <- m[parent_indices]
  z_par <- seq_along(m_par)
  for (i in seq_along(z_par)){
    sig_i <- (sig_mb * !parent_triggered[i]) + (sig_ma * parent_triggered[i])
    xi_i <- (xi_mb  * !parent_triggered[i]) + (xi_ma  * parent_triggered[i])
    print(c(i, sig_i, xi_i))
    z_par[i] <- corrETAS::gpd_to_gaus( x = m_par[i], sig = sig_i, xi = xi_i, mu = m_0 )
  }

  plot(x = z_par, y = z_trig, pch = 16, col = parent_triggered + 1, ...)
  text(x = text_loc[1],
       y =  text_loc[2],
       labels = paste0("Sample correlation: \n ", round(cor(z_par, z_trig),4)),
       cex = text_cex)
}

