
library(corrETAS)
# Set true parameter values for simulation and plotting
mu_value <- 0.01
K_value <- 0.2
alpha_value <- 0.01
nu_t_value <- 0.1
xi_t_value <- 0

m_0_value <- 1.5

mag_values <- c(0.4,0)
mag_names <- c('nu_m', 'xi_m')

# Simulate catalgoue
cat <-  simulateETAScorr(
  mu = mu_value,
  K = K_value,
  alpha = alpha_value,
  gpd_t = c(nu_t_value, xi_t_value),
  rho = 0,
  gpd_m0 = mag_values,
  gpd_m1 = mag_values,
  m_0 = m_0_value,
  t_max = 50000,
  displayOutput = TRUE
)
plot(x = cat$ts,y = cat$ms, col = (cat$b>0)+1)

saveRDS(object = cat, file = "../corrETAS_examples/single_mag_example/cat.RDS")

C_value <- K_value *exp( alpha_value * (mean(cat$ms)-m_0_value))
etas_values <- c(mu_value, C_value, alpha_value, nu_t_value, xi_t_value)
etas_names <- c('mu', 'C', 'alpha', 'nu_t', 'xi_t')


# Collect posterior samples
mcmc <- sampleETASposterior(
  ts = cat$ts,
  ms = cat$ms,
  m_0 = 1.5,
  t_max = 50000,
  init_ETAS = c(0.01,0.1, 0.1, 0.1, 0.1),
  init_mag = c(1, 0),
  sims = 5000,
  B_samples = FALSE,
  B_fixed = FALSE,
  etas_sds = rep(0.1, 4),
  mag_sds = rep(0.1, 2),
  mu_prior = c(0.1, 0.1)
)

saveRDS(object = mcmc, file = "../corrETAS_examples/single_mag_example/mcmc.RDS")

pdf(file = "../corrETAS_examples/single_mag_example/mcmc_trace_plots.pdf")
opar <- par()
par(mar = c(5.1,4.1,2.1,1.1))
for(i in 1:5){
plot(mcmc$etas[,i], ylab = etas_names[i],type= 'l', bty = 'n')
  abline(h = quantile(x = mcmc$etas[,i],probs = c(0.025,0.975)), col = 3)
  abline(h = etas_values[i],col = 2)
}
for(i in 1:2){
  plot(mcmc$mag[,i], ylab = mag_names[i],type= 'l', bty = 'n')
  abline(h = quantile(x = mcmc$mag[,i],probs = c(0.025,0.975)), col = 3)
  abline(h = mag_values[i],col = 2)
}
par(mar = opar$mar)
dev.off()
