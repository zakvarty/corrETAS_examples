## Example code for simulating and fitting correlated magnitude ETAS model.
##
##  Two catalogues are simulated, one with dual magnitudes and one with
##  correlated magnitudes. The dual and correlated magnitude ETAS models are
##  each fitted to both of these catalgoues.
##
##  Plots are made using the resulting MCMC chains to show that:
##    - the correlated magnitude code can recover the true parameter values from
##      a correlated catalgoue;
##    - the correlated magntiude code can recover the true parameter values from
##      a dual catalogue;
##    - the branching vector is equally well recovered by each ETAS model for
##      a dual catalgoue;
##    - the branching vector recovery is impaired by falsely assuming magnitudes
##      are independent.

source("./aux_functions/mcmc_plotting.R")
source("./aux_functions/ETAS_RMSE.R")
## 0.0: Simulate dual and correlated ETAS catalogues ---------------------------

## 0.1: Set parameter values ---------------------------------------------------

mu_value <- 0.02
K_value <- 0.2
alpha_value <- 0.1
nu_t_value <- 0.1
xi_t_value <- 0

m_0_value <- 1.5

mag_values_dual <- c(0.6 , 0, 0.1, 0, 0.0)
mag_names_dual <- c('nu_m0', 'xi_m0', 'nu_m1', 'xi_m1', 'rho')

mag_values_corr <- c(0.6 , 0, 0.1, 0, 0.6)
mag_names_corr <- c('nu_m0', 'xi_m0', 'nu_m1', 'xi_m1', 'rho')

# 0.2: Simulate a dual magntude catalogue --------------------------------------
set.seed(2212)
cat_dual <-  corrETAS::simulateETAScorr(
  mu = mu_value,
  K = K_value,
  alpha = alpha_value,
  gpd_t = c(nu_t_value, xi_t_value),
  rho = mag_values_dual[5],
  gpd_m0 = mag_values_dual[1:2],
  gpd_m1 = mag_values_dual[3:4],
  m_0 = m_0_value,
  t_max = 50000,
  displayOutput = TRUE
)

# 0.3: Simulate a correlated magnitude catalogue -------------------------------
set.seed(2212)
cat_corr <-  corrETAS::simulateETAScorr(
  mu = mu_value,
  K = K_value,
  alpha = alpha_value,
  gpd_t = c(nu_t_value, xi_t_value),
  rho = mag_values_corr[5],
  gpd_m0 = mag_values_corr[1:2],
  gpd_m1 = mag_values_corr[3:4],
  m_0 = m_0_value,
  t_max = 50000,
  displayOutput = TRUE
)

# 0.4: Calculate useful summaries and statistics -------------------------------

type_dual <- cat_dual$b > 0
type_corr <- cat_corr$b > 0

colouring_dual <- ifelse(!type_dual, "grey", "blue")
colouring_corr <- ifelse(!type_corr, "grey", "blue")

plot(x = cat_dual$ts,y = cat_dual$ms, col = colouring_dual, type = "h")
plot(x = cat_corr$ts,y = cat_corr$ms, col = colouring_corr, type = "h")

C_value_dual <- K_value * exp( alpha_value * (mean(cat_dual$ms) - m_0_value))
C_value_corr <- K_value * exp( alpha_value * (mean(cat_corr$ms) - m_0_value))

etas_values_dual <- c(mu_value, C_value_dual, alpha_value, nu_t_value, xi_t_value)
etas_values_corr <- c(mu_value, C_value_corr, alpha_value, nu_t_value, xi_t_value)
etas_names <- c('mu', 'C', 'alpha', 'nu_t', 'xi_t')

n_trig_dual <- sum(cat_dual$b > 0)
n_trig_corr <- sum(cat_corr$b > 0)

m_bar_dual <- mean(cat_dual$ms)
m_bar_corr <- mean(cat_corr$ms)

## 1.0: MCMC for each model-data combination ------------------------------------

## 1.1: Dual magnitude model, Dual magnitude data -----------------------------

# set.seed(1993)
# initial_B <- seq_along(cat_dual$ts) - 1
# initial_B <- initial_B * rbinom(n = length(cat_dual$ts), size = 1, prob = 0.4)
#
# mcmc_dual_mod_dual_cat <- corrETAS::sampleETASposterior_dualmag(
#   ts = cat_dual$ts,
#   ms = cat_dual$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = c(0.1,0.1,0.1,0.5,0),
#   init_B = initial_B,
#   init_mag = c(1,0,0.5,0),
#   sims = 5000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05),
#   mu_prior = c(0.001,0.001)
# )
# saveRDS(mcmc_dual_mod_dual_cat, file = "./output/data/mcmc_dual_mod_dual_cat.RDS")
mcmc_dual_mod_dual_cat <- readRDS("./output/data/mcmc_dual_mod_dual_cat.RDS")

## 1.2: Dual magnitude model, Correlated magnitude data -----------------------

# set.seed(1993)
# initial_B <- (seq_along(cat_dual$ts)-1) * rbinom(n = length(cat_dual$ts),
#                                              size = 1,
#                                              prob = 0.4)
#
# mcmc_dual_mod_corr_cat <- corrETAS::sampleETASposterior_dualmag(
#   ts = cat_corr$ts,
#   ms = cat_corr$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = c(0.1,0.1,0.1,0.5,0),
#   init_B = initial_B,
#   init_mag = c(1,0,0.5,0),
#   sims = 5000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05),
#   mu_prior = c(0.001,0.001)
# )
#saveRDS(mcmc_dual_mod_corr_cat, file = "./output/data/mcmc_dual_mod_corr_cat.RDS")
mcmc_dual_mod_corr_cat <- readRDS("./output/data/mcmc_dual_mod_corr_cat.RDS")

## 1.3: Correlated magnitude model, Dual magnitude data -----------------------

# init_prop_trig <- 0.4
# init_mu <- (1 - init_prop_trig) * length(cat_dual$b) / 50000
#
# set.seed(1993)
#
# init_etas <- c(init_mu, 0.4, 0.01, 50, 0)
# init_mag <- c(1,0,0.5,0,0)
# init_b <- (seq_along(cat_dual$b) - 1) * rbinom(length(cat_dual$b),
#                                                size = 1,
#                                                prob = 1 - init_prop_trig)
#
# plot(x = seq(0,10,length.out = 1001), y = dlnorm(x = seq(0,10,length.out = 1001),meanlog = log(0.2),sdlog = 1000))
#
# mcmc_corr_mod_dual_cat <- corrETAS::sampleETASposterior_corrmag(
#   ts = cat_dual$ts,
#   ms = cat_dual$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = init_etas,
#   init_B = init_b,
#   init_mag = init_mag,
#   sims = 5000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05,0.05),
#   mu_prior = c(0.0001,0.0001),
#   logC_prior = c(log(0.2), 10)
# )
# saveRDS(mcmc_corr_mod_dual_cat, file = "./output/mcmc_corr_mod_dual_cat.RDS")
mcmc_corr_mod_dual_cat <- readRDS("./output/data/mcmc_corr_mod_dual_cat.RDS")


## 1.4: Correlated magnitude model, Correlated magnitude data -----------------

# init_prop_trig <- 0.4
# init_mu <- (1 - init_prop_trig) * length(cat_dual$b) / 50000
#
# set.seed(1993)
#
# init_etas <- c(init_mu, 0.4, 0.01, 50, 0)
# init_mag <- c(1,0,0.5,0,0)
# init_b <- (seq_along(cat_corr$b) - 1) * rbinom(length(cat_corr$b),
#                                                size = 1,
#                                                prob = 1 - init_prop_trig)
#
# mcmc_corr_mod_corr_cat <- corrETAS::sampleETASposterior_corrmag(
#   ts = cat_corr$ts,
#   ms = cat_corr$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = init_etas,
#   init_B = init_b,
#   init_mag = init_mag,
#   sims = 5000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05,0.05),
#   mu_prior = c(0.0001,0.0001),
#   logC_prior = c(log(0.2), 10)
# )
# saveRDS(mcmc_corr_mod_corr_cat, file = "./output/data/mcmc_corr_mod_corr_cat.RDS")
mcmc_corr_mod_corr_cat <- readRDS("./output/data/mcmc_corr_mod_corr_cat.RDS")

## 2.0: Plotting chains --------------------------------------------------------

## 2.1: Traceplots for MCMC using dual model and dual catalogue ----------------
pdf(file = "./output/plots/traceplots/mcmc_dual_mod_dual_cat_traceplots.pdf", width = 7, height = 5)
  plot_ETAS_chains(mcmc_dual_mod_dual_cat, true_values = etas_values_dual)
  plot_mag_chains(mcmc_dual_mod_dual_cat, true_values = mag_values_dual, is_corr = FALSE)
  plot( rowSums(mcmc_dual_mod_dual_cat$b > 0), type = "l", col = 2, xlab = "Index", ylab = "n_trig")
  abline(h = n_trig_corr,lwd = 2)
dev.off()

## 2.2: Traceplots for MCMC using dual model and corr catalogue ----------------
pdf(file = "./output/plots/traceplots/mcmc_dual_mod_corr_cat_traceplots.pdf", width = 7, height = 5)
  plot_ETAS_chains(mcmc_dual_mod_corr_cat, true_values = etas_values_corr)
  plot_mag_chains(mcmc_dual_mod_corr_cat, true_values = mag_values_corr, is_corr = FALSE)
  plot( rowSums(mcmc_dual_mod_corr_cat$b > 0), type = "l", col = 2, xlab = "Index", ylab = "n_trig")
  abline(h = n_trig_corr,lwd = 2)
dev.off()

## 2.3: Traceplots for MCMC using corr model and dual catalogue ----------------
pdf(file = "./output/plots/traceplots/mcmc_corr_mod_dual_cat_traceplots.pdf", width = 7, height = 5)
  plot_ETAS_chains(mcmc_corr_mod_dual_cat, true_values = etas_values_dual)
  plot_mag_chains(mcmc_corr_mod_dual_cat, true_values = mag_values_dual, is_corr = FALSE)
  plot( rowSums(mcmc_corr_mod_dual_cat$b > 0), type = "l", col = 2, xlab = "Index", ylab = "n_trig")
  abline(h = n_trig_corr,lwd = 2)
dev.off()

## 2.4: Traceplots for MCMC using corr model and corr catalogue ----------------
pdf(file = "./output/plots/traceplots/mcmc_corr_mod_corr_cat_traceplots.pdf", width = 7, height = 5)
  plot_ETAS_chains(mcmc_corr_mod_corr_cat, true_values = etas_values_corr)
  plot_mag_chains(mcmc_corr_mod_corr_cat, true_values = mag_values_corr, is_corr = FALSE)
  plot( rowSums(mcmc_corr_mod_corr_cat$b > 0), type = "l", col = 2, xlab = "Index", ylab = "n_trig")
  abline(h = n_trig_corr,lwd = 2)
dev.off()

## 3.0: Comparing models -------------------------------------------------------

BURNIN <- 500
etas_expressions <- c(expression(mu), "C", "a", expression(nu[t]), expression(xi[t]))
magnitude_expressions = c(expression(nu[m0]), expression(xi[m0]), expression(nu[m1]), expression(xi[m1]), expression(rho))
source("./aux_functions/ETAS_RMSE.R")

## 3.1:  ETAS parameter recovery (dual catalogue) ------------------------------

pdf(file = "./output/plots/ETAS_posteriors/dual_catalogue_ETAS_posteriors.pdf",
    width = 4, height = 4)
par(mfrow = c(1,1), mar = c(4.5,4.5,2.1,2.1))
for(i in 1:5){
  d_dual <- density(mcmc_dual_mod_dual_cat$etas[[i]][-(1:BURNIN)])
  d_corr <- density(mcmc_corr_mod_dual_cat$etas[[i]][-(1:BURNIN)])
  ymax = max(c(d_dual$y, d_corr$y)) * 1.05
  plot(d_corr,
       col = "purple",
       lwd = 2,
       ylim = c(0,ymax),
       main = "",
       ylab = "Posterior density",
       xlab = etas_expressions[i],
       cex.lab = 1.25)
  lines(d_dual, col = "blue", lwd = 2)
  abline(v = etas_values_dual[i], lwd = 2)
}
dev.off()

## 3.2:  ETAS parameter recovery (corr catalogue) ------------------------------

pdf(file = "./output/plots/ETAS_posteriors/corr_catalogue_ETAS_posteriors.pdf",
    width = 4, height = 4)
par(mfrow = c(1,1), mar = c(4.5,4.5,2.1,2.1))
for(i in 1:5){
  d_dual <- density(mcmc_dual_mod_corr_cat$etas[[i]][-(1:BURNIN)])
  d_corr <- density(mcmc_corr_mod_corr_cat$etas[[i]][-(1:BURNIN)])
  ymax = max(c(d_dual$y, d_corr$y)) * 1.05
  plot(d_corr,
       col = "purple",
       lwd = 2,
       ylim = c(0,ymax),
       main = "",
       ylab = "Posterior density",
       xlab = etas_expressions[i],
       cex.lab = 1.25)
  lines(d_dual, col = "blue", lwd = 2)
  abline(v = etas_values_corr[i], lwd = 2)
}
dev.off()

## 3.3:  RMSE of ETAS parameters  ----------------------------------------------

# RMSEs for dual catalogue
ETAS_RMSE(mcmc_dual_mod_dual_cat, etas_values_dual, burnin = 500)
ETAS_RMSE(mcmc_corr_mod_dual_cat, etas_values_dual, burnin = 500)
# RRMSE
ETAS_RMSE(mcmc_corr_mod_dual_cat, etas_values_dual, burnin = 500) / ETAS_RMSE(mcmc_dual_mod_dual_cat, etas_values_dual, burnin = 500)

# RMSEs for correlated catalogue
ETAS_RMSE(mcmc_dual_mod_corr_cat, etas_values_corr, burnin = 500)
ETAS_RMSE(mcmc_corr_mod_corr_cat, etas_values_corr, burnin = 500)
# RRMSE
ETAS_RMSE(mcmc_corr_mod_corr_cat, etas_values_corr, burnin = 500) / ETAS_RMSE(mcmc_dual_mod_corr_cat, etas_values_corr, burnin = 500)

# Write to file for easy reference
fileConn<-file("./output/data/ETAS_RMSEs.txt")
writeLines(
  text = c(
    "RMSE dual catalogue dual model",
    as.character(ETAS_RMSE(mcmc_dual_mod_dual_cat, etas_values_dual, burnin = 500)),
    "RMSE dual catalogue corr model",
    as.character(ETAS_RMSE(mcmc_corr_mod_dual_cat, etas_values_dual, burnin = 500)),
    "RMSE corr catalogue dual model",
    as.character(ETAS_RMSE(mcmc_dual_mod_dual_cat, etas_values_corr, burnin = 500)),
    "RMSE corr catalogue corr model",
    as.character(ETAS_RMSE(mcmc_corr_mod_dual_cat, etas_values_corr, burnin = 500))
  ),
  sep = "\n",
  con = fileConn)
close(fileConn)

## 3.4: Magnitude parameter recovery (dual catalogue) --------------------------

pdf(file = "./output/plots/mag_posteriors/dual_catalogue_mag_posteriors.pdf",
    width = 4, height = 4)
par(mfrow = c(1,1), mar = c(4.5,4.5,2.1,2.1))
for(i in 1:5){
  if(i < 5){d_dual <- density(mcmc_dual_mod_dual_cat$mag[[i]][-(1:BURNIN)])}
  d_corr <- density(mcmc_corr_mod_dual_cat$mag[[i]][-(1:BURNIN)])
  ymax = max(c(d_dual$y, d_corr$y)) * 1.05
  plot(d_corr,
       col = "purple",
       lwd = 2,
       ylim = c(0,ymax),
       main = "",
       ylab = "Posterior density",
       xlab = magnitude_expressions[i],
       cex.lab = 1.25)
  if(i < 5){lines(d_dual, col = "blue", lwd = 2)}
  abline(v = mag_values_dual[i], lwd = 2)
  #if(i == 5){ abline(v = 0, col = "blue", lty = 2, lwd = 2)}
}
dev.off()


## 3.5: Magnitude parameter recovery (corr catalogue) -------------------------

pdf(file = "./output/plots/mag_posteriors/corr_catalogue_mag_posteriors.pdf",
    width = 4, height = 4)
par(mfrow = c(1,1), mar = c(4.5,4.5,2.1,2.1))
for(i in 1:4){
  d_dual <- density(mcmc_dual_mod_corr_cat$mag[[i]][-(1:BURNIN)])
  d_corr <- density(mcmc_corr_mod_corr_cat$mag[[i]][-(1:BURNIN)])
  ymax = max(c(d_dual$y, d_corr$y)) * 1.05
  plot(d_corr,
       col = "purple",
       lwd = 2,
       ylim = c(0,ymax),
       main = "",
       ylab = "Posterior density",
       xlab = magnitude_expressions[i],
       cex.lab = 1.25)
  if(i < 5){lines(d_dual, col = "blue", lwd = 2)}
  abline(v = mag_values_corr[i], lwd = 2)
}

for(i in c(5)){
  d_corr <- density(mcmc_corr_mod_corr_cat$mag[[i]][-(1:BURNIN)], from = -1, to = 1)
  ymax = max(c(d_dual$y, d_corr$y)) * 1.05
  plot(d_corr,
       col = "purple",
       lwd = 2,
       ylim = c(0,ymax),
       xlim = c(-1,1),
       main = "",
       ylab = "Posterior density",
       xlab = magnitude_expressions[i],
       cex.lab = 1.25)
  abline(v = mag_values_corr[i], lwd = 2)
  abline(v = 0, lwd = 2, col = "blue")
  #if(i == 5){ abline(v = 0, col = "blue", lty = 2, lwd = 2)}
}
dev.off()

## 4.0: Branching vector recovery ---------------------------------------------

## 4.1: Proportion of B correct at each iteration (dual catalogue) -------------

prop_dual <- apply(mcmc_dual_mod_dual_cat$b[-(1:BURNIN),], 1, FUN =  function(x){mean(x==cat_dual$b)})
prop_corr <- apply(mcmc_corr_mod_dual_cat$b[-(1:BURNIN),], 1, FUN =  function(x){mean(x==cat_dual$b)})

pdf(file = "./output/plots/B_posteriors/dual_catalogue_post_prop_B_correct.pdf",
    width = 7, height = 5)
plot(density(prop_dual, from = 0.95, to = 0.99),
     col = "blue",
     lwd = 2,
     main = "",
     xlab = "Proportion of B correct",
     ylab = "Posterior density",
     cex.lab = 1.25)
lines(density(prop_corr, from = 0.95, to = 0.99),
      col = "purple",
      lwd = 2)
dev.off()

c(mean(prop_dual), mean(prop_corr))

## 4.2: Proportion of B correct at each iteration (corr catalogue) -------------

prop_dual <- apply(mcmc_dual_mod_corr_cat$b[-(1:BURNIN),], 1, FUN =  function(x){mean(x==cat_corr$b)})
prop_corr <- apply(mcmc_corr_mod_corr_cat$b[-(1:BURNIN),], 1, FUN =  function(x){mean(x==cat_corr$b)})

pdf(file = "./output/plots/B_posteriors/corr_catalogue_post_prop_B_correct.pdf",
    width = 7, height = 5)
plot(density(prop_dual, from = 0.95, to = 0.99),
     col = "blue",
     lwd = 2,
     main = "",
     xlab = "Proportion of B correct",
     ylab = "Posterior density",
     cex.lab = 1.25)
lines(density(prop_corr, from = 0.95, to = 0.99),
      col = "purple",
      lwd = 2)
dev.off()

mean(prop_dual) ; mean(prop_corr)

## 4.3: Estimate the increase in Expected proportion of B correct by using corr model ----

n_sim <- 10000
corr_mean_is_greater <- rep(NA_real_, n_sim)
for(i in 1:n_sim){
  boot_prop_dual <- sample(prop_dual, length(prop_dual),replace = TRUE)
  boot_prop_corr <- sample(prop_corr, length(prop_corr),replace = TRUE)
  corr_mean_is_greater[i] <- mean(boot_prop_corr) - mean(boot_prop_dual)
}
plot(density(corr_mean_is_greater))
mean(corr_mean_is_greater)
quantile(corr_mean_is_greater, c(0.025, 0.975))
# Posterior mean of the Proporion  of B correct is increased by 0.43%  when using the corr model.
# 95% Credible interval for change is (0.41%,0.45%).

## EOF -------------------------------------------------------------------------




