# Diagnosing issues in correlated magnitude code. ------------------------------

## 0.0: Simulate some catalogues to experiment with ----------------------------

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

# Simulate a dual magntude catalogue ---------
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

# Simulate a correlated magnitude catalogue ---------
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

# Calculate useful summaries and statistics --------
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


## 01: What does using the dual mangitude code do? -----------------------------

# 01.1: Applying dual magntiude MCMC to dual catalgoue -------------------------
# set.seed(1234)
# mcmc_dual_mod_dual_cat <- corrETAS::sampleETASposterior_dualmag(
#   ts = cat_dual$ts,
#   ms = cat_dual$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = c(0.1,0.1,0.1,0.5,0),
#   init_B = (seq_along(cat_dual$ts)-1) * rbinom(n = length(cat_dual$ts), size = 1, prob = 0.4),
#   init_mag = c(1,0,0.5,0),
#   sims = 5000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05),
#   mu_prior = c(0.001,0.001)
# )
# saveRDS(mcmc_dual_mod_dual_cat, file = "./output/diagnosing/mcmc_dual_mod_dual_cat.RDS")
mcmc_dual_mod_dual_cat <- readRDS("./output/diagnosing/mcmc_dual_mod_dual_cat.RDS")

# 01.2: Applying dual magnitude MCMC to corr catalogue -------------------------
# set.seed(1234)
# mcmc_dual_mod_corr_cat <- corrETAS::sampleETASposterior_dualmag(
#   ts = cat_corr$ts,
#   ms = cat_corr$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = c(0.1,0.1,0.1,0.5,0),
#   init_B = (seq_along(cat_corr$ts)-1) * rbinom(n = length(cat_corr$ts), size = 1, prob = 0.4),
#   init_mag = c(1,0,0.5,0),
#   sims = 5000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05),
#   mu_prior = c(0.001,0.001)
# )
# saveRDS(mcmc_dual_mod_corr_cat, file = "./output/diagnosing/mcmc_dual_mod_corr_cat.RDS")
mcmc_dual_mod_corr_cat <- readRDS("./output/diagnosing/mcmc_dual_mod_corr_cat.RDS")

# 01.3 Inspect chains ---------------------------------------------------------- (COME BACK TO THIS)

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



plot_ETAS_chains(mcmc_dual_mod_dual_cat, true_values = etas_values_dual)
plot_ETAS_chains(mcmc_dual_mod_corr_cat, true_values = etas_values_corr)

plot_mag_chains(mcmc_dual_mod_dual_cat, true_values = mag_values_dual)
plot_mag_chains(mcmc_dual_mod_corr_cat, true_values = mag_values_corr)

for(i in 1:4){
  plot(density(mcmc_dual_mod_dual_cat$mag[[i]]),
       main = "dual magitude model",
       xlab = paste0(mag_names_dual[i]))
  lines(density(mcmc_dual_mod_corr_cat$mag[[i]]), col = 2)
  abline(v = mag_values_corr[i])
  legend("topright", legend = c("dual data", "corr data"), col = c(1,2), lwd = c(2,2), bty = "n")
}

## 02: Correlated magnitude fitting --------------------------------------------

## When rho is fixed to its true value can the other parameters be recovered?
# set.seed(1234)
# mcmc_2_1 <- corrETAS::sampleETASposterior_corrmag(
#   ts = cat_corr$ts,
#   ms = cat_corr$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = c(0.1,0.1,0.1,0.5,0),
#   init_B = (seq_along(cat_corr$ts)-1) * rbinom(n = length(cat_corr$ts), size = 1, prob = 0.4),
#   init_mag = c(1,0,0.4,0,0.5),
#   sims = 1000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.0), # rho is fixed to 0.5
#   mu_prior = c(0.001,0.001)
# )

# overestimating the number of triggered events. Try starting at truth, see if it wanders away.

# mcmc_2_2 <- corrETAS::sampleETASposterior_corrmag(
#   ts = cat_corr$ts,
#   ms = cat_corr$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = etas_values_corr,
#   init_B = cat_corr$b,
#   init_mag = mag_values_corr,
#   sims = 100,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05,0.0), # rho is fixed to 0.5
#   mu_prior = c(0.001,0.001)
# )
#
# plot_ETAS_chains(mcmc_2_2, true_values = etas_values_corr)
# plot_mag_chains(mcmc_2_2, true_values = mag_values_corr,is_corr = TRUE)


# mcmc_2_3 <- corrETAS::sampleETASposterior_corrmag(
#   ts = cat_corr$ts,
#   ms = cat_corr$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = etas_values_corr,
#   init_B = cat_corr$b,
#   init_mag = mag_values_corr,
#   sims = 1000,
#   B_samples = TRUE,
#   B_fixed = TRUE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05,0.0), # rho is fixed to 0.5
#   mu_prior = c(0.0001,0.0001),
#   logC_prior = c(log(0.2), 0.8)
# )
# plot_ETAS_chains(mcmc_2_3, true_values = etas_values_corr)
# plot_mag_chains(mcmc_2_3, true_values = mag_values_corr,is_corr = TRUE)
#
#
# # What about when corr is 0?
# mcmc_2_4 <- corrETAS::sampleETASposterior_corrmag(
#   ts = cat_dual$ts,
#   ms = cat_dual$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = etas_values_dual,
#   init_B = cat_dual$b,
#   init_mag = c(mag_values_dual,0),
#   sims = 1000,
#   B_samples = TRUE,
#   B_fixed = TRUE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05,0.0), # rho is fixed to 0.0
#   mu_prior = c(0.0001,0.0001),
#   logC_prior = c(log(0.2), 0.8)
# )
# plot_ETAS_chains(mcmc_2_4, true_values = etas_values_corr)
# plot_mag_chains(mcmc_2_4, true_values = mag_values_corr,is_corr = TRUE)
#
# saveRDS(mcmc_dual_mod_corr_cat, file = "./output/diagnosing/mcmc_dual_mod_corr_cat.RDS")
# mcmc_dual_mod_corr_cat <- readRDS("./output/diagnosing/mcmc_dual_mod_corr_cat.RDS")
#
# # When corr is 0 and B is known can I recover rho?
# mcmc_2_5 <- corrETAS::sampleETASposterior_corrmag(
#   ts = cat_dual$ts,
#   ms = cat_dual$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = etas_values_dual,
#   init_B = cat_dual$b,
#   init_mag = c(mag_values_dual,0.2),
#   sims = 1000,
#   B_samples = TRUE,
#   B_fixed = TRUE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05,0.05), # rho is fixed to 0.0
#   mu_prior = c(0.0001,0.0001),
#   logC_prior = c(log(0.2), 0.8)
# )
# plot_ETAS_chains(mcmc_2_5, true_values = etas_values_corr)
# plot_mag_chains(mcmc_2_5, true_values = mag_values_corr,is_corr = TRUE)
#
# # When corr is 0 and I start at the truth can I recover everything?
# mcmc_2_6 <- corrETAS::sampleETASposterior_corrmag(
#   ts = cat_dual$ts,
#   ms = cat_dual$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = etas_values_dual,
#   init_B = cat_dual$b,
#   init_mag = c(mag_values_dual,0.0),
#   sims = 1000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05,0.05), # rho is fixed to 0.0
#   mu_prior = c(0.0001,0.0001),
#   logC_prior = c(log(0.2), 0.8)
# )
# #saveRDS(mcmc_2_6, file = "./output/diagnosing/mcmc_2_6.RDS")
#
# plot_ETAS_chains(mcmc_2_6, true_values = etas_values_corr)
# plot_mag_chains(mcmc_2_6, true_values = mag_values_corr,is_corr = TRUE)
#
# # When corr is 0 and I start at a sensible-ish first guess can I recover everything?
# init_prop_trig <- 0.4
# init_mu <- (1 - init_prop_trig) * length(cat_dual$b) / 50000
#
# set.seed(2212)
#
# init_etas <- c(init_mu, 0.4, 0.01, 50, 0)
# init_mag <- c(1,0,0.5,0,0)
# init_b <- (seq_along(cat_dual$b) - 1) * rbinom(length(cat_dual$b), size = 1, prob = 1 - init_prop_trig)
#
# plot(x = seq(0,10,length.out = 1001), y = dlnorm(x = seq(0,10,length.out = 1001),meanlog = log(0.2),sdlog = 1000))
#
# mcmc_2_7 <- corrETAS::sampleETASposterior_corrmag(
#   ts = cat_dual$ts,
#   ms = cat_dual$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = init_etas,
#   init_B = init_b,
#   init_mag = init_mag,
#   sims = 1000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05,0.05),
#   mu_prior = c(0.0001,0.0001),
#   logC_prior = c(log(0.2), 10)
# )
# #saveRDS(mcmc_2_7, file = "./output/diagnosing/mcmc_2_7.RDS")
#
# #pdf(file = "./output/diagnosing/mcmc_2_7_plots/mcmc_2_7_traceplots.pdf", width = 7, height = 5)
# plot_ETAS_chains(mcmc_2_7, true_values = etas_values_dual)
# plot_mag_chains(mcmc_2_7, true_values = c(mag_values_dual,0),is_corr = TRUE)
# plot( rowSums(mcmc_2_7$b > 0), type = "l", col = 2, xlab = "Index", ylab = "n_trig")
# abline(h = n_trig_dual,lwd = 2)
# #dev.off()
#
# # Try the same but for the correlated catalogue now.
#
# init_prop_trig <- 0.4
# init_mu <- (1 - init_prop_trig) * length(cat_corr$b) / 50000
#
# set.seed(2212)
#
# init_etas <- c(init_mu, 0.4, 0.01, 50, 0)
# init_mag <- c(1,0,0.5,0,0)
# init_b <- (seq_along(cat_corr$b) - 1) * rbinom(length(cat_corr$b), size = 1, prob = 1 - init_prop_trig)
#
# plot(x = seq(0,10,length.out = 1001), y = dlnorm(x = seq(0,10,length.out = 1001),meanlog = log(0.2),sdlog = 1000))
#
# mcmc_2_8 <- corrETAS::sampleETASposterior_corrmag(
#   ts = cat_corr$ts,
#   ms = cat_corr$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = init_etas,
#   init_B = init_b,
#   init_mag = init_mag,
#   sims = 1000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   etas_sds = c(0.1,0.1,0.1,0.1),
#   mag_sds = c(0.05,0.05,0.05,0.05,0.05),
#   mu_prior = c(0.0001,0.0001),
#   logC_prior = c(log(0.2), 10)
# )
# #saveRDS(mcmc_2_8, file = "./output/diagnosing/mcmc_2_8.RDS")
#
# pdf(file = "./output/diagnosing/mcmc_2_8_plots/mcmc_2_8_traceplots.pdf", width = 7, height = 5)
# plot_ETAS_chains(mcmc_2_8, true_values = etas_values_corr)
# plot_mag_chains(mcmc_2_8, true_values = c(mag_values_corr),is_corr = TRUE)
# plot( rowSums(mcmc_2_8$b > 0), type = "l", col = 2, xlab = "Index", ylab = "n_trig")
# abline(h = n_trig_corr,lwd = 2)
# dev.off()


# DO LONG RUNS
# init_prop_trig <- 0.4
# init_mu <- (1 - init_prop_trig) * length(cat_dual$b) / 50000
#
# set.seed(1993)
#
# init_etas <- c(init_mu, 0.4, 0.01, 50, 0)
# init_mag <- c(1,0,0.5,0,0)
# init_b <- (seq_along(cat_dual$b) - 1) * rbinom(length(cat_dual$b), size = 1, prob = 1 - init_prop_trig)
#
# plot(x = seq(0,10,length.out = 1001), y = dlnorm(x = seq(0,10,length.out = 1001),meanlog = log(0.2),sdlog = 1000))
#
# mcmc_2_9 <- corrETAS::sampleETASposterior_corrmag(
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
# saveRDS(mcmc_2_9, file = "./output/diagnosing/mcmc_2_9.RDS")
mcmc_corr_mod_dual_cat <- readRDS("./output/diagnosing/mcmc_2_9.RDS")

# init_prop_trig <- 0.4
# init_mu <- (1 - init_prop_trig) * length(cat_dual$b) / 50000
#
# set.seed(1993)
#
# init_etas <- c(init_mu, 0.4, 0.01, 50, 0)
# init_mag <- c(1,0,0.5,0,0)
# init_b <- (seq_along(cat_corr$b) - 1) * rbinom(length(cat_corr$b), size = 1, prob = 1 - init_prop_trig)
#
# mcmc_2_10 <- corrETAS::sampleETASposterior_corrmag(
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
# saveRDS(mcmc_2_10, file = "./output/diagnosing/mcmc_2_10.RDS")

# pdf(file = "./output/diagnosing/mcmc_2_10_plots/mcmc_2_10_traceplots.pdf", width = 7, height = 5)
# plot_ETAS_chains(mcmc_2_10, true_values = etas_values_corr)
# plot_mag_chains(mcmc_2_10, true_values = c(mag_values_corr),is_corr = TRUE)
# plot( rowSums(mcmc_2_10$b > 0), type = "l", col = 2, xlab = "Index", ylab = "n_trig")
# abline(h = n_trig_corr,lwd = 2)
# dev.off()

mcmc_corr_mod_corr_cat <- readRDS("./output/diagnosing/mcmc_2_10.RDS")


## 03: Comparing models --------------------------------------------------------

## 03.1: ETAS parameter recovery -----------------------------------------------
BURNIN <- 500
etas_expressions <- c(expression(mu), "C", "a", expression(nu[t]), expression(xi[t]))

ETAS_RMSE <- function(mcmc, true_etas_values, burnin){
  squared_errors <- apply(mcmc$etas[-(1:burnin),], 1, function(x){(x - true_etas_values)^2})
  rmse <- apply(squared_errors, 1, function(x){sqrt(mean(x))})
  return(rmse)
}

# 0.3.11: Dual magnitude catalogue ---------------------------------------------

pdf(file = "./output/diagnosing/comparisons/dual_catalogue_ETAS_posteriors.pdf",
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

# Correlated magnitude catalogue --
pdf(file = "./output/diagnosing/comparisons/corr_catalogue_ETAS_posteriors.pdf",
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
       main = "corr cat",
       ylab = "Posterior density",
       xlab = etas_expressions[i],
       cex.lab = 1.25)
  lines(d_dual, col = "blue", lwd = 2)
  abline(v = etas_values_corr[i], lwd = 2)
}
dev.off()

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

fileConn<-file("./output/diagnosing/comparisons/ETAS_RMSEs.txt")
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

## 03.1: Magnitude parameter recovery ------------------------------------------
magnitude_expressions = c(
  expression(nu[m0]),
  expression(xi[m0]),
  expression(nu[m1]),
  expression(xi[m1]),
  expression(rho))

## Dual catalogue -----
pdf(file = "./output/diagnosing/comparisons/dual_catalogue_mag_posteriors.pdf",
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

## Correlated catalogue -----
pdf(file = "./output/diagnosing/comparisons/corr_catalogue_mag_posteriors.pdf",
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


## 03.3: Branching vector recovery ---------------------------------------------

# proportion of branching vector correct at each iteration
# dual magnitude catalogue -
prop_dual <- apply(mcmc_dual_mod_dual_cat$b[-(1:BURNIN),], 1, FUN =  function(x){mean(x==cat_dual$b)})
prop_corr <- apply(mcmc_corr_mod_dual_cat$b[-(1:BURNIN),], 1, FUN =  function(x){mean(x==cat_dual$b)})
plot(prop_dual, type = "l")
lines(prop_corr, col = "purple")

plot(density(prop_dual, from = 0.95, to = 0.99), col = "blue", lwd = 2, main = "", xlab = "Proportion of B correct", ylab = "Posterior density", cex.lab = 1.25)
lines(density(prop_corr, from = 0.95, to = 0.99), col = "purple", lwd = 2)
mean(prop_dual)
mean(prop_corr)

# correlated catalogue -
prop_dual <- apply(mcmc_dual_mod_corr_cat$b[-(1:BURNIN),], 1, FUN =  function(x){mean(x==cat_corr$b)})
prop_corr <- apply(mcmc_corr_mod_corr_cat$b[-(1:BURNIN),], 1, FUN =  function(x){mean(x==cat_corr$b)})
plot(prop_dual, type = "l")
lines(prop_corr, col = "purple")

plot(density(prop_dual, from = 0.95, to = 0.99), col = "blue", lwd = 2, main  = "", xlab = "Proportion of B correct", ylab = "Posterior density", cex.lab= 1.25)
lines(density(prop_corr, from = 0.95, to = 0.99), col = "purple", lwd = 2)
mean(prop_dual)
mean(prop_corr)

## Bootstrap probability that E prop using corr model is greater than E prop using dual model
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
    # Posterior mean of the Proporion  of B correct is increased by 0.42% (0.40%,0.44%) when using the corr model.

# probability that a randomly chosen element is allocated correctly
prob_dual <- rowMeans(apply(mcmc_dual_mod_dual_cat$b, 1, FUN = function(x){x==cat_dual$b}))
prob_corr <- rowMeans(apply(mcmc_corr_mod_dual_cat$b, 1, FUN = function(x){x==cat_dual$b}))
plot(density(prob_corr - prob_dual))
mean(prob_corr - prob_dual)

plot(density(c(prob_dual, 1 + (1 - prob_dual)), from = 0, to = 1, bw = 0.001), main = "dual cat", lwd = 2)
lines(density(c(prob_corr, 1 + (1 - prob_corr)), from = 0.9, to = 1,bw = 0.001), col = "purple", lwd = 2)


prob_dual <- rowMeans(apply(mcmc_dual_mod_corr_cat$b, 1, FUN = function(x){x==cat_corr$b}))
prob_corr <- rowMeans(apply(mcmc_corr_mod_corr_cat$b, 1, FUN = function(x){x==cat_corr$b}))
plot(density(prob_corr - prob_dual))
mean(prob_corr - prob_dual)

plot(density(c(prob_dual, 1 + (1 - prob_dual)), from = 0.9, to = 1, bw = 0.001), main = "corr cat", lwd = 2)
lines(density(c(prob_corr, 1 + (1 - prob_corr)), from = 0.9, to = 1,bw = 0.001), col = "purple", lwd = 2)







