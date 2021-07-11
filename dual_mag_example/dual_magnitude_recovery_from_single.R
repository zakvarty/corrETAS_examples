## 00: Set up ------------------------------------------------------------------
source("./aux_functions.R")

# Set true parameter values for simulation and plotting ------------------------
mu_value <- 0.02
K_value <- 0.2
alpha_value <- 0.1
nu_t_value <- 0.1
xi_t_value <- 0

m_0_value <- 1.5

mag_values <- c(0.4 , 0, 0.4, 0)
mag_names <- c('nu_m0', 'xi_m0', 'nu_m1', 'xi_m1')

# Simulate catalgoue -----------------------------------------------------------
set.seed(2212)
cat <-  corrETAS::simulateETAScorr(
  mu = mu_value,
  K = K_value,
  alpha = alpha_value,
  gpd_t = c(nu_t_value, xi_t_value),
  rho = 0,
  gpd_m0 = mag_values[1:2],
  gpd_m1 = mag_values[3:4],
  m_0 = m_0_value,
  t_max = 50000,
  displayOutput = TRUE
)

colouring <- ifelse(cat$b == 0, "grey", "darkorange")
plot(x = cat$ts,y = cat$ms, col = colouring, type = "h")

C_value <- K_value * exp( alpha_value * (mean(cat$ms) - m_0_value))
etas_values <- c(mu_value, C_value, alpha_value, nu_t_value, xi_t_value)
etas_names <- c('mu', 'C', 'alpha', 'nu_t', 'xi_t')

n_trig <- sum(cat$b > 0)
m_bar <- mean(cat$ms)


# 01: Posterior samples 1 - Branching vector known -----------------------------

# sample joint posterior -------------------------------------------------------

# set.seed(2212)
# mcmc_1 <- corrETAS::sampleETASposterior_dualmag(
#   ts = cat$ts,
#   ms = cat$ms,
#   m_0 = 1.5,
#   t_max = 50000,
#   init_ETAS = c(0.01,0.1, 0.1, 0.1, 0.1),
#   init_mag = c(0.5, 0, 0.01, 0),
#   sims = 5000,
#   B_samples = TRUE,
#   B_fixed = TRUE,
#   init_B =  cat$b,
#   etas_sds = rep(0.1, 4),
#   mag_sds = rep(0.01, 4),
#   mu_prior = c(0.1, 0.1)
# )
# saveRDS(object = mcmc_1, file = "./output/mcmc/single_magnitude_case/mcmc_1.RDS")
mcmc_1 <- readRDS("./output/mcmc/single_magnitude_case/mcmc_1.RDS")

# calculate magnitude MLEs and transform onto nuxi space -----------------------
mag_mles <- c(ismev::gpd.fit(xdat = cat$ms[cat$b == 0],threshold = 1.5)$mle,
              ismev::gpd.fit(xdat = cat$ms[cat$b > 0],threshold = 1.5)$mle)
mag_mles[1:2] <- corrETAS::sig2nu(mag_mles[1], mag_mles[2])
mag_mles[3:4] <- corrETAS::sig2nu(mag_mles[3], mag_mles[4])

# parameter trace plots --------------------------------------------------------
xrange = 1:5000
pdf(
  file = "./output/mcmc/single_magnitude_case/trace_plots/trace_plots_1.pdf",
  width = 7,
  height = 5)
opar <- par()
par(mar = c(4.1,4.1,2.1,1.1))
for(i in 1:5){
  plot(mcmc_1$etas[xrange,i],
       ylab = etas_names[i],
       type= 'l',
       bty = 'n')
  abline(h = quantile(x = mcmc_1$etas[xrange,i],probs = c(0.025,0.975)), col = 3, lwd = 2)
  abline(h = etas_values[i],col = 2, lwd = 2)
}

for(i in 1:4){
  plot(mcmc_1$mag[xrange,i], ylab = mag_names[i],type= 'l', bty = 'n')
  abline(h = quantile(x = mcmc_1$mag[xrange,i],probs = c(0.025,0.975)), col = 3, lwd = 2)
  abline(h = mag_values[i],col = "red", lwd = 2)
  abline(h = mag_mles[i], col = "red", lty = 2, lwd = 2)
}

par(mar = opar$mar)


# 02: Posterior samples 2 - B unknown, starting at true value -------------------

# sample joint posterior -------------------------------------------------------

# set.seed(2212)
# mcmc_2 <- corrETAS::sampleETASposterior_dualmag(
#   ts = cat$ts,
#   ms = cat$ms,
#   m_0 = 1.5,
#   t_max = 50000,
#   init_ETAS = c(0.01,0.1, 0.1, 0.1, 0.1),
#   init_mag = c(0.5, 0, 0.01, 0),
#   sims = 5000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   init_B =  cat$b,
#   etas_sds = rep(0.1, 4),
#   mag_sds = rep(0.01, 4),
#   mu_prior = c(0.1, 0.1)
# )
# saveRDS(object = mcmc_2, file = "./output/mcmc/single_magnitude_case/mcmc_2.RDS")
mcmc_2 <- readRDS("./output/mcmc/single_magnitude_case/mcmc_2.RDS")


# parameter trace plots --------------------------------------------------------
pdf(
  file = "./output/mcmc/single_magnitude_case/trace_plots/trace_plots_2.pdf",
  width = 7,
  height = 5)
opar <- par()
par(mar = c(4.1,4.1,2.1,1.1))
for(i in 1:5){
  plot(mcmc_2$etas[xrange,i],
       ylab = etas_names[i],
       type= 'l',
       bty = 'n')
  abline(h = quantile(x = mcmc_2$etas[xrange,i],probs = c(0.025,0.975)), col = 3, lwd = 2)
  abline(h = etas_values[i],col = 2, lwd = 2)
}
xrange <- 1:5000
for(i in 1:4){
  plot(mcmc_2$mag[xrange,i], ylab = mag_names[i],type= 'l', bty = 'n')
  abline(h = quantile(x = mcmc_2$mag[xrange,i],probs = c(0.025,0.975)), col = 3, lwd = 2)
  abline(h = mag_values[i],col = "red", lwd = 2)
  abline(h= mag_mles[i], col = "red", lty = 2, lwd = 2)
}
dev.off()
par(mar = opar$mar)

# 03: Posterior samples 3 -B unknown, random initialisation --------------------

# sample joint posterior ----
# set.seed(2212)
# branching_init <- (seq_along(cat$b) - 1) * rbinom(length(cat$b), size = 1, prob = 0.4)
# mcmc_3 <- corrETAS::sampleETASposterior_dualmag(
#   ts = cat$ts,
#   ms = cat$ms,
#   m_0 = 1.5,
#   t_max = 50000,
#   init_ETAS = c(0.01,0.1, 0.1, 0.1, 0.1),
#   init_mag = c(0.5, 0, 0.01, 0),
#   sims = 5000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   init_B =  branching_init,
#   etas_sds = rep(0.1, 4),
#   mag_sds = rep(0.01, 4),
#   mu_prior = c(0.1, 0.1)
# )
# saveRDS(object = mcmc_3, file = "./output/mcmc/single_magnitude_case/mcmc_3.RDS")
mcmc_3 <- readRDS("./output/mcmc/single_magnitude_case/mcmc_3.RDS")

# parameter trace plots --------------------------------------------------------
pdf(
  file = "./output/mcmc/trace_plots/trace_plots_3.pdf",
  width = 7,
  height = 5)
opar <- par()
par(mar = c(4.1,4.1,2.1,1.1))
for(i in 1:5){
  plot(mcmc_3$etas[xrange,i],
       ylab = etas_names[i],
       type= 'l',
       bty = 'n')
  abline(h = quantile(x = mcmc_3$etas[xrange,i],probs = c(0.025,0.975)), col = 3, lwd = 2)
  abline(h = etas_values[i],col = 2, lwd = 2)
}
xrange <- 1:5000
for(i in 1:4){
  plot(mcmc_3$mag[xrange,i], ylab = mag_names[i],type= 'l', bty = 'n')
  abline(h = quantile(x = mcmc_3$mag[xrange,i],probs = c(0.025,0.975)), col = 3, lwd = 2)
  abline(h = mag_values[i],col = "red", lwd = 2)
  #abline(h= mag_mles[i], col = "red", lty = 2, lwd = 2)
}
dev.off()
par(mar = opar$mar)


## 04: Posterior samples 4 - Correctly assume single magnitude distribtion -------
# set.seed(2212)
# branching_init <- (seq_along(cat$b) - 1) * rbinom(length(cat$b), size = 1, prob = 0.4)
# mcmc_4 <- corrETAS::sampleETASposterior(
#   ts = cat$ts,
#   ms = cat$ms,
#   m_0 = m_0_value,
#   t_max = 50000,
#   init_ETAS = c(0.01,0.1, 0.1, 0.1, 0.1),
#   init_mag = c(0.3,0.1),
#   sims = 5000,
#   B_samples = TRUE,
#   B_fixed = FALSE,
#   init_B = branching_init,
#   etas_sds = rep(01.,4),
#   mag_sds = rep(0.01,2),
#   mu_prior = c(0.1,0.1))
# saveRDS(object = mcmc_4, file = "./output/mcmc/single_magnitude_case/mcmc_4.RDS")
mcmc_4 <- readRDS("./output/mcmc/single_magnitude_case/mcmc_4.RDS")

# parameter trace plots --------------------------------------------------------
pdf(
  file = "./output/mcmc/single_magnitude_case/trace_plots/trace_plots_4.pdf",
  width = 7,
  height = 5)
opar <- par()
par(mar = c(4.1,4.1,2.1,1.1))
for(i in 1:5){
  plot(mcmc_4$etas[xrange,i],
       ylab = etas_names[i],
       type= 'l',
       bty = 'n')
  abline(h = quantile(x = mcmc_4$etas[xrange,i],probs = c(0.025,0.975)), col = 3, lwd = 2)
  abline(h = etas_values[i],col = 2, lwd = 2)
}
xrange <- 1:5000
single_mag_names <- c(expression(nu[m]), expression(xi[m]))
for(i in 1:2){
  plot(mcmc_4$mag[xrange,i], ylab = single_mag_names[i],type= 'l', bty = 'n')
  abline(h = quantile(x = mcmc_4$mag[xrange,i],probs = c(0.025,0.975)), col = 3, lwd = 2)
  abline(h = mag_values[i],col = "red", lwd = 2)
  abline(h= mag_mles[i], col = "red", lty = 2, lwd = 2)
}
dev.off()
par(mar = opar$mar)

## 05: Assessing Magnitude recovery --------------------------------------------

# Marginal magnitude paramter posterirors
pdf(file = "./output/parameter_recovery/single_magnitude_case/magnitude_parameter_posteriors.pdf",
    width = 4, height = 4)
par(mar = c(4.5,4.5,2.1,2.1))

x_low <- 0.25
x_high <- 0.45
plot(density(mcmc_4$mag$nu, from = x_low, to = x_high),
     xlab = expression(nu[m]),
     ylab = "posterior density",
     main = "",
     col = "darkorange",
     lwd = 2,
     cex.lab = 1.5)
lines(density(mcmc_3$mag$nu_m1, from = x_low, to = x_high),
      col = 2,
      lwd = 2)
lines(density(mcmc_3$mag$nu_m0,  from = x_low, to = x_high),
      col = "black",
      lwd = 2)
abline(v = mag_values[1], lty = 2 , lwd = 2, col = "blue")

x_low <- -0.3
x_high <- 0.3
plot(density(mcmc_4$mag$xi, from = x_low, to = x_high),
     xlab = expression(xi[m]),
     ylab = "posterior density",
     main = "",
     col = "darkorange",
     lwd = 2,
     cex.lab = 1.5)
lines(density(mcmc_3$mag$xi_m1, from = x_low, to = x_high),
      col = 2,
      lwd = 2)
lines(density(mcmc_3$mag$xi_m0,  from = x_low, to = x_high),
      col = "black",
      lwd = 2)
abline(v = mag_values[2], lty = 2 , lwd = 2, col = "blue")
# Joint marginal posteriors

contour(
  x = MASS::kde2d(x = mcmc_3$mag$nu_m0, y = mcmc_3$mag$xi_m0),
  levels = seq(0,350, by = 20),
  xlim = c(0.3,0.45),
  ylim = c(-0.2,0.1),
  col = "black",
  xlab = expression(nu[m]),
  ylab = expression(xi[m]),
  cex.lab = 1.5
)

contour(
  x = MASS::kde2d(
    x = mcmc_3$mag$nu_m1, y = mcmc_3$mag$xi_m1),
  levels = seq(0,350,by = 20),
  add= TRUE,
  col = "red")

contour(
  x = MASS::kde2d(
    x = mcmc_4$mag$nu_m, y = mcmc_4$mag$xi_m),
  levels = seq(0,600,by = 20),
  add= TRUE,
  col = "darkorange")
points(x = mag_values[1], y = mag_values[2], cex = 1.5, pch = 16, col = "blue")
dev.off()


# Posteriors of parameter differences
pdf("./output/parameter_recovery/single_magnitude_case/parameter_difference_posteriors.pdf", width = 4, height = 4)
# Background vs triggered
par(mar = c(4.5, 4.5, 2.1,2.1))
emdbook::HPDregionplot(
  coda::as.mcmc(x = cbind(mcmc_3$mag$nu_m0 - mcmc_3$mag$nu_m1, mcmc_3$mag$xi_m0 - mcmc_3$mag$xi_m1)),
  col = 2,
  xlim = c(-0.05,0.15),
  ylim = c(-0.15, 0.25),
  xlab = expression(nu[m0] - nu[m1]),
  ylab = expression(xi[m0] - xi[m1]),
  cex.lab = 1.5
  )
contour(
  x = MASS::kde2d(
    x = mcmc_3$mag$nu_m0 - mcmc_3$mag$nu_m1,
    y =  mcmc_3$mag$xi_m0 - mcmc_3$mag$xi_m1),
  add= TRUE)
points(x = 0 , y = 0, pch = 16, col = "blue", cex = 1.5)
#abline(h = 0)
#abline(v = 0)

# Background vs single
par(mar = c(4.5, 4.5, 2.1,2.1))
emdbook::HPDregionplot(
  coda::as.mcmc(x = cbind(mcmc_3$mag$nu_m0 - mcmc_4$mag$nu_m, mcmc_3$mag$xi_m0 - mcmc_4$mag$xi_m)),
  col = 2,
  xlim = c(-0.02,0.05),
  ylim = c(-0.05, 0.075),
  xlab = expression(nu[m0] - nu[m]),
  ylab = expression(xi[m0] - xi[m]),
  cex.lab = 1.5
)
contour(
  x = MASS::kde2d(
    x = mcmc_3$mag$nu_m0 - mcmc_4$mag$nu_m,
    y =  mcmc_3$mag$xi_m0 - mcmc_4$mag$xi_m),
  add= TRUE)
points(x = 0 , y = 0, pch = 16, col = "blue", cex = 1.5)
#abline(h = 0)
#abline(v = 0)

#Triggered vs single
par(mar = c(4.5, 4.5, 2.1,2.1))
emdbook::HPDregionplot(
  coda::as.mcmc(x = cbind(mcmc_3$mag$nu_m1 - mcmc_4$mag$nu_m, mcmc_3$mag$xi_m1 - mcmc_4$mag$xi_m)),
  col = 2,
  xlim = c(-0.1,0.05),
  ylim = c(-0.25, 0.2),
  xlab = expression(nu[m1] - nu[m]),
  ylab = expression(xi[m1] - xi[m]),
  cex.lab = 1.5
)
contour(
  x = MASS::kde2d(
    x = mcmc_3$mag$nu_m1 - mcmc_4$mag$nu_m,
    y =  mcmc_3$mag$xi_m1 - mcmc_4$mag$xi_m),
  add= TRUE)
points(x = 0 , y = 0, pch = 16, col = "blue", cex = 1.5)
#abline(h = 0)
#abline(v = 0)
dev.off()
# plot magnitude return levels and confidence bounds  --------------------------------------

# Dual magnitude model : per event type ----------------------------------------
background_rls <- get_return_levels(
  MCMC = mcmc_3,
  return_periods = 2^(1:15),
  m_0 = 1.5,
  event_type = 0)

triggered_rls <- get_return_levels(
  MCMC = mcmc_3,
  return_periods = 2^(1:15),
  m_0 = 1.5,
  event_type = 1)

true_return_levels_background <- corrETAS::qgpd(
  p = 1 - 1/2^(1:15),
  mu = m_0_value,
  nu = mag_values[1],
  shape = mag_values[2]
)

true_return_levels_triggered <- corrETAS::qgpd(
  p = 1 - 1/2^(1:15),
  mu = m_0_value,
  nu = mag_values[1],
  shape = mag_values[2]
)

pdf("./output/parameter_recovery/single_magnitude_case/dual_mag_return_level_plots.pdf", width = 4, height = 4)
par(mfrow = c(1,1), mar = c(4.5,4.5,1.1,1.1))
plot_sampled_rls(
  sampled_rls = background_rls,
  alpha = 0.05,
  catalogue = cat,
  event_type = 0,
  #main = "background events",
  ylim = c(1.5,10),
  cex.lab = 1.5)
lines(x = 2^(1:15), y = true_return_levels_background, lty = 3, lwd = 2)

plot_sampled_rls(
  sampled_rls = triggered_rls,
  alpha = 0.05,
  catalogue = cat,
  event_type = 1,
  #main = "triggered events",
  ylim = c(1.5,10),
  cex.lab = 1.5)
lines(x = 2^(1:15), y = true_return_levels_triggered, lty = 3, lwd = 2)
dev.off()

# Dual vs Single magnitude model : magnitude return levels  --------------------

single_mag_mcmc <- corrETAS::gpd_mcmc_nu(
  x = cat$ms,
  threshold = 1.5,
  init = c(1,0),
  step_sds = c(0.05,0.05),
  n_samples = 10000,
  verbose = TRUE)

plot(single_mag_mcmc$nu, type = "l")
plot(single_mag_mcmc$xi, type = "l")
plot(single_mag_mcmc$lpost, type = "l")

single_mag_rls <- get_return_levels_single(
  MAG_MCMC = single_mag_mcmc,
  m_0 = 1.5,
  return_periods = 2^(1:12))

# # Test using true parameters to show that mixture return levels are being
# # calculated correctly.
# true_return_levels <- rep(NA, length(return_periods))
# for (i in seq_along(return_periods)){
#   true_return_levels[i] <- qgpd_mix(
#     p = 1 - 1/return_periods[i],
#     nu_0 = mag_values[1],
#     xi_0 = mag_values[2],
#     nu_1 = mag_values[3],
#     xi_1 = mag_values[4],
#     m_0 = 1.5,
#     prop_trig = sum(cat$b > 1) / length(cat$b),
#     lower = 1.5,
#     upper = 50)
# }
#
# plot(x = return_periods,
#      y = quantile(cat$ms, probs = 1 - 1 / return_periods),
#      log = "xy")
# lines(x = return_periods, y = true_return_levels)

return_periods <- c(2,4,8,16,32,64,128, 256, 512, 1024,2048, 4096, 8196)

dual_return_levels <- data.frame(a = rep(NA, nrow(mcmc_3$b)))
for(i in seq_along(return_periods)){
  dual_return_levels[,i] <- mcmc_return_levels(mcmc_3, return_periods[i], m_0_value, upper = 10)
  print(return_periods[i])
}
names(dual_return_levels) <- paste("rl_",return_periods)

true_return_levels_combined <- purrr::map_dbl(
  .x = 1 - 1/return_periods,
  .f = qgpd_mix,
  nu_0 = mag_values[1],
  xi_0 = mag_values[2],
  nu_1 = mag_values[3],
  xi_1 = mag_values[4],
  m_0 = m_0_value,
  prop_trig = mean(cat$b > 0),
  lower = 0,
  upper = 10)

pdf("./output/parameter_recovery/single_magnitude_case/mag_return_level_plots_dual_vs_single.pdf",
    width = 4, height =  4)
  par(mar = c(4.5,4.5,1.1,1.1))
undebug(plot_sampled_rls_single)
  plot_sampled_rls_single(
    sampled_rls = single_mag_rls[-(1:100),],
    alpha = 0.05,
    catalogue = cat,
    line_colours = c("lightsalmon","darkorange","lightsalmon"),
    line_type = 2,
    cex.lab = 1.5,
    ylim = c(1,10))

  plot_sampled_rls(
    sampled_rls = dual_return_levels[-(1:100),],
    alpha = 0.05,
    catalogue = cat,
    event_type = 2,
    add = TRUE)

  plot_sampled_rls_single(
    sampled_rls = single_mag_rls[-(1:100),],
    alpha = 0.05,
    catalogue = cat,
    line_colours = c("lightsalmon","darkorange","lightsalmon"),
    line_type = 2,
    add = TRUE)

  lines(x = return_periods, y = true_return_levels_combined, lty = 3, lwd = 2 )
dev.off()

## 06: Assessing Branching recovery --------------------------------------------


# Posterior proportion triggered
pdf("./output/parameter_recovery/single_magnitude_case/B_posterior_proportion_triggered_dual_vs_single.pdf",
    width = 7,
    height = 5)
par(mar = c(4.5,4.5,2.1,2.1))
plot(density(rowMeans(mcmc_3$b >0),from = 0.1, to = 0.25),
     col = "blue",
     lwd = 3,
     xlim = c(0.18,0.22),
     xlab = "proportion of events triggered",
     ylab = "posterior density",
     main ="",
     cex.lab = 1.25)
lines(density(rowMeans(mcmc_4$b > 0)), col = "orange", lwd = 3)
abline(v = mean(cat$b >0), col = 1, lty = 2, lwd = 3)
dev.off()

# Posterior proportion correctly allocated
prop_correct_samples_3 <- apply(mcmc_3$b, MARGIN = 1, FUN = function(x){mean(x == cat$b)})
prop_correct_samples_4 <- apply(mcmc_4$b, MARGIN = 1, FUN = function(x){mean(x == cat$b)})

pdf("./output/parameter_recovery/single_magnitude_case/B_posterior_proportion_correct_dual_vs_single.pdf",
    width = 7,
    height = 5)
par(mar = c(4.5,4.5,2.1,2.1))
plot(density(prop_correct_samples_4,from = 0.9, to = 1),
     col = "darkorange",
     lwd = 3,
     xlim = c(0.9,1.0),
     xlab = "proportion of B correct",
     ylab = "posterior density",
     main ="",
     cex.lab = 1.25)
lines(density(prop_correct_samples_3, from = 0.9, to = 1), col = "blue", lwd = 3)
abline(v = 1, col = 1, lty = 2, lwd = 3)
dev.off()

# Element-wise change in performance
elementwise_proportion_correct_3 <- rowMeans(apply(mcmc_3$b, MARGIN =  1, function(x){x == cat$b}))
elementwise_proportion_correct_4 <- rowMeans(apply(mcmc_4$b, MARGIN =  1, function(x){x == cat$b}))
change_in_proportion_correct <- elementwise_proportion_correct_3 - elementwise_proportion_correct_4
colouring <- ifelse(change_in_proportion_correct>0, yes= "blue", no = "darkorange")

pdf("./output/parameter_recovery/single_magnitude_case/B_elementwise_change_in_posterior_proportion_correct_dual_vs_single.pdf",
    width = 7,
    height = 5)
par(mar = c(4.5,5.1,2.1,2.1))
plot(change_in_proportion_correct, col = colouring, type = "h", xlab = "Event index", ylab = "difference in posterior probability \n of correct allocation")
dev.off()

mean(change_in_proportion_correct)
mean(change_in_proportion_correct[change_in_proportion_correct != 0])

## 07 : ETAS parameter recovery ------------------------------------------------

pdf(file = "./output/parameter_recovery/single_magnitude_case/ETAS_parameter_recovery_dual_vs_single.pdf",
    width = 4, height = 4)

par(mfrow = c(1,1), mar = c(4.5,4.5,2.1,2.1))

d_3 <- density(mcmc_3$etas$mu, from = 0.017, to = 0.022)
d_4 <- density(mcmc_4$etas$mu, from = 0.017, to = 0.022)
plot(d_3, lwd = 2, xlab = expression(mu), ylab = "Posterior density", col = "blue", main ='')
lines(d_4, col = "darkorange", lwd = 2)
abline(v = etas_values[1], lwd = 2, lty = 2)

d_3 <- density(mcmc_3$etas$C, from =  0.15, to =  0.25)
d_4 <- density(mcmc_4$etas$C, from =  0.15, to =  0.25)
plot(d_4, lwd = 2, xlab = "C", ylab = "Posterior density", col = "darkorange", main ='')
lines(d_3, col = "blue", lwd = 2)
abline(v = etas_values[2], lwd = 2, lty = 2)

d_3 <- density(mcmc_3$etas$a, from = -0.1, to = 1)
d_4 <- density(mcmc_4$etas$a, from =  -0.1, to =  1)
plot(d_3, lwd = 2, xlab = "a", ylab = "Posterior density", col = "blue", main ='')
lines(d_4, col = "darkorange", lwd = 2)
abline(v = etas_values[3], lwd = 2, lty = 2)

d_3 <- density(mcmc_3$etas$nu_t, from = 0.07, to = 0.15)
d_4 <- density(mcmc_4$etas$nu_t, from = 0.07, to = 0.15)
plot(d_3, lwd = 2, xlab = expression(nu[t]), ylab = "Posterior density", col = "blue", main ='')
lines(d_4, col = "darkorange", lwd = 2)
abline(v = etas_values[4], lwd = 2, lty = 2)

d_3 <- density(mcmc_3$etas$xi_t, from = -0.3, to = 0.6)
d_4 <- density(mcmc_4$etas$xi_t, from = -0.3, to = 0.6)
plot(d_4, lwd = 2, xlab = expression(xi[t]), ylab = "Posterior density", col = "darkorange", main ='')
lines(d_3, col = "blue", lwd = 2)
abline(v = etas_values[4], lwd = 2, lty = 2)

dev.off()


## EOF -------------------------------------------------------------------------
