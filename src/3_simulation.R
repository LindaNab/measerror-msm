# -----------------------------------------------------------------------------
# This code entails the simulation study in 'Quantitative bias analysis for a 
# misclassified confounding variable in point-treatment studies: a comparison 
# between marginal structural models and conditional models'.
# Author: Linda Nab, l.nab@lumc.nl
# Date: 05122019
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Load required libraries
require(rsimsum)
require("xtable")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1. Load simulated datasets ---------------------------------------------------
simDataDir <- "../data/processed"
seeds <- readRDS(paste0(dataDir, "/seeds.rds"))
# read simulated sets-----------------------------------------------------------
simData <- list("scen0" = list("ss1000" = matrix(), "ss100" = matrix()),
                "scen1" = list("ss1000" = matrix(), "ss100" = matrix()),
                "scen2" = list("ss1000" = matrix(), "ss100" = matrix()),
                "scen3" = list("ss1000" = matrix(), "ss100" = matrix()),
                "scen4" = list("ss1000" = matrix(), "ss100" = matrix()))
simData$scen0$ss1000 <- readRDS(file = paste0(simDataDir, "/sim_s0_1000.rds"))
simData$scen1$ss1000 <- readRDS(file = paste0(simDataDir, "/sim_s1_1000.rds"))
simData$scen2$ss1000 <- readRDS(file = paste0(simDataDir, "/sim_s2_1000.rds"))
simData$scen3$ss1000 <- readRDS(file = paste0(simDataDir, "/sim_s3_1000.rds"))
simData$scen4$ss1000 <- readRDS(file = paste0(simDataDir, "/sim_s4_1000.rds"))
simData$scen0$ss100 <- readRDS(file = paste0(simDataDir, "/sim_s0_100.rds"))
simData$scen1$ss100 <- readRDS(file = paste0(simDataDir, "/sim_s1_1000.rds"))
simData$scen2$ss100 <- readRDS(file = paste0(simDataDir, "/sim_s2_1000.rds"))
simData$scen3$ss100 <- readRDS(file = paste0(simDataDir, "/sim_s3_1000.rds"))
simData$scen4$ss100 <- readRDS(file = paste0(simDataDir, "/sim_s4_1000.rds"))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 2. Load functions
# Read for the derivations of the bias expressions the appendix of the 
# original paper
# ------------------------------------------------------------------------------
# bias in conditional model
bias_cm <- function(par){
  alpha <- 1
  beta <- par[1]; gamma <- par[2]; lambda <- par[3]
  p_0 <- par[4]; p_1 <- par[5]; pi_0 <- par[6]; pi_1 <- par[7]
  #omega = P(A)
  omega <- pi_0 * (1-lambda) + pi_1 * lambda
  #ell = P(L^*)
  ell <- p_0 * (1-lambda) + p_1 * lambda
  #pi_star = P(A|L_star=l_star)
  pi_star <- function(l_star){
    term0 <- pi_0 * ((1-p_0)^(1-l_star)*p_0^l_star*(1-lambda)) / 
      ((1-ell)^(1-l_star)*ell^l_star)
    term1 <- pi_1 * ((1-p_1)^(1-l_star)*p_1^l_star*lambda) / 
      ((1-ell)^(1-l_star)*ell^l_star)
    return(term0 + term1)
  }
  #phi = P(L|A=a,L_star=l_star)
  phi <- function(a, l_star){
    out <- (lambda*(1-pi_1)^(1-a)*pi_1^a*(1-p_1)^(1-l_star)*p_1^l_star) / 
      ((1-pi_star(l_star))^(1-a)*pi_star(l_star)^a*
         (1-ell)^(1-l_star)*ell^l_star)
    return(out)
  }
  #bias
  temp <- (pi_star(1)*ell*(1-omega)-pi_star(1)*(pi_star(1) - 
                                                  pi_star(0))*ell*(1-ell)) / 
    (omega*(1-omega)-(pi_star(1)-pi_star(0))^2*ell*(1-ell))
  bias_cm <- gamma * (phi(1,0) - phi(0,0)) * (1 - temp) + 
    gamma * (phi(1,1) - phi(0,1)) * temp
  return(bias_cm)
}
# bias in msm estimated using ipw
bias_msm <- function(par){
  alpha <- 1
  beta <- par[1]; gamma <- par[2]; lambda <- par[3]
  p_0 <- par[4]; p_1 <- par[5]; pi_0 <- par[6]; pi_1 <- par[7]
  #omega = P(A)
  omega <- pi_0 * (1-lambda) + pi_1 * lambda
  #ell = P(L^*)
  ell <- p_0 * (1-lambda) + p_1 * lambda
  #pi_star = P(A|L_star=l_star)
  pi_star <- function(l_star){
    term0 <- pi_0 * ((1-p_0)^(1-l_star)*p_0^l_star*(1-lambda)) / 
      ((1-ell)^(1-l_star)*ell^l_star)
    term1 <- pi_1 * ((1-p_1)^(1-l_star)*p_1^l_star*lambda) / 
      ((1-ell)^(1-l_star)*ell^l_star)
    return(term0 + term1)
  }
  #phi = P(L|A=a,L_star=l_star)
  phi <- function(a, l_star){
    out <- (lambda*(1-pi_1)^(1-a)*pi_1^a*(1-p_1)^(1-l_star)*p_1^l_star) / 
      ((1-pi_star(l_star))^(1-a)*pi_star(l_star)^a*
         (1-ell)^(1-l_star)*ell^l_star)
    return(out)
  }
  #bias
  bias_msm <- gamma * (phi(1,0) - phi(0,0)) * (1 - ell) + 
    gamma * (phi(1,1) - phi(0,1)) * ell
  return(bias_msm)
}
# 5 different simulation scenarios
create_scen <- function(beta = 1, gamma = 2, lambda, p_0, p_1, pi_0, pi_1) {
  c(beta, gamma, lambda, p_0, p_1, pi_0, pi_1)
}
# statistics from summary
stat_from_summ <- function(summ, stat, method, digitsb, digitse){
  stats_from_summ <- summ[summ$stat == stat,]
  stat_from_summ <- stats_from_summ[stats_from_summ$method == method,]
  out <- paste0(round(stat_from_summ$est, digitsb),
                " (", round(stat_from_summ$mcse, digitse), ")")
  return(out)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 3. Create scenarios
# ------------------------------------------------------------------------------
scen0 <- create_scen(lambda = 0.5, p_0 = 0, p_1 = 1, pi_0 = 0.5, 
                     pi_1 = 0.75)
scen1 <- create_scen(lambda = 0.5, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.9, 
                     pi_1 = 0.45)
scen2 <- create_scen(lambda = 0.8, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.5, 
                     pi_1 = 0.75)
scen3 <- create_scen(lambda = 0.8, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.25, 
                     pi_1 = 0.75)
scen4 <- create_scen(lambda = 0.45, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.5, 
                     pi_1 = 0.75)

# create simsum objects of simulated sets --------------------------------------
simsum.sim_s0_1000 <- simsum(simData$scen0$ss1000, estvarname = "b", true = scen0[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s1_1000 <- simsum(simData$scen1$ss1000, estvarname = "b", true = scen1[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s2_1000 <- simsum(simData$scen2$ss1000, estvarname = "b", true = scen2[1],
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s3_1000 <- simsum(simData$scen3$ss1000, estvarname = "b", true = scen3[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s4_1000 <- simsum(simData$scen4$ss1000, estvarname = "b", true = scen4[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s0_100 <- simsum(simData$scen0$ss100, estvarname = "b", true = scen0[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s1_100 <- simsum(simData$scen1$ss100, estvarname = "b", true = scen1[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s2_100 <- simsum(simData$scen2$ss100, estvarname = "b", true = scen2[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s3_100 <- simsum(simData$scen3$ss100, estvarname = "b", true = scen3[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s4_100 <- simsum(simData$scen4$ss100, estvarname = "b", true = scen4[1], 
                             se = "se", methodvar = "method", x = TRUE)

# Create tables for in paper ---------------------------------------------------
summ_s0_1000 <- summary(simsum.sim_s0_1000)$summ
summ_s1_1000 <- summary(simsum.sim_s1_1000)$summ
summ_s2_1000 <- summary(simsum.sim_s2_1000)$summ
summ_s3_1000 <- summary(simsum.sim_s3_1000)$summ
summ_s4_1000 <- summary(simsum.sim_s4_1000)$summ
summ_s0_100 <- summary(simsum.sim_s0_100)$summ
summ_s1_100 <- summary(simsum.sim_s1_100)$summ
summ_s2_100 <- summary(simsum.sim_s2_100)$summ
summ_s3_100 <- summary(simsum.sim_s3_100)$summ
summ_s4_100 <- summary(simsum.sim_s4_100)$summ
#
summ_list <- list(summ_s0_1000, summ_s1_1000, summ_s2_1000, summ_s3_1000, 
                  summ_s4_1000, summ_s0_100, summ_s1_100, summ_s2_100, 
                  summ_s3_100, summ_s4_100)
# 
stat_from_summ <- function(summ, stat, method, digitsb, digitse){
  stats_from_summ <- summ[summ$stat == stat,]
  stat_from_summ <- stats_from_summ[stats_from_summ$method == method,]
  out <- paste0(round(stat_from_summ$est, digitsb),
                " (", round(stat_from_summ$mcse, digitse), ")")
  return(out)
}


results <- data.frame(method = numeric(20), 
                      sample_size = numeric(20), 
                      scenario = numeric(20), 
                      bias_form = numeric(20), 
                      bias = numeric(20), 
                      mse = numeric(20), 
                      coverage = numeric(20), 
                      empse = numeric(20),
                      modelse = numeric(20))
results$method <- c(rep("msm", 10), rep("cm", 10))
results$sample_size <- c(rep("1000", 5), rep("100", 5), 
                         rep("1000", 5), rep("100", 5))
results$scenario <- rep(c("0", "1", "2", "3", "4"), 4)
# results$scenario <- c(rep("s0", 4), rep("s1", 4), rep("s2", 4), 
# rep("s3", 4), rep("s4", 4))
# results$sample_size <- rep(c("1000", "1000", "100", "100"), 5)
# results$method <- rep(c("msm", "cm"), 10)
results$bias_form[1] <- round(bias_msm(scen0), 2)
results$bias_form[2] <- round(bias_msm(scen1), 2)
results$bias_form[3] <- round(bias_msm(scen2), 2)
results$bias_form[4] <- round(bias_msm(scen3), 2)
results$bias_form[5] <- round(bias_msm(scen4), 2)
results$bias_form[6:10] <- results$bias_form[1:5]
#
results$bias_form[11] <- round(bias_cm(scen0), 2)
results$bias_form[12] <- round(bias_cm(scen1), 2)
results$bias_form[13] <- round(bias_cm(scen2), 2)
results$bias_form[14] <- round(bias_cm(scen3), 2)
results$bias_form[15] <- round(bias_cm(scen4), 2)
results$bias_form[16:20] <- results$bias_form[11:15]
#bias
for(j in 1:2){
  for(i in 1:10){
    results$bias[10*(j-1)+i] <- stat_from_summ(summ_list[[i]], "bias", 
                                               results$method[10*(j-1)+i], 2, 3)
  }
}
#mse
for(j in 1:2){
  for(i in 1:10){
    results$mse[10*(j-1)+i] <- stat_from_summ(summ_list[[i]], "mse", 
                                              results$method[10*(j-1)+i], 2, 3)
  }
}
#coverage
for(j in 1:2){
  for(i in 1:10){
    results$coverage[10*(j-1)+i] <- stat_from_summ(summ_list[[i]], "cover", 
                                                   results$method[10*(j-1)+i], 
                                                   2, 3)
  }
}
#empse
for(j in 1:2){
  for(i in 1:10){
    results$empse[10*(j-1)+i] <- stat_from_summ(summ_list[[i]], "empse", 
                                                results$method[10*(j-1)+i], 
                                                2, 3)
  }
}
#modelse
for(j in 1:2){
  for(i in 1:10){
    results$modelse[10*(j-1)+i] <- stat_from_summ(summ_list[[i]], "modelse", 
                                                  results$method[10*(j-1)+i], 
                                                  2, 3)
  }
}
xtable(results, include.rownames = F)

