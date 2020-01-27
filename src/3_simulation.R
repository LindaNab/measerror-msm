# -----------------------------------------------------------------------------
# This code entails the simulation study in 'Quantitative bias analysis for a 
# misclassified confounding variable in point-treatment studies: a comparison 
# between marginal structural models and conditional models'.
# Author: Linda Nab, l.nab@lumc.nl
# Date: 05122019
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load required libraries
# require(ipw)
# require(survey)
require(rsimsum)
require("xtable")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Bias formulae -------------------------------------------------------------
#-------------------------------------------------------------------------------
# read for the derivations of these bias expressions the appendix of the 
# original paper
#-------------------------------------------------------------------------------
#bias in conditional model
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

#bias in msm estimated using ipw
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

# ------------------------------------------------------------------------------
# 1. Five scenarios 's'
# s is a matrix with entries: 1) beta = 1, 2) gamma = 2, 3) lambda, 4) p_0, 
#                             5) p_1,      6) pi_0,     7) pi_1
# ------------------------------------------------------------------------------
# Scenario 0 ('s0'): no measurement error---------------------------------------
# lambda = 0.5; p_0 = 0; p_1 = 1, pi_0 = 0.5; pi_1 = 0.75
s0 <- matrix(nrow = 1, ncol = 7)
s0[,1] <- 1; s0[,2] <- 2; s0[,3] <- 0.5; s0[,4] <- 0; s0[,5] <- 1
s0[,6] <- 0.5; s0[,7] <- 0.75 
# bias_cm(s0)
# bias_msm(s0)
# Scenario 1 ('s1')-------------------------------------------------------------
# lambda = 0.5; p_0 = 0.05; p_1 = 0.9; pi_0 = 0.9; pi_1 = 0.45, pi_1/pi_0 = 0.5
s1 <- matrix(nrow = 1, ncol = 7)
s1[,1] <- 1; s1[,2] <- 2; s1[,3] <- 0.5; s1[,4] <- 0.05; s1[,5] <- 0.9
s1[,6] <- 0.9; s1[,7] <- 0.45 
# bias_cm(s1)
# bias_msm(s1)
# Scenario 2 ('s2')-------------------------------------------------------------
# lambda = 0.8, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.5, pi_1 = 0.75, (pi_1/pi_0 = 3)
s2 <- matrix(nrow = 1, ncol = 7)
s2[,1] <- 1; s2[,2] <- 2; s2[,3] <- 0.8; s2[,4] <- 0.05; s2[,5] <- 0.9 
s2[,6] <- 0.5; s2[,7] <- 0.75 
# bias_cm(s2)
# bias_msm(s2)
# Scenario 3 ('s3')-------------------------------------------------------------
# lambda = 0.8, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.25, pi_1 = 0.75, (pi_1/pi_0 = 3)
s3 <- matrix(nrow = 1, ncol = 7)
s3[,1] <- 1; s3[,2] <- 2; s3[,3] <- 0.8; s3[,4] <- 0.05; s3[,5] <- 0.9
s3[,6] <- 0.25; s3[,7] <- 0.75 
# bias_cm(s3)
# bias_msm(s3)
# Scenario 4 ('s4'): equal bias-------------------------------------------------
# lambda = 0.45, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.5, pi_1 = 0.75
# (pi_1/pi_0 = 3)
s4 <- matrix(nrow = 1, ncol = 7)
s4[,1] <- 1; s4[,2] <- 2; s4[,3] <- 0.45; s4[,4] <- 0.05; s4[,5] <- 0.9
s4[,6] <- 0.5; s4[,7] <- 0.75
# bias_cm(s4)
# bias_msm(s4)

# Function that generates data--------------------------------------------------
# data_gen <- function(m, nobs, seed){
#   set.seed(seed)
#   alpha <- 1
#   beta  <- m[1]
#   gamma <- m[2]
#   lambda <- m[3]
#   p_0 <- m[4]
#   p_1 <- m[5]
#   pi_0 <- m[6]
#   pi_1 <- m[7]
#   L <- rbinom(nobs, 1, lambda) 
#   A <- ifelse(L==1, rbinom(nobs, 1, pi_1), rbinom(nobs, 1, pi_0)) 
#   Y <- alpha + beta * A + gamma * L + rnorm(nobs, 0, 1) 
#   Ls <- ifelse(L==1, rbinom(nobs, 1, p_1), rbinom(nobs, 1, p_0)) 
#   data <- data.frame(cbind(L, A, Y, Ls))
#   return(data)
# }

# Function that performs simulation---------------------------------------------
# sim <- function(m, S, nobs, seeds){
#   out <- data.frame(dataset = numeric(2*S), method = numeric(2*S),
#                     b = numeric(2*S), se = numeric(2*S))
#   for(i in 1:S){
#       d <- data_gen(m, nobs, seeds[i])
#       #sim msm
#       out[(2*i - 1), 1] <- i; out[i + (i - 1), 2] <- "msm"
#       temp <- ipwpoint(exposure = A, 
#                        family = "binomial", 
#                        link = "logit", 
#                        denominator = ~ Ls,
#                        numerator = ~ 1,
#                        data = d)
#       msm <- svyglm(Y ~ A, 
#                     design = svydesign(~ 1, 
#                                        weights = ~ temp$ipw.weights, 
#                                        data = d))
#       out[i + (i - 1), 3] <- msm$coefficients[2]
#       out[i + (i - 1), 4] <- summary(msm)$coefficients[2,2]
#       #sim cond m
#       fit <- lm(Y ~ A + Ls, d)
#       out[(2*i), 1] <- i; out[(2*i), 2] <- "cm"
#       out[(2*i), 3] <- fit$coef[2]
#       out[(2*i), 4] <- summary(fit)$coefficients[2,2]}
#   return(out)
# }

# Simulate data ----------------------------------------------------------------
# generate seeds
# seeds <- sample(1:1e8, 10000, replace = F)
simsum_dir <- "../data/processed"#directory where output of sim study is located
# saveRDS(seeds, paste0(simsum_dir, "/seeds.rds"))
seeds <- readRDS(paste0(simsum_dir, "/seeds.rds"))
# sample size 1000
# sim_s0_1000 <- sim(s0, 5000, 1e3, seeds[1:5000])
# sim_s1_1000 <- sim(s1, 5000, 1e3, seeds[1:5000])
# sim_s2_1000 <- sim(s2, 5000, 1e3, seeds[1:5000])
# sim_s3_1000 <- sim(s3, 5000, 1e3, seeds[1:5000])
# sim_s4_1000 <- sim(s4, 5000, 1e3, seeds[1:5000])
# sample size 100
# sim_s0_100 <- sim(s0, 5000, 1e2, seeds[5001:1e4])
# sim_s1_100 <- sim(s1, 5000, 1e2, seeds[5001:1e4])
# sim_s2_100 <- sim(s2, 5000, 1e2, seeds[5001:1e4])
# sim_s3_100 <- sim(s3, 5000, 1e2, seeds[5001:1e4])
# sim_s4_100 <- sim(s4, 5000, 1e2, seeds[5001:1e4])

# save simulated sets-----------------------------------------------------------
# saveRDS(sim_s0_1000, paste0(simsum_dir, "/sim_s0_1000.rds"))
# saveRDS(sim_s1_1000, paste0(simsum_dir, "/sim_s1_1000.rds"))
# saveRDS(sim_s2_1000, paste0(simsum_dir, "/sim_s2_1000.rds"))
# saveRDS(sim_s3_1000, paste0(simsum_dir, "/sim_s3_1000.rds"))
# saveRDS(sim_s4_1000, paste0(simsum_dir, "/sim_s4_1000.rds"))
# saveRDS(sim_s0_100, paste0(simsum_dir, "/sim_s0_100.rds"))
# saveRDS(sim_s1_100, paste0(simsum_dir, "/sim_s1_1000.rds"))
# saveRDS(sim_s2_100, paste0(simsum_dir, "/sim_s2_1000.rds"))
# saveRDS(sim_s3_100, paste0(simsum_dir, "/sim_s3_1000.rds"))
# saveRDS(sim_s4_100, paste0(simsum_dir, "/sim_s4_1000.rds"))

# read simulated sets-----------------------------------------------------------
sim_s0_1000 <- readRDS(file = paste0(simsum_dir, "/sim_s0_1000.rds"))
sim_s1_1000 <- readRDS(file = paste0(simsum_dir, "/sim_s1_1000.rds"))
sim_s2_1000 <- readRDS(file = paste0(simsum_dir, "/sim_s2_1000.rds"))
sim_s3_1000 <- readRDS(file = paste0(simsum_dir, "/sim_s3_1000.rds"))
sim_s4_1000 <- readRDS(file = paste0(simsum_dir, "/sim_s4_1000.rds"))
sim_s0_100 <- readRDS(file = paste0(simsum_dir, "/sim_s0_100.rds"))
sim_s1_100 <- readRDS(file = paste0(simsum_dir, "/sim_s1_1000.rds"))
sim_s2_100 <- readRDS(file = paste0(simsum_dir, "/sim_s2_1000.rds"))
sim_s3_100 <- readRDS(file = paste0(simsum_dir, "/sim_s3_1000.rds"))
sim_s4_100 <- readRDS(file = paste0(simsum_dir, "/sim_s4_1000.rds"))

# create simsum objects of simulated sets --------------------------------------
simsum.sim_s0_1000 <- simsum(sim_s0_1000, estvarname = "b", true = s0[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s1_1000 <- simsum(sim_s1_1000, estvarname = "b", true = s1[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s2_1000 <- simsum(sim_s2_1000, estvarname = "b", true = s2[1],
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s3_1000 <- simsum(sim_s3_1000, estvarname = "b", true = s3[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s4_1000 <- simsum(sim_s4_1000, estvarname = "b", true = s4[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s0_100 <- simsum(sim_s0_100, estvarname = "b", true = s0[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s1_100 <- simsum(sim_s1_100, estvarname = "b", true = s1[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s2_100 <- simsum(sim_s2_100, estvarname = "b", true = s2[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s3_100 <- simsum(sim_s3_100, estvarname = "b", true = s3[1], 
                             se = "se", methodvar = "method", x = TRUE)
simsum.sim_s4_100 <- simsum(sim_s4_100, estvarname = "b", true = s4[1], 
                             se = "se", methodvar = "method", x = TRUE)

# Create tables for in paper ---------------------------------------------------
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
results$bias_form[1] <- round(bias_msm(s0), 2)
results$bias_form[2] <- round(bias_msm(s1), 2)
results$bias_form[3] <- round(bias_msm(s2), 2)
results$bias_form[4] <- round(bias_msm(s3), 2)
results$bias_form[5] <- round(bias_msm(s4), 2)
results$bias_form[6:10] <- results$bias_form[1:5]
#
results$bias_form[11] <- round(bias_cm(s0), 2)
results$bias_form[12] <- round(bias_cm(s1), 2)
results$bias_form[13] <- round(bias_cm(s2), 2)
results$bias_form[14] <- round(bias_cm(s3), 2)
results$bias_form[15] <- round(bias_cm(s4), 2)
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

