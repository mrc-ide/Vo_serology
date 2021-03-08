################################################################################
# This code estimates the true prevalence and the sensitivity, specificity,    #
# positive predictive value (PPV) and negative predictive value (NPV) for each #
# assay, using MCMC and a multinomial likelihood (see Supplementary Methods).  #                                                      #
#                                                                              #
# Cite as:                                                                     #
# Dorigatti I et al, SARS-CoV-2 antibody dynamics, within-household            #
# transmission and the impact of contact tracing from community-wide           #
# serological testing in the Italian municipality of Vo'.                      #
# Pre-print available at: http://dx.doi.org/10.2139/ssrn.3792892               #
#                                                                              #
# Description:                                                                 #
# See README.md                                                                #
################################################################################

cat("\n\n### Running 'scripts/Estimate_true_seroprevalence_script_2.R'\n\n")


# source scripts
source("R/functions_multinomial_likelihood.R")

# read data
ds <- read.csv("data/Vo_serology_data.csv",
               header = TRUE)

# output dir
path_figure <- file.path(getwd(), "figures")
path_table <- file.path(getwd(), "tables")

# Frequency assay results -----------------------------------------------------#

# select Vo residents
ds <- ds[which(ds$town == "VO"), ]

tab <- table(ds$Abbot_qualitative,
      ds$Diasorin_IgG_qualitative,
      ds$Roche_Total_qualitative)

# Estimate true prevalence ----------------------------------------------------#

data <- c(tab["Positive", "Positive", "Positive"],
          tab["Positive", "Negative", "Positive"],
          tab["Positive", "Positive", "Negative"],
          tab["Negative", "Positive", "Positive"],
          tab["Positive", "Negative", "Negative"],
          tab["Negative", "Positive", "Negative"],
          tab["Negative", "Negative", "Positive"],
          tab["Negative", "Negative", "Negative"])

# observed sensitivity
se_A_n <- 31 # positives Abbott
se_A_d <- 33 # tested Abbott

se_D_n <- 86 # positives DiaSorin
se_D_d <- 101 # tested DiaSorin

se_R_n <- 31 # positives Roche
se_R_d <- 33 # tested Roche

# observed specificity
sp_A_n <- 54 # negatives Abbott
sp_A_d <- 54 # tested Abbott

sp_D_n <- 20 # negatives DiaSorin
sp_D_d <- 20 # tested DiaSorin

sp_R_n <- 41 # negatives Roche
sp_R_d <- 42 # tested Roche

# sensitivites
sens_A <- se_A_n/se_A_d
sens_D <- se_D_n/se_D_d
sens_R <- se_R_n/se_R_d

# specificities
spec_A <- sp_A_n/sp_A_d
spec_D <- sp_D_n/sp_D_d
spec_R <- sp_R_n/sp_R_d

# prevalence (initial value)
theta <- 0.04

# observed numerators and denominators
obs <- c(se_A_n, se_A_d, se_D_n, se_D_d, se_R_n, se_R_d,
         sp_A_n, sp_A_d, sp_D_n, sp_D_d, sp_R_n, sp_R_d)

# parameters to estimate
param <- c(theta, sens_A, sens_D, sens_R, spec_A, spec_D, spec_R)
parameterNames <- c("theta", "sens_A", "sens_D", "sens_R", "spec_A", "spec_D",
                    "spec_R")
parametersToEstimate <- rep(1, 7)

# parameters for MCMC
# niter, thinning, lowe_bound, upper_bound, jump_rate
mcmc_par <- c(100000, 10, 0.000001, 1, 0.05)

# run MCM
res <- run_MCMC(data, mcmc_par, param, obs)

# plot posterior estimates
burnin <- 100
plot_posterior_estimates(burnin,
                         parametersToEstimate,
                         parameterNames,
                         res)

# output chains
write.csv(res$ParameterValues,
          file.path(path_table,"Posterior_estimates_sens_spec.csv"))

# true prevalence
round(c(mean(res$ParameterValues[,1]),
        quantile(res$ParameterValues[,1], c(0.025, 0.975))), 3)

#sens Abbott
round(c(mean(res$ParameterValues[,2]),
        quantile(res$ParameterValues[,2], c(0.025, 0.975))),3)

#sens DiaSorin
round(c(mean(res$ParameterValues[,3]),
        quantile(res$ParameterValues[,3], c(0.025, 0.975))), 3)

#sens Roche
round(c(mean(res$ParameterValues[,4]),
      quantile(res$ParameterValues[,4], c(0.025, 0.975))),3)

#spec Abbott
round(c(mean(res$ParameterValues[,5]),
        quantile(res$ParameterValues[,5], c(0.025, 0.975))), 3)

#spec DiaSorin
round(c(mean(res$ParameterValues[,6]),
        quantile(res$ParameterValues[,6], c(0.025, 0.975))),3)

#spec Roche
round(c(mean(res$ParameterValues[,7]),
        quantile(res$ParameterValues[,7], c(0.025, 0.975))),3)

# PPV Abbott
PPV_A <- res$ParameterValues[,1]*res$ParameterValues[,2]/
  (res$ParameterValues[,1]*res$ParameterValues[,2] +
     (1-res$ParameterValues[,5])*(1-res$ParameterValues[,1]))
mean(PPV_A)
quantile(PPV_A, c(0.025, 0.975))

# PPV DiaSorin
PPV_D <- res$ParameterValues[,1]*res$ParameterValues[,3]/
  (res$ParameterValues[,1]*res$ParameterValues[,3] +
     (1-res$ParameterValues[,6])*(1-res$ParameterValues[,1]))
mean(PPV_D)
quantile(PPV_D, c(0.025, 0.975))

# PPV Roche
PPV_R <- res$ParameterValues[,1]*res$ParameterValues[,4]/
  (res$ParameterValues[,1]*res$ParameterValues[,4] +
     (1-res$ParameterValues[,7])*(1-res$ParameterValues[,1]))
mean(PPV_R)
quantile(PPV_R, c(0.025, 0.975))

# NPV Abbott
NPV_A <- (1-res$ParameterValues[,1])*res$ParameterValues[,5]/
  ((1-res$ParameterValues[,1])*res$ParameterValues[,5] +
     (1-res$ParameterValues[,2])*(res$ParameterValues[,1]))
mean(NPV_A)
quantile(NPV_A, c(0.025, 0.975))

# NPV DiaSorin
NPV_D <- (1-res$ParameterValues[,1])*res$ParameterValues[,6]/
  ((1-res$ParameterValues[,1])*res$ParameterValues[,6] +
     (1-res$ParameterValues[,3])*res$ParameterValues[,1])
mean(NPV_D)
quantile(NPV_D, c(0.025, 0.975))

# NPV Roche
NPV_R <- (1-res$ParameterValues[,1])*res$ParameterValues[,7]/
  ((1-res$ParameterValues[,1])*res$ParameterValues[,7] +
     (1-res$ParameterValues[,4])*res$ParameterValues[,1])
mean(NPV_R)
quantile(NPV_R, c(0.025, 0.975))

################################################################################
################################################################################



