################################################################################
# Cite as:                                                                     #
# Dorigatti I et al, SARS-CoV-2 antibody dynamics, within-household            #
# transmission and the impact of contact tracing from community-wide           #
# serological testing in the Italian municipality of Vo'.                      #
# Pre-print available at: http://dx.doi.org/10.2139/ssrn.3792892               #
#                                                                              #
# Code to estimate the impact of contact tracing on the spread of COVID-19 in  #
# Vo' Euganeo in spring 2020                                                   #
#                                                                              #
# Produces Figure 1, Figure 6, and Table S9.                                   #
#                                                                              #
# code adapted from: Lavezzo et. al., 2020, Suppression of COVID-19 outbreak   #
# in the municipality of Vo, Italy, medRxiv, doi: 10.1101/2020.04.17.20053157  #
#                                                                              #
# data input:                                                                  #
# - observed number of negative, pre-symptomatic, symptomatic and asymptomatic #
#   individuals after two screening runs of large parts of the population of   #
#   Vo' Euganeo                                                                #
# - the number of positive cases among all traced contacts                     #
# - the number of traced subjects among all positive cases                     #
#                                                                              #
################################################################################


# Source functions

source("R/functions_model_11.R")
source("R/functions_clean_11.R")
source("R/functions_tables_11.R")
source("R/functions_figures_11.R")


dir_figures <- "figures" # Figure 1, Figure 6
dir_output  <- "11_impact_of_contact_tracing"  # all other output (including Table S9)


# MCMC parameters -------------------------------------------------------------#

id_chain   <- 1:3
mcmc_iterations <- 200000 # number MCMC iterations
sample_spacing  <- 100    # every how many iterations do we save the MCMC output
nr_burnin <- 200 # number stored parameter combinations we condiser to be burnin
nr_sample <- 100 # number stored parameter combinations we sample for plotting

# mcmc_iterations <- 10000
# sample_spacing  <- 50
# nr_burnin <- mcmc_iterations/sample_spacing/10
# nr_sample <- mcmc_iterations/sample_spacing - nr_burnin


# Observed data ---------------------------------------------------------------#

# Observed data at first and second sampling
# weighted average time of first sampling = 24 Feb 2020
# weighted average time of second sampling = 07 Mar 2020
data <- data.frame(Tested         = c(2812, 2343),
                   Asymptomatic   = c(  29,   13),
                   Symptomatic    = c(  34,   15),
                   Presymptomatic = c(  10,    1))

data2 <- data.frame(traced_contacts = 190,          # nr of traced contacts
                    negative_traced_contacts = 146) # nr of traced contacts that never tested positive

data3 <- data.frame(PCR_positive_subjects = 100,  # all the individuals who had at least on positive PCR test
                    PCR_positive_not_traced = 56) # traced contacts that had a positive PCR test


# Prepare model ---------------------------------------------------------------#

# Parameters whose value is fixed in the model for all analyses
fixed_parameters <- data.frame(
  tSeed   = 0L,   # first day of model (4th February 2020)
  time1   = 22L,  # weighted average day of first sampling (26 February 2020)
  time2   = 32L,  # weighted average day of second sampling (8 March 2020)
  tQ1     = 20L,  # day lockdown and contact tracing started (24 February 2020)
  tQ2     = 57L,  # day contact tracing ended (1 April 2020)
  tf      = 240L, # final day when running the simulation after MCMC
  N       = 3275, # Vo' cluster size
  q_TPpre = 1,    # relative infectiousness of pre-symptomatics w.r.t. symptomatics
  q_A     = 1,    # relative infectiousness of asymptomatics w.r.t. symptomatics
  q_Q     = 0,    # relative infectiousness of quarantined w.r.t. symptomatics
  sigma   = 1/4,
  generation_time = 7
)

# Parameters for which we test different values
looped_parameters <- data.frame(
  R0_1 = seq(2.1, 2.7, by = 0.3)
)

# Parameters to fit
lims <- data.frame(seed      = c(1, 3), # initialisation (uniform sample)
                   p         = c(0.35, 0.45),
                   w         = c(0.05, 0.15),
                   inv_nu    = c(1.5,  2.5),
                   inv_delta = c(1,    2),
                   p_traced  = c(0.01, 0.1),
                   p_traced_S= c(0.01, 0.1))
limits <- data.frame(seed      = c(1,    10), # MCMC interval boundaries
                     p         = c(0.01,  1),
                     w         = c(0.001, 1),
                     inv_nu    = c(0.01, fixed_parameters$generation_time),
                     inv_delta = c(0.01, fixed_parameters$generation_time),
                     p_traced  = c(0.01,  1),
                     p_traced_S= c(0.001, 1))
random_walk_rate <- data.frame(seed      = 0.05,
                               p         = 0.01,
                               w         = 0.05,
                               inv_nu    = 0.01,
                               inv_delta = 0.01,
                               p_traced  = 0.01,
                               p_traced_S= 0.001)


# Run MCMC (loop) -------------------------------------------------------------#

side_effects <- lapply(seq_len(nrow(looped_parameters)),
                       wrapper_model2, # main model function
                       data              = data,
                       fixed_parameters  = fixed_parameters,
                       looped_parameters = looped_parameters,
                       lims              = lims,
                       limits            = limits,
                       random_walk_rate  = random_walk_rate,
                       mcmc_iterations   = mcmc_iterations,
                       sample_spacing    = sample_spacing,
                       dir_output        = dir_output,
                       id_chain          = id_chain)


# Process results -----------------------------------------------------------#

print("Clean simulation output")

# Clean results
clean_posteriors(nr_burnin, nr_sample)


# Create figures and tables ---------------------------------------------------#

dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)
print("Create figures and tables")

# Make test figures
fig_chains()
fig_acceptance_rates()

# Tables
table_fitted(data)
table_final_size()

# Figure 1
pA <- fig_flowchart()
pB <- fig_timeline()

p <- arrangeGrob(pA, pB, ncol = 1, heights = c(1, 0.3))
ggsave(filename = file.path(dir_figures, "Fig1.tiff"),
       plot = p,
       device = "tiff",
       width = 183, height = 200,
       units = "mm", dpi = 300,
       type = "cairo", compression = "lzw")
ggsave(filename = file.path(dir_figures, "Fig1.pdf"),
       plot = p,
       device = "pdf",
       width = 183, height = 200,
       units = "mm", dpi = 300)

# Figure 6
pA <- fig_venn(panel = "a")
pB <- fig_contact_tracing(panel = "b")
pC <- fig_PCR_testing(panel = "c")
pD <- fig_prevalence(panel = "d")
pE <- fig_incidence(panel = "e")
pF <- fig_relative_final_size(panel = "f")
pG <- fig_relative_final_size(panel = "g")
pH <- fig_relative_final_size(panel = "h")

lay <- matrix(c(1,1,1,1,1,1,2,2,2,3,3,3,
                4,4,4,4,4,4,5,5,5,5,5,5,
                6,6,6,6,7,7,7,7,8,8,8,8), byrow = TRUE, ncol = 12)
p <- arrangeGrob(pA, pB, pC, pD, pE, pF, pG, pH,
                 layout_matrix = lay,
                 widths = c(1,1,1,0.5,1.5,1,1,1,1,1,1,1))
ggsave(filename = file.path(dir_figures, "Fig6.tiff"),
       plot = p,
       device = "tiff",
       width = 183, height = 118,
       units = "mm", dpi = 300,
       type = "cairo", compression = "lzw")
ggsave(filename = file.path(dir_figures, "Fig6.pdf"),
       plot = p,
       device = "pdf",
       width = 183, height = 118,
       units = "mm", dpi = 300)
