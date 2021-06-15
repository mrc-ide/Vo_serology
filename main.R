################################################################################
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


# Load packages ---------------------------------------------------------------#

library(binom)
library(car)
library(cowplot)
library(data.table)
library(dplyr)
library(gdata)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(grid)
library(gridExtra)
library(gtable)
library(latex2exp)
library(LaplacesDemon)
library(odin)
library(stringr)
library(tidyr)
library(wesanderson)
library(rmdformats)
library(viridis)

# Run scripts -----------------------------------------------------------------#

source("scripts/Descriptive_analysis_script_1.R")

source("scripts/Estimate_true_seroprevalence_script_2.R")

source("scripts/Seroprevalence_script_3.R")

source("scripts/Association_antibody_titres_script_4.R")

# N.B. make sure to run step 4 before running step 5, as script 5 uses the
# output of script 4
source("scripts/Association_antibody_decay_script_5.R")

# N.B. this script takes time to run,
# uncomment lines 121-122 for a quicker test run
source("scripts/Fit_original_household_model_script_6.R")

# N.B. this script takes time to run,
# uncomment lines 126-127 for a quicker test run
source("scripts/Fit_extended_household_model_script_7.R")

# N.B. this script takes time to run,
# uncomment lines 123-124 for a quicker test run
source("scripts/Fit_2_groups_household_model_script_8.R")

# N.B. make sure to run step 6, 7 and 8 before this step, as this script uses
# the output from these script
source("scripts/Plot_DIC_script_9.R")

# N.B. make sure to run step 6 with the original parameterisation,
# i.e. with lines 121-122 commented
source("scripts/Plot_SITP_overdisp_fit_original_model_script_10.R")

# Produces Figure 1, Figure 6, and Table S9
# N.B. uncomment lines 46-49 for a quicker test run
source("scripts/Impact_of_contact_tracing_11.R")

# Produces Figure S1, Tables S5 and S6
rmarkdown::render("scripts/Complementary_SI_info_12.Rmd",
                  knit_root_dir = getwd(),
                  output_dir = "figures")

source("scripts/Neutralization_decay_script_13.R")
