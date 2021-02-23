################################################################################
#                                                                              #
# Cite as:                                                                     #
# Dorigatti I et al, SARS-CoV-2 antibody dynamics, within-household            #
# transmission and the impact of contact tracing from community-wide           #
# serological testing in the Italian municipality of Vo'.                      #
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
library(LaplacesDemon)
library(odin)
library(stringr)
library(tidyr)
library(wesanderson)


# Run scripts -----------------------------------------------------------------#

source("scripts/Descriptive_analysis_script_1.R")

source("scripts/Estimate_true_seroprevalence_script_2.R")

source("scripts/Seroprevalence_script_3.R")

source("scripts/Association_antibody_titres_script_4.R")

# N.B. make sure to run step 4 before this step, as this script uses an output
# from the previous script
source("scripts/Association_antibody_decay_script_5.R")

# N.B. this script takes time to run,
# uncomment lines 118-119 for a quicker test run
source("scripts/Fit_original_household_model_script_6.R")

# N.B. this script takes time to run,
# uncomment lines 124-125 for a quicker test run
source("scripts/Fit_extended_household_model_script_7.R")

# N.B. this script takes time to run,
# uncomment lines 121-122 for a quicker test run
source("scripts/Fit_2_groups_household_model_script_8.R")

source("scripts/Plot_DIC_script_9.R")

source("scripts/Plot_SITP_overdisp_fit_original_model_script_10.R")

source("scripts/Impact_of_contact_tracing_11.R")
