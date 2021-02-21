################################################################################
# This code processes the posterior estimates and generates Figure S4.         #
#                                                                              # 
# Cite as:                                                                     #
# Dorigatti I et al, SARS-CoV-2 antibody dynamics, within-household            #
# transmission and the impact of contact tracing from community-wide           # 
# serological testing in the Italian municipality of Vo'.                      #   
#                                                                              #
# Description:                                                                 #
# See README.md                                                                # 
################################################################################

# load packages
library(ggplot2)
library(gtable)
library(cowplot)
library(grid)
library(binom)
library(LaplacesDemon)
library(stringr)
library(grid)
library(gridExtra)


# source script ---------------------------------------------------------------#
source("R/functions_plot_fitted_model.R")
source("R/functions_plot_deviance_DIC.R")

# read data -------------------------------------------------------------------#

res <- read.csv(file.path("mcmc_posterior_chains_original",
                "mcmc_posterior_estimates_original.csv"), 
                header = TRUE)

rese <- read.csv(file.path("mcmc_posterior_chains_extended",
                 "mcmc_posterior_estimates_extended.csv"), 
                header = TRUE)

resg <- read.csv(file.path("mcmc_posterior_chains_2_groups",
                 "mcmc_posterior_estimates_2_groups.csv"), 
                 header = TRUE)

models <- read.csv("data/household_model/Model_variants_original.csv", 
                   header = TRUE)

modelse <- read.csv("data/household_model/Model_variants_extended.csv", 
                   header = TRUE)

modelsg <- read.csv("data/household_model/Model_variants_2_groups.csv", 
                   header = TRUE)

default_par <- read.csv("data/household_model/Default_parameters_original.csv", 
                        header = TRUE)

default_parg <- read.csv("data/household_model/Default_parameters_2_groups.csv", 
                         header = TRUE)

data <- read.csv("tables/Household_final_size_Vo_baseline.csv", 
                 header = TRUE)

data <- data[, -1]
colnames(data) <- seq(1,7,1)
rownames(data) <- seq(0,4,1)

colnames(models)[7] <- "p_i"
colnames(models)[9] <- "p_sd"
colnames(models)[5] <- "k"

colnames(modelse)[7] <- "p_i"
colnames(modelse)[9] <- "p_sd"
colnames(modelse)[5] <- "k"

colnames(modelsg)[6] <- "p_i"
colnames(modelsg)[8] <- "p_sd"

summary_res <- res 
summary_res_e <- rese
summary_res_g <- resg

# set directories --------------------------------------------------------------#

dir_output  <- file.path("summary_results_household_model")
dir_figures <- file.path(dir_output, "figures")
dir.create(dir_output,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

#---------------------------- Original model ----------------------------------#

pout <- plot_deviance_dic_params(res, models, "original model", NA)

r_dev <- pout[[1]]
r_dic <- pout[[2]]
q <- pout[[3]]
mindic <- pout[[4]]
selected_chains <- pout[[5]]

column1 <- plot_grid(r_dic, q, 
                    labels = c('a','b'),
                    rel_heights = c(1,1),
                    label_size = 8,
                    ncol = 1,
                    label_fontface = "bold")

#---------------------------- Extended model ----------------------------------#

pout1 <- plot_deviance_dic_params(rese, modelse, "extended model", mindic)

r_dev_e <- pout1[[1]]
r_dic_e <- pout1[[2]]
q_e <- pout1[[3]]
selected_chains_e <- pout1[[5]]

column2 <- plot_grid(r_dic_e, q_e, 
                     labels = c('c','d'),
                     rel_heights = c(1,1),
                     label_size = 8,
                     ncol = 1,
                     label_fontface = "bold")

#---------------------------- 2-groups model ----------------------------------#

pout2 <- plot_deviance_dic_params(resg, modelsg, "2-groups model", mindic)

r_dev_2 <- pout2[[1]]
r_dic_2 <- pout2[[2]]
q_2 <- pout2[[3]]
selected_chains_2 <- pout2[[5]]

column3 <- plot_grid(r_dic_2, q_2, 
                     labels = c('e','f'),
                     rel_heights = c(1,1),
                     label_size = 8,
                     ncol = 1,
                     label_fontface = "bold")

# save plot -------------------------------------------------------------------#

ggsave(filename =
         file.path(file.path(dir_output, "figures"),
                   "Figure_S4.tiff"),
       plot = plot_grid(column1, column2, column3, 
                        rel_widths = c(1.5,1,1), 
                        ncol = 3),
       device = "tiff",
       width = 240, height = 120,
       units = "mm", dpi = 300, limitsize = TRUE)
# -----------------------------------------------------------------------------#

################################################################################
################################################################################
