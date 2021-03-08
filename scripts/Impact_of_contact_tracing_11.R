################################################################################
# Cite as:                                                                     #
# Dorigatti I et al, SARS-CoV-2 antibody dynamics, within-household            #
# transmission and the impact of contact tracing from community-wide           #
# serological testing in the Italian municipality of Vo'.                      #
# Pre-print available at: https://ssrn.com/abstract=3792892                    #
#                                                                              #
# Code to estimate the impact of contact tracing on the spread of COVID-19 in  #
# Vo' Euganeo in spring 2020                                                   #
#                                                                              #
# Produces Figure 5.                                                           #
#                                                                              #
# adapted from: Lavezzo et. al., 2020, Suppression of COVID-19 outbreak in the #
# municipality of Vo, Italy, medRxiv, doi: 10.1101/2020.04.17.20053157         #
#                                                                              #
# data input:                                                                  #
# - observed number of negative, symptomatic and asymptomatic individuals      #
#   after two screening runs of large parts of the population of Vo' Euganeo   #
# - output file from Lavezzo et. al. 2020: parameter estimates for SEIR model  #
#   fitted to spring outbreak in Vo' Euganeo                                   #
#                                                                              #
################################################################################


# Source functions

source("R/functions_model_11.R")
source("R/functions_analysis_11.R")
source("R/functions_figures_11.R")


# Observed data at first and second sampling
# weighted average time of first sampling = 24 Feb 2020
# weighted average time of second sampling = 07 Mar 2020

data <- data.frame(Tested         = c(2812, 2343),
                   Asymptomatic   = c(  29,   13),
                   Symptomatic    = c(  34,   15),
                   Presymptomatic = c(  10,    1))


for (do.sensitivity_analysis in c(FALSE, TRUE)) {


  # define folder paths

  dir_in  <- "data"
  dir_fig <- "figures" # these figures will be published as part of the main
  dir_out <- "impact_of_contact_tracing" # all other output goes here
  dir_out <- if (do.sensitivity_analysis) file.path(dir_out, "sensitivity") else  dir_out
  dir.create(dir_out)


  # Compute and save SEIR for all interventions using posterior from Nature paper

  compute_interventions(do.sensitivity_analysis)


  # Table comparing epidemic final size across interventions

  tab_final_size()


  # for paper styling
  update_geom_defaults("text", list(size = 8))


  # create figures

  for (v in 1:3) {

    if (v == 1L) {
      do_conttrac <- 0.3
    } else if (v == 2L) {
      do_conttrac <- 0.06
    } else {
      do_conttrac <- 0.5
    }

    dt.filter <- data.frame(
      Intervention = c("No intervention", "Mass testing & lockdown", rep("Contact tracing", 3)),
      Quarantined  = c(rep("no transmission assumed", 3),
                       "30% transmission assumed", "50% transmission assumed"),
      v            = c(0, 0, v, v, v),
      do_conttrac  = c(0, 0, do_conttrac, do_conttrac, do_conttrac)
    )

    # Figure (not used for paper output)
    # grid of prevalence trends across types of interventions (rows) and
    # transmission scenarios (columns)

    p <- fig_prevalence()
    ggsave(filename = file.path(dir_out, paste0("Prevalence_v", v, ".tiff")),
           plot = p,
           device = "tiff",
           width = 130, height = 150,
           units = "mm", dpi = 300,
           type = "cairo", compression = "lzw")


    # Figure 5
    # compare incidence trends (a) and final size (b) for the interventions:
    # - No intervention
    # - Mass testing & lockdown
    # - Contact tracing at
    #   - 0%
    #   - 30%
    #   - 50%

    pA <- fig_venn()
    pB <- fig_incidence(dt.filter)
    pC <- fig_final_size(dt.filter)

    lay <- matrix(c(1, 3,
                    2, 3), ncol = 2, byrow = TRUE)
    p <- arrangeGrob(pA, pB, pC, layout_matrix = lay,
                     heights = c(1, 1.6), widths = c(1, 1.2))
    ggsave(filename = file.path(dir_out, paste0("FigX_v", v, ".tiff")),
           plot = p,
           device = "tiff",
           width = 183, height = 90,
           units = "mm", dpi = 300,
           type = "cairo", compression = "lzw")

    if (! do.sensitivity_analysis & v == 2) {
      ggsave(filename = file.path(dir_fig, "Figure_5.tiff"),
             plot = p,
             device = "tiff",
             width = 183, height = 90,
             units = "mm", dpi = 300,
             type = "cairo", compression = "lzw")
    }

  }


  # Figure 1

  if (! do.sensitivity_analysis) {
    pA <- fig_flowchart()
    pB <- fig_timeline()

    p <- arrangeGrob(pA, pB, ncol = 1, heights = c(1, 0.3))
    ggsave(filename = file.path(dir_fig, "Figure_1.tiff"),
           plot = p,
           device = "tiff",
           width = 183, height = 200,
           units = "mm", dpi = 300,
           type = "cairo", compression = "lzw")
  }

}
