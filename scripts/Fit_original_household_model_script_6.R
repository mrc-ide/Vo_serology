################################################################################
# This code implements the method described in Fraser et al (2011), Am J       #
# Epidemiol (https://doi.org/10.1093/aje/kwr122). Here we fit the model to the #
# SARS-CoV-2 household data collected in Vo'.                                  #
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

cat("\n\n### Running 'scripts/Fit_original_household_model_script_6.R' ###\n\n")


# source scripts
source("R/functions_fit_model.R")
source("R/functions_plot_fitted_model.R")

# # connect to the cluster ----------------------------------------------------#
#
# options(
#   didehpc.username = "idorigat",
#   didehpc.cluster = "fi--didemrchnb")
#
# root <- "contexts"
#
# # load packages
# list_packages <- c("gdata", "data.table")
#
# # load sources
# list_sources <- c("R/functions_fit_model.R",
#                   "R/functions_plot_fitted_model.R")
#
# ctx <- context::context_save(root,
#                              packages = list_packages,
#                              sources = list_sources)
#
# # Create a queue within the context (/environment)
# obj <- didehpc::queue_didehpc(ctx)

# read data -------------------------------------------------------------------#
data <- read.csv("tables/Household_final_size_Vo_baseline.csv", header = TRUE)
data <- data[, -1]
colnames(data) <- seq(1,7,1)
rownames(data) <- seq(0,4,1)

# read defult parameters
parameters <- read.csv("data/household_model/Default_parameters_original.csv",
                       header = TRUE)

# read model variants
par <- read.csv("data/household_model/Model_variants_original.csv",
                header = TRUE)

# read MCMC parameters
mcmc_par <- read.csv("data/household_model/MCMC_parameters_original.csv",
                     header = TRUE)

# Create output folders -------------------------------------------------------#
dir_output  <- file.path("mcmc_posterior_chains_original")
dir_figures <- file.path(dir_output, "figures")
dir.create(dir_output,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

# Create output file ----------------------------------------------------------#

df <- data.frame(model = character(0),
                 chain = numeric(0),
                 estimate = character(0),
                 beta = numeric(0),
                 Q = numeric(0),
                 alpha = numeric(0),
                 shape = numeric(0),
                 Qprior = numeric(0),
                 p_pr = numeric(0),
                 p_srev = numeric(0),
                 p_com = numeric(0),
                 DIC = numeric(0),
                 deviance = numeric(0),
                 llike = numeric(0))

file_res <- "mcmc_posterior_estimates_original.csv"
write.csv(df, file.path(dir_output, file_res),
          row.names = FALSE)

# output plot of single chains
single_chains_plots <- FALSE

# Run MCMC --------------------------------------------------------------------#

# set number of chains to run for each model and burnin
id_chain <- 3
burnin <- 100

sizes <- seq_len(nrow(par))

# # runs on the cluster
# t <- obj$enqueue(sessionInfo())
# grp <- obj$lapply(sizes,
#                   wrapper_model,
#                   data                 = data,
#                   par                  = par,
#                   mcmc_par             = mcmc_par,
#                   parameters           = unlist(parameters),
#                   parameterNames       = colnames(parameters),
#                   dir_output           = dir_output,
#                   filename             = file_res,
#                   id_chain             = id_chain,
#                   burnin               = burnin)
#
# grp$status()
# res <- grp$wait(Inf)

# runs locally
burnin <- 100

#mcmc_par$nbIter <- 1000 # uncomment for a quicker test run
#burnin <- 20 # uncomment for a quicker test run

res <- lapply(sizes,
              wrapper_model,
              data                 = data,
              par                  = par,
              mcmc_par             = mcmc_par,
              parameters           = parameters,
              parameterNames       = colnames(parameters),
              dir_output           = dir_output,
              filename             = file_res,
              id_chain             = id_chain,
              burnin               = burnin)

# save outputs
height_plot <- c(120,120,180,180,240,240,300,300)/120
parameters_default <- parameters

for(j in 1:length(res)){

  for(i in 1:length(res[[j]][[1]])){

    parametersToEstimate <- unlist(par[j, 2:9])

    ## output posterior chains
    save_posterior(burnin,
                   parametersToEstimate,
                   colnames(parameters),
                   res[[j]][[1]][[i]],
                   par[j,1],
                   i,
                   dir_output)

    if(single_chains_plots == TRUE){

      ## plot posterior chains (single chains)
      v <- length(which(parametersToEstimate == 1))
      ggsave(filename =
              file.path(file.path(dir_output, "figures"),
                        paste(paste(par[j, 1], i, sep= "_"),
                              ".tiff", sep="")),
            plot = plot_chains(burnin,
                               parametersToEstimate,
                               colnames(parameters),
                               res[[j]][[1]][[i]]),
            device = "tiff",
            width = 240, height = height_plot[v]*120,
            units = "mm", dpi = 300, limitsize = TRUE)

    }else{

      ## plot posterior chains (overlaid, all in a plot)
      v <- length(which(parametersToEstimate == 1))
      ggsave(filename =
               file.path(file.path(dir_output, "figures"),
                         paste(par[j, 1],".tiff", sep="")),
             plot = plot_multiple_chains(burnin,
                                         parametersToEstimate,
                                         colnames(parameters),
                                         res[[j]][[1]]),
             device = "tiff",
             width = 240, height = height_plot[v]*120,
             units = "mm", dpi = 300, limitsize = TRUE)

    }

    # compute posterior estimates
    df <- compute_post_estimates(burnin,
                                  as.numeric(unlist(parameters_default)),
                                  parametersToEstimate,
                                  colnames(parameters),
                                  res[[j]][[1]][[i]])

    # reformat output
    vec <- vector()
    for(k in 1:5){

      v1 <- rep("-", 8)
      v1[which(parametersToEstimate == 1)] <- unlist(df[[1]][k])

      vec <- c(vec,
               as.character(par[j, 1]),
               i,
               colnames(df[[1]][k]),
               v1,
               df[[length(df)-2]]) # DIC
    }

    dev <- c(df[[length(df)-3]], "-") # deviance
    llik <- c(df[[length(df)]], "-") # llike

    df_res <- matrix(vec, nrow = 5, byrow = TRUE)

    df_res1 <- matrix(0, nrow = 5, ncol = dim(df_res)[2]+2)
    df_res1[,1:dim(df_res)[2]] <- df_res
    df_res1[, dim(df_res)[2]+1] <- dev
    df_res1[, dim(df_res)[2]+2] <- llik

    fwrite(df_res1, file.path(dir_output, file_res), append = TRUE,
           row.names = FALSE, quote = FALSE)

  }

}

################################################################################
################################################################################
