get_paths <- function (pattern, x) {

  list.files(path = dir_output,
             pattern = paste0(x, "-", pattern),
             full.names = TRUE)
}

get_tables <- function (pattern, x) {

  files <- get_paths(pattern = pattern, x = x)

  df <- lapply(seq_along(files), function(i) {
    read.csv(files[i]) %>%
      mutate(iter  = row_number() - 1,
             chain = i)
  })

  return(do.call(bind_rows, df))

}

get_SEIR <- function (intervention,
                      nr_sample, times, iter,
                      post, looped_parameters, fixed_parameters) {

  res <-
    lapply(seq_len(nr_sample), function(i) {
      params <- c(post[i, ], looped_parameters, fixed_parameters)
      sol <- suppressWarnings(model_gen(user = as.list(params))$run(times))
      sol <- data.frame(sol, p = params$p, id_sample = i)
    }) %>%
    do.call(bind_rows, .) %>%
    tibble::add_column(looped_parameters, fixed_parameters, iter, intervention) %>%
    relocate(intervention, R0_1)

  # check disease has died out
  disease_ongoing <- res %>%
    select(E, t, id_sample) %>%
    filter(t %in% c(max(t), max(t) - 1)) %>%
    group_by(id_sample) %>%
    summarise(diff_E = -diff(E)) %>%
    ungroup() %>%
    summarise(upper_quantile = unname(quantile(diff_E, 0.975)) > 0.5) %>%
    pull()

  if (disease_ongoing)
    stop("Disease is ongoing ", iter)

  return(res)
}

# Clean the output of the model run
clean_posteriors <- function (nr_burnin, nr_sample) {

  set.seed(1)

  # prepare for loop read
  SEIR             <- list()
  acceptance_rate  <- list()
  random_walk_rate <- list()

  iter_string <-
    list.files(path = dir_output, pattern = "^fixed_parameters-") %>%
    str_remove("fixed_parameters-") %>%
    str_remove(".csv")

  # loop over results of separate jobs
  for (iter in seq_along(iter_string)) {

    i_string <- iter_string[iter]

    # read in -----------------------------------------------------------------#
    fixed_parameters  <- read.csv(get_paths(paste0(i_string, ".csv"), "fixed_parameters"))
    looped_parameters <- read.csv(get_paths(paste0(i_string, ".csv"), "looped_parameters"))
    posterior        <- get_tables(pattern = paste0(i_string, "-"), x = "posterior")
    acceptance_rate[[iter]] <- get_tables(pattern = paste0(i_string, "-"), x = "acceptance_rate") %>%
      tibble::add_column(looped_parameters)
    random_walk_rate[[iter]] <- get_tables(pattern = paste0(i_string, "-"), x = "random_walk_rate") %>%
      tibble::add_column(looped_parameters)


    # clean posterior ---------------------------------------------------------#

    # remove burnin
    posterior <- posterior %>%
      filter(iter > nr_burnin + 0.5)

    # save all chains
    posterior_all_chains <- list(posterior         = posterior,
                                 fixed_parameters  = fixed_parameters,
                                 looped_parameters = looped_parameters)
    saveRDS(posterior_all_chains, file.path(dir_output, paste0(
      "posterior_all_chains_", i_string, ".rds")))

    # find best chain and sample
    posterior <- posterior %>%
      group_by(chain) %>%
      mutate(mean_ll = mean(ll)) %>%
      ungroup() %>%
      filter(mean_ll == max(mean_ll)) %>%
      select(-iter, -chain, -mean_ll)

    # save best chain
    posterior_best_chain <- list(posterior         = posterior,
                                 fixed_parameters  = fixed_parameters,
                                 looped_parameters = looped_parameters)
    saveRDS(posterior_best_chain, file.path(dir_output, paste0(
      "posterior_best_chain_", i_string, ".rds")))

    # sample best posterior and compute SEIR models ---------------------------#

    times <- seq(fixed_parameters$tSeed, fixed_parameters$tf, 1L)
    id_sample <- sample(seq_len(nrow(posterior)), nr_sample)
    SEIR_ <- list()

    # baseline: MT + CT (100%)
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    SEIR_[[1]] <- get_SEIR(intervention = "MT + CT 100% (baseline)",
                           nr_sample, times, iter,
                           post, looped_parameters, fp)

    # No intervention
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    post$w <- 1
    post$p_traced <- 0
    post$p_traced_S <- 0
    SEIR_[[2]] <- get_SEIR(intervention = "No intervention",
                           nr_sample, times, iter,
                           post, looped_parameters, fp)

    # MT
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    post$p_traced <- 0
    post$p_traced_S <- 0
    SEIR_[[3]] <- get_SEIR(intervention = "MT",
                           nr_sample, times, iter,
                           post, looped_parameters, fp)

    # CT 100%
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    post$w <- 1
    SEIR_[[4]] <- get_SEIR(intervention = "CT 100%",
                           nr_sample, times, iter,
                           post, looped_parameters, fp)

    # CT 50%
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    post$w <- 1
    fp$q_Q <- 0.5
    SEIR_[[5]] <- get_SEIR(intervention = "CT 50%",
                           nr_sample, times, iter,
                           post, looped_parameters, fp)

    # CT 30%
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    post$w <- 1
    fp$q_Q <- 0.7
    SEIR_[[6]] <- get_SEIR(intervention = "CT 30%",
                           nr_sample, times, iter,
                           post, looped_parameters, fp)

    # Sensitivity analysis: Baseline + 0.5 * p_traced
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    post$p_traced <- 0.5 * post$p_traced
    post$p_traced_S <- 0.5 * post$p_traced_S
    SEIR_[[7]] <- get_SEIR(intervention = "Baseline + 0.5 * p_traced",
                           nr_sample, times, iter,
                           post, looped_parameters, fp)

    # Sensitivity analysis: Baseline + 2 * p_traced
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    post$p_traced <- 2 * post$p_traced
    post$p_traced_S <- 2 * post$p_traced_S
    SEIR_[[8]] <- get_SEIR(intervention = "Baseline + 2 * p_traced",
                           nr_sample, times, iter,
                           post, looped_parameters, fp)

    # Sensitivity analysis: CT 100% + 2 * p_traced
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    post$w <- 1
    post$p_traced <- 2 * post$p_traced
    post$p_traced_S <- 2 * post$p_traced_S
    SEIR_[[9]] <- get_SEIR(intervention = "CT 100% + 2 * p_traced",
                           nr_sample, times, iter,
                           post, looped_parameters, fp)

    # Sensitivity analysis: CT 100% + 4 * p_traced
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    post$w <- 1
    post$p_traced <- 4 * post$p_traced
    post$p_traced_S <- 4 * post$p_traced_S
    SEIR_[[10]] <- get_SEIR(intervention = "CT 100% + 4 * p_traced",
                            nr_sample, times, iter,
                            post, looped_parameters, fp)

    # Sensitivity analysis: CT 100% + 16 * p_traced
    post <- posterior[id_sample, ]
    fp <- fixed_parameters
    post$w <- 1
    post$p_traced <- 16 * post$p_traced
    post$p_traced_S <- 16 * post$p_traced_S
    SEIR_[[11]] <- get_SEIR(intervention = "CT 100% + 16 * p_traced",
                            nr_sample, times, iter,
                            post, looped_parameters, fp)


    SEIR[[iter]] <- SEIR_ %>% bind_rows()

  }

  # clean lists and save to file
  SEIR             <- do.call(bind_rows, SEIR)
  acceptance_rate  <- do.call(bind_rows, acceptance_rate)
  random_walk_rate <- do.call(bind_rows, random_walk_rate)

  saveRDS(SEIR,             file.path(dir_output, "SEIR.rds"))
  saveRDS(acceptance_rate,  file.path(dir_output, "acceptance_rate.rds"))
  saveRDS(random_walk_rate, file.path(dir_output, "random_walk_rate.rds"))

}
