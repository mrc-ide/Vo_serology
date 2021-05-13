# compute CrI in print format
compute_CrI <- function(X, what = "string") {


  digits <- if (mean(X) < 0.01) 4 else 2

  m  <- round(mean(X), digits)
  q1 <- unname(round(quantile(X, 0.025), digits))
  q2 <- unname(round(quantile(X, 0.975), digits))
  string <- paste0(m, " (", q1, ", ", q2, ")")

  if (what == "m") {
    return(m)
  } else if (what == "q1") {
    return(q1)
  } else if (what == "q2") {
    return(q2)
  } else if (what == "string") {
    return(string)
  }

}

table_fitted <- function (data) {

  # initialise table
  table <- list()

  # loop over posterior files
  files <- list.files(path = dir_output, pattern = "^posterior_best_chain_",
                      full.names = TRUE)
  for (iter in seq_along(files)) {

    # extract data
    dt <- readRDS(files[iter])
    posterior         <- dt$posterior
    looped_parameters <- dt$looped_parameters
    fixed_parameters  <- dt$fixed_parameters

    # compute DIC
    times <- seq(fixed_parameters$tSeed, fixed_parameters$tQ2, 1)
    params <- c(apply(posterior, 2, mean), unlist(looped_parameters), unlist(fixed_parameters))
    sol <- suppressWarnings(model_gen(user = as.list(params))$run(times))
    llike_of_mean <- compute_loglikelihood(sol, data, params)
    DIC <- - 4 * mean(posterior$ll) + 2 * unname(llike_of_mean)

    # compute DIC.median
    times <- seq(fixed_parameters$tSeed, fixed_parameters$tQ2, 1)
    params <- c(apply(posterior, 2, median), unlist(looped_parameters), unlist(fixed_parameters))
    sol <- suppressWarnings(model_gen(user = as.list(params))$run(times))
    llike_of_median <- compute_loglikelihood(sol, data, params)
    DIC.median <- - 4 * mean(posterior$ll) + 2 * unname(llike_of_median)

    # compute dependent parameters
    posterior["inv_gamma"] <- with(as.list(c(fixed_parameters, posterior)),
                                   (generation_time - inv_nu - inv_delta))
    posterior["1-w"] <- 1 - posterior["w"]
    posterior["incubation_period"] <- posterior["inv_nu"] + posterior["inv_delta"]

    # compute 95% credible intervals
    CrI <- apply(posterior, 2, compute_CrI)

    # better presentation
    names(CrI) <- sub("^inv_", "1/", names(CrI))

    # combine
    table[[iter]] <- bind_cols(c(
      looped_parameters, CrI, DIC = round(DIC, 2), DIC.median = round(DIC.median, 2)
    ))

  }

  # clean lists and save to file
  table <- do.call(bind_rows, table) %>%
    arrange(R0_1) %>%
    select(-ll, -w)

  saveRDS(table, file.path(dir_output, "table.rds"))

}


table_final_size <- function () {

  # compute final size (%)
  dt <- readRDS(file.path(dir_output, "SEIR.rds")) %>%
    filter(t == max(t)) %>%
    group_by(intervention, R0_1, id_sample) %>%
    summarise(final_size = 100 * (N - S - S_Q) / N) %>%
    ungroup() %>%
    group_by(intervention, R0_1) %>%
    summarise(final_size = compute_CrI(final_size)) %>%
    ungroup() %>%
    mutate(intervention = factor(
      intervention,
      c("No intervention",
        "MT + CT 100% (baseline)", "MT", "CT 100%", "CT 50%", "CT 30%",
        "Baseline + 0.5 * p_traced", "Baseline + 2 * p_traced",
        "CT 100% + 2 * p_traced", "CT 100% + 4 * p_traced", "CT 100% + 16 * p_traced"))) %>%
    arrange(R0_1, intervention)

  saveRDS(dt, file.path(dir_output, "table_final_size.rds"))

}
