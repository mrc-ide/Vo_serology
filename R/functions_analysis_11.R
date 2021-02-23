compute_interventions <- function (do.sensitivity_analysis) {
  
  set.seed(1)
  
  # read in parameters from Nature paper (some were fixed, some fitted)
  params <- readRDS(file.path(dir_in, "posterior_best_chain_14.rds")) %>%
    bind_cols() %>%
    # remove superfluous parameters
    select(-X, -ll, -inv_gamma, -new_beta, -`1/sigma`) %>%
    # sample 100 combinations
    slice_sample(n = 100) 
  
  
  # sensitivity analysis:
  # anticipate start of contact tracing from 24th to 22nd February 2020
  if (do.sensitivity_analysis)
    params <- params %>% mutate(tQ = 18)
  
  
  SEIR <- list()
  times <- seq(unique(params$tSeed), 290, 1)
  
  
  # compute different scenarios -----------------------------------------------#
  
  # compute SEIR without interventions
  SEIR[[1]] <- params %>%
    mutate(do_lockdown = 0,
           do_conttrac = 0,
           v           = 0L,
           do_q        = 0) %>%
    get_SEIR(times) %>%
    is_disease_ongoing() %>%
    bind_cols(Intervention = "No intervention",
              Quarantined = "no transmission assumed")
  
  # compute SEIR with lockdown
  SEIR[[2]] <- params %>%
    mutate(do_lockdown = 1,
           do_conttrac = 0,
           v           = 0L,
           do_q        = 0) %>%
    get_SEIR(times) %>%
    is_disease_ongoing() %>%
    bind_cols(Intervention = "Mass testing & lockdown",
              Quarantined = "no transmission assumed")
  
  # compute SEIR with contact tracing, quarantined don't transmit
  SEIR[[3]] <- params %>%
    crossing(do_lockdown = 0,
             do_conttrac = c(0.06, 0.07, 0.08, 0.09, seq(0.05, 0.5, by = 0.05)),
             v           = c(1L, 2L, 3L),
             do_q        = 0) %>%
    get_SEIR(times) %>%
    is_disease_ongoing() %>%
    bind_cols(Intervention = "Contact tracing",
              Quarantined = "no transmission assumed")
  
  # compute SEIR with contact tracing, quarantined transmit at 30%
  SEIR[[4]] <- params %>%
    crossing(do_lockdown = 0,
             do_conttrac = c(0.06, 0.07, 0.08, 0.09, seq(0.05, 0.5, by = 0.05)),
             v           = c(1L, 2L, 3L),
             do_q        = 0.30) %>%
    get_SEIR(times) %>%
    is_disease_ongoing() %>%
    bind_cols(Intervention = "Contact tracing",
              Quarantined = "30% transmission assumed")
  
  # compute SEIR with contact tracing, quarantined transmit at 50%
  SEIR[[5]] <- params %>%
    crossing(do_lockdown = 0,
             do_conttrac = c(0.06, 0.07, 0.08, 0.09, seq(0.05, 0.5, by = 0.05)),
             v           = c(1L, 2L, 3L),
             do_q        = 0.5) %>%
    get_SEIR(times) %>%
    is_disease_ongoing() %>%
    bind_cols(Intervention = "Contact tracing",
              Quarantined = "50% transmission assumed")
  
  
  # save to file --------------------------------------------------------------#
  
  SEIR <- bind_rows(SEIR) %>%
    mutate(Intervention = factor(
      Intervention, c("No intervention", "Mass testing & lockdown", "Contact tracing")))
  
  SEIR %>% saveRDS(file.path(dir_out, "SEIR.rds"))
  
}

# compute SEIR over time from params
get_SEIR <- function (df, times) {
  
  lapply(seq_len(nrow(df)), 
         function(i) {
           
           model_gen(user = as.list(df[i, ]))$run(times) %>%
             suppressWarnings() %>%
             bind_cols(id_sample = i, df[i, ])
           
         }) %>%
    bind_rows()
  
}

# check disease has effectively died out
is_disease_ongoing <- function (df) {
  
  disease_ongoing <- df %>%
    select(S, t, id_sample) %>%
    filter(t %in% c(max(t), max(t) - 1)) %>%
    group_by(id_sample) %>%
    summarise(diff_S = -diff(S)) %>%
    ungroup() %>%
    summarise(upper_quantile = unname(quantile(diff_S, 0.975)) > 0.5) %>%
    pull()
  
  if (disease_ongoing)
    stop("Disease is ongoing!")
  
  return(df)
  
}
