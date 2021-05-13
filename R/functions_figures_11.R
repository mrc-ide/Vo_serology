theme_set(
  theme_bw() +
    theme(
      # Nature requirement:
      text = element_text(size = 8, family = "sans"),
      plot.tag = element_text(size = 8, face = "bold", family = "sans"),
      legend.text = element_text(size = 6))
)

fig_chains <- function () {

  # loop over posterior files
  files <- list.files(path = dir_output, pattern = "^posterior_all_chains_",
                      full.names = TRUE)
  for (iter in seq_along(files)) {

    # extract data
    dt <- readRDS(files[iter])
    posterior         <- dt$posterior
    looped_parameters <- dt$looped_parameters

    # melt data
    posterior <- gather(posterior, "variable", "value", -iter, -chain) %>%
      mutate(chain    = as.factor(chain),
             variable = as.factor(variable) %>% relevel("ll"))

    # title
    title <- paste("Intervention = MT + CT 100%, R0_1 =", looped_parameters["R0_1"])

    p <- ggplot(data = posterior, aes(x = iter, y = value,
                                      group = chain, colour = chain)) +
      geom_line() +
      facet_wrap(vars(variable), scales = "free_y", ncol = 2) +
      ggtitle(title)

    ggsave(filename = file.path(dir_output, paste0("mcmc_chains_", iter, ".png")),
           plot = p, device = "png",
           width = 17, height = 20, units = "cm", dpi = 300, limitsize = TRUE)

  }

}

fig_acceptance_rates <- function () {

  # read in
  dt <- readRDS(file.path(dir_output, "acceptance_rate.rds"))

  # melt
  dt <- gather(dt, "variable", "value", -R0_1, -chain)
  dt$chain <- as.factor(as.integer(dt$chain))

  p <- ggplot(data = dt, aes(x = value, fill = chain)) +
    geom_histogram(alpha = 0.5, position = "identity") +
    facet_wrap(vars(variable))

  ggsave(filename = file.path(dir_output, "acceptance_rate.png"),
         plot = p, device = "png",
         width = 17, height = 20, units = "cm", dpi = 300, limitsize = TRUE)

}


fig_flowchart <- function () {

  x0 <- 0
  x1 <- 2.5
  x2 <- 5
  x3 <- 9
  x4 <- 13
  x5 <- 15
  x6 <- 19

  y0 <- 34.3
  y1 <- 28.5
  y2 <- 23.5
  y3 <- 18
  y4 <- 16
  y5 <- 12
  y6 <- 5.5
  y7 <- 2.2
  y8 <- -1.5

  dt.text <- data.frame(X = c(x4-1, x3+0.7),
                        Y = c(y0-2.8, y6+3.2),
                        text = c('First serosurvey 1-3 May 2020',
                                 'Second serosurvey 28-29 November 2020'))

  dt.labs <- data.frame(X = c(x0, x1, x3+0.5, x3, x3, x5, x4, x1, x1, x3, x5, x3, x1, x5),
                        Y = c(y0, y2, y1, y2, y3, y3, y5, y4, y6, y6, y6, y8, y8, y8),
                        text = c("Vo' cluster\nn = 3,329",
                                 'Tested\nn = 2,602\n(78.2%)',
                                 'Tested with a single assay,\nnegative result\nn = 137 (5.3%)',
                                 'Negative to\nall three assays*\nn = 2,303 (88.5%)',
                                 'Positive to at\nleast one assay\nn = 162 (6.2%)',
                                 'Neutralising antibody\ntitres > 1:40\nn = 44 (1.7%)',
                                 'Drop out\nn = 24 (0.9%)',
                                 'Not tested\nn = 727 (21.8%)',
                                 'Positive swab\nn = 8 (0.3%)',
                                 'Tested\nn = 156 (6.0%)',
                                 'Positive swab\nor T-cell-positive\nn = 10 (0.4%)',
                                 'Positive to at\nleast one assay\nn = 129 (5.0%)',
                                 'Neutralising antibody\ntitres > 1:40\nn = 26 (1.0%)',
                                 'Negative to all three assays\nn = 25 (1.0%)'))

  dt.rect <- data.frame(X1 = c(x0-1, x0-1), Y1 = c(y8-3, y5+1.5), # lower left corners
                        X2 = c(x6-1, x6-1), Y2 = c(y6+4, y0-2)) # upper right corners

  dt.segments <- data.frame(X1 = c(x0,     x2, x3+3.5, x3+2.5, x6, x3),
                            X2 = c(x0,     x2, x6,     x6,     x6, x5),
                            Y1 = c(y0-1.4, y1, y1,     y2,     y1, y7),
                            Y2 = c(y4,     y3, y1,     y2,     y6, y7))

  dt.arrows.topdown <- data.frame(X1 = c(x3, x1, x3, x5),
                                  X2 = c(x3, x1, x3, x5),
                                  Y1 = c(y3-2, y4-2, y6-2, y7),
                                  Y2 = c(y6+2, y6+2, y8+2, y8+2))
  dt.arrows.left2right <- data.frame(X1 = c(x0,     x0,     x1+1.5, x2,     x2,     x3+2, x3,   x1+2),
                                     X2 = c(x1-1.5, x1-1.5, x3-2.5, x3-2.5, x3-2.5, x5-2, x4-2, x3-1.5),
                                     Y1 = c(y2,     y4,     y2,     y1,     y3,     y3,   y5,   y6),
                                     Y2 = c(y2,     y4,     y2,     y1,     y3,     y3,   y5,   y6))
  dt.arrows.right2left <- data.frame(X1 = c(x6,     x5-2,   x3-2),
                                     X2 = c(x5+2, x3+1.5, x1+2),
                                     Y1 = c(y6,     y6,     y8),
                                     Y2 = c(y6,     y6,     y8))

  p <- ggplot() +
    geom_rect(data = dt.rect, aes(xmin = X1, ymin = Y1, xmax = X2, ymax = Y2),
              linetype = "dashed", colour = "black", fill = NA, size = 0.3) +
    geom_segment(data = dt.segments, aes(x = X1, xend = X2, y = Y1, yend = Y2)) +
    geom_segment(data = dt.arrows.topdown, aes(x = X1, xend = X2, y = Y1, yend = Y2),
                 lineend = "round", linejoin = "round",
                 arrow = arrow(length = unit(0.08, "inches"), angle = 20)) +
    geom_segment(data = dt.arrows.left2right, aes(x = X1, xend = X2, y = Y1, yend = Y2),
                 lineend = "round", linejoin = "round",
                 arrow = arrow(length = unit(0.08, "inches"), angle = 20)) +
    geom_segment(data = dt.arrows.right2left, aes(x = X1, xend = X2, y = Y1, yend = Y2),
                 lineend = "round", linejoin = "round",
                 arrow = arrow(length = unit(0.08, "inches"), angle = 20)) +
    geom_text(data = dt.text, aes(x = X, y = Y, label = text), size = 3.5, hjust = 0) +
    geom_text(data = dt.labs, aes(x = X, y = Y, label = text), size = 3) +
    xlim(x0-1, x6) +
    ylim(y8-3, y0) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0),
          text = element_text(size = 10, family = "sans"),
          plot.tag = element_text(size = 14, face = "bold", family = "sans")) +
    labs(tag = "a")

  return(p)

}

fig_timeline <- function () {

  dt <- data.frame(
    position = c(5, 16, 22, 28.5),
    start    = c(1, NA, 21, 28),
    end      = c(9, NA, 23, 29),
    event = c("First oro-nasopharyngeal survey", "Second oro-nasopharyngeal survey",
              "First serosurvey and\noro-nasopharyngeal survey",
              "Second serosurvey and\noro-nasopharyngeal survey")
  )

  dt.axis <- data.table::fread(
    "position, date
     1,        21 February
     9,        29 February
     16,       7 March
     21,       1 May
     23,       3 May
     28,       28 November
     29,       29 November"
  )

  p <- ggplot(data = dt) +
    # labels
    geom_text_repel(data = dt[c(2), ],
                    aes(x = position, label = event),
                    y = 0, ylim = c(1.2, 2.4), size = 3.3) +
    geom_text_repel(data = dt[c(1,3), ],
                    aes(x = position, label = event),
                    y = 0, ylim = c(3, NA), size = 3.3) +
    geom_text_repel(data = dt[c(4), ],
                    aes(x = position, label = event),
                    y = 0, ylim = c(0.8, 1.4), size = 3.3) +
    # dotted line throughout
    geom_segment(data = dt.axis,
                 aes(x = min(position) - 1, xend = max(position) + 1.8),
                 y = 0, yend = 0, linetype = "dotted") +
    # first continuous segment (February-March)
    geom_segment(data = dt.axis,
                 aes(x = min(position) - 1, xend = position[3] + 1.5),
                 y = 0, yend = 0) +
    # second continuous segment (May)
    geom_segment(data = dt.axis,
                 aes(x = position[4] - 1.5, xend = position[5] + 1.5),
                 y = 0, yend = 0) +
    # third continuous segment (November)
    geom_segment(data = dt.axis,
                 aes(x = position[6] - 1.5, xend = max(position) + 1.8),
                 y = 0, yend = 0) +
    # data points
    geom_point(data = dt %>% filter(is.na(start)),
               aes(x = position),
               y = 0) +
    geom_point(data = dt.axis,
               aes(x = position),
               y = 0, pch = 3) +
    # coloured segments
    geom_segment(aes(x = start, xend = end, colour = as.factor(start)),
                 y = 0.04, yend = 0.04, size = 5,
                 show.legend = FALSE) +
    # axes
    ylim(-0.5, 5) +
    scale_x_continuous(breaks = dt.axis$position,
                       labels = dt.axis$date) +
    theme_void() +
    theme(plot.margin = margin(-10,0,35,0),
          axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t=-20)),
          legend.position = "none",
          text = element_text(size = 10, family = "sans"),
          plot.tag = element_text(size = 14, face = "bold", family = "sans")) +
    labs(tag = "b")

  return(p)

}


fig_venn <- function (panel) {

  dt.venn <- tibble(
    value = c(   56,   146,   44),
    Set1  = c( TRUE, FALSE, TRUE),
    Set2  = c(FALSE,  TRUE, TRUE)
  )

  dt.labels <- data.frame(
    # x = c(-2,   0, 2),
    # y = c(-1.3, 1.45, -1.3),
    x = c(-2.5,   0, 2.3),
    y = c(0, 1.45, 0),
    label = c("PCR positive\nsubjects", "PCR positive named contacts", "Traced\ncontacts")
  )

  p <- ggplot() +
    geom_venn(data = dt.venn,
              mapping = aes(A = Set1, B = Set2, label = value),
              set_names = "", text_size = 2.5,
              show_percentage = FALSE, stroke_color = NA, fill_color = c("red", "blue")) +
    xlim(-3.1, 3) +
    ylim(-1, 1.5) +
    geom_text(data = dt.labels,
              mapping = aes(x = x, y = y, label = label),
              size = 2.5) +
    geom_segment(aes(x = 0, xend = 0, y = 1.2, yend = 0.2)) +
    coord_fixed() +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0),
          text = element_text(size = 3.5, family = "sans"),
          # text = element_text(size = 8, family = "sans"),
          plot.tag = element_text(size = 8, face = "bold", family = "sans"),
    ) +
    labs(tag = panel)

  p

  return(p)

}

fig_contact_tracing <- function (panel) {

  # get compartment counts
  dt <- readRDS(file.path(dir_output, "SEIR.rds"))

  # get relevant times
  times <- dt %>%
    select(tQ2, N) %>%
    unique()

  # clean data
  dt <- dt %>%
    filter(intervention == "MT + CT 100% (baseline)" & t == times$tQ2 &
             R0_1 > 2.3 & R0_1 < 2.5) %>%
    # compute proportion of traced contacts that have always been negative
    mutate(p_negative_traced = S_Q / (S_Q + E_Q + TPpre_Q + I_Q + TPpost_Q + TN_Q)) %>%
    mutate(p_positive_traced = 100 * (1 - p_negative_traced)) %>%
    # compute mean and 95% CrI
    summarise(column = "Model",
              m = mean(p_positive_traced),
              q1 = quantile(p_positive_traced, probs = 0.025),
              q2 = quantile(p_positive_traced, probs = 0.975))

  # Contact tracing data
  dt.data <- data2 %>%
    transmute(column = "Data",
              m  = 100 * binom.test(traced_contacts - negative_traced_contacts, traced_contacts)$estimate,
              q1 = 100 * binom.test(traced_contacts - negative_traced_contacts, traced_contacts)$conf.int[1],
              q2 = 100 * binom.test(traced_contacts - negative_traced_contacts, traced_contacts)$conf.int[2])

  # combine data.frames
  dt <- bind_rows(dt, dt.data)

  # Plot
  p <- ggplot(data = dt, aes(x = column,
                             y = m, ymin = q1, ymax = q2,
                             colour = column)) +
    # model mean and 95% CrI
    geom_point(size = 1) +
    geom_errorbar(width = 0.08) +
    # axes
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Traced contacts\ntesting positive (%)",
                       limits = c(2, 98),
                       breaks = c(0, 20, 40, 60, 80, 100)) +
    # colour palette and legend
    scale_colour_brewer("", palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.margin = margin(0, 5, 0, 0)) +
    labs(tag = panel)

  return(p)

}

fig_PCR_testing <- function (panel) {

  # get compartment counts
  dt <- readRDS(file.path(dir_output, "SEIR.rds"))

  # get relevant times
  times <- dt %>%
    select(tQ2, N) %>%
    unique()

  # clean data
  dt <- dt %>%
    filter(intervention == "MT + CT 100% (baseline)" & t == times$tQ2 &
             R0_1 > 2.3 & R0_1 < 2.5) %>%
    # compute proportion of traced contacts that have always been negative
    mutate(temp1 = E_Q + TPpre_Q + I_Q       + TPpost_Q            + TN_Q,
           temp2 = E   + TPpre   + I_A + I_S + TPpost_A + TPpost_S + TN) %>%
    mutate(p_PCR_not_traced = 100 * (1 - (temp1 / (temp1 + temp2)))) %>%
    # compute mean and 95% CrI
    summarise(column = "Model",
              m = mean(p_PCR_not_traced),
              q1 = quantile(p_PCR_not_traced, probs = 0.025),
              q2 = quantile(p_PCR_not_traced, probs = 0.975))

  # Contact tracing data
  dt.data <- data3 %>%
    transmute(column = "Data",
              m  = 100 * binom.test(PCR_positive_subjects - PCR_positive_not_traced, PCR_positive_subjects)$estimate,
              q1 = 100 * binom.test(PCR_positive_subjects - PCR_positive_not_traced, PCR_positive_subjects)$conf.int[1],
              q2 = 100 * binom.test(PCR_positive_subjects - PCR_positive_not_traced, PCR_positive_subjects)$conf.int[2])

  # combine data.frames
  dt <- bind_rows(dt, dt.data)

  # Plot
  p <- ggplot(data = dt, aes(x = column,
                             y = m, ymin = q1, ymax = q2,
                             colour = column)) +
    # model mean and 95% CrI
    geom_point(size = 1) +
    geom_errorbar(width = 0.08) +
    # axes
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "PCR positive subjects\ndetected by CT (%)",
                       limits = c(2, 98),
                       breaks = c(0, 20, 40, 60, 80, 100)) +
    # colour palette and legend
    scale_colour_brewer("", palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.margin = margin(0, 5, 0, 0)) +
    labs(tag = panel)

  return(p)

}

fig_prevalence <- function (panel) {


  # get compartment counts
  dt <- readRDS(file.path(dir_output, "SEIR.rds"))

  # get relevant times
  times <- dt %>%
    select(tSeed, time1, time2, tQ1, tQ2, tf, N) %>%
    unique()

  dt <- dt %>%
    filter(intervention == "MT + CT 100% (baseline)" & t <= times$tQ2 &
             R0_1 > 2.3 & R0_1 < 2.5) %>%
    # Compute prevalence by pre-symptomatic, symptomatic and asymptomatic
    transmute(t, R0_1,
              Presymptomatic = (1 - p) * TPpre,
              Symptomatic    = I_S + TPpost_S,
              Asymptomatic   = p * TPpre + I_A + TPpost_A) %>%
    pivot_longer(c(Presymptomatic, Symptomatic, Asymptomatic),
                 names_to = "Compartment", values_to = "value") %>%
    # compute mean and 95% CrI
    group_by(t, R0_1, Compartment) %>%
    mutate(value = 100 * value / times$N) %>%
    summarise(mean = mean(value),
              low  = quantile(value, probs = 0.025),
              high = quantile(value, probs = 0.975)) %>%
    ungroup() %>%
    # column names for plotting
    mutate(R0_1  = paste("R0 =", R0_1),
           Compartment = sub("Presym", "Pre-sym", Compartment) %>%
             factor(c("Pre-symptomatic", "Symptomatic", "Asymptomatic")))

  # Screening data
  data$Time <- c(times$time1, times$time2)
  my.data <- data %>%
    gather("Compartment", "value", -Time, -Tested) %>%
    group_by(Time, Tested, Compartment) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>%                           # don't remove this and following line
    group_by(Time, Tested, Compartment) %>% # else result won't be correct!
    mutate(mean  = 100 * binom.test(value, Tested)$estimate,
           lower = 100 * binom.test(value, Tested)$conf.int[1],
           upper = 100 * binom.test(value, Tested)$conf.int[2]) %>%
    ungroup() %>%
    # column names for plotting
    mutate(Compartment = sub("Presym", "Pre-sym", Compartment) %>%
             factor(c("Pre-symptomatic", "Symptomatic", "Asymptomatic")))

  # Dates for x axis
  my.times <- times[c("tSeed", "tQ1", "time2", "tQ2")]
  my.dates <- get_dates(my.times)

  # Plot
  dotsize <- 1
  errorsize <- 1.2
  dashsize <- .3

  p <- ggplot(data = dt,
              aes(x = t, colour = Compartment, fill = Compartment)) +
    # model mean and 95% CrI
    geom_line(aes(y = mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.4, colour = NA) +
    # screening data
    geom_errorbar(data = my.data, aes(x = Time, ymin = lower, ymax = upper),
                  width = errorsize) +
    geom_point(data = my.data, aes(x = Time, y = mean), size = dotsize) +
    # lockdown date
    geom_vline(xintercept = as.numeric(my.times["tQ1"]),
               linetype = "dashed", size = dashsize) +
    # colour palette and legend
    scale_colour_brewer("", palette = "Dark2") +
    scale_fill_brewer("", palette = "Dark2") +
    scale_x_continuous(name = "Date",
                       breaks = unlist(my.times),
                       minor_breaks = seq(min(my.times), max(my.times), 7),
                       labels = my.dates[which(my.dates != "")]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-8, 0, 0, 0),
          # Nature requirement:
          text = element_text(size = 8, family = "sans"),
          plot.tag = element_text(size = 8, face = "bold", family = "sans")) +
    labs(y = "Prevalence (%)", tag = panel)

  return(p)

}

fig_incidence <- function (panel) {

  # get compartment counts
  dt <- readRDS(file.path(dir_output, "SEIR.rds"))

  # get relevant times
  times <- dt %>%
    select(tSeed, time1, time2, tQ1, tQ2, tf, N) %>%
    unique()

  # clean data
  dt <- dt %>%
    filter(intervention %in% c("MT + CT 100% (baseline)",
                               "MT", "CT 100%", "No intervention") &
             R0_1 > 2.3 & R0_1 < 2.5) %>%
    # compute incidence
    mutate(Scombo = S + S_Q) %>%
    arrange(t) %>%
    group_by(R0_1, id_sample, intervention) %>%
    transmute( # (what left S(t-1)) - (what entered S_Q(t))
      t, Scombo,
      incid = (lag(Scombo, default = Scombo[1]) - Scombo)
    ) %>%
    ungroup() %>%
    # compute mean and 95% CrI
    group_by(R0_1, intervention, t) %>%
    mutate(incid = 100 * incid / times$N) %>%
    summarise(mean = mean(incid),
              low  = quantile(incid, probs = 0.025),
              high = quantile(incid, probs = 0.975)) %>%
    ungroup() %>%
    mutate(intervention = if_else(intervention == "No intervention", NA_character_, intervention)) %>%
    # clean
    mutate(intervention = factor(intervention,
                                 c("MT + CT 100% (baseline)",
                                   "MT", "CT 100%", "No intervention")))

  # find time of last infection
  time_final <- dt %>%
    filter(mean > 0.01) %>%
    summarise(max(t)) %>%
    pull()
  time_final <- min(time_final, 70)
  time_final <- max(time_final, times$time2)
  dt <- dt %>% filter(t < time_final + 0.5)

  # legend labels
  legend.intv <- c("MT + CT", "MT", "CT", "No intervention")

  # Dates for x axis
  if (time_final > 70) {
    my.times <- c(unlist(times[c("tSeed", "tQ1", "time2")]), 42, 56, 70)
  } else if (time_final > times$time2) {
    my.times <- c(unlist(times[c("tSeed", "tQ1", "time2")]), time_final)
    ## add a date in the last long interval
    # my.times <- sort(c(my.times, as.integer(mean(rev(my.times)[1:2]))))
    my.times <- sort(c(my.times, 42, 56, 70))
  } else {
    my.times <- times[c("tSeed", "tQ1", "time2")]
  }
  my.dates <- get_dates(unname(my.times))

  # Plot
  p <- ggplot() +
    # mean and 95% CrI
    geom_line(data = dt,
              aes(x = t, y = mean, colour = intervention)) +
    geom_ribbon(data = dt,
                aes(x = t, ymin = low, ymax = high, fill = intervention),
                alpha = 0.4, colour = NA) +
    # lockdown date
    geom_vline(xintercept = as.numeric(my.times["tQ1"]),
               linetype = "dashed", size = .3) +
    # colour palette and legend
    scale_colour_brewer("", palette = "Set1", labels = legend.intv, na.value = "black") +
    scale_fill_brewer("", palette = "Set1", labels = legend.intv, na.value = "black") +
    scale_y_continuous(name = "Incidence (%)   ",
                       breaks = 0:5) +
    scale_x_continuous(name = "Date",
                       breaks = my.times,
                       minor_breaks = seq(min(my.times), max(my.times), 7),
                       labels = my.dates[which(my.dates != "")]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-8, 0, 0, 0),
          # Nature requirement:
          text = element_text(size = 8, family = "sans"),
          plot.tag = element_text(size = 8, face = "bold", family = "sans")) +
    labs(tag = panel)

  return(p)

}

fig_relative_final_size <- function (panel) {

  if (panel == "f") {
    intv <- c("MT + CT 100% (baseline)", "MT", "CT 100%")
    legintv <- c("MT + CT", "MT", "CT")
    legtitle <- ""
    pd <- 0.8
    wd <- 0.5
  } else if (panel == "g") {
    intv <- c("MT + CT 100% (baseline)", "MT", "Baseline + 0.5 * p_traced", "Baseline + 2 * p_traced")
    legintv <- c("MT + CT", "MT", "MT + CT x 0.5", "MT + CT x 2")
    legtitle <- ""
    pd <- 0.7
    wd <- 0.4
  } else if (panel == "h") {
    intv <- c("MT", "CT 100%", "CT 100% + 2 * p_traced", "CT 100% + 4 * p_traced")
    legintv <- c("MT", "CT", "CT x 2", "CT x 4")
    legtitle <- ""
    pd <- 0.7
    wd <- 0.5
  }

  # compute final size of each single run
  dt <- readRDS(file.path(dir_output, "SEIR.rds")) %>%
    filter(t == max(t)) %>%
    group_by(intervention, R0_1, id_sample) %>%
    summarise(final_size = 100 * (N - S - S_Q) / N) %>%
    ungroup()

  # split into no intervention and all other interventions
  dt.nointv <- dt %>%
    filter(intervention == "No intervention") %>%
    rename(final_size_nointv = final_size) %>%
    select(-intervention)

  dt.intv <- dt %>%
    filter(intervention %in% intv)

  # join no intervention values to intervention table
  dt <- dt.intv %>%
    left_join(dt.nointv) %>%
    # compute reduction
    mutate(reduction = 100 * (1 - final_size / final_size_nointv)) %>%
    # compute CrI
    group_by(intervention, R0_1) %>%
    summarise(m  = compute_CrI(reduction, "m"),
              q1 = compute_CrI(reduction, "q1"),
              q2 = compute_CrI(reduction, "q2")) %>%
    ungroup() %>%
    # clean
    mutate(R0_1 = factor(R0_1),
           intervention = factor(
             intervention,
             c("No intervention",
               "MT + CT 100% (baseline)", "MT", "CT 100%", "CT 50%", "CT 30%",
               "Baseline + 0.5 * p_traced", "Baseline + 2 * p_traced",
               "CT 100% + 2 * p_traced", "CT 100% + 4 * p_traced", "CT 100% + 16 * p_traced")))

  # Plot
  p <- ggplot(data = dt, aes(x = R0_1,
                             y = m, ymin = q1, ymax = q2,
                             colour = intervention)) +
    # model mean and 95% CrI
    geom_point(size = 1, position = position_dodge(pd)) +
    geom_errorbar(width = wd, position = position_dodge(pd)) +
    # axes
    xlab(bquote(~R[0])) +
    scale_y_continuous(name = "Relative reduction in\nepidemic final size (%)",
                       limits = c(2, 100),
                       breaks = seq(0L, 100L, by = 25L)) +
    # colour palette and legend
    scale_colour_brewer(palette = "Set1",
                        name = legtitle,
                        labels = legintv) +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          plot.margin = margin(0, 2, 0, 0),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, -5),
          # Nature requirement:
          text = element_text(size = 8, family = "sans"),
          plot.tag = element_text(size = 8, face = "bold", family = "sans")) +
    labs(tag = panel)

  return(p)

}


get_dates <- function (my.times) {

  # day 0 is the 4th of Fabruary 2020

  my.times.pos    <- my.times[my.times > -0.5]
  my.times.nonpos <- my.times[my.times < -0.5]

  dates.pos <- c(paste0(sprintf("%02d", 4:29), "/02"),
                 paste0(sprintf("%02d", 1:31), "/03"),
                 paste0(sprintf("%02d", 1:30), "/04"),
                 paste0(sprintf("%02d", 1:31), "/05"),
                 paste0(sprintf("%02d", 1:30), "/06"))
  dates.nonpos <- c(paste0(sprintf("%02d",  3:1), "/02"),
                    paste0(sprintf("%02d", 31:1), "/01"))

  dates.pos    <- dates.pos[my.times.pos + 1]
  dates.nonpos <- dates.nonpos[-my.times.nonpos]

  my.dates <- c(rev(dates.nonpos), dates.pos)

  return(my.dates)

}
