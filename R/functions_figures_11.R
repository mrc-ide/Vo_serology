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


fig_prevalence <- function () {
  
  # get compartment counts
  dt <- readRDS(file.path(dir_out, "SEIR.rds")) %>%
    # filter
    right_join(dt.filter)
  
  fixed_parameters <- dt %>%
    select(time1, time2, tQ, tSeed, N) %>%
    unique() %>%
    unlist()
  
  # Compute prevalence by pre-symptomatic, symptomatic and asymptomatic
  dt <- dt %>%
    filter(t < fixed_parameters["time2"] + 0.5) %>%
    # summarise compartments of interest
    transmute(Intervention, Quarantined, t,
              `Pre-symptomatic` = (1 - p) * TPp,
              Symptomatic       = I_S + TP_S,
              Asymptomatic      = p * TPp + I_A + TP_A) %>%
    gather("Compartment", "value", -t, -Intervention, -Quarantined) %>%
    # compute mean and 95% CrI
    group_by(Intervention, Quarantined, t, Compartment) %>%
    mutate(value = 100 * value / fixed_parameters["N"]) %>%
    summarise(mean = mean(value),
              low  = quantile(value, probs = 0.025),
              high = quantile(value, probs = 0.975)) %>%
    ungroup() %>%
    # column names for plotting
    mutate(Compartment = factor(Compartment, 
                                c("Pre-symptomatic", "Symptomatic", "Asymptomatic")))
  
  # Screening data
  data$Time <- fixed_parameters[c("time1", "time2")]
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
  my.times <- fixed_parameters[c("tSeed", "tQ", "time2")]
  my.dates <- get_dates(my.times)
  
  # Plot
  p <- ggplot(data = dt,
              aes(x = t, colour = Compartment, fill = Compartment)) +
    # facets
    facet_wrap(Intervention ~ Quarantined, scales = "free") +
    # model mean and 95% CrI
    geom_line(aes(y = mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.4, colour = NA) +
    # screening data
    geom_errorbar(data = my.data, aes(x = Time, ymin = lower, ymax = upper),
                  width = 1.2) +
    geom_point(data = my.data, aes(x = Time, y = mean), size = 1) +
    # lockdown date
    geom_vline(xintercept = as.numeric(my.times["tQ"]),
               linetype = "dashed", size = 0.3) +
    # axes
    labs(y = "Prevalence (%)") +
    # colour palette and legend
    scale_colour_brewer("", palette = "Dark2") +
    scale_fill_brewer("", palette = "Dark2") +
    scale_x_continuous(name = "Date",
                       breaks = my.times,
                       minor_breaks = seq(min(my.times), max(my.times), 7),
                       labels = my.dates[which(my.dates != "")]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-8, 0, 0, 0),
          text = element_text(size = 8, family = "sans")) +
    # legend
    theme(legend.position = "bottom",
          plot.tag = element_text(size = 8, face = "bold", family = "sans"))
  
  return(p)
  
}

fig_venn <- function () {
  
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
    theme(text = element_text(size = 3.5, family = "sans"),
          # text = element_text(size = 8, family = "sans"),
          plot.tag = element_text(size = 8, face = "bold", family = "sans"),
          # plot.margin = margin(0, 10, -8, 10)
    ) +
    labs(tag = "a")
  
  p
  
  return(p)
  
}

fig_incidence <- function (dt.filter) {
  
  # get compartment counts
  dt <- readRDS(file.path(dir_out, "SEIR.rds"))
  
  fixed_parameters <- dt %>%
    select(time1, time2, tQ, tSeed, N) %>%
    unique() %>%
    unlist()
  
  # filter SEIR
  dt <- dt %>%
    right_join(dt.filter) %>%
    mutate(Quarantined = if_else(Intervention != "Contact tracing", "", Quarantined),
           Intervention = paste(Intervention, Quarantined, sep = "\n") %>% 
             factor() %>%
             relevel("Contact tracing\nno transmission assumed") %>%
             relevel("Mass testing & lockdown\n") %>%
             relevel("No intervention\n"))
  
  # clean data
  dt <- dt %>%
    arrange(t) %>%
    group_by(id_sample, Intervention) %>%
    # summarise compartments of interest
    transmute(t, incid = -(S - lag(S, default = S[1]))) %>%
    ungroup() %>%
    # compute mean and 95% CrI
    group_by(Intervention, t) %>%
    mutate(incid = 100 * incid / fixed_parameters["N"]) %>%
    summarise(mean = mean(incid),
              low  = quantile(incid, probs = 0.025),
              high = quantile(incid, probs = 0.975)) %>%
    ungroup()
  
  # find time of disease termination
  time_final <- dt %>%
    filter(mean > 0.01) %>%
    summarise(max(t)) %>%
    pull()
  time_final <- max(time_final, fixed_parameters["time2"])
  dt <- dt %>%
    filter(t < time_final + 0.5)
  
  # Dates for x axis
  if (time_final > fixed_parameters["time2"]) {
    my.times <- c(fixed_parameters[c("tSeed", "tQ")],
                  35, 49, 63, 77, 91, 105)
  } else {
    my.times <- fixed_parameters[c("tSeed", "tQ", "time2")]
  }
  my.dates <- get_dates(unname(my.times))
  
  # Plot
  p <- ggplot(data = dt,
              aes(x = t, colour = Intervention, fill = Intervention)) +
    # model mean and 95% CrI
    geom_line(aes(y = mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.4, colour = NA) +
    # lockdown date
    geom_vline(xintercept = as.numeric(my.times["tQ"]),
               linetype = "dashed", size = .3) +
    # axes
    labs(y = "Incidence (%)") +
    # colour palette and legend
    scale_colour_brewer("", palette = "Set1") +
    scale_fill_brewer("", palette = "Set1") +
    scale_y_continuous(breaks = 0:5) +
    scale_x_continuous(name = "Date",
                       breaks = my.times,
                       minor_breaks = seq(min(my.times), max(my.times), 7),
                       labels = my.dates[which(my.dates != "")]) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          text = element_text(size = 8, family = "sans"),
          plot.tag = element_text(size = 8, face = "bold", family = "sans")) +
    labs(tag = "b")
  
  return(p)
  
}

filter_by_DIC <- function () {
  
  # get best in terms of DIC
  readRDS(file.path(dir_out, "table.rds")) %>%
    mutate(`1/sigma` = as.integer(`1/sigma`)) %>%
    filter(DIC.median < 36.42) %>%
    select(R0_1, `1/sigma`)
  
}

tab_final_size <- function () {
  
  # compute influx into Q 
  dt1 <- readRDS(file.path(dir_out, "SEIR.rds")) %>%
    mutate(
      p_3  = if_else(v == 3L, do_conttrac, 0),
      p_12 = if_else(v != 3L, do_conttrac, 0),
      q_Q = if_else(do_q > 0.01, do_q, 0),
      ww = if_else(do_lockdown > 0.5, w,  1),
      ww = if_else(t > tQ - 0.5,      ww, 1),
      beta = ww * R0_1 / (generation_time - inv_nu),
      SoverN = S/N,
      lambda = beta * (q_TPp * TPp + q_A * I_A + I_S + q_Q * (TP_Q1 + I_Q)) * SoverN) %>%
    group_by(v, Intervention, Quarantined, do_conttrac, id_sample) %>%
    summarise(influx_into_Q = sum(p_3 * lambda + p_12 * E)) %>%
    ungroup() %>%
    group_by(v, Intervention, Quarantined, do_conttrac) %>%
    summarise(m             = round(mean(influx_into_Q), 2),
              q1            = round(quantile(influx_into_Q, 0.025), 2),
              q2            = round(quantile(influx_into_Q, 0.975), 2),
              influx_into_Q = paste0(m, " (", q1, ", ", q2, ")")) %>%
    ungroup() %>%
    select(-m, -q1, -q2)
  
  
  # compute final size (%)
  dt2 <- readRDS(file.path(dir_out, "SEIR.rds")) %>%
    filter(t == max(t)) %>%
    group_by(v, Intervention, Quarantined, do_conttrac, id_sample) %>%
    summarise(final_size = 100 * (TN + TN_Q) / N,
              prop_Qarantined = 100 * TN_Q / (TN + TN_Q)) %>%
    ungroup() %>%
    group_by(v, Intervention, Quarantined, do_conttrac) %>%
    summarise(m          = round(mean(final_size), 2),
              q1         = round(quantile(final_size, 0.025), 2),
              q2         = round(quantile(final_size, 0.975), 2),
              final_size = paste0(m, " (", q1, ", ", q2, ")"),
              m_Q             = round(mean(prop_Qarantined), 2),
              q1_Q            = round(quantile(prop_Qarantined, 0.025), 2),
              q2_Q            = round(quantile(prop_Qarantined, 0.975), 2),
              prop_Qarantined = paste0(m_Q, " (", q1_Q, ", ", q2_Q, ")")) %>%
    ungroup()
  
  # join tables
  dt <- left_join(dt1, dt2)
  
  # table
  saveRDS(dt, file.path(dir_out, "table_final_size.rds"))
  
}

fig_final_size <- function (dt.filter) {
  
  # read in final size (%)
  dt <- readRDS(file.path(dir_out, "table_final_size.rds"))
  
  # filter SEIR
  dt <- dt %>%
    right_join(dt.filter) %>%
    mutate(Quarantined = if_else(Intervention != "Contact tracing", "", Quarantined),
           Intervention = paste(Intervention, Quarantined, sep = "\n") %>% 
             factor() %>%
             relevel("Contact tracing\nno transmission assumed") %>%
             relevel("Mass testing & lockdown\n") %>%
             relevel("No intervention\n"))
  
  # figure
  p <- ggplot(data = dt, aes(x = Intervention, y = m,
                             ymin = q1, ymax = q2, 
                             fill = Intervention)) +
    # model mean and 95% CrI
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(width = 0.2, position = position_dodge(0.8)) +
    # axes
    scale_y_continuous(name = "Epidemic final size (%)    ",
                       limits = c(0, 98),
                       breaks = seq(0L, 100L, by = 25L)) +
    # colour palette and legend
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 8, family = "sans"),
          plot.tag = element_text(size = 8, face = "bold", family = "sans")) +
    labs(tag = "c")
  
  return(p)
  
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
    event = c("First nasopharyngeal survey", "Second nasopharyngeal survey",
              "First serosurvey and\nnasopharyngeal survey",
              "Second serosurvey and\nnasopharyngeal survey")
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
