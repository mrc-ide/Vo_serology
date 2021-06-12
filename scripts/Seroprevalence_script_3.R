################################################################################
# This code produced the seroprevalence estimates shown in Figure 2 and the    #
# estimates summarised in Supplementary Table S1.                              #
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

cat("\n\n### Running 'scripts/Seroprevalence_script_3.R'\n\n")


# fix seed
set.seed(1)

# source scripts
source("R/functions_seroprevalence.R")

# read data
ds <- read.csv("data/Vo_serology_data.csv",
               header = TRUE)

# output dir
path_figure <- file.path(getwd(), "figures")
path_table <- file.path(getwd(), "tables")

# colors for plotting
cols <- gg_colour_hue(10)

# full dataset
ds_cluster <- ds

# select Vo residents
ds <- ds[which(ds$town == "VO"), ]

# Seroprevalence by assay -----------------------------------------------------#

# raw seroprevalence Abbott
raw_serop_Abbott <- raw_seroprevalence(ds$Abbot_qualitative)

# raw seroprevalence DiaSorin
raw_serop_DiaSorin <- raw_seroprevalence(ds$Diasorin_IgG_qualitative)

# raw seroprevalence Roche
raw_serop_Roche <- raw_seroprevalence(ds$Roche_Total_qualitative)

# Adjusted seroprevalence by assay-specific sensitivity and specificity -------#

# observed (Table 5)
se_A <- 0.939
sp_A <- 1

se_D <- 0.852
sp_D <- 1

se_R <- 0.939
sp_R <- 0.976

# adjusted seroprevalence Abbott
adj_serop_Abbott <- adjusted_seroprevalence(ds$Abbot_qualitative,
                                            se_A, sp_A)

# adjusted seroprevalence DiaSorin
adj_serop_DiaSorin <- adjusted_seroprevalence(ds$Diasorin_IgG_qualitative,
                                              se_D, sp_D)

# adjusted seroprevalence Roche
adj_serop_Roche <- adjusted_seroprevalence(ds$Roche_Total_qualitative,
                                           se_R, sp_R)

# estimated (Table 5)
se_A_est <- 0.969
sp_A_est <- 0.997

se_D_est <- 0.856
sp_D_est <- 0.982

se_R_est <- 0.968
sp_R_est <- 0.997

# adjusted seroprevalence Abbott
adj_est_serop_Abbott <- adjusted_seroprevalence(ds$Abbot_qualitative,
                                                se_A_est, sp_A_est)

# adjusted seroprevalence DiaSorin
adj_est_serop_DiaSorin <- adjusted_seroprevalence(ds$Diasorin_IgG_qualitative,
                                              se_D_est, sp_D_est)

# adjusted seroprevalence Roche
adj_est_serop_Roche <- adjusted_seroprevalence(ds$Roche_Total_qualitative,
                                           se_R_est, sp_R_est)

# Plot seroprevalence by assay ------------------------------------------------#
df_serop <- data.frame(seroprevalence = c(raw_serop_Abbott[1],
                                        raw_serop_DiaSorin[1],
                                        raw_serop_Roche[1],
                                        adj_serop_Abbott[1],
                                        adj_serop_DiaSorin[1],
                                        adj_serop_Roche[1],
                                        adj_est_serop_Abbott[1],
                                        adj_est_serop_DiaSorin[1],
                                        adj_est_serop_Roche[1]),

                     assay = rep(c("Abbott","DiaSorin","Roche"), 3),

                     lower_CI = c(raw_serop_Abbott[2],
                                  raw_serop_DiaSorin[2],
                                  raw_serop_Roche[2],
                                  adj_serop_Abbott[2],
                                  adj_serop_DiaSorin[2],
                                  adj_serop_Roche[2],
                                  adj_est_serop_Abbott[2],
                                  adj_est_serop_DiaSorin[2],
                                  adj_est_serop_Roche[2]),

                     upper_CI = c(raw_serop_Abbott[3],
                                  raw_serop_DiaSorin[3],
                                  raw_serop_Roche[3],
                                  adj_serop_Abbott[3],
                                  adj_serop_DiaSorin[3],
                                  adj_serop_Roche[3],
                                  adj_est_serop_Abbott[3],
                                  adj_est_serop_DiaSorin[3],
                                  adj_est_serop_Roche[3]),

                     type = as.factor(c(rep("raw", 3),
                                      rep("adjusted",3),#in-house
                                      rep("adjusted",3))),#model

                     line = as.factor(c(rep("combined", 3),
                                        rep("VEO",3),
                                        rep("combined", 3))))

p_adj <- ggplot(df_serop)+ylim(0.0095, 0.085)+
  geom_point(aes(x = assay, y = seroprevalence, colour = type, linetype = line),
             position = position_dodge(0.5))+
  geom_errorbar(aes(x = assay, y = seroprevalence,
                    ymin = lower_CI, ymax = upper_CI, colour = type,
                    linetype = line),
                width = 0.2,
                position = position_dodge(0.5))+
  geom_rect(aes(xmin = 0.3, xmax = 3.8, ymin=0.028, ymax=0.043),
            fill = "#FDE725FF", alpha = 0.01)+
  geom_hline(aes(yintercept = 0.035), colour = "#FDE725FF")+
  geom_point(aes(x = assay, y = seroprevalence, colour = type, linetype = line),
             position = position_dodge(0.5))+
  geom_errorbar(aes(x = assay, y = seroprevalence,
                    ymin = lower_CI, ymax = upper_CI,
                    colour = type, linetype = line),
                width = 0.2,
                position = position_dodge(0.5))+
  xlab("")+
  theme_bw()+
  theme(text = element_text(size=9),
        legend.position = c(0.45, 0.90),
        legend.key.height=unit(0.5,"line"),
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        legend.background = element_rect(colour ="transparent",
                                         fill = "transparent"))+
  scale_color_manual(values = cols[c(10,6,8)])+
  scale_linetype_manual(values = c("solid","dashed"))

# Plot seroprevalence by defintion  -------------------------------------------#

colnames_gt <- c("Diasorin_IgG_qualitative",
                 "Roche_Total_qualitative",
                 "Abbot_qualitative",
                 "Neutralization",
                 "positive_swab")

# discard subjects with no PCR or serological test taken
discard <- vector()
for(i in 1:dim(ds)[1]){

  if(all(is.na(ds[i, colnames_gt])) == TRUE) discard <- c(discard, i)

}

# number of subjects with at least one PCR or serological result
denom_GT <- dim(ds[-discard, ])[1]

num_GT1 <- length(which(ds$Groundtruth == 1))
num_GT2 <- length(which(ds$Groundtruth_direct_contacts == 1))
num_GT3 <- length(which(ds$Groundtruth_indirect_contacts == 1))

confint_GT1  <- binom.confint(num_GT1,
                              denom_GT,
                              conf.level = 0.95,
                              methods = "exact")

confint_GT2  <- binom.confint(num_GT2,
                              denom_GT,
                              conf.level = 0.95,
                              methods = "exact")

confint_GT3  <- binom.confint(num_GT3,
                              denom_GT,
                              conf.level = 0.95,
                              methods = "exact")

df_GT <- data.frame(seroprevalence = c(confint_GT1$mean,
                                       confint_GT2$mean,
                                       confint_GT3$mean),

                    definition = c("GT\nbaseline",
                                   "GT\ndirect",
                                   "GT\nindirect"),

                     lower_CI = c(confint_GT1$lower,
                                  confint_GT2$lower,
                                  confint_GT3$lower),

                     upper_CI = c(confint_GT1$upper,
                                  confint_GT2$upper,
                                  confint_GT3$upper))

q <- ggplot(df_GT)+ylim(0.0095, 0.085)+
  geom_point(aes(x = definition, y = seroprevalence, colour = definition))+
  geom_errorbar(aes(x = definition, y = seroprevalence,
                    ymin = lower_CI, ymax = upper_CI, colour = definition),
                width = 0.2,
                position = position_dodge(0.9))+
  geom_rect(aes(xmin = 0.3, xmax = 3.8, ymin=0.028, ymax=0.043),
            fill = "#FDE725FF", alpha = 0.03)+
  geom_hline(aes(yintercept = 0.035), colour = "#FDE725FF")+
  geom_point(aes(x = definition, y = seroprevalence, colour = definition))+
  geom_errorbar(aes(x = definition, y = seroprevalence,
                    ymin = lower_CI, ymax = upper_CI, colour = definition),
                width = 0.2,
                position = position_dodge(0.9))+
  theme_bw() +
  theme(text = element_text(size=9),
        legend.position = "none",
        )+xlab("")+ylab("seroprevalence")

# Estimate assay performance with respect to the different GT defintions
# Supplementary Table S1

# Vo residents
# Abbott baseline
(est <- estimate_sens_spec_ppv_npv(ds$Abbot_qualitative,ds$Groundtruth))
seroprev_Abbott_baseline <- est[[2]]

# DiaSorin baseline
(est <- estimate_sens_spec_ppv_npv(ds$Diasorin_IgG_qualitative,
                                   ds$Groundtruth))
seroprev_DiaSorin_baseline <- est[[2]]

# Roche baseline
(est <- estimate_sens_spec_ppv_npv(ds$Roche_Total_qualitative,ds$Groundtruth))
seroprev_Roche_baseline <- est[[2]]

# Abbott direct
(est <- estimate_sens_spec_ppv_npv(ds$Abbot_qualitative,
                                   ds$Groundtruth_direct_contacts))
seroprev_Abbott_direct <- est[[2]]

# DiaSorin direct
(est <- estimate_sens_spec_ppv_npv(ds$Diasorin_IgG_qualitative,
                                   ds$Groundtruth_direct_contacts))
seroprev_DiaSorin_direct <- est[[2]]

# Roche direct
(est <- estimate_sens_spec_ppv_npv(ds$Roche_Total_qualitative,
                                   ds$Groundtruth_direct_contacts))
seroprev_Roche_direct <- est[[2]]

# Abbott indirect
(est <- estimate_sens_spec_ppv_npv(ds$Abbot_qualitative,
                                   ds$Groundtruth_indirect_contacts))
seroprev_Abbott_indirect <- est[[2]]

# DiaSorin indirect
(est <- estimate_sens_spec_ppv_npv(ds$Diasorin_IgG_qualitative,
                                   ds$Groundtruth_indirect_contacts))
seroprev_DiaSorin_indirect <- est[[2]]

# Roche indirect
(est <- estimate_sens_spec_ppv_npv(ds$Roche_Total_qualitative,
                                   ds$Groundtruth_indirect_contacts))
seroprev_Roche_indirect <- est[[2]]

# Vo cluster
# Abbott baseline
(est_cluster <- estimate_sens_spec_ppv_npv(ds_cluster$Abbot_qualitative,
                                           ds_cluster$Groundtruth))

# DiaSorin baseline
(est_cluster <- estimate_sens_spec_ppv_npv(ds_cluster$Diasorin_IgG_qualitative,
                                           ds_cluster$Groundtruth))

# Roche baseline
(est_cluster <- estimate_sens_spec_ppv_npv(ds_cluster$Roche_Total_qualitative,
                                           ds_cluster$Groundtruth))

# Abbott direct
(est_cluster <- estimate_sens_spec_ppv_npv(ds_cluster$Abbot_qualitative,
                                   ds_cluster$Groundtruth_direct_contacts))

# DiaSorin direct
(est_cluster <- estimate_sens_spec_ppv_npv(ds_cluster$Diasorin_IgG_qualitative,
                                   ds_cluster$Groundtruth_direct_contacts))

# Roche direct
(est_cluster <- estimate_sens_spec_ppv_npv(ds_cluster$Roche_Total_qualitative,
                                   ds_cluster$Groundtruth_direct_contacts))

# Abbott indirect
(est_cluster <- estimate_sens_spec_ppv_npv(ds_cluster$Abbot_qualitative,
                                   ds_cluster$Groundtruth_indirect_contacts))

# DiaSorin indirect
(est_cluster <- estimate_sens_spec_ppv_npv(ds_cluster$Diasorin_IgG_qualitative,
                                   ds_cluster$Groundtruth_indirect_contacts))

# Roche indirect
(est_cluster <- estimate_sens_spec_ppv_npv(ds_cluster$Roche_Total_qualitative,
                                   ds_cluster$Groundtruth_indirect_contacts))


# plot seroprevalence adjusted by se and sp by definition  --------------------#
df_serop_adj_def <- data.frame(serop = c(seroprev_Abbott_baseline[1],
                                          seroprev_Abbott_direct[1],
                                          seroprev_Abbott_indirect[1],
                                          seroprev_DiaSorin_baseline[1],
                                          seroprev_DiaSorin_direct[1],
                                          seroprev_DiaSorin_indirect[1],
                                          seroprev_Roche_baseline[1],
                                          seroprev_Roche_direct[1],
                                          seroprev_Roche_indirect[1]),

                               lower = c(seroprev_Abbott_baseline[2],
                                         seroprev_Abbott_direct[2],
                                         seroprev_Abbott_indirect[2],
                                         seroprev_DiaSorin_baseline[2],
                                         seroprev_DiaSorin_direct[2],
                                         seroprev_DiaSorin_indirect[2],
                                         seroprev_Roche_baseline[2],
                                         seroprev_Roche_direct[2],
                                         seroprev_Roche_indirect[2]),

                               upper = c(seroprev_Abbott_baseline[3],
                                         seroprev_Abbott_direct[3],
                                         seroprev_Abbott_indirect[3],
                                         seroprev_DiaSorin_baseline[3],
                                         seroprev_DiaSorin_direct[3],
                                         seroprev_DiaSorin_indirect[3],
                                         seroprev_Roche_baseline[3],
                                         seroprev_Roche_direct[3],
                                         seroprev_Roche_indirect[3]),

                                definition = rep(c("GT baseline",
                                                   "GT direct contacts",
                                                   "GT indirect contacts"), 3),

                                assay = c(rep("Abbott", 3),
                                          rep("DiaSorin", 3),
                                          rep("Roche", 3)))


p_serop_adj_def <- ggplot(df_serop_adj_def)+ylim(0.0095, 0.085)+
  geom_point(aes(x = assay, y = serop, colour = definition),
             position = position_dodge(0.5))+
  geom_errorbar(aes(x = assay, y = serop,
                    ymin = lower, ymax = upper, colour = definition),
                width = 0.2,
                position = position_dodge(0.5))+
  geom_rect(aes(xmin = 0.3, xmax = 3.8, ymin=0.028, ymax=0.043),
            fill = "#FDE725FF", alpha = 0.01)+
  geom_hline(aes(yintercept = 0.035), colour = "#FDE725FF")+
  geom_point(aes(x = assay, y = serop, colour = definition),
             position = position_dodge(0.5))+
  geom_errorbar(aes(x = assay, y = serop,
                    ymin = lower, ymax = upper,
                    colour = definition),
                width = 0.2,
                position = position_dodge(0.5))+
  xlab("")+ ylab("seroprevalence")+
  theme_bw()+
  theme(text = element_text(size=9),
        legend.key.height=unit(0.5,"line"),
        legend.position = c(0.30, 0.90),
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        legend.background = element_rect(colour ="transparent",
                                         fill = "transparent"))

# Figure 2 -------------------------------------------------------------------#
row1 <- plot_grid(p_adj, q, p_serop_adj_def,
                  labels = c('a', 'b', 'c'),
                  label_size = 8, ncol = 3,
                  label_fontface = "bold")

 ggsave(filename = file.path(path_figure, "Figure_2.eps"),
        plot = row1,
        device = "eps",
        width = 180, height = 60,
        units = "mm", dpi = 300, limitsize = TRUE)
 
 ggsave(filename = file.path(path_figure, "Figure_2.tiff"),
        plot = row1,
        device = "tiff",
        width = 180, height = 60,
        units = "mm", dpi = 300, limitsize = TRUE)
# -----------------------------------------------------------------------------#

################################################################################
################################################################################
