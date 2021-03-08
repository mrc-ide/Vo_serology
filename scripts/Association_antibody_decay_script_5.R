################################################################################
# This code runs the association analysis on antibody decay rates observed     #
# between May and November vs age and BMI.                                     #
#                                                                              #
# Cite as:                                                                     #
# Dorigatti I et al, SARS-CoV-2 antibody dynamics, within-household            #
# transmission and the impact of contact tracing from community-wide           #
# serological testing in the Italian municipality of Vo'.                      #
# Pre-print available at: https://ssrn.com/abstract=3792892                    #
#                                                                              #
# Description:                                                                 #
# See README.md                                                                #
################################################################################

cat("\n\n### Running 'scripts/Association_antibody_decay_script_5.R'\n\n")


# set seed for bootstrapping
set.seed(1)

# source scripts
source("R/functions_association_antibody_slopes.R")

# read data
ds <- read.csv("data/Vo_serology_data.csv",
               header = TRUE)

# output dir
path_figure <- file.path(getwd(), "figures")
path_table <- file.path(getwd(), "tables")

# formatting ------------------------------------------------------------------#

ds$first_symptoms_date <- as.character(ds$first_symptoms_date)

ds$Abbott_quantitative_november_2020 <-
   as.numeric(as.character(ds$Abbott_quantitative_november_2020))

ds$Diasorin_quantitative_november_2020 <-
  as.numeric(as.character(ds$Diasorin_quantitative_november_2020))

ds$Roche_quantitative_november_2020 <-
   as.numeric(as.character(ds$Roche_quantitative_november_2020))

levels(ds$hospitalized)[levels(ds$hospitalized)==""] <- "no"

# BMI
ds$BMI <- ds$Weight_kg/(ds$Height_cm/100)^2

# BMI categories
ds$BMI_category <- NA
ds$BMI_category[which(ds$BMI < 18.5)] <- "Underweight"
ds$BMI_category[which(ds$BMI >= 18.5 & ds$BMI< 25)] <- "Normal"
ds$BMI_category[which(ds$BMI >= 25 & ds$BMI< 30)] <- "Overweight"
ds$BMI_category[which(ds$BMI >= 30)] <- "Obese"

ds$BMI_category <- as.factor(ds$BMI_category)

ds$BMI_category <- factor(ds$BMI_category, levels= c("Underweight",
                                                     "Normal",
                                                     "Overweight",
                                                     "Obese",
                                                     "NA"))

levels(ds$Symptomatic_any_time_from_1_january_to_1_may)[
  levels(ds$Symptomatic_any_time_from_1_january_to_1_may) == ""] <- "no"


# select baseline groud truth
ds <- ds[which(ds$Groundtruth == 1), ]

# number of days between serosurveys
days_surveys <- as.numeric(as.Date("11/29/2020", "%m/%d/%Y") -
  as.Date("05/01/2020", "%m/%d/%Y"))

################################################################################
# Decay rate                                                                   #
################################################################################

ds$slope_Abbott <- log(ds$Abbott_quantitative_november_2020/
                          ds$Abbot_semiquantitative)/days_surveys

ds$slope_Diasorin <- log(ds$Diasorin_quantitative_november_2020/
                           ds$Diasorin_IgG_semiquantitative)/days_surveys

ds$slope_Roche <- log(ds$Roche_quantitative_november_2020/
                        ds$Roche_Total_ICO)/days_surveys

################################################################################
# Half life                                                                    #
################################################################################

ds$hl_Abbott <- log(0.5)/ds$slope_Abbott
ds$hl_Diasorin <- log(0.5)/ds$slope_Diasorin
ds$hl_Roche <- log(0.5)/ds$slope_Roche

# Abbott
compute_half_life(ds, "Abbott")
bootstrap_hl(ds, "Abbott")

# DiaSorin
compute_half_life(ds, "DiaSorin")
bootstrap_hl(ds, "DiaSorin")

# Roche
compute_half_life(ds, "Roche")
bootstrap_hl(ds, "Roche")

################################################################################
# Decay rate vs symptoms                                                       #
################################################################################

# Abbott
var1 <- "Abbot_qualitative"
var2 <- "slope_Abbott"
type <- "symptom occurrence"
type2 <- "decay rate"

(res <- test_decay_rate_half_life(var1, var2, type, type2))

# Diasorin
var1 <- "Diasorin_IgG_qualitative"
var2 <- "slope_Diasorin"
type <- "symptom occurrence"
type2 <- "decay rate"

(res <- test_decay_rate_half_life(var1, var2, type, type2))

# Roche
var1 <- "Roche_Total_qualitative"
var2 <- "slope_Roche"
type <- "symptom occurrence"
type2 <- "decay rate"

(res <- test_decay_rate_half_life(var1, var2, type, type2))

################################################################################
# Decay rate vs hospitalisation                                                #
################################################################################

# Abbott
var1 <- "Abbot_qualitative"
var2 <- "slope_Abbott"
type <- "hospitalized"
type2 <- "decay rate"

(res <- test_decay_rate_half_life(var1, var2, type, type2))

# Diasorin
var1 <- "Diasorin_IgG_qualitative"
var2 <- "slope_Diasorin"
type <- "hospitalized"
type2 <- "decay rate"

(res <- test_decay_rate_half_life(var1, var2, type, type2))

# Roche
var1 <- "Roche_Total_qualitative"
var2 <- "slope_Roche"
type <- "hospitalized"
type2 <- "decay rate"

(res <- test_decay_rate_half_life(var1, var2, type, type2))

################################################################################
# Decay rate vs sex                                                            #
################################################################################

# Abbott
var1 <- "Abbot_qualitative"
var2 <- "slope_Abbott"
type <- "sex"
type2 <- "decay rate"

(res <- test_decay_rate_half_life(var1, var2, type, type2))

# Diasorin
var1 <- "Diasorin_IgG_qualitative"
var2 <- "slope_Diasorin"
type <- "sex"
type2 <- "decay rate"

(res <- test_decay_rate_half_life(var1, var2, type, type2))

# Roche
var1 <- "Roche_Total_qualitative"
var2 <- "slope_Roche"
type <- "sex"
type2 <- "decay rate"

(res <- test_decay_rate_half_life(var1, var2, type, type2))

################################################################################
# Decay rate vs age group                                                      #
################################################################################

# Abbott
var1 <- "Abbot_qualitative"
var2 <- "slope_Abbott"
type <- "age group"
type2 <- "decay rate"

(res <- Anova_Kuskal_decay_rate_half_life(var1, var2, type, type2))
p1 <- res[[4]]


# Diasorin
var1 <- "Diasorin_IgG_qualitative"
var2 <- "slope_Diasorin"
type <- "age group"
type2 <- "decay rate"

(res <- Anova_Kuskal_decay_rate_half_life(var1, var2, type, type2))
p2 <- res[[4]]

# Roche
var1 <- "Roche_Total_qualitative"
var2 <- "slope_Roche"
type <- "age group"
type2 <- "decay rate"

(res <- Anova_Kuskal_decay_rate_half_life(var1, var2, type, type2))
p3 <- res[[4]]

################################################################################
# Decay rate vs BMI category                                                   #
################################################################################

# Abbott
var1 <- "Abbot_qualitative"
var2 <- "slope_Abbott"
type <- "BMI category"
type2 <- "decay rate"

(res <- Anova_Kuskal_decay_rate_half_life(var1, var2, type, type2))
q1 <- res[[4]]

# Diasorin
var1 <- "Diasorin_IgG_qualitative"
var2 <- "slope_Diasorin"
type <- "BMI category"
type2 <- "decay rate"

(res <- Anova_Kuskal_decay_rate_half_life(var1, var2, type, type2))
q2 <- res[[4]]

# Roche
var1 <- "Roche_Total_qualitative"
var2 <- "slope_Roche"
type <- "BMI category"
type2 <- "decay rate"

(res <- Anova_Kuskal_decay_rate_half_life(var1, var2, type, type2))
q3 <- res[[4]]

################################################################################
# Multiple linear regression: decay ~ symptoms + age + BMI                  #
################################################################################

# Abbott
idx <- which(ds$Abbot_qualitative == "Positive")

df <- data.frame(
  antibody_slope = ds$slope_Abbott[idx],
  symptomatic = as.factor(ds$Symptomatic_any_time_from_1_january_to_1_may[idx]),
  age_group = as.factor(ds$age_group[idx]),
  BMI = ds$BMI[idx])

fit <- lm(antibody_slope ~ BMI*symptomatic, data = df)
summary(fit)

fit <- lm(antibody_slope ~ BMI + symptomatic, data = df)
summary(fit)

fit <- lm(antibody_slope ~ BMI , data = df[df$symptomatic == "no", ])
summary(fit)

fit <- lm(antibody_slope ~ BMI , data = df[df$symptomatic == "yes", ])
summary(fit)

# Diasorin
idx <- which(ds$Diasorin_IgG_qualitative == "Positive")

df <- data.frame(
  antibody_slope = ds$slope_Diasorin[idx],
  symptomatic = as.factor(ds$Symptomatic_any_time_from_1_january_to_1_may[idx]),
  age_group = as.factor(ds$age_group[idx]),
  BMI = ds$BMI[idx])

fit <- lm(antibody_slope ~ BMI*symptomatic, data = df)
summary(fit)

fit <- lm(antibody_slope ~ BMI + symptomatic, data = df)
summary(fit)

fit <- lm(antibody_slope ~ BMI , data = df[df$symptomatic == "no", ])
summary(fit)

fit <- lm(antibody_slope ~ BMI , data = df[df$symptomatic == "yes", ])
summary(fit)

# Roche
idx <- which(ds$Roche_Total_qualitative == "Positive")

df <- data.frame(
  antibody_slope = ds$slope_Roche[idx],
  symptomatic = as.factor(ds$Symptomatic_any_time_from_1_january_to_1_may[idx]),
  age_group = as.factor(ds$age_group[idx]),
  BMI = ds$BMI[idx])

fit <- lm(antibody_slope ~ BMI*symptomatic, data = df)
summary(fit)

fit <- lm(antibody_slope ~ BMI + symptomatic, data = df)
summary(fit)

fit <- lm(antibody_slope ~ BMI , data = df[df$symptomatic == "no", ])
summary(fit)

fit <- lm(antibody_slope ~ BMI , data = df[df$symptomatic == "yes", ])
summary(fit)

r3 <- ggplot(df[df$symptomatic == "yes", ],aes(y=antibody_slope,
                                               x=BMI,
                                               col=symptomatic))+
  geom_point()+
  ylab("Antibody decay rate (Roche)")+ xlab("BMI")+
  stat_smooth(method="lm",se=TRUE)+
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_blank())

# Figure S2 -------------------------------------------------------------------#

cols12 <- readRDS(file.path(path_figure, "temp.rds"))
col3 <- plot_grid(p1, p2, p3,
          labels = c('g', 'h', 'i'),
          label_size = 8,
          ncol = 3,
          nrow = 1,
          label_fontface = "bold")

ggsave(filename = file.path(path_figure, "Figure_S2.tiff"),
       plot = plot_grid(cols12, col3,
                        labels = c('', ''),
                        label_size = 8,
                        ncol = 1,
                        rel_heights = c(2,1),
                        label_fontface = "bold"),
       device = "tiff",
       width = 300, height = 270,
       units = "mm", dpi = 300, limitsize = TRUE)
# -----------------------------------------------------------------------------#

################################################################################
################################################################################
