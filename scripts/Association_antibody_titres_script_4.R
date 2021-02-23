################################################################################
# This code runs the association analysis on antibody titres observed in  May  #
# and November vs days since symptom onset, age and BMI.                       #
#                                                                              #
# Cite as:                                                                     #
# Dorigatti I et al, SARS-CoV-2 antibody dynamics, within-household            #
# transmission and the impact of contact tracing from community-wide           #
# serological testing in the Italian municipality of Vo'.                      #
#                                                                              #
# Description:                                                                 #
# See README.md                                                                #
################################################################################

cat("\n\n### Running 'scripts/Association_antibody_titres_script_4.R'\n\n")


# source scripts
source("R/functions_association_antibody_titres.R")

# read data
ds <- read.csv("data/Vo_serology_data.csv",
               header = TRUE)

# output dir
path_figure <- file.path(getwd(), "figures")
path_table <- file.path(getwd(), "tables")

# formatting ------------------------------------------------------------------#
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

################################################################################
# Logistic regression symptoms vs age + BMI                                    #
################################################################################

# symptoms vs age
tab <- table(ds$Symptomatic_any_time_from_1_january_to_1_may, ds$age_group)

confint <- binom.confint(tab[2,],
                         colSums(tab),
                         conf.level = 0.95,
                         methods = "exact")

df <- data.frame(prob = tab[2,]/colSums(tab),
                 age = colnames(tab),
                 low = confint$lower,
                 up = confint$upper)

plot1 <- ggplot(data = df)+
  geom_point(aes(y = prob, x = age))+
  geom_errorbar(aes(x = age, ymin = low, ymax = up),
                width = 0.2)+
  ylab("probability of being symptomatic")+
  xlab("age groups (years)")

# symptoms vs BMI
tab2 <- table(ds$Symptomatic_any_time_from_1_january_to_1_may, ds$BMI_category)

confint2 <- binom.confint(tab2[2,],
                          colSums(tab2),
                          conf.level = 0.95,
                          methods = "exact")

df2 <- data.frame(prob = tab2[2,]/colSums(tab2),
                 bmi = colnames(tab2),
                 low = confint2$lower,
                 up = confint2$upper)

df2 <- df2[-5, ]

df2$bmi <- factor(df2$bmi,
                  levels = c("Underweight",
                             "Normal",
                             "Overweight",
                             "Obese"))

plot2 <- ggplot(data = df2)+
  geom_point(aes(y = prob, x = bmi))+
  geom_errorbar(aes(x = bmi, ymin = low, ymax = up),
                width = 0.2)+
  ylab("probability of being symptomatic")+
  xlab("BMI category")


ds$age_group <- relevel(ds$age_group, ref = "51-60")
ds$BMI_category <- relevel(ds$BMI_category, ref = "Normal")

# logistic regression symptoms ~ age
res <- glm(Symptomatic_any_time_from_1_january_to_1_may ~ age_group,
           data = ds, family = "binomial")
summary(res)
confint(res)

exp(coef(res))
exp(cbind(OR = coef(res), confint(res)))

# logistic regression symptoms ~ BMI
res <- glm(Symptomatic_any_time_from_1_january_to_1_may ~ BMI_category,
           data = ds, family = "binomial")
summary(res)
confint(res)

exp(coef(res))
exp(cbind(OR = coef(res), confint(res)))

# logistic regression symptoms ~ age + BMI
res <- glm(Symptomatic_any_time_from_1_january_to_1_may ~
             age_group + BMI_category,
           data = ds, family = "binomial")
summary(res)
confint(res)

exp(coef(res))
exp(cbind(OR = coef(res), confint(res)))

################################################################################
# Antibody titres in time since symptom onset                                  #
################################################################################

ds$first_symptoms_date <- as.character(ds$first_symptoms_date)

ds$Abbott_quantitative_november_2020 <-
  as.numeric(as.character(ds$Abbott_quantitative_november_2020))

ds$Diasorin_quantitative_november_2020 <-
  as.numeric(as.character(ds$Diasorin_quantitative_november_2020))

ds$Roche_quantitative_november_2020 <-
  as.numeric(as.character(ds$Roche_quantitative_november_2020))

# Abbott
# May
var1 <- "Abbot_qualitative"
var2 <- "Abbot_semiquantitative"

res <- association_titres_days_symptom_onset(var1, var2)
res[[1]]
res[[2]]

# Nov
var1 <- "Abbott_qualitative_november_2020"
var2 <- "Abbott_quantitative_november_2020"

res <- association_titres_days_symptom_onset(var1, var2)
res[[1]]
res[[2]]
res[[3]]
res[[4]]

# Diasorin
# May
var1 <- "Diasorin_IgG_qualitative"
var2 <- "Diasorin_IgG_semiquantitative"

res <- association_titres_days_symptom_onset(var1, var2)
res[[1]]
res[[2]]
res[[3]]
res[[4]]

# Nov
var1 <- "Diasorin_qualitative_november_2020"
var2 <- "Diasorin_quantitative_november_2020"

res <- association_titres_days_symptom_onset(var1, var2)
res[[1]]
res[[2]]
res[[3]]
res[[4]]

# Roche
# May
var1 <- "Roche_Total_qualitative"
var2 <- "Roche_Total_ICO"

res <- association_titres_days_symptom_onset(var1, var2)
res[[1]]
res[[2]]

# Nov
var1 <- "Roche_qualitative_november_2020"
var2 <- "Roche_quantitative_november_2020"

res <- association_titres_days_symptom_onset(var1, var2)
res[[1]]
res[[2]]

################################################################################
# Antibody titres by symptoms occurrence                                       #
################################################################################

# Abbott
# May
var1 <- "Abbot_qualitative"
var2 <- "Abbot_semiquantitative"
type <- "symptom occurrence"

(res <- test_antibody_titres(var1, var2, type))

# Nov
var1 <- "Abbott_qualitative_november_2020"
var2 <- "Abbott_quantitative_november_2020"
type <- "symptom occurrence"

(res <- test_antibody_titres(var1, var2, type))

# Diasorin
# May
var1 <- "Diasorin_IgG_qualitative"
var2 <- "Diasorin_IgG_semiquantitative"
type <- "symptom occurrence"

(res <- test_antibody_titres(var1, var2, type))

# Nov
var1 <- "Diasorin_qualitative_november_2020"
var2 <- "Diasorin_quantitative_november_2020"
type <- "symptom occurrence"

(res <- test_antibody_titres(var1, var2, type))

# Roche
# May
var1 <- "Roche_Total_qualitative"
var2 <- "Roche_Total_ICO"
type <- "symptom occurrence"

(res <- test_antibody_titres(var1, var2, type))

# Nov
var1 <- "Roche_qualitative_november_2020"
var2 <- "Roche_quantitative_november_2020"
type <- "symptom occurrence"

(res <- test_antibody_titres(var1, var2, type))

################################################################################
# Antibody titres by hospitalisation status                                    #
################################################################################

# Abbott
# May
var1 <- "Abbot_qualitative"
var2 <- "Abbot_semiquantitative"
type <- "hospitalized"

(res <- test_antibody_titres(var1, var2, type))

# Nov
var1 <- "Abbott_qualitative_november_2020"
var2 <- "Abbott_quantitative_november_2020"
type <- "hospitalized"

(res <- test_antibody_titres(var1, var2, type))

# Diasorin
# May
var1 <- "Diasorin_IgG_qualitative"
var2 <- "Diasorin_IgG_semiquantitative"
type <- "hospitalized"

(res <- test_antibody_titres(var1, var2, type))

# Nov
var1 <- "Diasorin_qualitative_november_2020"
var2 <- "Diasorin_quantitative_november_2020"
type <- "hospitalized"

(res <- test_antibody_titres(var1, var2, type))

# Roche
# May
var1 <- "Roche_Total_qualitative"
var2 <- "Roche_Total_ICO"
type <- "hospitalized"

(res <- test_antibody_titres(var1, var2, type))

# Nov
var1 <- "Roche_qualitative_november_2020"
var2 <- "Roche_quantitative_november_2020"
type <- "hospitalized"

(res <- test_antibody_titres(var1, var2, type))

################################################################################
# Antibody titres by sex                                                       #
################################################################################

# Abbott
# May
var1 <- "Abbot_qualitative"
var2 <- "Abbot_semiquantitative"
type <- "sex"

(res <- test_antibody_titres(var1, var2, type))

# Nov
var1 <- "Abbott_qualitative_november_2020"
var2 <- "Abbott_quantitative_november_2020"
type <- "sex"

(res <- test_antibody_titres(var1, var2, type))

# Diasorin
# May
var1 <- "Diasorin_IgG_qualitative"
var2 <- "Diasorin_IgG_semiquantitative"
type <- "sex"

(res <- test_antibody_titres(var1, var2, type))

# Nov
var1 <- "Diasorin_qualitative_november_2020"
var2 <- "Diasorin_quantitative_november_2020"
type <- "sex"

(res <- test_antibody_titres(var1, var2, type))

# Roche
# May
var1 <- "Roche_Total_qualitative"
var2 <- "Roche_Total_ICO"
type <- "sex"

(res <- test_antibody_titres(var1, var2, type))

# Nov
var1 <- "Roche_qualitative_november_2020"
var2 <- "Roche_quantitative_november_2020"
type <- "sex"

(res <- test_antibody_titres(var1, var2, type))

################################################################################
# Antibody titres by age group                                                 #
################################################################################

# Abbott
# May
var1 <- "Abbot_qualitative"
var2 <- "Abbot_semiquantitative"
type <- "age group"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))
p1 <- res[[4]]

# Nov
var1 <- "Abbott_qualitative_november_2020"
var2 <- "Abbott_quantitative_november_2020"
type <- "age group"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))

# Diasorin
# May
var1 <- "Diasorin_IgG_qualitative"
var2 <- "Diasorin_IgG_semiquantitative"
type <- "age group"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))
p2 <- res[[4]]

# Nov
var1 <- "Diasorin_qualitative_november_2020"
var2 <- "Diasorin_quantitative_november_2020"
type <- "age group"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))

# Roche
# May
var1 <- "Roche_Total_qualitative"
var2 <- "Roche_Total_ICO"
type <- "age group"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))
p3 <- res[[4]]

# Nov
var1 <- "Roche_qualitative_november_2020"
var2 <- "Roche_quantitative_november_2020"
type <- "age group"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))
p4 <- res[[4]]

################################################################################
# Antibody titres by BMI                                                       #
################################################################################

# Abbott
# May
var1 <- "Abbot_qualitative"
var2 <- "Abbot_semiquantitative"
type <- "BMI category"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))
q1 <- res[[4]]

# Nov
var1 <- "Abbott_qualitative_november_2020"
var2 <- "Abbott_quantitative_november_2020"
type <- "BMI category"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))

# Diasorin
# May
var1 <- "Diasorin_IgG_qualitative"
var2 <- "Diasorin_IgG_semiquantitative"
type <- "BMI category"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))
q2 <- res[[4]]

# Nov
var1 <- "Diasorin_qualitative_november_2020"
var2 <- "Diasorin_quantitative_november_2020"
type <- "BMI category"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))

# Roche
# May
var1 <- "Roche_Total_qualitative"
var2 <- "Roche_Total_ICO"
type <- "BMI category"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))
q3 <- res[[4]]

# Nov
var1 <- "Roche_qualitative_november_2020"
var2 <- "Roche_quantitative_november_2020"
type <- "BMI category"

(res <- Anova_Kruskal_antibody_titres(var1, var2, type))

################################################################################
# Multiple linear regression: antibody ~ symptoms + age + BMI                  #
################################################################################

# Abbott
# May
idx <- which(ds$Abbot_qualitative == "Positive")

df <- data.frame(
  antibody_titre = ds$Abbot_semiquantitative[idx],
  symptomatic = as.factor(ds$Symptomatic_any_time_from_1_january_to_1_may[idx]),
  age_group = as.factor(ds$age_group[idx]),
  BMI = ds$BMI[idx])

fit <- lm(antibody_titre ~ BMI*symptomatic, data = df)
summary(fit)

fit <- lm(antibody_titre ~ BMI + symptomatic, data = df)
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "no", ])
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "yes", ])
summary(fit)

r1 <- ggplot(df[df$symptomatic == "yes", ],aes(y=antibody_titre,
                    x=BMI,
                    col=symptomatic))+
  geom_point()+
  ylab("Antibody titre May (Abbott)")+ xlab("BMI")+
  stat_smooth(method="lm",se=TRUE)+
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_blank())

# Nov
idx <- which(ds$Abbott_qualitative_november_2020 == "Positive")

df <- data.frame(
  antibody_titre = ds$Abbott_quantitative_november_2020[idx],
  symptomatic = as.factor(ds$Symptomatic_any_time_from_1_january_to_1_may[idx]),
  age_group = as.factor(ds$age_group[idx]),
  BMI = ds$BMI[idx])

fit <- lm(antibody_titre ~ BMI*symptomatic, data = df)
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "no", ])
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "yes", ])
summary(fit)

r2 <- ggplot(df[df$symptomatic == "yes",],
             aes(y=antibody_titre,
                    x=BMI,
                    col=symptomatic))+
  geom_point()+
  ylab("Antibody titre Nov 2020 (Abbott)")+ xlab("BMI")+
  stat_smooth(method="lm",se=TRUE)+
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_blank())

# DiaSorin
# May
idx <- which(ds$Diasorin_IgG_qualitative == "Positive")

df <- data.frame(
  antibody_titre = ds$Diasorin_IgG_semiquantitative[idx],
  symptomatic = as.factor(ds$Symptomatic_any_time_from_1_january_to_1_may[idx]),
  age_group = as.factor(ds$age_group[idx]),
  BMI = ds$BMI[idx])

fit <- lm(antibody_titre ~ BMI*symptomatic, data = df)
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "no", ])
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "yes", ])
summary(fit)

r3 <- ggplot(df[df$symptomatic == "yes", ],aes(y=antibody_titre,
                    x=BMI,
                    col=symptomatic))+
  geom_point()+
  ylab("Antibody titre May (Diasorin)")+ xlab("BMI")+
  stat_smooth(method="lm",se=TRUE)+
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_blank())

# Nov
idx <- which(ds$Diasorin_qualitative_november_2020 == "Positive")

df <- data.frame(
  antibody_titre = ds$Diasorin_quantitative_november_2020[idx],
  symptomatic = as.factor(ds$Symptomatic_any_time_from_1_january_to_1_may[idx]),
  age_group = as.factor(ds$age_group[idx]),
  BMI = ds$BMI[idx])

fit <- lm(antibody_titre ~ BMI*symptomatic, data = df)
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "no", ])
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "yes", ])
summary(fit)

# Roche
# May
idx <- which(ds$Roche_Total_qualitative == "Positive")

df <- data.frame(
  antibody_titre = ds$Roche_Total_ICO[idx],
  symptomatic = as.factor(ds$Symptomatic_any_time_from_1_january_to_1_may[idx]),
  age_group = as.factor(ds$age_group[idx]),
  BMI = ds$BMI[idx])

fit <- lm(antibody_titre ~ BMI*symptomatic, data = df)
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "no", ])
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "yes", ])
summary(fit)

r4 <- ggplot(df[df$symptomatic == "yes", ],aes(y=antibody_titre,
                    x=BMI,
                    color=symptomatic))+
  geom_point()+
  ylab("Antibody titre May (Roche)")+ xlab("BMI")+
  stat_smooth(method="lm",se=TRUE)+
  theme(legend.position = "bottom")

# Nov
idx <- which(ds$Roche_qualitative_november_2020 == "Positive")

df <- data.frame(
  antibody_titre = ds$Roche_quantitative_november_2020[idx],
  symptomatic = as.factor(ds$Symptomatic_any_time_from_1_january_to_1_may[idx]),
  age_group = as.factor(ds$age_group[idx]),
  BMI = ds$BMI[idx])

fit <- lm(antibody_titre ~ BMI*symptomatic, data = df)
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "no", ])
summary(fit)

fit <- lm(antibody_titre ~ BMI , data = df[df$symptomatic == "yes", ])
summary(fit)

r5 <- ggplot(df[df$symptomatic == "yes", ],aes(y=antibody_titre,
                    x=BMI,
                    color=symptomatic))+
  geom_point()+
  ylab("Antibody titre Nov (Roche)")+ xlab("BMI")+
  stat_smooth(method="lm",se=TRUE)+
  theme(legend.position = "bottom")

# Rows 1 and 2 of Figure S2----------------------------------------------------#

cols12 <- plot_grid(p1, p2, p3, q1, q2, q3,
                    labels = c('a', 'b', 'c', 'd', 'e','f'),
                    label_size = 8,
                    ncol = 3,
                    nrow = 2,
                    label_fontface = "bold")

saveRDS(cols12, file = file.path(path_figure, "temp.rds"))


# -----------------------------------------------------------------------------#


################################################################################
################################################################################

