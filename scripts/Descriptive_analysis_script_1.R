################################################################################
# This code outputs the case counts in Figure 1 and Table 4. It also outputs a #
# contigency table containing the number of SARS-CoV-2 infections by household #
# size (Table 6), Figure S3 and the attack rates by household size             #
#                                                                              #
# Cite as:                                                                     #
# Dorigatti I et al, SARS-CoV-2 antibody dynamics, within-household            #
# transmission and the impact of contact tracing from community-wide           #
# serological testing in the Italian municipality of Vo'.                      #
#                                                                              #
# Description:                                                                 #
# See README.md                                                                #
################################################################################

cat("\n\n### Running 'scripts/Descriptive_analysis_script_1.R'\n\n")


# source scripts
source("R/functions_descriptive_analysis.R")

# read data
ds <- read.csv("data/Vo_serology_data.csv", header = TRUE)

# Create output folders -------------------------------------------------------#
dir_output  <- file.path("tables")
dir_figures <- file.path("figures")
dir.create(dir_output,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

# Data counts for Figure 1 ----------------------------------------------------#

output_case_counts(ds)
output_case_counts(ds[which(ds$town == "VO"),])

# Analyse data ----------------------------------------------------------------#

# select Vo residents
ds <- ds[which(ds$town == "VO"), ]

# identify subjects of households with all members screened for SARS-CoV-2
ds$included <- 0
for(i in 1:dim(ds)[1]){

  if(!is.na(ds$household_id[i])){

    if(length(which(ds$household_id == ds$household_id[i])) ==
       ds$household_size[i]) ds$included[i] <- 1
  }

}

# select subjects included in household transmission analysis
ds1 <- ds[which(ds$included == 1), ]

# households size distribution
df <- data.frame(size = names(table(table(ds1$household_id))),
                 frequency = c(table(table(ds1$household_id))))

# plot household size distribution for households included in the survey
p <- ggplot(df, aes(x = size, y = frequency)) +
    geom_bar(stat="identity", position=position_dodge())+
    theme(legend.position = "bottom")+
    xlab("household size")+
    ylab("frequency")+
  theme_bw()

# household size distribution by outcome
house_id <- vector()
age_group <- vector()
size <- vector()
gender <- vector()
outcome <- vector()

freq <- table(ds1$household_id)

for(i in 1:max(as.numeric(names(table(freq))))){

  house <- names(which(freq == i))

  house_id <- c(house_id,
                ds1$household_id[which(ds1$household_id %in% house)])

  age_group <- c(age_group,
                 as.character(ds1$age_group[which(ds1$household_id %in% house)]))

  size <- c(size, rep(i,length(which(ds1$household_id %in% house))))

  gender <- c(gender, ds1$Gender[which(ds1$household_id %in% house)])

  outcome <- c(outcome,
                    ds1$Groundtruth[
                      which(ds1$household_id %in% house)])
}

# redefine outcome levels
outcome[which(is.na(outcome))] <- "Negative"
outcome[which(outcome == 1)] <- "Positive"

# redefine gender levels
gender[which(gender == 1)] <- "Female"
gender[which(gender == 2)] <- "Male"

# create linelist
df1 <- data.frame(household_id = house_id,
                  size = paste("households of size", size, sep = " "),
                  age_group = as.factor(age_group),
                  gender = as.factor(gender),
                  size_plot = size,
                  outcome = outcome)

# plot age distribution of households by size
q <- ggplot(df1, aes(x = age_group, fill = outcome)) +
  geom_bar(aes(fill = outcome))+
  facet_wrap(~ size,  ncol=3, scales="free_y")+
  theme(legend.position = "none")+
  xlab("age group (years)")+
  ylab("frequency")+
  labs(fill = "outcome")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30))+
  scale_fill_manual(values = rev(gg_colour_hue(2)),
                    na.value = "grey50")

df1$size <- unlist(lapply(df1$size, str_remove, pattern = "households of size "))

# define household dataset
df2 <- data.frame(household_id = unique(ds1$household_id),
                  size = 0,
                  pos = 0,
                  neg = 0,
                  check = 0)

v <- unique(ds1$household_id)

for(i in 1:length(v)){

  df2$size[i] <- length(which(ds1$household_id == v[i]))

  df2$pos[i] <- length(which(ds1$Groundtruth[
    which(ds1$household_id == v[i])] == 1))

  df2$neg[i] <- length(which(is.na(ds1$Groundtruth[
    which(ds1$household_id == v[i])])))

  if(df2$pos[i] + df2$neg[i] == df2$size[i]) df2$check[i] <- 1
}

# ---------------------- Contingency table  -----------------------------------#

# output final outcome contingency table
cont_tab <- matrix(rep(0, max(df2$size)*(max(df2$pos)+1)), ncol = max(df2$size))
size <- seq(1,max(df2$size),1)
pos <- seq(0,max(df2$pos),1)

colnames(cont_tab) <- size
rownames(cont_tab) <- pos

for(i in 1:length(size)){

  v <- table(df2$pos[which(df2$size == size[i])])
  cont_tab[1:length(v),i] <- v

}

write.csv(cont_tab, file.path(dir_output, "Household_final_size_Vo_baseline.csv"))

# -----------------------------------------------------------------------------#

# attack rate
AR <- vector()
for(i in 1:length(size)){
  AR[i] <- sum(as.numeric(rownames(cont_tab))*cont_tab[,i])/
    (as.numeric(colnames(cont_tab))[i]*sum(cont_tab[,i]))
}

# secondary attack rate
SAR <- vector()
for(i in 2:length(size)){
  tmp <- (as.numeric(rownames(cont_tab))-1)*cont_tab[,i]
  SAR[i] <- sum(tmp[2:length(tmp)])/
    ((as.numeric(colnames(cont_tab))-1)[i]*sum(cont_tab[2:length(tmp),i]))
}

# plot final size distribution and attack rates
s <- ggplot() +
  geom_bar(data = df2,
           aes(x = as.factor(size), fill = as.factor(pos)),
           position = position_fill(reverse = TRUE))+
  scale_fill_manual(values = wes_palette("Zissou1", n = 5))+
  theme(legend.position = "bottom")+
  xlab("household size")+
  ylab("proportion")+
  labs(fill = "positive")+
  theme_bw()

t <-s +
  geom_point(aes(x = size[1:7], y = 2*SAR[1:7]), shape = 21, size = 3,
             fill = "white") +
  geom_line(aes(x = size[1:7], y = 2*SAR[1:7]), linetype = "longdash",
            col="dodgerblue4")+
  geom_point(aes(x = size[1:7], y = 2*AR[1:7]), shape = 22, size = 3,
             fill = "yellow") +
  geom_line(aes(x = size[1:7], y = 2*AR[1:7]), linetype = "dashed",
            col="dodgerblue4")+
  scale_y_continuous(sec.axis = sec_axis(~./2,
                                         name = "Attack rate"))+
  theme(axis.title.y = element_text(color = "black", size=13),
        axis.title.y.right = element_text(color = "dodgerblue4", size=13),
        axis.text.y.right = element_text(color = "dodgerblue4"))

# ------------------------------ Figure S3 ------------------------------------#
ggsave(filename = file.path(dir_figures, "Figure_S3.tiff"),
        plot = plot_grid(plot_grid(p, t,
                 labels = c('a', 'b'),
                 label_size = 8,
                 ncol = 2,
                 rel_widths = c(1,2),
                 label_fontface = "bold"),
                 q,
                 labels = c('', 'c'),
                 label_size = 8, ncol = 1,
                 rel_heights = c(1,2),
                 label_fontface = "bold"),
       device = "tiff",
       width = 270, height = 270,
       units = "mm", dpi = 300, limitsize = TRUE)
# -----------------------------------------------------------------------------#

# Overall household attack rate
binom.confint(sum(cont_tab[-1,]), sum(cont_tab),
              conf.level = 0.95, methods = "exact")

################################################################################
################################################################################

