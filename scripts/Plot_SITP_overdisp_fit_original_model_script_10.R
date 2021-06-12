################################################################################
# This code processes the posterior estimates obtained with the original       #
# model and generates Figure 4.                                                #
#                                                                              #                                                                             #
# Cite as:                                                                     #
# Dorigatti I et al, SARS-CoV-2 antibody dynamics, within-household            #
# transmission and the impact of contact tracing from community-wide           #
# serological testing in the Italian municipality of Vo'.                      #
# Pre-print available at: http://dx.doi.org/10.2139/ssrn.3792892               #
#                                                                              #
# Description:                                                                 #
# See README.md                                                                #
################################################################################

cat("\n\n### Running 'scripts/Plot_SITP_overdisp_fit_original_model_script_10.R' 
    ###\n\n")

# fix seed for reproducibility ------------------------------------------------#
set.seed(123)

# source script ---------------------------------------------------------------#
source("R/functions_plot_fitted_model.R")
source("R/functions_fit_model.R")

# set directories -------------------------------------------------------------#

dir_output  <- file.path("mcmc_posterior_chains_original")
dir_figures <- file.path("figures")

# read data -------------------------------------------------------------------#

res <- read.csv(file.path(dir_output, "mcmc_posterior_estimates_original.csv"),
                header = TRUE)

models <- read.csv("data/household_model/Model_variants_original.csv",
                   header = TRUE)

default_par <- read.csv("data/household_model/Default_parameters_original.csv",
                       header = TRUE)

data <- read.csv("tables/Household_final_size_Vo_baseline.csv",
                 header = TRUE)
data <- data[, -1]
colnames(data) <- seq(1,7,1)
rownames(data) <- seq(0,4,1)

colnames(models)[7] <- "p_i"
colnames(models)[9] <- "p_sd"

summary_res <- res

# Select chains with min deviance for each model ------------------------------#

selected_chains <- vector()
dev_score <- vector()
dev_l <- vector()
dev_u <- vector()
dic <- vector()

res1 <- res[which(res$estimate == "mean"), ]
res2 <- res[which(res$estimate == "lower_CrI"), ]
res3 <- res[which(res$estimate == "upper_CrI"), ]

for(i in 1:dim(models)[1]){
  
  idx1 <- which(res1[,1] == models$model[i])
  
  min_mean_dev <- min(as.numeric(as.character(res1$deviance[idx1])))
  idx <- which(res1$deviance[idx1] == min_mean_dev)
  selected_chains <- c(selected_chains, res1$chain[idx1[idx]]) 
  dev_score <- c(dev_score, min_mean_dev)
  
  dev_l <- c(dev_l, as.numeric(as.character(res2$deviance[idx1[idx]])))
  dev_u <- c(dev_u, as.numeric(as.character(res3$deviance[idx1[idx]])))
  dic <- c(dic, res1$DIC[idx1[idx]])
  
}

# selected chains
file_chains <- paste(paste(as.character(models$model),
                           paste("chain", selected_chains, sep = ""),
                           sep = "_"),
                     ".csv",
                     sep="")

#write.csv(file_chains, file.path(dir_output, "selected_chains.csv"))

df <- data.frame(model = models[,1],
                 deviance_mean = dev_score, 
                 deviance_lower = dev_l, 
                 deviance_upper = dev_u,
                 dic = dic)

# Plot fit to data as proportions ---------------------------------------------# 

# selected models 
min_dic <- which(df$dic < 16)

selected_mod <- df$model[which(df$dic < 16)]

# generate predictions 
j <- 1 # initialise
N <- 1000 # number of draws from posterior distribution 
res <- list()
estimatedParameters <- list()
expected <- list(matrix(0, ncol = dim(data)[2], nrow = dim(data)[2]))
expected_vec <- matrix(0, ncol = N, nrow = dim(data)[1]*dim(data)[2])
pred_mean <- list()
pred_CI_low <- list()
pred_CI_up <- list()
pre <- list()
pre_low <- list()
pre_up <- list()

for(i in which(df$dic < 16)){
  
  res[[j]] <- read.csv(file.path(dir_output, file_chains[i]), header = TRUE)
  
  res[[j]] <- res[[j]][, -c(1,2)]
  estimatedParameters[[j]] <- which(models[i, 2:9] == 1)
  
  for(n in 1:N){
    
    r <- sample(seq(1, dim(res[[j]])[1], 1), 1) 
    
    param <- default_par
    param[estimatedParameters[[j]]] <- res[[j]][r, ]
    
    expected[[n]] <- model_predictions(data, unlist(param)) 
    expected_vec[, n] <- as.vector(expected[[n]])
    
  }
  
  pred_mean[[j]] <- matrix(apply(expected_vec, 1, mean), 
                           ncol = dim(data)[2], 
                           byrow = FALSE) 
  
  pred_CI_low[[j]] <- matrix(apply(expected_vec, 1, quantile, c(0.025)), 
                             ncol = dim(data)[2], 
                             byrow = FALSE) 
  
  pred_CI_up[[j]] <- matrix(apply(expected_vec, 1, quantile, c(0.975)), 
                            ncol = dim(data)[2], 
                            byrow = FALSE) 
  
  pre[[j]] <- sweep(pred_mean[[j]], 2, colSums(pred_mean[[j]]), FUN = '/')
  pre_low[[j]] <- sweep(pred_CI_low[[j]], 2,colSums(pred_CI_low[[j]]),FUN = '/')
  pre_up[[j]] <- sweep(pred_CI_up[[j]], 2, colSums(pred_CI_up[[j]]), FUN = '/')
  
  j <- j+1
}

# data 
dat <- sweep(data, 2, colSums(data), FUN = '/')
dat_ci_low <- matrix(NA, ncol = dim(data)[2], nrow = dim(data)[1])
dat_ci_up <- matrix(NA, ncol = dim(data)[2], nrow = dim(data)[1])

for(j in 1:dim(data)[2]){
  
  for(i in 1:dim(data)[1]){
    
    ci <- binom.confint(data[i,j], sum(data[,j]), 
                        conf.level = 0.95, methods = "exact")
    
    dat_ci_low[i,j] <- ci$lower
    dat_ci_up[i,j] <- ci$upper
    
  }
}

# create data frame 
df_prop <- data.frame(household_size = rep(paste( "Households of size", 
                                                  rep(colnames(data), 
                                                      each = dim(data)[1]),
                                                  sep = " "), 9), 
                      
                      frequency = c(unlist(dat), 
                                    unlist(pre)), 
                      
                      infections = rep(rep(rownames(data), dim(data)[2]), 9),
                      
                      type = c(rep("observed", dim(data)[1]*dim(data)[2]),
                               rep(str_split_fixed(file_chains[which(df$dic<17)], 
                                                   "_", n = 2)[,1], 
                                   each = dim(data)[1]*dim(data)[2])), 
                      
                      lower = c(unlist(dat_ci_low), 
                                unlist(pre_low)), 
                      
                      upper = c(unlist(dat_ci_up), 
                                unlist(pre_up)))


p_model_fit <- ggplot(df_prop, aes(x = infections, y = frequency)) +
  geom_point(aes(x = infections, y = frequency, colour= type), 
             position = position_dodge(0.7))+
  geom_errorbar(aes(x = infections, ymin=lower, ymax=upper, colour = type), 
                position=position_dodge(0.7), 
                width = 0) +
  facet_wrap(~ household_size,  ncol=2)+
  theme(legend.position = "bottom")+
  xlab("SARS-CoV-2 positive household members")+
  ylab("proportion of households")+
  theme(legend.title=element_blank())

p_new2 <- p_model_fit +
  guides(fill = guide_legend(title.position = "top",
                             label.position = "right",
                             nrow = 3)) +
  theme(legend.direction = "horizontal", 
        legend.text = element_text(size = 6))

# save plot -------------------------------------------------------------------#
# ggsave(filename =
#          file.path(file.path(dir_output, "figures"),
#                    "Model_fit_proportions.tiff"),
#        plot = grid.draw(shift_legend(p_new2)),
#        device = "tiff",
#        width = 180, height = 120,
#        units = "mm", dpi = 300, limitsize = TRUE)
# -----------------------------------------------------------------------------#

# Plot household AR -----------------------------------------------------------# 

house_AR <- vector() 
house_AR_ci_low <- vector() 
house_AR_ci_up <- vector() 

for(i in 1:dim(data)[2]){
  
  n <- colSums(data[-1,])[i]
  d <- colSums(data)[i]
  
  bc <- binom.confint(n, d, conf.level = 0.95, methods = "exact")
  
  house_AR[i] <- bc$mean
  house_AR_ci_low[i] <- bc$lower 
  house_AR_ci_up[i] <- bc$upper
  
}

# best DIC model 
i <- 2 
models[i, ]
res <- read.csv(file.path("mcmc_posterior_chains_original", 
                          file_chains[i]), 
                header = TRUE)

res <- res[, -c(1,2)]
estimatedParameters <- which(models[i, 2:9] == 1)

for(n in 1:N){
  
  r <- sample(seq(1, dim(res)[1], 1), 1) 
  
  param <- default_par
  param[estimatedParameters] <- res[r, ]
  
  expected[[n]] <- model_predictions(data, unlist(param)) 
  expected_vec[, n] <- as.vector(expected[[n]])
  
}

house_AR_exp <- matrix(0, nrow = N, ncol = dim(data)[2])
for(j in 1:length(expected)){
  
  for(i in 1:dim(expected[[j]])[2]){
    
    n <- colSums(expected[[j]][-1,])[i]
    d <- colSums(expected[[j]])[i]
    
    house_AR_exp[j,i] <- n/d 
  }
}

house_AR_exp_ci <- apply(house_AR_exp, 2, FUN = quantile, c(0.025, 0.975))
house_AR_exp_mean <- apply(house_AR_exp, 2, FUN = mean)


df_house_AR <- data.frame(mean = c(house_AR,
                                   house_AR_exp_mean),
                          lower = c(house_AR_ci_low, 
                                    house_AR_exp_ci[1,]), 
                          upper = c(house_AR_ci_up, 
                                    house_AR_exp_ci[2,]),
                          size = as.factor(rep(seq(1,7,1),2)),
                          type = c(rep("observed",7), 
                                   rep("expected",7)))

plot_house_AR <- ggplot(df_house_AR, aes(x = size, y = mean)) +
  geom_point(aes(x = size, y = mean, colour= type), 
             position = position_dodge(0.1))+
  geom_errorbar(aes(x = size, ymin=lower, ymax=upper, colour = type), 
                position=position_dodge(0.1), 
                width = 0) +
  theme(legend.position = "bottom")+
  xlab("household size")+
  ylab("household attack rate (%)")+
  theme(legend.title=element_blank(),
        legend.position = c(0.1, 0.9))

# Plot individual-level attack rates ------------------------------------------#

# observed attack rate 
AR <- vector()
AR_low <- vector()
AR_up <- vector()
for(i in 1:dim(data)[2]){
  
  n <- sum(as.numeric(rownames(data))*data[,i])
  d <- (as.numeric(colnames(data))[i]*sum(data[,i]))
  
  ar <- binom.confint(n, d, methods = "exact")
  
  AR[i] <- ar$mean
  AR_low[i] <- ar$lower
  AR_up[i] <- ar$upper
  
}

# observed secondary attack rate 
SAR <- vector()
SAR_low <- vector()
SAR_up <- vector()
for(i in 2:dim(data)[2]){
  
  tmp <- (as.numeric(rownames(data))-1)*data[,i]
  
  n <- sum(tmp[2:length(tmp)])
  d <- ((as.numeric(colnames(data))-1)[i]*sum(data[2:length(tmp),i]))
  
  ar <- binom.confint(n, d, methods = "exact")
  
  SAR[i] <- ar$mean
  SAR_low[i] <- ar$lower
  SAR_up[i] <- ar$upper
  
}

# expected attack rates
AR_mat <- matrix(0, ncol = dim(data)[2], nrow = N)
SAR_mat <- matrix(0, ncol = dim(data)[2], nrow = N)
SAR_mat_num <- matrix(0, ncol = dim(data)[2], nrow = N)
SAR_mat_denom <- matrix(0, ncol = dim(data)[2], nrow = N)
num_inf_house <- matrix(0, ncol = dim(data)[2], nrow = N)

for(j in 1:N){
  
  dataexp <- expected[[j]]
  
  for(i in 1:dim(data)[2]){
    
    AR_mat[j,i] <- sum(as.numeric(rownames(data))*dataexp[,i])/
      (as.numeric(colnames(data))[i]*sum(dataexp[,i]))
    
    tmp <- (as.numeric(rownames(data))-1)*dataexp[,i]
    
    SAR_mat[j,i] <- sum(tmp[2:length(tmp)])/
      ((as.numeric(colnames(data))-1)[i]*sum(dataexp[2:length(tmp),i]))
    
    SAR_mat_denom[j,i] <- ((as.numeric(colnames(data))-1)[i]*
                             sum(dataexp[2:length(tmp),i]))
    
    SAR_mat_num[j,i] <- sum(tmp[2:length(tmp)])
    num_inf_house[j,i] <- sum(dataexp[2:dim(dataexp)[1],i])
  }
}

SAR_mat_num <- SAR_mat_num/num_inf_house

df_AR_final <- data.frame(value = c(AR,
                                    apply(AR_mat, 2, mean),
                                    SAR,
                                    c(NA, apply(SAR_mat[,2:7], 2, mean)),
                                    house_AR,
                                    house_AR_exp_mean),
                          
                          low = c(AR_low,
                                  apply(AR_mat, 2, quantile, 0.025),
                                  SAR_low,
                                  c(NA, apply(SAR_mat[,2:7], 2, quantile,0.025)),
                                  house_AR_ci_low,
                                  house_AR_exp_ci[1,]),
                          
                          up = c(AR_up,
                                 apply(AR_mat, 2, quantile, 0.975),
                                 SAR_up,
                                 c(NA, apply(SAR_mat[, 2:7], 2,quantile,0.975)),
                                 house_AR_ci_up,
                                 house_AR_exp_ci[2,]),
                          
                          size = as.factor(rep(seq(1,7,1), 6)),
                          
                          type = c(rep("individual attack rate", 14),
                                   rep("individual secondary attack rate", 14),
                                   rep("household attack rate", 14)),
                          
                          shape = rep(c(rep("observed", 7), 
                                        rep("expected", 7)), 3))

# plot for final figure 
plot_AR_final <- ggplot(df_AR_final)+
  geom_point(aes(x = size, 
                 y = value,
                 shape = type,
                 colour = shape,
                 fill = type), 
             size = 1, 
             position = position_dodge(0.7))+
  geom_errorbar(aes(x = size, 
                    ymin = low,
                    ymax = up,
                    shape = type, 
                    colour = shape,
                    fill = type), 
                width= 0, 
                position = position_dodge(0.7))+
  xlab("household size")+
  ylab("attack rate (%)")+
  theme_bw()+
  theme(axis.title=element_text(size=8),
        legend.position = c(0.28, 0.85),
        legend.spacing.y = unit(-0.2, "cm"),
        legend.key.height=unit(0.5,"line"), 
        legend.title = element_blank(),
        legend.key = element_rect(colour = "transparent"),
        legend.text = element_text(size=6),
        legend.background = element_blank())


df_AR_final_num <- data.frame(value = apply(SAR_mat_num[,2:7], 2, mean, na.rm = TRUE),
                              low = apply(SAR_mat_num[,2:7], 2, quantile, 0.025, na.rm = TRUE),
                              up = apply(SAR_mat_num[,2:7], 2, quantile, 0.975, na.rm = TRUE),
                              size = as.factor(rep(seq(2,7,1), 6)))

plot_AR_final_num <- ggplot(df_AR_final_num)+
  ylim(c(0,max(df_AR_final_num$up)))+
  geom_point(aes(x = size,
                 y = value),
             shape = 15)+
  geom_errorbar(aes(x = size,
                    ymin = low,
                    ymax = up),
                width= 0)+
  xlab("household size")+
  ylab("non-primary\ninfections")+
  theme_bw()+
  theme(axis.title=element_text(size=8),
        legend.title = element_blank(),
        legend.key = element_rect(colour = "transparent"),
        legend.background = element_blank())

# Plot SITP for best models ---------------------------------------------------#

# selected model 
min_dic <- df$dic[order(df$dic)][1:8]
selected_models <- 
  df$model[which(df$dic %in% min_dic)][order(df$dic[which(df$dic %in% min_dic)])]

pred_mean <- matrix(NA, ncol = 6, nrow = length(selected_models)) 
pred_CI_low <- matrix(NA, ncol = 6, nrow = length(selected_models)) 
pred_CI_up <- matrix(NA, ncol = 6, nrow = length(selected_models)) 
expected_vec <- matrix(NA, ncol = N, nrow = 6) 
overall_SITP <- matrix(NA, ncol = N, nrow = 7) 
SITP_pop <- vector()
SITP_mean <- vector()
SITP_CI_low <- vector()
SITP_CI_up <- vector()

# household size distribution 
prn <- colSums(sweep(data, 1, sum(data), FUN = '/'))

# mean household size
mean_house_size <- sum(seq(1,7,1)*prn)

fn <- (seq(1,7,1)*prn)/mean_house_size 

for(i in 1:length(selected_models)){
  
  selected_mod <- selected_models[i]
  
  res <- read.csv(file.path("mcmc_posterior_chains_original", 
                            file_chains[which(df$model == selected_mod)]), 
                  header = TRUE)
  
  res <- res[, -c(1,2)]
  estimatedParameters <- which(models[which(df$model == selected_mod),2:9] == 1)
  
  for(n in 1:N){
    
    r <- sample(seq(1, dim(res)[1], 1), 1) 
    
    param <- default_par
    param[estimatedParameters] <- res[r, ]
    param <- unlist(param)
    
    for(t in 2:7){
      expected_vec[t-1, n] <- 1 - MGF_gamma(1, t, param[1], param[3], param[4])
    }
    
    for(t in 1:7){
      
      overall_SITP[t, n] <- MGF_gamma(1, t, param[1], param[3], param[4])
      
    }
    
    SITP_pop[n] <- 1 - sum(overall_SITP[,n]*fn)   
    
  }
  
  pred_mean[i,] <- apply(expected_vec, 1, mean)
  pred_CI_low[i,] <- apply(expected_vec, 1, quantile, c(0.025))
  pred_CI_up[i,] <- apply(expected_vec, 1, quantile, c(0.975))
  
  SITP_mean[i] <- mean(SITP_pop) 
  SITP_CI_low[i] <- quantile(SITP_pop, 0.025) 
  SITP_CI_up[i] <- quantile(SITP_pop, 0.975) 
  
}

# overall SITP estimate -------------------------------------------------------#
df_SITP_overall <- data.frame(model = selected_models,
                              SITP_overall_mean = SITP_mean, 
                              SITP_overall_CI_low = SITP_CI_low,
                              SITP_overall_CI_up = SITP_CI_up)

print(df_SITP_overall)
#------------------------------------------------------------------------------#

df_SITP <- data.frame(model = rep(selected_models, each = 6), 
                      size = rep(seq(2,7,1), length(selected_models)), 
                      mean = as.vector(t(pred_mean)), 
                      low = as.vector(t(pred_CI_low)), 
                      up = as.vector(t(pred_CI_up)))

df_SITP$size <- as.factor(df_SITP$size)

#select models estimating alpha
Pmodels <- lapply(selected_models, gregexpr, pattern ='P')
idx <- which(unlist(Pmodels) == TRUE)

df_SITP$withP <- NA 
df_SITP$withP[which(df_SITP$model %in% selected_models[idx])] <- 
  "estimated alpha"
df_SITP$withP[which(!df_SITP$model %in% selected_models[idx])] <- 
  "fixed alpha"

df_SITP$withP <- as.factor(df_SITP$withP)

p_SITP1 <- ggplot(df_SITP[which(as.numeric(as.character(df_SITP$size)) > 1),], 
                  aes(x = size, y = mean, 
                      colour = model, group = model)) + 
  geom_point()+ 
  geom_line()+
  geom_ribbon(aes(x= size, ymin=low, ymax=up, 
                  fill = model, group = model), alpha=0.2, colour = NA)+
  facet_wrap(~ model,  ncol=4)+
  xlab("household size")+
  ylab("SITP")+
  theme(legend.position="none",
        axis.text=element_text(size=8), 
        strip.text.x = element_text(size = 8))

# save plot -------------------------------------------------------------------#
# ggsave(filename =
#          file.path(file.path(dir_output, "figures"),
#                    "SITP_best_models_indirect.tiff"),
#        plot = p_SITP1,
#        device = "tiff",
#        width = 180, height = 120,
#        units = "mm", dpi = 300, limitsize = TRUE)
# -----------------------------------------------------------------------------#

# plot overdispersion ---------------------------------------------------------#

shape_low <- matrix(NA, ncol = 1, nrow = length(selected_models)) 
shape_up <- matrix(NA, ncol = 1, nrow = length(selected_models)) 

cols <- gg_colour_hue(length(levels(ordered(df_SITP$model))))
ncol <- which(levels(ordered(df_SITP$model)) %in% selected_models)

df_overdisp <- list()
for(i in 1:length(selected_models)){
  
  selected_mod <- selected_models[i]
  
  res <- read.csv(file.path("mcmc_posterior_chains_original", 
                            file_chains[which(df$model == selected_mod)]), 
                  header = TRUE)
  
  res <- res[, -c(1,2)]
  estimatedParameters <- which(models[which(df$model == selected_mod),2:9] == 1)
  
  x <- seq(0.0001, 50, 0.1)
  xaxis <- matrix(NA, ncol = N, nrow = length(x))
  yaxis <- matrix(NA, ncol = N, nrow = length(x))
  nlab <- vector()
  
  for(n in 1:N){
    
    r <- sample(seq(1, dim(res)[1], 1), 1) 
    
    param <- default_par
    param[estimatedParameters] <- res[r, ]
    param <- unlist(param)
    
    R0 <- 2.4
    shape <- param[4]
    xaxis[,n] <- 1-pgamma(x, shape = shape, scale = R0/shape)
    yaxis[,n] <- 1-(1/R0) * 
      cumsum(x*dgamma(x, shape = shape, scale = R0/shape)*c(diff(x),diff(x)[1]))
    
    nlab <- c(nlab, rep(n, length(x)))
    
  }
  
  df_overdisp[[i]] <- data.frame(x = as.vector(xaxis), 
                                 y = as.vector(yaxis), 
                                 n = nlab,
                                 model = selected_mod,
                                 col = cols[ncol[i]])
  
}

# select trajectories within 95% CI (according to sum y values)
df_overdisp_low <- list()
df_overdisp_up <- list()
df_mu <- list()
for(j in 1:length(df_overdisp)){
  df_mu[[j]] <- data.frame(x = rep(0, 1000), y = rep(0, 1000))
}

for(j in 1:length(df_overdisp)){
  
  ds <- df_overdisp[[j]]
  chunks <- seq(500, dim(ds)[1], 500)
  for(i in 1:length(chunks)){
    
    if(i == 1){
      df_overdisp[[j]][1:chunks[1],6] <- sum(ds[1:chunks[1],2])
    }else{
      df_overdisp[[j]][(chunks[i-1]+1):chunks[i],6] <- 
        sum(ds[(chunks[i-1]+1):chunks[i],2])
    }
  }
  
  # median 
  lim <- quantile(df_overdisp[[j]][,6], c(0.025, 0.975))
  mu <- quantile(df_overdisp[[j]][,6], c(0.5))
  
  indexes <- which(abs(df_overdisp[[j]][,6] - mu) == 
                     min(abs(df_overdisp[[j]][,6] - mu) ))
  
  index_ordered <- indexes[order(df_overdisp[[j]][indexes,1])]
  
  if(length(index_ordered) > 1000){
    idx <- c(min(index_ordered),
             sample(index_ordered, 998, replace = FALSE),
             max(index_ordered))
  }else{
    idx <- index_ordered
  }
  
  df_mu[[j]][,1] <- df_overdisp[[j]][idx,1]
  df_mu[[j]][,2] <- df_overdisp[[j]][idx,2]
  
  # lower limit 
  indexes_low <- which(abs(df_overdisp[[j]][,6] - lim[1]) == 
                         min(abs(df_overdisp[[j]][,6] - lim[1]) ))
  
  indexes_low_ordered <- indexes_low[order(df_overdisp[[j]][indexes_low,6])]
  
  if(length(indexes_low_ordered) > 1000){
    idx <- c(min(indexes_low_ordered),
             sample(indexes_low_ordered, 998, replace = FALSE),
             max(indexes_low_ordered))
  }else{
    idx <- indexes_low_ordered
  }
  
  df_overdisp_low[[j]] <- df_overdisp[[j]][idx,]
  
  # upper limit 
  indexes_up <- which(abs(df_overdisp[[j]][,6] - lim[2]) == 
                        min(abs(df_overdisp[[j]][,6] - lim[2]) ))
  
  indexes_up_ordered <- indexes_up[order(df_overdisp[[j]][indexes_up,6])]
  
  if(length(indexes_up_ordered) > 1000){
    idx <- c(min(indexes_up_ordered),
             sample(indexes_up_ordered, 998, replace = FALSE),
             max(indexes_up_ordered))
  }else{
    idx <- indexes_up_ordered
  }
  
  df_overdisp_up[[j]] <- df_overdisp[[j]][idx,]
}

plot_overdisp <- list()

for(i in 1:length(selected_models)){
  
  h1 <- df_overdisp_low[[i]]$y[which(abs(df_overdisp_low[[i]]$x - 0.20) == 
                                       min(abs(df_overdisp_low[[i]]$x - 0.20)))]
  
  h2 <- df_overdisp_up[[i]]$y[which(abs(df_overdisp_up[[i]]$x - 0.20) == 
                                      min(abs(df_overdisp_up[[i]]$x - 0.20)))]
  
  mu <- df_mu[[i]]$y[which(abs(df_mu[[i]]$x - 0.20) == 
                             min(abs(df_mu[[i]]$x - 0.20)))]
  
  plot_overdisp[[i]] <- ggplot(df_overdisp_low[[i]])+
    geom_line(aes(x = x, y = y, group = model), size = 1.5, 
              colour = vir_lite(cols[ncol[i]]))+
    geom_line(data = df_overdisp_up[[i]],
              aes(x = x, y = y, group = model),  size = 1.5, 
              colour = vir_lite(cols[ncol[i]]))+
    geom_line(data = df_mu[[i]], 
              aes(x = x, y = y), size = 1.5, colour = cols[ncol[i]])+
    geom_vline(xintercept=0.20, linetype="dashed", size = 0.5, col= "grey")+
    geom_hline(yintercept= h1[1], 
               linetype="dashed", size = 0.5, col= "grey")+
    geom_hline(yintercept= h2[1], 
               linetype="dashed", size = 0.5, col= "grey")+
    ylab("")+
    xlab("")
  
  if(i == 1){
    message("prop transmission due to 20% most infectious individuals")
    message("mean, lower CI and upper CI")
  }
  print(c(mu[1], h1[1],h2[1]))
  
}

plot_over <- plot_grid(plot_overdisp[[1]],
                       plot_overdisp[[2]],
                       plot_overdisp[[3]],
                       plot_overdisp[[4]],
                       plot_overdisp[[5]],
                       plot_overdisp[[6]],
                       plot_overdisp[[7]],
                       plot_overdisp[[8]],
                       labels = c('a', 'b', 'c', 'd',
                                  'e','f', 'g', 'h'),
                       label_size = 8,
                       ncol = 4,
                       label_fontface = "bold")

y.grob <- textGrob("proportion of transmission due to X%  most infectious subjects",
                   gp=gpar(fontsize=12), rot=90)

x.grob <- textGrob("proportion X%",
                   gp=gpar(fontsize=12))

# save plot -------------------------------------------------------------------#
# ggsave(filename =
#         file.path(file.path(dir_output, "figures"),
#                   "Oversipersion_best_models.tiff"),
#       plot = grid.arrange(arrangeGrob(plot_over,
#                                       left = y.grob,
#                                       bottom = x.grob)),
#       device = "tiff",
#       width = 360, height = 240,
#       units = "mm", dpi = 150, limitsize = TRUE)
# -----------------------------------------------------------------------------#

# Selected plot for main figure -----------------------------------------------#
# V
i <- which(selected_models == "V")

h1 <- df_overdisp_low[[i]]$y[which(abs(df_overdisp_low[[i]]$x - 0.20) == 
                                     min(abs(df_overdisp_low[[i]]$x - 0.20),
                                         na.rm = TRUE))]

h2 <- df_overdisp_up[[i]]$y[which(abs(df_overdisp_up[[i]]$x - 0.20) == 
                                    min(abs(df_overdisp_up[[i]]$x - 0.20), 
                                        na.rm = TRUE))]


plot_overdisp_final <- ggplot(df_overdisp_low[[i]])+
  geom_line(aes(x = x, y = y, group = model), size = 1, 
            colour = vir_lite(cols[ncol[1]]))+
  geom_line(data = df_overdisp_up[[i]],
            aes(x = x, y = y, group = model),  size = 1, 
            colour = vir_lite(cols[ncol[1]]))+
  geom_line(data = df_mu[[i]], 
            aes(x = x, y = y), size = 1, colour = cols[ncol[1]])+
  geom_vline(xintercept=0.20, linetype="dashed", size = 0.5, col= "grey33")+
  geom_hline(yintercept= h1[1], 
             linetype="dashed", size = 0.5, col= "grey33")+
  geom_hline(yintercept= h2[1], 
             linetype="dashed", size = 0.5, col= "grey33")+
  ylab("prop transm due to X%")+
  xlab("proportion X%")+
  theme_bw()+
  theme(axis.title=element_text(size=8))

# Plot Figure 4 ---------------------------------------------------------------#

# SITP for models V, PV
df_SITP_selected <- df_SITP[which(df_SITP$model %in% c("V", "PV")),]
df_SITP_selected$cols <- NA
df_SITP_selected$cols[which(df_SITP_selected$model == "V")] <- cols[ncol[1]]
df_SITP_selected$cols[which(df_SITP_selected$model == "PV")] <- cols[ncol[3]]

plot_SITP_final <- ggplot(df_SITP_selected) +
  geom_errorbar(aes(x= size, 
                    ymin=low, 
                    ymax=up,
                    group = model), 
                position = position_dodge(0.5),
                width= 0, 
                size = 0.7, 
                colour = c(rep(cols[ncol[1]], 6),
                           rep(cols[ncol[3]], 6)))+ 
  geom_line(aes(x = size, 
                y = mean, 
                group = model),
            linetype = "longdash",
            position = position_dodge(0.5),
            size = 0.5, 
            colour = c(rep(cols[ncol[3]], 6),
                       rep(cols[ncol[1]], 6)))+
  geom_point(aes(x = size, 
                 y = mean,
                 shape = model,
                 colour = model),
             position = position_dodge(0.5), 
             size = 2,
             colour = c(rep(cols[ncol[1]], 6),
                        rep(cols[ncol[3]], 6)))+
  xlab("household size")+
  ylab("SITP")+
  theme_bw()+
  theme(legend.position="bottomleft",
        legend.title=element_blank(),
        axis.title=element_text(size=8))

# combine plots 
top_row <- plot_grid(plot_SITP_final,
                     plot_overdisp_final,
                     labels = c('a', 'c'), 
                     label_size = 8,
                     ncol=1,
                     label_fontface = "bold")

col2 <- plot_grid(plot_AR_final,
                  plot_AR_final_num,
                  labels = c('b', 'd'),
                  rel_heights = c(2.5,1),
                  label_size = 8,
                  ncol = 1,
                  label_fontface = "bold")

p_final <- plot_grid(top_row, 
                     col2, 
                     labels = c('', ''), 
                     label_size = 8, 
                     ncol = 2,
                     label_fontface = "bold")


# save plot -------------------------------------------------------------------#
ggsave(filename =
         file.path(dir_figures, "Figure_5.tiff"),
       plot = p_final,
       device = "tiff",
       width = 180, height = 100,
       units = "mm", dpi = 300, limitsize = TRUE)

ggsave(filename =
         file.path(dir_figures, "Figure_5.eps"),
       plot = p_final,
       device = "eps",
       width = 180, height = 100,
       units = "mm", dpi = 300, limitsize = TRUE)
# -----------------------------------------------------------------------------#

################################################################################
################################################################################