################################################################################
# Reproducing ggplot standard colours for factors                              #
################################################################################

gg_colour_hue <- function(n){
  
  hues= seq(15, 375, length = n + 1)
  
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# Calculate raw seroprevalence [= P(test positive)] and Binomial exact  95% CI #
################################################################################

raw_seroprevalence <- function(x){
  
  tab <- table(x)
  
  num <- tab[which(rownames(tab) == "Positive")]
  denom <- sum(tab[which(rownames(tab) %in% c("Positive", "Negative"))])
  
  raw_serop  <- binom.confint(num, 
                              denom, 
                              conf.level = 0.95,
                              methods = "exact")
  
  return(c(raw_serop$mean, raw_serop$lower, raw_serop$upper))
  
}

################################################################################
# Calculate adjusted seroprevalence (by sensitivity and specificity)           #
# with Binomial exact  95% CI                                                  #
################################################################################

adjusted_seroprevalence <- function(x, se, sp){
  
  tab <- table(x)
  
  num <- tab[which(rownames(tab) == "Positive")]
  denom <- sum(tab[which(rownames(tab) %in% c("Positive", "Negative"))])
  
  raw_serop  <- binom.confint(num, 
                              denom, 
                              conf.level = 0.95,
                              methods = "exact")
  
  
  num <- ((raw_serop$mean + sp -1)/(se + sp -1))*denom
  
  adj_serop  <- binom.confint(num, 
                              denom, 
                              conf.level = 0.95,
                              methods = "exact")
  
  return(c(adj_serop$mean, adj_serop$lower, adj_serop$upper))
  
}

################################################################################
# Calculate sensitivity and specificity, mean and Binomial exact  95% CI       #
################################################################################

estimate_performance <- function(x, y){
  
  tab <- table(x)
  num <- tab[which(rownames(tab) == "Positive")]
  denom <- sum(tab[which(rownames(tab) %in% c("Positive", "Negative"))])
  
  tp <- table(x, y, useNA = "always")
  tp <- tp[-which(!rownames(tp) %in% c("Positive", "Negative")),]
  
  sens <- binom.confint(tp["Positive",1], 
                        as.numeric(colSums(tp)[1]), 
                        methods = "exact")
  
  spec <- binom.confint(tp["Negative",2], 
                        as.numeric(colSums(tp)[2]), 
                        methods = "exact")
  
  performance <- data.frame(sens =c(sens$mean, sens$lower, sens$upper), 
                            spec =c(spec$mean, spec$lower, spec$upper))
  
  return(list(round(performance, 3)))
  
}

################################################################################
# Calculate adjusted seroprevalence, PPV and NPV, mean and bootstrap 95% CI    #
################################################################################

estimate_serop_ppv_npv <- function(x, y){
  
  tab <- table(x)
  num <- tab[which(rownames(tab) == "Positive")]
  denom <- sum(tab[which(rownames(tab) %in% c("Positive", "Negative"))])
  
  tp <- table(x, y, useNA = "always")
  tp <- tp[-which(!rownames(tp) %in% c("Positive", "Negative")),]
  
  sens <- binom.confint(tp["Positive",1], 
                        as.numeric(colSums(tp)[1]), 
                        methods = "exact")
  
  spec <- binom.confint(tp["Negative",2], 
                        as.numeric(colSums(tp)[2]), 
                        methods = "exact")
  
  falp <- binom.confint(tp["Positive",2], 
                        as.numeric(colSums(tp)[2]), 
                        methods = "exact")
  
  faln <- binom.confint(tp["Negative",1], 
                        as.numeric(colSums(tp)[1]), 
                        methods = "exact")
  
  # mean estimates 
  seroprev_corr_mean <- 
    ((num/denom) - falp$mean)/(sens$mean - falp$mean)
  
  ppv_mean <- (sens$mean/(num/denom))*((num/denom) - falp$mean)/
    (sens$mean - falp$mean)
  
  npv_mean <- (spec$mean/((denom-num)/denom))*(1-((num/denom) - falp$mean)/
                           (sens$mean - falp$mean))

  # bootstrap
  N <- 1000
  seroprev_corr <- vector()
  ppv <- vector()
  npv <- vector()

  for(n in 1:N){
    
    idx <- sample(seq(1,length(x),1), length(x), replace = TRUE)
    x1 <- x[idx]
    y1 <- y[idx]
    
    tp <- table(x1, y1, useNA = "always")
    tp <- tp[-which(!rownames(tp) %in% c("Positive", "Negative")),]
    
    num <- sum(tp["Positive",])
    denom <- sum(tp)
    
    sens <- binom.confint(tp["Positive",1], 
                          as.numeric(colSums(tp)[1]), 
                          methods = "exact")
    
    spec <- binom.confint(tp["Negative",2], 
                          as.numeric(colSums(tp)[2]), 
                          methods = "exact")
    
    falp <- binom.confint(tp["Positive",2], 
                          as.numeric(colSums(tp)[2]), 
                          methods = "exact")
    
    faln <- binom.confint(tp["Negative",1], 
                          as.numeric(colSums(tp)[1]), 
                          methods = "exact")
    
    # mean estimates 
    seroprev_corr[n] <- 
      ((num/denom) - falp$mean)/(sens$mean - falp$mean)
    
    ppv[n] <- (sens$mean/(num/denom))*((num/denom) - falp$mean)/
      (sens$mean - falp$mean)
    
    npv[n] <- (spec$mean/((denom-num)/denom))*(1-((num/denom) - falp$mean)/
                                                 (sens$mean - falp$mean))
    
  }
  
  # extract 2.5 - 97.5 percentiles
  df <- data.frame(
    serop = round(c(seroprev_corr_mean, 
                    quantile(seroprev_corr, c(0.025, 0.975))), 3), 
    ppv = round(c(ppv_mean, quantile(ppv, c(0.025, 0.975))), 3), 
    npv = round(c(npv_mean, quantile(npv, c(0.025, 0.975))), 3))
  
  return(df)
  
}

################################################################################
# Create output Table S1                                                       #
################################################################################

estimate_sens_spec_ppv_npv <- function(x,y){
  
  sesp <- estimate_performance(x,y)
  est <- estimate_serop_ppv_npv(x,y)
  
  return(list(rbind(t(sesp[[1]]), t(est[,2:3])), est[,1]))
  
}

################################################################################
################################################################################
