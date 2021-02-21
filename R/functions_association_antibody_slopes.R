################################################################################
# Compute mean half-life antibody decay rate (i) over all positive subjects in # 
# May; (ii) over all positive subjects in May and with decaying antibodies     # 
# between May and November; (iii) over all positive subjects in May and with   #
# antibodies not doubling between May and November.                            # 
################################################################################
compute_half_life <- function(ds, type){
  
  days <- as.numeric(as.Date("11/29/2020", "%m/%d/%Y") - 
                       as.Date("05/01/2020", "%m/%d/%Y"))
  
  if(type == "Abbott"){
    
    # positive in May 
    mean_p <- mean(ds$hl_Abbott[which(ds$Abbot_qualitative == "Positive")], 
                   na.rm = TRUE)
    median_p <- median(ds$hl_Abbott[which(ds$Abbot_qualitative == "Positive")], 
                       na.rm = TRUE)
    
    # positive in May and decreasing 
    mean_pn <- mean(ds$hl_Abbott[which(ds$Abbot_qualitative == "Positive" & 
                              ds$slope_Abbott <= 0)], na.rm = TRUE)
    median_pn <- median(ds$hl_Abbott[which(ds$Abbot_qualitative == "Positive"& 
                                  ds$slope_Abbott <= 0)], na.rm = TRUE)
    
    # positive in May and not doubling 
    mean_pd <- mean(ds$hl_Abbott[which(ds$Abbot_qualitative == "Positive" & 
                              ds$slope_Abbott <= log(2)/days)], na.rm = TRUE)
    median_pd <- median(ds$hl_Abbott[which(ds$Abbot_qualitative == "Positive"& 
                                  ds$slope_Abbott <= log(2)/days)], na.rm = TRUE)
    
  }else if(type == "DiaSorin"){
    
    # positive in May 
    mean_p <- mean(ds$hl_Diasorin[which(ds$Diasorin_IgG_qualitative == "Positive" & 
                                          is.infinite(ds$hl_Diasorin) == FALSE) ], 
                   na.rm = TRUE)
    median_p <- median(ds$hl_Diasorin[which(ds$Diasorin_IgG_qualitative == "Positive"& 
                                              is.infinite(ds$hl_Diasorin) == FALSE)], 
                       na.rm = TRUE)
    
    # positive in May and decreasing 
    mean_pn <- mean(ds$hl_Diasorin[which(ds$Diasorin_IgG_qualitative == "Positive" & 
                              ds$slope_Diasorin <= 0 & 
                                is.infinite(ds$hl_Diasorin) == FALSE)], na.rm = TRUE)
    median_pn <- median(ds$hl_Diasorin[which(ds$Diasorin_IgG_qualitative == "Positive"& 
                                  ds$slope_Diasorin <= 0 & 
                                    is.infinite(ds$hl_Diasorin) == FALSE)], na.rm = TRUE)
    
    # positive in May and not doubling 
    mean_pd <- mean(ds$hl_Diasorin[which(ds$Diasorin_IgG_qualitative == "Positive" & 
                              ds$slope_Diasorin <= log(2)/days & 
                                is.infinite(ds$hl_Diasorin) == FALSE)], na.rm = TRUE)
    median_pd <- median(ds$hl_Diasorin[which(ds$Diasorin_IgG_qualitative == "Positive"& 
                                  ds$slope_Diasorin <= log(2)/days & 
                                    is.infinite(ds$hl_Diasorin) == FALSE)], na.rm = TRUE)
    
    
  }else if(type == "Roche"){
    
    # positive in May 
    mean_p <- mean(ds$hl_Roche[which(ds$Roche_Total_qualitative == "Positive")], 
                   na.rm = TRUE)
    median_p <- median(ds$hl_Roche[which(ds$Roche_Total_qualitative == "Positive")], 
                       na.rm = TRUE)
    
    # positive in May and decreasing 
    mean_pn <- mean(ds$hl_Roche[which(ds$Roche_Total_qualitative == "Positive" & 
                                           ds$slope_Roche <= 0)], na.rm = TRUE)
    median_pn <- median(ds$hl_Roche[which(ds$Roche_Total_qualitative == "Positive"& 
                                               ds$slope_Roche <= 0)], na.rm = TRUE)
    
    # positive in May and not doubling 
    mean_pd <- mean(ds$hl_Roche[which(ds$Roche_Total_qualitative == "Positive" & 
                                           ds$slope_Roche <= log(2)/days)], na.rm = TRUE)
    median_pd <- median(ds$hl_Roche[which(ds$Roche_Total_qualitative == "Positive"& 
                                               ds$slope_Roche <= log(2)/days)], na.rm = TRUE)

  }else{
    
    message("type not valid")
    return()
    
  }
  
  return(rbind(c(mean_p, median_p), 
               c(mean_pn, median_pn),
               c(mean_pd, median_pd)))

}

################################################################################
# Compute mean half-life and 95% CI by bootstrapping                           #
################################################################################
bootstrap_hl <- function(ds, type){
  
  ds_orig <- ds
  
  hl <- vector()
  hl_db <- vector()
  
  for(i in 1:1000){
    
    idx <- sample(seq(1, dim(ds_orig)[1], 1), dim(ds_orig)[1], replace = TRUE)
    ds <- ds_orig[idx,]
    res <- compute_half_life(ds, type)
    
    hl[i] <- res[1,2]
    hl_db[i] <- res[3,2]
    
  }
  
  return(rbind(quantile(hl, c(0.025, 0.975)),
               quantile(hl_db, c(0.025, 0.975))))
  
}

################################################################################
# Test for difference in response, two categories (t-test/Wilcoxon test)       #
################################################################################
test_decay_rate_half_life <- function(var_qual, var_quant, type, type2){
  
  ncl <- which(colnames(ds) == var_qual)
  ncl1 <- which(colnames(ds) == var_quant)
  
  if(str_detect(var_qual, "Diasorin", negate = FALSE) & 
     str_detect(var_quant, "Diasorin", negate = FALSE)){
    
    out <- "DiaSorin"
    
  }else if(str_detect(var_qual, "Roche", negate = FALSE) & 
           str_detect(var_quant, "Roche", negate = FALSE)){
    
    out <- "Roche"
    
  }else if(str_detect(var_qual, "Abbot", negate = FALSE) & 
           str_detect(var_quant, "Abbot", negate = FALSE)){
    
    out <- "Abbott"
    
  }else{
    
    message("Check variable names, they are not matching!")
    return()
    
  }
  
  idx <- which(ds[,ncl] == "Positive")
  
  if(type == "pcr positive"){
    
    idx1 <- which(ds$first_sampling[idx] == "Positive")
    idx2 <- which(ds$second_sampling[idx] == "Positive")
    
    df <- data.frame(
      
      antibody_titres = c(ds[idx[idx1], ncl1], 
                          ds[idx[idx2], ncl1]), 
      
      type = c(rep("positive at\nfirst survey", length(ds[idx[idx1], ncl1])), 
               rep("positive at\nsecond survey", length(ds[idx[idx2], ncl1]))))
    
  }else if(type == "symptom occurrence"){
    
    idx1 <- which(ds$Symptomatic_any_time_from_1_january_to_1_may[idx] == "yes")
    idx2 <- which(ds$Symptomatic_any_time_from_1_january_to_1_may[idx] != "yes")
    
    df <- data.frame(
      
      antibody_titres = c(ds[idx[idx1], ncl1], 
                          ds[idx[idx2], ncl1]), 
      
      type = c(rep("symptomatic", length(ds[idx[idx1], ncl1])), 
               rep("asymptomatic", length(ds[idx[idx2], ncl1]))))
    
    rmidx <- which(is.infinite(df$antibody_titres))
    
    if(length(rmidx) > 0){
      df <- df[-rmidx, ]
    }
    
  }else if(type == "hospitalized"){
    
    idx1 <- which(ds$hospitalized[idx] == "yes")
    idx2 <- which(ds$hospitalized[idx] != "yes")
    
    df <- data.frame(
      
      antibody_titres = c(ds[idx[idx1], ncl1], 
                          ds[idx[idx2], ncl1]), 
      
      type = c(rep("hospitalized", length(ds[idx[idx1], ncl1])), 
               rep("non-hospitalized", length(ds[idx[idx2], ncl1]))))
    
    rmidx <- which(is.infinite(df$antibody_titres))
    
    if(length(rmidx) > 0){
      df <- df[-rmidx, ]
    }
    
  }else if(type == "sex"){
    
    idx1 <- which(ds$Gender[idx] == "F")
    idx2 <- which(ds$Gender[idx] == "M")
    
    df <- data.frame(
      
      antibody_titres = c(ds[idx[idx1], ncl1], 
                          ds[idx[idx2], ncl1]), 
      
      type = c(rep("Female", length(ds[idx[idx1], ncl1])), 
               rep("Male", length(ds[idx[idx2], ncl1]))))
    
    rmidx <- which(is.infinite(df$antibody_titres))
    
    if(length(rmidx) > 0){
      df <- df[-rmidx, ]
    }
    
  }else{
    
    message("check type, entry is not valid!")
    return()
  }
  
  # test normality of the data 
  normality <- 1
  for(i in seq_along(levels(df$type))){
    
    data <- df$antibody_titres[which(df$type == levels(df$type)[i])]
    data <- data[!is.na(data)]
    
    if(length(data) > 2){
      
      sh <- shapiro.test(data)
      
      if(sh$p.value < 0.05){
        normality <- 0 # FALSE
        message(paste("normality is not supported for category",
                      levels(df$type)[i], 
                      sep = " "))
      }else{
        message("data are normally distributed (univariate test)")
      }
    }
  }
  
  ttest <- t.test(antibody_titres ~ type, data = df)
  wtest <- wilcox.test(antibody_titres ~ type, data = df) 
  
  r <- ggplot(df, aes(x = type, 
                      y = antibody_titres, 
                      fill = type))+
    geom_boxplot()+
    ylab(paste(paste(paste(paste("Antibody", type2, sep = " "), "(",sep = " "),
               out, sep= ""), ")", sep=""))+ 
    xlab("")+
    theme_bw()+
    theme(legend.position = "none")
  
  return(list(ttest, r, wtest))
  
}

################################################################################
# Adding case counts in boxplot                                                #
################################################################################
stat_box_data <- function(y, upper_limit = max(df$antibody_titres) * 1.15) {
  return( 
    data.frame(
      y = -0.018, 
      label = length(y)
    )
  )
}

################################################################################
# Test for difference in response, multiple categories (Anova/Kruskal-Wallis)  #
################################################################################
Anova_Kuskal_decay_rate_half_life <- function(var_qual, var_quant, type, type2){
  
  ncl <- which(colnames(ds) == var_qual)
  ncl1 <- which(colnames(ds) == var_quant)
  
  if(str_detect(var_qual, "Diasorin", negate = FALSE) & 
     str_detect(var_quant, "Diasorin", negate = FALSE)){
    
    out <- "DiaSorin"
    
  }else if(str_detect(var_qual, "Roche", negate = FALSE) & 
           str_detect(var_quant, "Roche", negate = FALSE)){
    
    out <- "Roche"
    
  }else if(str_detect(var_qual, "Abbot", negate = FALSE) & 
           str_detect(var_quant, "Abbot", negate = FALSE)){
    
    out <- "Abbott"
    
  }else{
    
    message("Check variable names, they are not matching!")
    return()
    
  }
  
  idx <- which(ds[,ncl] == "Positive")
  
  if(type == "age group"){
    
    df <- data.frame(antibody_titres = ds[idx, ncl1], 
                     type = ds$age_group[idx])
    
    df$type <- factor(df$type, 
                      levels = c("00-10", 
                                 "11-20", 
                                 "21-30",
                                 "31-40",
                                 "41-50",
                                 "51-60",
                                 "61-70",
                                 "71-80",
                                 "81-90",
                                 "91+"))
    
    xlab <- "age group (years)"
    
  }else if(type == "BMI category"){
    
    df <- data.frame(antibody_titres = ds[idx, ncl1], 
                     type =ds$BMI_category[idx])
    
    df$type <- factor(df$type, 
                      levels = c("Underweight",
                                 "Normal", 
                                 "Overweight",
                                 "Obese"))
    
    xlab <- "BMI category"
    
  }else{
    
    message("check type, entry not valid!")
    return()
  }
  
  rmidx <- which(is.infinite(df$antibody_titres))
  
  if(length(rmidx) > 0){
    df <- df[-rmidx, ]
  }
  
  t <- ggplot(df, aes(x = type, 
                      y = antibody_titres,
                      fill = type))+
    geom_boxplot()+
    ylab(paste(paste(paste(paste("Antibody", type2, sep = " "),
                           "(", sep = " "), 
                     out, sep= ""), ")", sep=""))+ 
    xlab(xlab)+
    stat_summary(
      fun.data = stat_box_data, 
      geom = "text", 
      hjust = 0.5,
      vjust = 0.9,
      size = 3) +
    theme_bw()+
    theme(legend.position = "none")
  
  if(type == "age group"){
    t <- t +
      theme(axis.text.x = element_text(angle = 30))
  }
  
  
  # test normality of the data 
  normality <- 1
  for(i in seq_along(levels(df$type))){
    
    data <- df$antibody_titres[which(df$type == levels(df$type)[i])]
    data <- data[!is.na(data)]
    
    if(length(data) > 2){
      
      sh <- shapiro.test(data)
      if(sh$p.value < 0.05){
        normality <- 0 # FALSE
        message(paste("normality is not supported for category",
                      levels(df$type)[i], 
                      sep = " "))
      }
    }
  }
  
  equal_variance <- 0
  if(normality == 1){
    # test equal variance 
    equal_var <- leveneTest(antibody_titres ~ type, data = df)
    if(equal_var$`Pr(>F)`[1] < 0.05){
      equal_variance <- 0
      message("equal within group variance is not supported")
    }else{
      equal_variance <- 1
      message("equal within group variance is supported")
    }
  }
  
  if(normality == 1 & equal_variance == 1){

      # anova with equal variances
      res.aov <- aov(antibody_titres ~ type, data = df)
      sumaov <- summary(res.aov)
      
      if(sumaov[[1]][5][[1]][1] < 0.05) {
        message("statistically significant difference (equal variance)")
      }
      
      res <- sumaov[[1]][5][[1]][1]
      test <- "one-way Anova equal within group variance"
      
      sumtur <- TukeyHSD(res.aov)
      
      if(length(which(sumtur$type[,4] < 0.05)) > 0) {
        message("statistically significant difference in:")
        message(rownames(sumtur$type)[which(sumtur$type[,4] < 0.05)])
      }
  }else{
    
    # non-paramteric Kruskal-Wallis test 
    kruskal <- kruskal.test(antibody_titres ~ type, data = df)
    
    if(kruskal$p.value< 0.05) {
      message("statistically significant difference (non-parametric)")
    }else{
      message("no significant difference (non-parametric)")
    }
    
    res <- kruskal$p.value
    test <- "Kruskal-Wallis test"
    sumtur <- NA
    
    # remove classes with less than 3 observations and rerun test
    tb <- table(df$type)
    
    idx <- which(df$type %in% names(tb)[which(tb < 3)])
    
    if(length(idx) > 0 ){
      df1 <- df[-idx, ]
      
      kruskal1 <- kruskal.test(antibody_titres ~ type, data = df1)
      
      if(kruskal1$p.value< 0.05) {
        message("statistically significant difference (non-parametric) having 
                excluded the classes with < 3 observations")
        print(kruskal1)
      }else{
        message("no significant difference (non-parametric) having excluded the 
                classes with < 3 observations")
      }
    }
  }
  
  return(list(res, test, sumtur, t)) 
  
}

################################################################################
################################################################################








