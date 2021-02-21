################################################################################
# Test for association antibody titres in time (days from symptoms onset)      #
################################################################################
association_titres_days_symptom_onset <- function(var_qual, var_quant){
  
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
  
  idx1 <- which(ds[, ncl] == "Positive")
  idx <- which(ds$first_symptoms_date[idx1] != "")
  
  onset_dates <- ds$first_symptoms_date[idx1[idx]]
  
  if(str_detect(var_qual, "2020", negate = TRUE) & 
     str_detect(var_quant, "2020", negate = TRUE)){
    
          time_span <- as.Date("05/02/2020", "%m/%d/%Y") - 
            as.Date(onset_dates, "%m/%d/%Y")
          
  }else if(str_detect(var_qual, "2020", negate = FALSE) & 
           str_detect(var_quant, "2020", negate = FALSE)){
    
    time_span <- as.Date("11/29/2020", "%m/%d/%Y") - 
      as.Date(onset_dates, "%m/%d/%Y")
    
  }else{
    
    message("Check variable names, they are not matching!")
    return()
  }
  
  ab <- ds[idx1[idx], ncl1]
  
  df_ab_time <- data.frame(days_since_symptom_onset = as.numeric(time_span), 
                           antibody_titres = ab)
  
  # test trend 
  fit <- lm(antibody_titres ~ days_since_symptom_onset, data = df_ab_time)
  sumfit <- summary(fit)
  
  # plot  
  q <- ggplot(df_ab_time,
              aes(x = days_since_symptom_onset, 
                  y = antibody_titres))+
    geom_point()+ 
    ylab(paste(paste("Antibody titre (", out, sep= ""), ")", sep=""))+ 
    xlab("days since symptom onset")+
    geom_smooth(method='lm', formula= y~x)
  
  if(sumfit$coefficients[2,4] < 0.05){
    
    message("significant association")
  }
  
  # check outliers 
  bx <- boxplot(df_ab_time$antibody_titres)
  
  if(length(bx$out) == 0){
    
    message("no outliers")
    
    return(list(sumfit, q))
    
  }else{
    
    df_ab_time1 <- df_ab_time[which(df_ab_time$antibody_titres < min(bx$out)), ]
    
    # test trend 
    fit1 <- lm(antibody_titres ~ days_since_symptom_onset, data = df_ab_time1)
    sumfit1 <- summary(fit1)
    
    # plot  
    q1 <- ggplot(df_ab_time1,
                aes(x = days_since_symptom_onset, 
                    y = antibody_titres))+
      geom_point()+ 
      ylab(paste(paste("Antibody titre (", out, sep= ""), ")", sep=""))+ 
      xlab("days since symptom onset")+
      geom_smooth(method='lm', formula= y~x)
    
    if(sumfit1$coefficients[2,4] < 0.05){
      
      message("significant association without outliers")
    }
    
    return(list(sumfit, q, sumfit1, q1))
    
  }

}

################################################################################
# Test for difference in response, two categories (t-test/Wilcoxon test)       #
################################################################################
test_antibody_titres <- function(var_qual, var_quant, type){
  
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
    
  }else if(type == "hospitalized"){
    
    idx1 <- which(ds$hospitalized[idx] == "yes")
    idx2 <- which(ds$hospitalized[idx] != "yes")
    
    df <- data.frame(
      
      antibody_titres = c(ds[idx[idx1], ncl1], 
                          ds[idx[idx2], ncl1]), 
      
      type = c(rep("hospitalized", length(ds[idx[idx1], ncl1])), 
               rep("not hospitalized", length(ds[idx[idx2], ncl1]))))
    
  }else if(type == "sex"){
    
    idx1 <- which(ds$Gender[idx] == "F")
    idx2 <- which(ds$Gender[idx] == "M")
    
    df <- data.frame(
      
      antibody_titres = c(ds[idx[idx1], ncl1], 
                          ds[idx[idx2], ncl1]), 
      
      type = c(rep("Female", length(ds[idx[idx1], ncl1])), 
               rep("Male", length(ds[idx[idx2], ncl1]))))
    
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
        message(paste("normality is supported for category",
                      levels(df$type)[i], 
                      sep = " "))      }
    }
  }
  
  ttest <- t.test(antibody_titres ~ type, data = df)
  wtest <- wilcox.test(antibody_titres ~ type, data = df) 
  
  r <- ggplot(df, aes(x = type, 
                      y = antibody_titres, 
                      fill = type))+
    geom_boxplot()+
    ylab(paste(paste("Antibody titre (", out, sep= ""), ")", sep=""))+ 
    xlab("")+
    theme(legend.position = "none")
  
  return(list(ttest, r, wtest))
  
}

################################################################################
# Adding case counts in boxplot                                                #
################################################################################
stat_box_data <- function(y, upper_limit = max(df$antibody_titres) * 1.15) {
  return( 
    data.frame(
      y = 0, 
      label = length(y)
    )
  )
}

################################################################################
# Test for difference in response, multiple categories (Anova/Kruskal-Wallis)  #
################################################################################
Anova_Kruskal_antibody_titres <- function(var_qual, var_quant, type){
  
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
  
  if(str_detect(var_qual, "2020", negate = TRUE) & 
     str_detect(var_quant, "2020", negate = TRUE)){
    
    month <- "May"
    
  }else if(str_detect(var_qual, "2020", negate = FALSE) & 
           str_detect(var_quant, "2020", negate = FALSE)){
    
    month <- "Nov"
    
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
  
  if(class(df$antibody_titres) == "factor")
    df$antibody_titres <- as.numeric(as.character(df$antibody_titres))
  
  t <- ggplot(df, aes(x = type, 
                      y = antibody_titres,
                      fill = type))+
    geom_boxplot()+
    ylab(paste(paste(paste(paste("Antibody titres", month, sep = " "),
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
  for(i in 1:length((levels(df$type)))){
               
    data <- df$antibody_titres[which(df$type == levels(df$type)[i])]
    
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

  if(normality == 1){
    
    if(equal_variance == 1){
      
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
    } else{
      
      # anova with unequal variances 
      res.aov_unequal <- oneway.test(antibody_titres ~ type, var.equal = FALSE,
                                     data = df)
      
      if(res.aov_unequal$p.value < 0.05) {
        message("statistically significant difference (unequal variance)")
      }
      
      res <- res.aov_unequal$p.value
      test <- "one-way Anova unequal within group variance"
      sumtur <- NA
      
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
