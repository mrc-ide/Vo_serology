################################################################################
# Function plotting the DIC and the estimated parameters                       #
################################################################################
plot_deviance_dic_params <- function(res, models, type, mindic){
  
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
  
  df <- data.frame(model = models[,1],
                   deviance_mean = dev_score,
                   deviance_lower = dev_l,
                   deviance_upper = dev_u,
                   dic = dic)
  
  # select models with best deviance
  bm <- which(df$deviance_mean == min(df$deviance_mean))
  idx_bm <- which(df$deviance_mean <= df$deviance_upper[bm])
  
  if(type == "original model"){
    
    xmax <- 16.5
    
  }else if(type %in% c("extended model", "2-groups model")){
    
    xmax <- 8.5
    
  }else{
    
    message("type not valid")
    return()
  }
  
  r_dev <- ggplot(df)+
    ylim(c(10,22))+
    geom_point(aes(x = model, y = deviance_mean))+
    geom_errorbar(aes(x = model, ymin=deviance_lower, ymax=deviance_upper), 
                  width=0) + 
    geom_rect(aes(xmin = 0.5, xmax = xmax, 
                  ymin=df$deviance_lower[bm], ymax=df$deviance_upper[bm]),
              fill = "#FDE725FF", alpha = 0.01)+
    geom_point(aes(x = model, y = deviance_mean))+
    geom_point(aes(x = model[bm], y = deviance_mean[bm]),
               colour = "red")+
    geom_errorbar(aes(x = model, ymin=deviance_lower, ymax=deviance_upper), 
                  width=0)+
    geom_errorbar(aes(x = model[bm], ymin=deviance_lower[bm], ymax=deviance_upper[bm]), 
                  colour = "red", 
                  width=0) +
    ylab("deviance")+xlab("")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(type)
  
  # ------------------------------------------------------------------------- #
  
  if(type == "original model") 
    mindic <- df$dic[which(df$dic == min(df$dic))]
  
  r_dic <- ggplot(df)+ylim(c(14,22))+ 
    geom_hline(aes(yintercept = mindic),  colour = "#FDE725FF")+
    geom_point(aes(x = model, y = dic))+
    ylab("DIC")+
    xlab("")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # ------------------------------------------------------------------------- #
  
  models_alphabetical <- models$model[order(models$model)]
  
  if(type == "original model"){
    
    estimated <- matrix(rep("gray", 6*dim(models)[1]), ncol = dim(models)[1]) 
    ncl <- c(2:5,7:8)
    add <- 0
    idx <- c(1,6,9)

  }else if(type == "extended model"){
    
    estimated <- matrix(rep("gray", 7*dim(models)[1]), ncol = dim(models)[1]) 
    ncl <- c(2:5,7:8,10)
    add <- 1
    idx <- c(1,6,9)

  }else if(type == "2-groups model"){
    
    estimated <- matrix(rep("gray", 7*dim(models)[1]), ncol = dim(models)[1])
    ncl <- c(2:4,6:7,9:10)
    add <- 1
    idx <- c(1,5,8)
  }
  
  for(i in 1:dim(models)[1]){
    
    estimated[which(models[which(models$model == 
                                   models_alphabetical[i]), ncl] == 1), i] <- 
      "green"
  }
  
  fake <- seq(0, 8, 1/8)[-1]
  fake <- fake
  
  if(type == "original model"){
    
    fake <- fake[1:16]
    
  }else{
    
    fake <- fake[1:8]
    
  }
  
  xmin <- 0.6
  xmax <- 1.4
  
  p <- ggplot(df)+
    geom_point(aes(x = model, y = fake))+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    geom_rect(xmin = xmin, xmax = xmax,  ymin = 0, ymax = 1, color = "black",  
              fill = estimated[1,1]) +
    geom_rect(xmin = xmin, xmax = xmax,  ymin = 1, ymax = 2, color = "black",
              fill = estimated[2,1]) +
    geom_rect(xmin = xmin, xmax = xmax,  ymin = 2, ymax = 3, color = "black",
              fill = estimated[3,1]) +
    geom_rect(xmin = xmin, xmax = xmax,  ymin = 3, ymax = 4, color = "black",
              fill = estimated[4,1]) +
    geom_rect(xmin = xmin, xmax = xmax,  ymin = 4, ymax = 5, color = "black", 
              fill = estimated[5,1]) +
    geom_rect(xmin = xmin, xmax = xmax,  ymin = 5, ymax = 6, color = "black",
              fill = estimated[6,1]) +
    scale_y_continuous("parameters", 
                       breaks = seq(0.5, (5.5+add), 1), 
                       labels = colnames(models)[-idx], 
                       limits = c(0.2, (5.8+add)))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank())
  
  if(type %in% c("extended model", "2-groups model")){
    
    p <- p +
      geom_rect(xmin = xmin, xmax = xmax,  ymin = 6, ymax = 7, color = "black",
                fill = estimated[7,1]) 
    
  } 
  
  for(i in 2:dim(models)[1]){
    
    xmin <- xmin + 1
    xmax <- xmax + 1
    
    q <- p +
      geom_rect(xmin = xmin, xmax = xmax,  ymin = 0, ymax = 1, color = "black",  
                fill = estimated[1,i]) +
      geom_rect(xmin = xmin, xmax = xmax,  ymin = 1, ymax = 2, color = "black",
                fill = estimated[2,i]) +
      geom_rect(xmin = xmin, xmax = xmax,  ymin = 2, ymax = 3, color = "black",
                fill = estimated[3,i]) +
      geom_rect(xmin = xmin, xmax = xmax,  ymin = 3, ymax = 4, color = "black",
                fill = estimated[4,i]) +
      geom_rect(xmin = xmin, xmax = xmax,  ymin = 4, ymax = 5, color = "black", 
                fill = estimated[5,i]) +
      geom_rect(xmin = xmin, xmax = xmax,  ymin = 5, ymax = 6, color = "black",
                fill = estimated[6,i]) 
    
    if(type %in% c("extended model", "2-groups model")){
      
      q <- q +
        geom_rect(xmin = xmin, xmax = xmax,  ymin = 6, ymax = 7, color = "black",
                  fill = estimated[7,i]) 
      
    }
    
    p <- q + xlab("")
  }
  
  return(list(r_dev, r_dic, p, mindic, selected_chains))
  
}

################################################################################
################################################################################


