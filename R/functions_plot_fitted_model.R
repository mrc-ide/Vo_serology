################################################################################
# function to set integer x axis                                               #
# https://www.r-bloggers.com/2019/11/setting-axes-to-integer-values-in-ggplot2/#
################################################################################

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

################################################################################
# Function to display single chains from MCMC                                  #
################################################################################

plot_chains <- function(burnin, parametersToEstimate, parameterNames, results){
  
  v <- which(parametersToEstimate == 1)
  nr <- ceiling((length(v)+2)/2)
  
  par(mfrow = c(nr, 2))
  
  plot(results$storedLnL[-seq(1,burnin,1), 1], type = "l", xlab="", ylab="", 
       main = "loglikelihood", cex.axis = 1.5, cex.lab = 1.5, 
       ylim=c(min(results$storedLnL[-seq(1,burnin,1), 1]), 
              max(results$storedLnL[-seq(1,burnin,1), 1])))
  
  plot(results$storedDev[-seq(1,burnin,1), 1], type = "l", xlab="", ylab="", 
       main = "deviance", cex.axis = 1.5, cex.lab = 1.5, 
       ylim=c(min(results$storedDev[-seq(1,burnin,1), 1]), 
              max(results$storedDev[-seq(1,burnin,1), 1])))
  
  for(i in 1:length(v)){
     
      mean <- mean(results$ParameterValues[-seq(1,burnin,1), i])
      median <- median(results$ParameterValues[-seq(1,burnin,1), i])

      plot(results$ParameterValues[-seq(1,burnin,1), i], 
           type = "l", xlab="",ylab="", 
           main = parameterNames[v[i]],
           cex.axis = 1.5, cex.lab=1.5)
    
      abline(h = mean, 
             col="red")
      
      abline(h = median, 
             col="blue")
      
      legend("topright", 
             legend = c(paste("mean =", round(mean, 2)), 
                        paste("median = ", round(median, 2))), 
             text.col=c("red", "blue"), 
             bty = "o")
      
  }

}

################################################################################
# Function to display multiple chains from MCMC                                # 
################################################################################

plot_multiple_chains <- function(burnin, parametersToEstimate, parameterNames, 
                                 results){
  
  v <- which(parametersToEstimate == 1)
  nr <- ceiling((length(v)+2)/2)
  
  par(mfrow = c(nr, 2), mar = c(2, 2.2, 2, 2))
  
  colours <- c("blue", "red", "forestgreen")
  legend_pos <- c("topright", "bottomright")

  for(k in 1:length(results)){
    
    res1 <- results[[k]]
    
    if(k == 1) {
      
      plot(res1$storedLnL[-seq(1,burnin,1), 1], type = "l", xlab="", ylab="", 
         main = "loglikelihood", cex.axis = 1.5, cex.lab = 1.5, 
         col = colours[k], 
         ylim=c(min(res1$storedLnL[-seq(1,burnin,1), 1]), 
                max(res1$storedLnL[-seq(1,burnin,1), 1])))
      } else{
        
        points(res1$storedLnL[-seq(1,burnin,1), 1], type = "l", 
               col = colours[k])
      }
  }
  for(k in 1:length(results)){
    
    res1 <- results[[k]]
    
    if(k == 1) {
      
      plot(res1$storedDev[-seq(1,burnin,1), 1], type = "l", xlab="", ylab="", 
         main = "deviance", cex.axis = 1.5, cex.lab = 1.5, col = colours[k], 
         ylim=c(min(res1$storedDev[-seq(1,burnin,1), 1]), 
                max(res1$storedDev[-seq(1,burnin,1), 1])))
    } else{
      
      points(res1$storedDev[-seq(1,burnin,1), 1], type = "l", 
             col = colours[k])
    }
  }
    
  for(i in 1:length(v)){
      
    mu <- vector()
    med <- vector()
    
    for(k in 1:length(results)){
      
      res1 <- results[[k]]
      
      mean <- mean(res1$ParameterValues[-seq(1,burnin,1), i])
      median <- median(res1$ParameterValues[-seq(1,burnin,1), i])
      
      if(k == 1) {
      
      plot(res1$ParameterValues[-seq(1,burnin,1), i], 
           type = "l", xlab="",ylab="", col = colours[k], 
           main = parameterNames[v[i]],
           cex.axis = 1.5, cex.lab=1.5)
        
      }else{
          
        points(res1$ParameterValues[-seq(1,burnin,1), i], type = "l", 
               col = colours[k])
        
      }
      
      abline(h = mean, 
             col=colours[k],
             lty = 1)
      
      abline(h = median, 
             col=colours[k], 
             lty = 2)
        
      mu <- c(mu, mean)
      med <- c(med, median)
      
    }
    
    legend(legend_pos[1], 
           legend = c(paste("mean", round(mu,2), sep = " ")),  
           text.col = colours,
           col = colours, 
           lty = c(1),
           bty = "o")
    
    legend(legend_pos[2], 
           legend = c(paste("median", round(med,2), sep = " ")),  
           text.col = colours,
           col = colours, 
           lty = c(2),
           bty = "o")
  }
}

################################################################################
# code to place the legend in an empty facet                                   #
# https://stackoverflow.com/questions/54438495/                                # 
# shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2                  #                                                                            #
################################################################################

shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from 
              ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% 
                                 class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. 
            Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), 
                             min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), 
                             max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}

################################################################################
# Reproducing ggplot standard colours for factors                              #
################################################################################

gg_colour_hue <- function(n){
  
  hues= seq(15, 375, length = n + 1)
  
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# Setting lighter colours                                                      #
################################################################################

vir_lite = function(cols, ds=0.4, dv=0.7) {
  cols = rgb2hsv(col2rgb(cols))
  cols["v", ] = cols["v", ] + dv*(1 - cols["v", ])
  cols["s", ] = ds*cols["s", ]
  apply(cols, 2, function(x) hsv(x[1], x[2], x[3]))
}

################################################################################
################################################################################

