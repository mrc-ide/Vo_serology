################################################################################
# Probability to obtain test result combinations (A+D+R+ = Ap_Dp_Rp,           #
# A+D-R+ = Ap_Dn_Rp etc.)                                                      #
################################################################################

p_Ap_Dp_Rp <- function(param, data){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  p <- data[1]*log(sens_A*sens_D*sens_R*theta + 
                     (1-spec_A)*(1-spec_D)*(1-spec_R)*(1-theta))
  
  return(p)
  
  }

p_Ap_Dn_Rp <- function(param, data){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  p <- data[2]*log(sens_A*(1-sens_D)*sens_R*theta + 
                     (1-spec_A)*spec_D*(1-spec_R)*(1-theta))
  
  return(p)
}

p_Ap_Dp_Rn <- function(param, data){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  p <- data[3]*log(sens_A*sens_D*(1-sens_R)*theta + 
                     (1-spec_A)*(1-spec_D)*spec_R*(1-theta))
  
  return(p)
}

p_An_Dp_Rp <- function(param, data){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  p <- data[4]*log((1-sens_A)*sens_D*sens_R*theta + 
                     spec_A*(1-spec_D)*(1-spec_R)*(1-theta))
  
  return(p)
  
}

p_Ap_Dn_Rn <- function(param, data){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  p <- data[5]*log(sens_A*(1-sens_D)*(1-sens_R)*theta + 
                     (1-spec_A)*spec_D*spec_R*(1-theta))
  
  return(p)
  
}

p_An_Dp_Rn <- function(param, data){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  p <- data[6]*log((1-sens_A)*sens_D*(1-sens_R)*theta + 
                     spec_A*(1-spec_D)*spec_R*(1-theta))
  
  return(p)
  
}

p_An_Dn_Rp <- function(param, data){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  p <- data[7]*log((1-sens_A)*(1-sens_D)*sens_R*theta + 
                     spec_A*spec_D*(1-spec_R)*(1-theta))
  
  return(p)
  
}

p_An_Dn_Rn <- function(param, data){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  p <- data[8]*log((1-sens_A)*(1-sens_D)*(1-sens_R)*theta + 
             spec_A*spec_D*spec_R*(1-theta))
  
  return(p)
  
}

################################################################################
# Probability to obtain observed sens and spec from in-house experiments       #
################################################################################

p_se_A_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- se_A_n*log(sens_A)
  
  return(p)
}

p_se_A_d_se_A_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- (se_A_d - se_A_n)*log(1-sens_A)
  
  return(p)
}

p_sp_A_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- sp_A_n*log(spec_A)
  
  return(p)
}

p_sp_A_d_sp_A_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- (sp_A_d - sp_A_n)*log(1-spec_A)
  
  return(p)
}

p_se_D_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- se_D_n*log(sens_D)
  
  return(p)
}

p_se_D_d_se_D_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- (se_D_d - se_D_n)*log(1-sens_D)
  
  return(p)
}

p_sp_D_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- sp_D_n*log(spec_D)
  
  return(p)
}

p_sp_D_d_sp_D_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- (sp_D_d - sp_D_n)*log(1-spec_D)
  
  return(p)
}

p_se_R_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- se_R_n*log(sens_R)
  
  return(p)
}

p_se_R_d_se_R_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- (se_R_d - se_R_n)*log(1-sens_R)
  
  return(p)
}

p_sp_R_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- sp_R_n*log(spec_R)
  
  return(p)
}

p_sp_R_d_sp_R_n <- function(param, obs){
  
  theta <- param[1]
  sens_A <- param[2]
  sens_D <- param[3]
  sens_R <- param[4]
  spec_A <- param[5]
  spec_D <- param[6]
  spec_R <- param[7]
  
  se_A_n <- obs[1]
  se_A_d <- obs[2]
  se_D_n <- obs[3]
  se_D_d <- obs[4]
  se_R_n <- obs[5]
  se_R_d <- obs[6]
  sp_A_n <- obs[7]
  sp_A_d <- obs[8]
  sp_D_n <- obs[9]
  sp_D_d <- obs[10]
  sp_R_n <- obs[11]
  sp_R_d <- obs[12]
  
  p <- (sp_R_d - sp_R_n)*log(1-spec_R)
  
  return(p)
}

################################################################################
# Log-likelihood function (see Supplementary Methods)                          #
################################################################################

loglikelihood <- function(param, data, obs){ 
  
  llike <- vector()
  
  llike[1] <- p_se_A_n(param, obs)
  
  llike[2] <- p_se_A_d_se_A_n(param, obs)
  
  llike[3] <- p_se_D_n(param, obs)
  
  llike[4] <- p_se_D_d_se_D_n(param, obs)
  
  llike[5] <- p_se_R_n(param, obs)
  
  llike[6] <- p_se_R_d_se_R_n(param, obs)
  
  llike[7] <- p_sp_A_n(param, obs)
  
  llike[8] <- p_sp_A_d_sp_A_n(param, obs)
  
  llike[9] <- p_sp_D_n(param, obs)
  
  llike[10] <- p_sp_D_d_sp_D_n(param, obs)
  
  llike[11] <- p_sp_R_n(param, obs)
  
  llike[12] <- p_sp_R_d_sp_R_n(param, obs)
  
  llike[13] <- p_Ap_Dp_Rp(param, data)
  
  llike[14] <- p_Ap_Dn_Rp(param, data)
  
  llike[15] <- p_Ap_Dp_Rn(param, data)
  
  llike[16] <- p_An_Dp_Rp(param, data)

  llike[17] <- p_Ap_Dn_Rn(param, data)
  
  llike[18] <- p_An_Dp_Rn(param, data)
  
  llike[19] <- p_An_Dn_Rp(param, data)
  
  llike[20] <- p_An_Dn_Rn(param, data)
  
  return(sum(llike, na.rm = TRUE)) 

}

################################################################################
# Run MCMC function 
################################################################################

run_MCMC <- function(data, mcmc_par, parameters, obs){
  
  # set up variables 
  nbIter <- mcmc_par[1]
  stepForStorage <- mcmc_par[2]
  
  lowerLimit <- rep(mcmc_par[3], 7)
  upperLimit <- rep(mcmc_par[4], 7)
  randomWalkRate <- rep(mcmc_par[5], 7)
  
  currentParVector <- unlist(parameters)
  
  storedPar <- matrix(0, ncol = (length(currentParVector)), 
                      nrow = nbIter%/%stepForStorage)  
  
  proposalSD <- matrix(0, ncol = length(currentParVector), 
                       nrow = nbIter%/%stepForStorage)
  
  logLik <- matrix(0, ncol = 1, nrow = nbIter%/%stepForStorage) 
  
  nbAccepted <- rep(0, length(currentParVector)) 
  
  currentLogLik <- loglikelihood(param = parameters, 
                                 data = data, 
                                 obs = obs)
  # MCMC
  for(iter in 2:nbIter){ 
    
    for(parID in 1:length(currentParVector)){ 
      
      oldValue <- currentParVector[parID]
      
      repeat{      
        newValue <- oldValue * exp(randomWalkRate[parID] * rnorm(1))   
        if(newValue > lowerLimit[parID] & newValue < upperLimit[parID]) break
      }
      
      currentParVector[parID] <- newValue
      
      parameters <- currentParVector
      parameters <- unlist(parameters)
      
      newLogLik <- loglikelihood(param = parameters, 
                                 data = data, 
                                 obs = obs)

      # maximise log-likelihood 
      Nu <- newLogLik - currentLogLik

      if(log(runif(1)) < Nu){  # acceptance rule
        
        currentLogLik <- newLogLik                              

        currentParVector[parID] <- newValue
        currentParVector <- unlist(currentParVector)
        nbAccepted[parID] <- nbAccepted[parID] + 1   
        
      }else{                                      
        
        currentParVector[parID] <- oldValue
        currentParVector <- unlist(currentParVector)
        nbAccepted[parID] <- nbAccepted[parID]      

      }
      
      # Adaptive MCMC - Tuning sample proposal variance
      Nu0 <- 0.234 #ideal acceptance probability
      
      if (iter < nbIter/(stepForStorage*2)){
        
        #tuning randomWalkRate
        temp <- randomWalkRate[parID]*
          exp(0.4*(exp(min(Nu,0)) - Nu0)/(35*iter/nbIter + 1))  
        
        #bounded RW between 0 - 10  
        randomWalkRate[parID] <- ifelse(temp > 1e-10, 
                                        ifelse(temp < 10, temp, 10), 1e-10)    
      }
      
      if(iter%%stepForStorage == 0){
        
        storedPar[iter%/%stepForStorage, ] <- c(currentParVector)
        logLik[iter%/%stepForStorage] <- c(currentLogLik)  
        proposalSD[iter%/%stepForStorage, ] <- unlist(randomWalkRate)
        
      }
    }
    
    AcceptRate <- nbAccepted/nbIter 
    
    results <- list(storedLnL=logLik, 
                    Acceptance=AcceptRate, 
                    ParameterValues=storedPar)
  }
  
  return(results)
  
}

################################################################################
# Plot MCMC chains                                                             #
################################################################################

plot_posterior_estimates <- function(burnin, 
                                     parametersToEstimate, 
                                     parameterNames, 
                                     res){
  
  v <- which(parametersToEstimate == 1)
  nr <- ceiling((length(v)+2)/2)
  
  par(mfrow = c(nr, 2), mar = c(2, 2.2, 2, 2))
  
  colours <- c("blue", "red", "forestgreen")
  legend_pos <- c("topright", "bottomright")
  
  plot(res$storedLnL[-seq(1,burnin,1), 1], type = "l", xlab="", ylab="", 
       main = "loglikelihood", cex.axis = 1.5, cex.lab = 1.5, 
       col = colours[1], 
       ylim=c(min(res$storedLnL[-seq(1,burnin,1), 1]), 
              max(res$storedLnL[-seq(1,burnin,1), 1])))
  
  for(i in 1:dim(res$ParameterValues)[2]){
    
    mean <- mean(res$ParameterValues[-seq(1,burnin,1), i])
    median <- median(res$ParameterValues[-seq(1,burnin,1), i])

    plot(res$ParameterValues[-seq(1,burnin,1), i], 
           type = "l", xlab="",ylab="", col = colours[1], 
           main = parameterNames[v[i]],
           cex.axis = 1.5, cex.lab=1.5) 
    
    abline(h = mean, 
           col=colours[1],
           lty = 1)
    
    abline(h = median, 
           col=colours[2], 
           lty = 2)
  }
  
}

################################################################################
################################################################################


