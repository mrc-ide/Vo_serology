################################################################################
# Moment generating function of infectiousness for generalised model,          #
# assuming gamma distributed hazard. See Fraser et al, AJE, 2011.              #
################################################################################
MGF_gamma <- function(x, n, beta, alpha, shape){
  
  if(is.infinite(shape)){ 
    
    return(exp((-beta/(n^alpha))*(x)))
    
  }else{
      
      return((shape/(shape+x*(beta/(n^alpha))))^shape)
      
  }
}

################################################################################
# Shift matrix elements upwards - implementing equation [10]                   #
################################################################################
shift_matrix <- function(mat, data, t, p) {
  
  nrow <- dim(mat)[1]
  ncol <- dim(mat)[2]
  
  mat_final <- matrix(0, nrow = nrow, ncol = ncol)
  
  if(t < nrow){
    
    mat_final[1:(nrow-t), 1:ncol] <- mat[(1+t):nrow, 1:ncol]
    
  }
  
  size <- as.numeric(rownames(data)) + t
  bin <- dbinom(t, size, p)
  
  mat_final <- t(t(mat_final) %*% diag(bin))
    
  return(mat_final)
  
}

################################################################################
# Define linear equation system needed to estimate the probability of observing#
# m cases in household of n people with s0 susceptibles. See equation [4] in   #
# SI of Fraser et al, AJE, 2011.                                               # 
################################################################################
compute_F <- function(n, s0, Q, parameters){
  
  beta <- parameters[1]
  alpha <- parameters[3]
  shape <- parameters[4]
  
  mat <- matrix(rep(0, (s0+1)*(s0+1)), ncol = (s0+1))
  y <- rep(0, (s0+1))

  for(i in 1:(s0+1)){

    k <- s0 - i + 1

    y[i] <- choose(n = s0, k = k)

    for(j in 1:(s0+1)){

      if(i <= j){

        m <- s0 - j + 1

        den <- ((MGF_gamma((s0-k), n, beta, alpha, shape))^m)*Q^(s0-k)
        mat[i, j] <- choose(n = (s0 - m), k = (k - m))/den

      }

    }

  }

  x <- backsolve (mat, y)

  return(rev(x))
  
}

################################################################################
# Solves equation [8] in SI of Fraser et al, AJE, 2011.                        #                       
################################################################################
compute_R <- function(data, parameters, r){
  
  Q <- parameters[2]
  Q_prior <- parameters[5]
  p_pr <- parameters[6]
  
  mat <- matrix(0, 
                ncol = dim(data)[2], 
                nrow = dim(data)[1])
  i <- 1
  for(n in as.numeric(colnames(data))){
    
    s0 <- (n - r)
    
    if(s0 >= 0){
      
      if(Q_prior != 1){
          
          matA <- matrix(0, ncol = (s0+1), nrow = (s0+1))
          
          Fq_prior <- compute_F(n, s0, Q_prior, parameters)
          
          for(l in 0:s0){
            
            v <- compute_F(n, (s0-l), Q, parameters)
            
            matA[1:length(v),(l+1)] <- v
            
          }
          
          x <- matA %*% Fq_prior
          
      }else{
          
          x <- compute_F(n, s0, Q, parameters)
          
      }
      if(length(x) > dim(data)[1]){
          
        v <- x[1:dim(data)[1]]
          
      }else{
          
          v <- x
          
      }
        
      mat[1:length(v), i] <- v*dbinom(r, n, p_pr)
      
    }
    
    i <- i+1
  }
  
  return(mat)
  
}

################################################################################
# Solves equation [9] in SI of Fraser et al, AJE, 2011.                        #                       
################################################################################
compute_S <- function(data, parameters){
  
  r <- seq(from = 0, to = dim(data)[1], by = 1) 

  S <- Reduce('+', lapply(r, compute_R, 
                         data = data, 
                         parameters = parameters))
  
  return(S)

}

################################################################################
# Solves equation [10] in SI of Fraser et al, AJE, 2011.                       #                       
################################################################################
compute_T <- function(data, parameters){
  
  mat <- compute_S(data, parameters)
  t <- seq(from = 0, to = dim(data)[1], by = 1)
  
  T <- Reduce('+', lapply(t, shift_matrix, 
                          mat = mat, 
                          data = data, 
                          p = parameters[7]))
  return(T)
}

################################################################################
# Calculate deviance, see equation [14] in SI of Fraser et al, AJE, 2011.      #
################################################################################
deviance <- function(data, parameters){
  
  mod <- compute_T(data, parameters)
  
  # implementing equation [11] (U_m)
  one_minus_pcom_mat <- matrix(0, nrow = dim(mod)[1], ncol = dim(mod)[2])
  one_minus_pcom_mat[1,] <- 1-parameters[8]
  model <- one_minus_pcom_mat + mod*parameters[8]
  
  # model must be a positive matrix (taking the log)
  model[which(model <= 0)] <- NA

  dev <- 2*sum(data * (log(t(t(data) / colSums(data)), base = exp(1)) - 
                         log(model, base = exp(1))), na.rm = TRUE)
  
  return(dev)
  
}

################################################################################
# Calculate loglikelihood, see equation [12] in SI of Fraser et al, AJE, 2011. #
################################################################################
loglikelihood <- function(data, parameters){
  
  mod <- compute_T(data, parameters)
  
  # implementing equation [11] (U_m)
  one_minus_pcom_mat <- matrix(0, nrow = dim(mod)[1], ncol = dim(mod)[2])
  one_minus_pcom_mat[1,] <- 1-parameters[8]
  model <- one_minus_pcom_mat + mod*parameters[8]
  
  # model must be a positive matrix (taking the log)
  model[which(model <= 0)] <- NA
  
  llike <- sum(data * (log(model, base = exp(1))), na.rm = TRUE)

  return(llike)
  
}

################################################################################
# Return model predictions                                                     #
################################################################################
model_predictions <- function(data, parameters){
  
  mod <- compute_T(data, parameters)
  
  # implementing equation [11] (U_m)
  one_minus_pcom_mat <- matrix(0, nrow = dim(mod)[1], ncol = dim(mod)[2])
  one_minus_pcom_mat[1,] <- 1-parameters[8]
  model <- one_minus_pcom_mat + mod*parameters[8]
  
  pred <- sweep(model, 2, colSums(data), FUN = '*')
  
  return(pred)
  
}

################################################################################
# Run MCMC                                                                     #
################################################################################
run_MCMC <- function(data, mcmc_par, parameters, parametersToEstimate){
  
  # set up variables 
  nbIter <- mcmc_par[1,1]
  stepForStorage <- mcmc_par[1,2]
  
  currentParVector <- unlist(parameters[which(parametersToEstimate > 0)])
  
  if(any(currentParVector == 0) | any(is.infinite(currentParVector))){
    
    currentParVector[which(currentParVector == 0)] <- 
      runif(length(which(currentParVector == 0)), 
            min = 0,
            max = 1)
    
    currentParVector[which(is.infinite(currentParVector))] <- 
      runif(length(which(is.infinite(currentParVector))), 
            min = 0,
            max = 1)
  }
  
  parameters[which(parametersToEstimate > 0)] <- currentParVector
  parameters <- unlist(parameters)
  
  currentLogLik <- loglikelihood(data, parameters)
  currentDev <- deviance(data, parameters)
  
  idx <- seq(3, 10, 1)
  lowerLimit <- mcmc_par[1, idx[which(parametersToEstimate == 1)]]
  upperLimit <- mcmc_par[2, idx[which(parametersToEstimate == 1)]]
  randomWalkRate <- mcmc_par[3, idx[which(parametersToEstimate == 1)]]
  
  storedPar <- matrix(0, ncol = (length(currentParVector)), 
                      nrow = nbIter%/%stepForStorage)  
  
  proposalSD <- matrix(0, ncol = length(currentParVector), 
                       nrow = nbIter%/%stepForStorage)
  
  logLik <- matrix(0, ncol = 1, nrow = nbIter%/%stepForStorage) 
  dev <- matrix(0, ncol = 1, nrow = nbIter%/%stepForStorage)             
  
  nbAccepted <- rep(0, length(currentParVector)) 
  
  # MCMC
  for(iter in 2:nbIter){ 
    
    for(parID in 1:length(currentParVector)){ 
      
      oldValue <- currentParVector[parID]
      
      repeat{      
        newValue <- oldValue * exp(randomWalkRate[parID] * rnorm(1))   
        if(newValue > lowerLimit[parID] & newValue < upperLimit[parID]) break
      }
      
      currentParVector[parID] <- newValue

      parameters[which(parametersToEstimate == 1)] <- currentParVector
      parameters <- unlist(parameters)
      
      newLogLik <- loglikelihood(data, parameters)
      newDev <- deviance(data, parameters)
      
      # minimise the deviance (same as maximise the loglikelihood)
      Nu <- - 0.5*(newDev - currentDev)   
      
      if(log(runif(1)) < Nu){  # acceptance rule
        
        currentLogLik <- newLogLik                              
        currentDev <- newDev                  
        
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
        
        #storing parameters and iterations
        if(parID == 1) message(paste("MCMC iteration ", iter, sep = " "))   
        
        storedPar[iter%/%stepForStorage, ] <- c(currentParVector)
        logLik[iter%/%stepForStorage] <- c(currentLogLik)  
        dev[iter%/%stepForStorage] <- c(currentDev)  
        proposalSD[iter%/%stepForStorage, ] <- unlist(randomWalkRate)
        
      }
    }
    
    AcceptRate <- nbAccepted/nbIter 
    
    results <- list(storedLnL=logLik, 
                    storedDev=dev, 
                    Acceptance=AcceptRate, 
                    ParameterValues=storedPar)
  }
  
  return(results)
  
}

################################################################################
# compute posterior estimates                                                  #
################################################################################
compute_post_estimates <- function(burnin, 
                                   default_parameters, 
                                   parametersToEstimate, 
                                   parameterNames, 
                                   results){
  
  v <- which(parametersToEstimate == 1)
  
  mean <- vector()
  median <- vector()
  lower_CrI <- vector()
  upper_CrI <- vector()

  for(i in 1:length(v)){
    
    mean[i] <- mean(results$ParameterValues[-seq(1,burnin,1), i])
    median[i] <- median(results$ParameterValues[-seq(1,burnin,1), i])
    
    lower_CrI[i] <- quantile(results$ParameterValues[-seq(1,burnin,1), i], 
                                 c(0.025))
      
    upper_CrI[i] <- quantile(results$ParameterValues[-seq(1,burnin,1), i], 
                                 c(0.975))
  }
  
  dev_mean <- mean(results$storedDev[-seq(1,burnin,1)])
  dev_median <- median(results$storedDev[-seq(1,burnin,1)])
  dev_lower <- quantile(results$storedDev[-seq(1,burnin,1)], c(0.025))
  dev_upper <- quantile(results$storedDev[-seq(1,burnin,1)], c(0.975)) 
  
  deviance <- c(dev_mean, dev_median, dev_lower, dev_upper)
  
  llike_mean <- mean(results$storedLnL[-seq(1,burnin,1)])
  llike_median <- median(results$storedLnL[-seq(1,burnin,1)])
  llike_lower <- quantile(results$storedLnL[-seq(1,burnin,1)], c(0.025))
  llike_upper <- quantile(results$storedLnL[-seq(1,burnin,1)], c(0.975)) 
  
  llike <- c(llike_mean, llike_median, llike_lower, llike_upper)
  
  df <- data.frame(mean = mean, 
             median = median, 
             lower_CrI = lower_CrI, 
             upper_CrI = upper_CrI, 
             acceptance = results$Acceptance)
  
  rownames(df) <- parameterNames[v]
  
  DIC <- compute_DIC(burnin, 
                     default_parameters, 
                     parametersToEstimate, 
                     parameterNames, 
                     results)

  return(list(df, deviance, DIC, NA, llike))

}

################################################################################
# save posterior samples                                                       # 
################################################################################
save_posterior <- function(burnin, 
                           parametersToEstimate, 
                           parameterNames, 
                           results, 
                           select_model_version,
                           idc, 
                           path){
  
  v <- which(parametersToEstimate == 1)
  
  df <- cbind(results$storedLnL[-seq(1,burnin,1), 1],
              results$storedDev[-seq(1,burnin,1), 1])
  
  for(i in 1:length(v)){
    
    df <- cbind(df, 
                results$ParameterValues[-seq(1,burnin,1), i])
  }
    
  colnames(df) <- c("logLike", "deviance", parameterNames[v])
  
  
  write.csv(df, 
            file.path(dir_output, 
                      paste(paste(paste(select_model_version, "chain", sep="_"), 
                            idc, sep=""), ".csv", sep="")), 
            row.names = FALSE)
  
}

################################################################################
# compute DIC (Deviance Information Criterion)                                 #                                                          
################################################################################
compute_DIC <- function(burnin, 
                        default_parameters, 
                        parametersToEstimate, 
                        parameterNames, 
                        results){
  
  
  v <- which(parametersToEstimate == 1)
  
  mean <- vector()
  
  for(i in 1:length(v)){
    
    mean[i] <- mean(results$ParameterValues[-seq(1,burnin,1), i])
    
  }
  
  # mean of parameters 
  parameters <- unlist(default_parameters)
  parameters[v] <- mean 

  DIC <- 2*mean(results$storedDev[-seq(1,burnin,1)]) - 
    deviance(data, parameters)
  
  return(DIC)

}

################################################################################
# Fit models to data                                                           #
################################################################################
wrapper_model <- function(model_variant, 
                          data, 
                          par, 
                          mcmc_par,
                          parameters,
                          parameterNames, 
                          burnin, 
                          dir_output,
                          filename, 
                          id_chain){
  
  
  parametersToEstimate <- unlist(par[model_variant, 2:9])
  
  results <- list()
  df_res <- list()
  
  # loop over chains
  for (idc in seq(1,id_chain,1)) {
    
    if(idc > 1){
      
      idx <- which(parametersToEstimate == 1)
      parameters[idx] <- runif(length(idx), 0, 1)
    }
    
    results[[idc]] <- run_MCMC(data, 
                      mcmc_par, 
                      parameters, 
                      parametersToEstimate)
  }
  
  return(list(results))
  
}
  
################################################################################
################################################################################

