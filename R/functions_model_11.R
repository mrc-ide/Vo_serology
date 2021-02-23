model_gen <- odin::odin({
  
  # define dependent parameters
  gamma <- 1 / (generation_time - inv_nu - inv_delta)
  nu    <- 1 / inv_nu
  delta <- 1 / inv_delta
  
  ww  <- if (do_lockdown > 0.5  && t > tQ - 0.5) w           else 1
  q_Q <- if (do_q > 0.01) do_q else 0
  p_1 <- if (do_conttrac > 0.01 && t > tQ - 0.5 && v == 1L) do_conttrac else 0 # p_traced
  p_2 <- if (do_conttrac > 0.01 && t > tQ - 0.5 && v == 2L) do_conttrac else 0 # p_traced
  p_3 <- if (do_conttrac > 0.01 && t > tQ - 0.5 && v == 3L) do_conttrac else 0 # p_traced
  
  
  # define time-dependent parameters
  beta <- ww * R0_1 / (inv_delta + 1/gamma)
  lambda <- beta * (q_TPp * TPp + q_A * I_A + I_S + q_Q * (TP_Q1 + I_Q)) * S/N
  
  
  # ODE system
  deriv(S)    <- - lambda
  deriv(E)    <- (1 - p_3) * lambda - (nu + p_1 + p_2) * E
  deriv(TPp)  <- nu * E - (delta + p_2) * TPp
  deriv(I_A)  <- p * delta * TPp - (gamma + p_2) * I_A
  deriv(I_S)  <- (1 - p) * delta * TPp - (gamma + p_2) * I_S
  deriv(TP_A) <- gamma * I_A - (sigma + p_2) * TP_A
  deriv(TP_S) <- gamma * I_S - (sigma + p_2) * TP_S
  deriv(TN)   <- sigma * (TP_A + TP_S)
  
  
  # quarantined classes
  deriv(Q)     <- p_3 * lambda + (p_1 + p_2) * E - nu * Q
  deriv(TP_Q1) <- p_2 * TPp + nu * Q - delta * TP_Q1
  deriv(I_Q)   <- p_2 * (I_A + I_S)  + delta * TP_Q1 - gamma * I_Q
  deriv(TP_Q2) <- p_2 * (TP_A + TP_S) + gamma * I_Q - sigma * TP_Q2
  deriv(TN_Q)  <- sigma * TP_Q2
  
  
  # initial conditions
  initial(S)     <- N - seed
  initial(E)     <- 0
  initial(TPp)   <- 0
  initial(I_A)   <- p * seed
  initial(I_S)   <- (1 - p) * seed
  initial(TP_A)  <- 0
  initial(TP_S)  <- 0
  initial(TN)    <- 0
  initial(Q)     <- 0
  initial(TP_Q1) <- 0
  initial(I_Q)   <- 0
  initial(TP_Q2) <- 0
  initial(TN_Q)  <- 0
  
  # input parameters
  seed      <- user()
  p         <- user()
  w         <- user()
  inv_nu    <- user()
  inv_delta <- user()
  R0_1      <- user()
  sigma     <- user()
  tSeed     <- user()
  time1     <- user()
  time2     <- user()
  tQ        <- user()
  N         <- user()
  q_TPp     <- user()
  q_A       <- user()
  do_q      <- user()
  generation_time <- user()
  do_lockdown <- user()
  do_conttrac <- user()
  v           <- user()
  
})
