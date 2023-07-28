# Evaluating the effects of nest management on a recovering raptor #
# using integrated population modeling #

# Code provided for peer review and is subject to change

# version: v17.1.3
# code author: Caroline D. Cappello
# last updated: June 2023

model <- nimbleCode({
  
  #### * Mark-resight-recovery * ###############################################
  
  # ------------------------------------------------.-
  # Parameters:
  # phi1: survival probability state 1 (Y1)
  # phi2: survival probability state 2 (Y2) 
  # phi3: survival probability state 3 (Y3)
  # phiAf: survival probability state 5 (AY3f adult floater)
  # phiAb: survival probability state 6 (AY3b adult breeder)
  # p2: resight probability Y2 (resight probability all ages carrying a GPS = 1)
  # p3: resight probability Y3
  # pAf: resight probability AY3 floater
  # pAb: resight probability AY3 breeder
  # psiFB: transition probability from floater to breeder (i.e. recruitment)
  # psiBF: transition probability from breeder to floater (i.e. end tenure)
  # psiGPS: transition probability, Y3 to YAf with GPS (only once instance of a previously tagged bird gaining a tag)
  # psiFAIL: transition probability from telemetered to non-telemetered (i.e. tag failure)
  # rBAND: probability of dead recovery of banded individual
  # rGPS: probability of dead recovery of telemetered individual
  
  # ------------------------------------------------.-
  # States (z):
  # 1 alive Y1
  # 2 alive Y2
  # 3 alive Y3
  # 4 alive Y4+ float
  # 5 alive Y4+ breed
  # 6 recently dead
  # 7 dead
  # 8 alive Y1 with GPS
  # 9 alive Y2 with GPS
  # 10 alive Y3 with GPS
  # 11 alive YAf with GPS
  # 12 alive YAb with GPS
  # 13 recently dead with GPS
  # 14 dead with GPS
  
  # Observations (y):  
  # 1 not seen
  # 2 seen Y1
  # 3 seen Y2
  # 4 seen Y3
  # 5 seen Y4+ floater
  # 6 seen Y4+ breeder 
  # 7 dead
  # 8 seen Y1 with GPS
  # 9 seen Y2 with GPS
  # 10 seen Y3 with GPS
  # 11 seen YAf with GPS
  # 12 seen YAb with GPS
  # 13 dead with GPS
  
  
  # ------------------------------------------------.-
  
  # priors
  beta[1] ~ dnorm(mean = 0, sd = 1.5) 
  beta[2] ~ dnorm(mean = 0, sd = 1.5) 
  beta[3] ~ dnorm(mean = 0, sd = 1.5) 
  beta[4] ~ dnorm(mean = 0, sd = 1.5) 
  beta[5] ~ dnorm(mean = 0, sd = 1.5)
  beta[6] ~ dnorm(mean = 0, sd = 1.5) 
  
  alpha_mr[1] ~ dnorm(mean = 0, sd = 1.5) 
  alpha_mr[2] ~ dnorm(mean = 0, sd = 1.5) 
  alpha_mr[3] ~ dnorm(mean = 0, sd = 1.5) 
  alpha_mr[4] ~ dnorm(mean = 0, sd = 1.5) 
  
  psiFB ~ dunif(0, 1)
  psiBF ~ dunif(0, 1)
  psiGPS ~ dunif(0, 1)
  psiFAIL ~ dunif(0, 1)
  
  rBAND ~ dunif(0, 1)
  rGPS ~ dunif(0, 1)
  
  sdeps_mr[1] ~ dunif(0,10)
  sdeps_mr[2] ~ dunif(0,10)
  sdeps_mr[3] ~ dunif(0,10)
  
  
  for(t in 1:(nyears_mr-1)+nyears_proj){
    phi1_m_t[t] <- 1/(1+exp(-(beta[1] + eps[t]))) 
    phi2_m_t[t] <- 1/(1+exp(-(beta[2] + eps[t]))) 
    phi3_m_t[t] <- 1/(1+exp(-(beta[3] + eps[t]))) 
    phiAf_m_t[t] <- 1/(1+exp(-(beta[4] + eps[t]))) 
    phiAb_m_t[t] <- 1/(1+exp(-(beta[5] + eps[t]))) 
  }
  
  for(t in 1:(nyears_mr-1)+nyears_proj){
    phi1_f_t[t] <- 1/(1+exp(-(beta[1] + beta[6] + eps[t]))) 
    phi2_f_t[t] <- 1/(1+exp(-(beta[2] + beta[6] + eps[t]))) 
    phi3_f_t[t] <- 1/(1+exp(-(beta[3] + beta[6] + eps[t]))) 
    phiAf_f_t[t] <- 1/(1+exp(-(beta[4] + beta[6] + eps[t]))) 
    phiAb_f_t[t] <- 1/(1+exp(-(beta[5] + beta[6] + eps[t]))) 
  }
  
  mean_phi1_m <- 1/(1+exp(-(beta[1]))) 
  mean_phi2_m <- 1/(1+exp(-(beta[2]))) 
  mean_phi3_m <- 1/(1+exp(-(beta[3]))) 
  mean_phiAf_m <- 1/(1+exp(-(beta[4]))) 
  mean_phiAb_m <- 1/(1+exp(-(beta[5]))) 
  
  mean_phi1_f <- 1/(1+exp(-(beta[1] + beta[6]))) 
  mean_phi2_f <- 1/(1+exp(-(beta[2] + beta[6]))) 
  mean_phi3_f <- 1/(1+exp(-(beta[3] + beta[6]))) 
  mean_phiAf_f <- 1/(1+exp(-(beta[4] + beta[6]))) 
  mean_phiAb_f <- 1/(1+exp(-(beta[5] + beta[6]))) 
  
  for(t in 1:(nyears_mr-1)+nyears_proj){
    p2_t[t] <- 1/(1+exp(-(alpha_mr[1] + zeta_f[t]))) 
    p3_t[t] <- 1/(1+exp(-(alpha_mr[2] + zeta_f[t]))) 
    pAf_t[t] <- 1/(1+exp(-(alpha_mr[3] + zeta_f[t]))) 
    pAb_t[t] <- 1/(1+exp(-(alpha_mr[4] + zeta_b[t]))) 
  }
  
  mean_p2 <- 1/(1+exp(-(alpha_mr[1]))) 
  mean_p3 <- 1/(1+exp(-(alpha_mr[2]))) 
  mean_pAf <- 1/(1+exp(-(alpha_mr[3]))) 
  mean_pAb <- 1/(1+exp(-(alpha_mr[4]))) 
  
  #Random effect for phi timestep
  for (t in 1:(nyears_mr-1)+nyears_proj) {
    eps[t] ~ dnorm(0, sd = sdeps_mr[1]) # residual variation not explained by the covariates
  }

  #Random effect for timestep p NOT BREEDING
  for (t in 1:(nyears_mr-1)+nyears_proj) {
    zeta_f[t] ~ dnorm(0, sd = sdeps_mr[2]) 
  }

  #Random effect for timestep p BREEDING
  for (t in 1:(nyears_mr-1)+nyears_proj) {
    zeta_b[t] ~ dnorm(0, sd = sdeps_mr[3]) 
  }

  
  # probabilities of state z(t+1) given z(t)
  for (i in 1:nind){ # loop over individuals
    for (t in 1:(nyears_mr-1)) {
      logit(phi1[i,t])  <- beta[1]  + beta[6]*sex[i] + eps[t]
      logit(phi2[i,t])  <- beta[2]  + beta[6]*sex[i] + eps[t]
      logit(phi3[i,t])  <- beta[3]  + beta[6]*sex[i] + eps[t]
      logit(phiAf[i,t]) <- beta[4] +  beta[6]*sex[i] + eps[t]
      logit(phiAb[i,t]) <- beta[5] +  beta[6]*sex[i] + eps[t]
      gamma[1,1,i,t] <- 0                                         # Pr(Y1 at time t+1, given Y1 at time t)
      gamma[1,2,i,t] <- phi1[i,t]                                 # Pr(Y1 -> Y2)
      gamma[1,3,i,t] <- 0                                         # Pr(Y1 -> Y3)
      gamma[1,4,i,t] <- 0                                         # Pr(Y1 -> Adult floater)
      gamma[1,5,i,t] <- 0                                         # Pr(Y1 -> Adult breeder)
      gamma[1,6,i,t] <- (1-phi1[i,t]) * rBAND                     # Pr(Y1 -> Recently dead)
      gamma[1,7,i,t] <- (1-phi1[i,t]) * (1-rBAND)                 # Pr(Y1 -> Long dead)
      gamma[1,8,i,t] <- 0                                         # Pr(Y1 -> Y1 GPS)
      gamma[1,9,i,t] <- 0                                         # Pr(Y1 -> Y2 GPS)
      gamma[1,10,i,t] <- 0                                        # Pr(Y1 -> Y3 GPS)
      gamma[1,11,i,t] <- 0                                        # Pr(Y1 -> Adult floater GPS)
      gamma[1,12,i,t] <- 0                                        # Pr(Y1 -> Adult breeder GPS)
      gamma[1,13,i,t] <- 0                                        # Pr(Y1 -> Recently dead GPS)
      gamma[1,14,i,t] <- 0                                        # Pr(Y1 -> Long dead GPS)
      
      gamma[2,1,i,t] <- 0                                         # Pr(Y2 -> Y1)
      gamma[2,2,i,t] <- 0                                         # Pr(Y2 -> Y2)
      gamma[2,3,i,t] <- phi2[i,t]                                 # Pr(Y2 -> Y3)
      gamma[2,4,i,t] <- 0                                         # Pr(Y2 -> Adult floater)
      gamma[2,5,i,t] <- 0                                         # Pr(Y2 -> Adult breeder)
      gamma[2,6,i,t] <- (1-phi2[i,t]) * rBAND                     # Pr(Y2 -> Recently dead)
      gamma[2,7,i,t] <- (1-phi2[i,t]) * (1-rBAND)                 # Pr(Y2 -> Long dead)
      gamma[2,8,i,t] <- 0                                         # Pr(Y2 -> Y1 GPS)
      gamma[2,9,i,t] <- 0                                         # Pr(Y2 -> Y2 GPS)
      gamma[2,10,i,t] <- 0                                        # Pr(Y2 -> Y3 GPS)
      gamma[2,11,i,t] <- 0                                        # Pr(Y2 -> Adult floater GPS)
      gamma[2,12,i,t] <- 0                                        # Pr(Y2 -> Adult breeder GPS)
      gamma[2,13,i,t] <- 0                                        # Pr(Y2 -> Recently dead GPS)
      gamma[2,14,i,t] <- 0                                        # Pr(Y2 -> Long dead GPS)
      
      gamma[3,1,i,t] <- 0                                         # Pr(Y3 -> Y1)
      gamma[3,2,i,t] <- 0                                         # Pr(Y3 -> Y2)
      gamma[3,3,i,t] <- 0                                         # Pr(Y3 -> Y3)
      gamma[3,4,i,t] <- phi3[i,t] * (1-psiFB) * (1-psiGPS)        # Pr(Y3 -> Adult floater)
      gamma[3,5,i,t] <- phi3[i,t] * psiFB                         # Pr(Y3 -> Adult breeder)
      gamma[3,6,i,t] <- (1-phi3[i,t]) * rBAND                     # Pr(Y3 -> Recently dead)
      gamma[3,7,i,t] <- (1-phi3[i,t]) * (1-rBAND)                 # Pr(Y3 -> Dead)
      gamma[3,8,i,t] <- 0                                         # Pr(Y3 -> Y1 GPS)
      gamma[3,9,i,t] <- 0                                         # Pr(Y3 -> Y2 GPS)
      gamma[3,10,i,t] <- 0                                        # Pr(Y3 -> Y3 GPS)
      gamma[3,11,i,t] <- phi3[i,t] * (1-psiFB) * psiGPS           # Pr(Y3 -> Adult floater GPS)
      gamma[3,12,i,t] <- 0                                        # Pr(Y3 -> Adult breeder GPS)
      gamma[3,13,i,t] <- 0                                        # Pr(Y3 -> Recently dead GPS)
      gamma[3,14,i,t] <- 0                                        # Pr(Y3 -> Long dead GPS)
      
      gamma[4,1,i,t] <- 0                                         # Pr(Adult floater -> Y1)
      gamma[4,2,i,t] <- 0                                         # Pr(Adult floater -> Y2)
      gamma[4,3,i,t] <- 0                                         # Pr(Adult floater -> Y3)
      gamma[4,4,i,t] <- phiAf[i,t] * (1-psiFB)                    # Pr(Adult floater -> Adult floater)
      gamma[4,5,i,t] <- phiAf[i,t] * psiFB                        # Pr(Adult floater -> Adult breeder)
      gamma[4,6,i,t] <- (1-phiAf[i,t]) * rBAND                    # Pr(Adult floater -> Recently dead)
      gamma[4,7,i,t] <- (1-phiAf[i,t]) * (1-rBAND)                # Pr(Adult floater -> Dead)
      gamma[4,8,i,t] <- 0                                         # Pr(Adult floater -> Y1 GPS)
      gamma[4,9,i,t] <- 0                                         # Pr(Adult floater -> Y2 GPS)
      gamma[4,10,i,t] <- 0                                        # Pr(Adult floater -> Y3 GPS)
      gamma[4,11,i,t] <- 0                                        # Pr(Adult floater -> Adult floater GPS)
      gamma[4,12,i,t] <- 0                                        # Pr(Adult floater -> Adult breeder GPS)
      gamma[4,13,i,t] <- 0                                        # Pr(Adult floater -> Recently dead GPS)
      gamma[4,14,i,t] <- 0                                        # Pr(Adult floater -> Long dead GPS)
      
      gamma[5,1,i,t] <- 0                                         # Pr(Adult breeder -> Y1)
      gamma[5,2,i,t] <- 0                                         # Pr(Adult breeder -> Y2)
      gamma[5,3,i,t] <- 0                                         # Pr(Adult breeder -> Y3)
      gamma[5,4,i,t] <- phiAb[i,t] * psiBF                        # Pr(Adult breeder -> Adult floater)
      gamma[5,5,i,t] <- phiAb[i,t] * (1-psiBF)                    # Pr(Adult breeder -> Adult breeder)
      gamma[5,6,i,t] <- (1-phiAb[i,t]) * rBAND                    # Pr(Adult breeder -> Recently dead)
      gamma[5,7,i,t] <- (1-phiAb[i,t]) * (1-rBAND)                # Pr(Adult breeder -> Dead)
      gamma[5,8,i,t] <- 0                                         # Pr(Adult breeder -> Y1 GPS
      gamma[5,9,i,t] <- 0                                         # Pr(Adult breeder -> Y2 GPS)
      gamma[5,10,i,t] <- 0                                        # Pr(Adult breeder -> Y3 GPS)
      gamma[5,11,i,t] <- 0                                        # Pr(Adult breeder -> Adult floater GPS)
      gamma[5,12,i,t] <- 0                                        # Pr(Adult breeder -> Adult breeder GPS)
      gamma[5,13,i,t] <- 0                                        # Pr(Adult breeder -> Recently dead GPS)
      gamma[5,14,i,t] <- 0                                        # Pr(Adult breeder -> Long dead GPS)
      
      gamma[6,1,i,t] <- 0                                         # Pr(Recently dead -> Y1)
      gamma[6,2,i,t] <- 0                                         # Pr(Recently dead -> Y2)
      gamma[6,3,i,t] <- 0                                         # Pr(Recently dead -> Y3)
      gamma[6,4,i,t] <- 0                                         # Pr(Recently dead -> Adult floater)
      gamma[6,5,i,t] <- 0                                         # Pr(Recently dead -> Adult breeder)
      gamma[6,6,i,t] <- 0                                         # Pr(Recently dead -> Recently dead)
      gamma[6,7,i,t] <- 1                                         # Pr(Recently dead -> Dead not recovered)
      gamma[6,8,i,t] <- 0                                         # Pr(Recently dead -> Y1 GPS)
      gamma[6,9,i,t] <- 0                                         # Pr(Recently dead -> Y2 GPS)
      gamma[6,10,i,t] <- 0                                        # Pr(Recently dead -> Y3 GPS)
      gamma[6,11,i,t] <- 0                                        # Pr(Recently dead -> Adult floater GPS)
      gamma[6,12,i,t] <- 0                                        # Pr(Recently dead -> Adult breeder GPS)
      gamma[6,13,i,t] <- 0                                        # Pr(Recently dead -> Recently dead GPS)
      gamma[6,14,i,t] <- 0                                        # Pr(Recently dead -> Long dead GPS)
      
      gamma[7,1,i,t] <- 0                                         # Pr(Long dead -> Y1)
      gamma[7,2,i,t] <- 0                                         # Pr(Long dead -> Y2)
      gamma[7,3,i,t] <- 0                                         # Pr(Long dead -> Y3)
      gamma[7,4,i,t] <- 0                                         # Pr(Long dead -> Adult floater)
      gamma[7,5,i,t] <- 0                                         # Pr(Long dead -> Adult breeder)
      gamma[7,6,i,t] <- 0                                         # Pr(Long dead -> Recently dead)
      gamma[7,7,i,t] <- 1                                         # Pr(Long dead -> Long dead)
      gamma[7,8,i,t] <- 0                                         # Pr(Long dead -> Y1 GPS)
      gamma[7,9,i,t] <- 0                                         # Pr(Long dead -> Y2 GPS)
      gamma[7,10,i,t] <- 0                                        # Pr(Long dead -> Y3 GPS)
      gamma[7,11,i,t] <- 0                                        # Pr(Long dead -> Adult floater GPS)
      gamma[7,12,i,t] <- 0                                        # Pr(Long dead -> Adult breeder GPS)
      gamma[7,13,i,t] <- 0                                        # Pr(Long dead -> Recently dead GPS)
      gamma[7,14,i,t] <- 0                                        # Pr(Long dead -> Long dead GPS)
      
      gamma[8,1,i,t] <- 0                                         # Pr(Y1 GPS alive at time t+1 -> Y1 at time t)
      gamma[8,2,i,t] <- phi1[i,t] * psiFAIL                       # Pr(Y1 GPS -> Y2)
      gamma[8,3,i,t] <- 0                                         # Pr(Y1 GPS -> Y3)
      gamma[8,4,i,t] <- 0                                         # Pr(Y1 GPS -> Adult floater)
      gamma[8,5,i,t] <- 0                                         # Pr(Y1 GPS -> Adult breeder)
      gamma[8,6,i,t] <- (1-phi1[i,t]) * rBAND * psiFAIL           # Pr(Y1 GPS -> Recently dead)
      gamma[8,7,i,t] <- (1-phi1[i,t]) * (1-rBAND) * psiFAIL       # Pr(Y1 GPS -> Dead)
      gamma[8,8,i,t] <- 0                                         # Pr(Y1 GPS -> Y1 GPS)
      gamma[8,9,i,t] <- phi1[i,t] * (1-psiFAIL)                   # Pr(Y1 GPS -> Y2 GPS)
      gamma[8,10,i,t] <- 0                                        # Pr(Y1 GPS -> Y3 GPS)
      gamma[8,11,i,t] <- 0                                        # Pr(Y1 GPS -> Adult floater GPS)
      gamma[8,12,i,t] <- 0                                        # Pr(Y1 GPS -> Adult breeder GPS)
      gamma[8,13,i,t] <- (1-phi1[i,t]) * rGPS * (1-psiFAIL)       # Pr(Y1 GPS -> Recently dead GPS)
      gamma[8,14,i,t] <- (1-phi1[i,t]) * (1-rGPS) * (1-psiFAIL)   # Pr(Y1 GPS -> Long dead GPS)
      
      gamma[9,1,i,t] <- 0                                         # Pr(Y2 GPS -> Y1)
      gamma[9,2,i,t] <- 0                                         # Pr(Y2 GPS -> Y2)
      gamma[9,3,i,t] <- phi2[i,t] * psiFAIL                       # Pr(Y2 GPS -> Y3)
      gamma[9,4,i,t] <- 0                                         # Pr(Y2 GPS -> Adult floater)
      gamma[9,5,i,t] <- 0                                         # Pr(Y2 GPS -> Adult breeder)
      gamma[9,6,i,t] <- (1-phi2[i,t]) * rBAND * psiFAIL           # Pr(Y2 GPS -> Recently dead)
      gamma[9,7,i,t] <- (1-phi2[i,t]) * (1-rBAND) * psiFAIL       # Pr(Y2 GPS -> Dead)
      gamma[9,8,i,t] <- 0                                         # Pr(Y2 GPS -> Y1 GPS)
      gamma[9,9,i,t] <- 0                                         # Pr(Y2 GPS -> Y2 GPS)
      gamma[9,10,i,t] <- phi2[i,t] * (1-psiFAIL)                  # Pr(Y2 GPS -> Y3 GPS)
      gamma[9,11,i,t] <- 0                                        # Pr(Y2 GPS -> Adult floater GPS)
      gamma[9,12,i,t] <- 0                                        # Pr(Y2 GPS -> Adult breeder GPS)
      gamma[9,13,i,t] <- (1-phi2[i,t]) * rGPS * (1-psiFAIL)       # Pr(Y2 GPS -> Recently dead GPS)
      gamma[9,14,i,t] <- (1-phi2[i,t]) * (1-rGPS) * (1-psiFAIL)   # Pr(Y2 GPS -> Long dead GPS)
      
      gamma[10,1,i,t] <- 0                                        # Pr(Y3 GPS -> Y1)
      gamma[10,2,i,t] <- 0                                        # Pr(Y3 GPS -> Y2)
      gamma[10,3,i,t] <- 0                                        # Pr(Y3 GPS -> Y3)
      gamma[10,4,i,t] <- phi3[i,t] * (1-psiFB) * psiFAIL          # Pr(Y3 GPS -> Adult floater)
      gamma[10,5,i,t] <- phi3[i,t] * psiFB * psiFAIL              # Pr(Y3 GPS -> Adult breeder)
      gamma[10,6,i,t] <- (1-phi3[i,t]) * rBAND * psiFAIL          # Pr(Y3 GPS -> Recently dead)
      gamma[10,7,i,t] <- (1-phi3[i,t]) * (1-rBAND) * psiFAIL      # Pr(Y3 GPS -> Dead)
      gamma[10,8,i,t] <-  0                                       # Pr(Y3 GPS -> Y1 GPS)
      gamma[10,9,i,t] <-  0                                       # Pr(Y3 GPS -> Y2 GPS)
      gamma[10,10,i,t] <- 0                                       # Pr(Y3 GPS -> Y3 GPS)
      gamma[10,11,i,t] <- phi3[i,t] * (1-psiFB) * (1-psiFAIL)     # Pr(Y3 GPS -> Adult floater GPS)
      gamma[10,12,i,t] <- phi3[i,t] * psiFB * (1-psiFAIL)         # Pr(Y3 GPS -> Adult breeder GPS)
      gamma[10,13,i,t] <- (1-phi3[i,t]) * rGPS * (1-psiFAIL)      # Pr(Y3 GPS -> Recently dead GPS)
      gamma[10,14,i,t] <- (1-phi3[i,t]) * (1-rGPS) * (1-psiFAIL)  # Pr(Y3 GPS -> Long dead GPS)
      
      gamma[11,1,i,t] <- 0                                        # Pr(Adult floater GPS -> Y1)
      gamma[11,2,i,t] <- 0                                        # Pr(Adult floater GPS -> Y2)
      gamma[11,3,i,t] <- 0                                        # Pr(Adult floater GPS -> Y3)
      gamma[11,4,i,t] <- phiAf[i,t] * (1-psiFB) * psiFAIL         # Pr(Adult floater GPS -> Adult floater)
      gamma[11,5,i,t] <- phiAf[i,t] * psiFB * psiFAIL             # Pr(Adult floater GPS -> Adult breeder)
      gamma[11,6,i,t] <- (1-phiAf[i,t]) * rBAND * psiFAIL         # Pr(Adult floater GPS -> Recently dead)
      gamma[11,7,i,t] <- (1-phiAf[i,t]) * (1-rBAND) * psiFAIL     # Pr(Adult floater GPS -> Dead)
      gamma[11,8,i,t] <-  0                                       # Pr(Adult floater GPS -> Y1 GPS)
      gamma[11,9,i,t] <-  0                                       # Pr(Adult floater GPS -> Y2 GPS)
      gamma[11,10,i,t] <- 0                                       # Pr(Adult floater GPS -> Y3 GPS)
      gamma[11,11,i,t] <- phiAf[i,t] * (1-psiFB) * (1-psiFAIL)    # Pr(Adult floater GPS -> Adult floater GPS)
      gamma[11,12,i,t] <- phiAf[i,t] * psiFB * (1-psiFAIL)        # Pr(Adult floater GPS -> Adult breeder GPS)
      gamma[11,13,i,t] <- (1-phiAf[i,t]) * rGPS * (1-psiFAIL)     # Pr(Adult floater GPS -> Recently dead GPS)
      gamma[11,14,i,t] <- (1-phiAf[i,t]) * (1-rGPS) * (1-psiFAIL) # Pr(Adult floater GPS -> Long dead GPS)
      
      gamma[12,1,i,t] <- 0                                        # Pr(Adult breeder GPS -> Y1)
      gamma[12,2,i,t] <- 0                                        # Pr(Adult breeder GPS -> Y2)
      gamma[12,3,i,t] <- 0                                        # Pr(Adult breeder GPS -> Y3)
      gamma[12,4,i,t] <- phiAb[i,t] * psiBF * psiFAIL             # Pr(Adult breeder GPS -> Adult floater)
      gamma[12,5,i,t] <- phiAb[i,t] * (1-psiBF) * psiFAIL         # Pr(Adult breeder GPS -> Adult breeder)
      gamma[12,6,i,t] <- (1-phiAb[i,t]) * rBAND * psiFAIL         # Pr(Adult breeder GPS -> Recently dead)
      gamma[12,7,i,t] <- (1-phiAb[i,t]) * (1-rBAND) * psiFAIL     # Pr(Adult breeder GPS -> Dead)
      gamma[12,8,i,t] <-  0                                       # Pr(Adult floater GPS -> Y1 GPS)
      gamma[12,9,i,t] <-  0                                       # Pr(Adult floater GPS -> Y2 GPS)
      gamma[12,10,i,t] <- 0                                       # Pr(Adult floater GPS -> Y3 GPS)
      gamma[12,11,i,t] <- phiAb[i,t] * (1-psiFB) * (1-psiFAIL)    # Pr(Adult floater GPS -> Adult floater GPS)
      gamma[12,12,i,t] <- phiAb[i,t] * psiFB * (1-psiFAIL)        # Pr(Adult floater GPS -> Adult breeder GPS)
      gamma[12,13,i,t] <- (1-phiAb[i,t]) * rGPS * (1-psiFAIL)     # Pr(Adult floater GPS -> Recently dead GPS)
      gamma[12,14,i,t] <- (1-phiAb[i,t]) * (1-rGPS) * (1-psiFAIL) # Pr(Adult floater GPS -> Long dead GPS)
      
      gamma[13,1,i,t] <- 0                                        # Pr(Recently dead GPS -> Y1)
      gamma[13,2,i,t] <- 0                                        # Pr(Recently dead GPS -> Y2)
      gamma[13,3,i,t] <- 0                                        # Pr(Recently dead GPS -> Y3)
      gamma[13,4,i,t] <- 0                                        # Pr(Recently dead GPS -> Adult floater)
      gamma[13,5,i,t] <- 0                                        # Pr(Recently dead GPS -> Adult breeder)
      gamma[13,6,i,t] <- 0                                        # Pr(Recently dead GPS -> Recently dead)
      gamma[13,7,i,t] <- 1                                        # Pr(Recently dead GPS -> Long dead)
      gamma[13,8,i,t] <-  0                                       # Pr(Recently dead GPS -> Y1 GPS)
      gamma[13,9,i,t] <-  0                                       # Pr(Recently dead GPS -> Y2 GPS)
      gamma[13,10,i,t] <- 0                                       # Pr(Recently dead GPS -> Y3 GPS)
      gamma[13,11,i,t] <- 0                                       # Pr(Recently dead GPS -> Adult floater GPS)
      gamma[13,12,i,t] <- 0                                       # Pr(Recently dead GPS -> Adult breeder GPS)
      gamma[13,13,i,t] <- 0                                       # Pr(Recently dead GPS -> Recently dead GPS)
      gamma[13,14,i,t] <- 0                                       # Pr(Recently dead GPS -> Long dead GPS)
      
      gamma[14,1,i,t] <- 0                                        # Pr(Long Dead GPS -> Y1)
      gamma[14,2,i,t] <- 0                                        # Pr(Long Dead GPS -> Y2)
      gamma[14,3,i,t] <- 0                                        # Pr(Long Dead GPS -> Y3)
      gamma[14,4,i,t] <- 0                                        # Pr(Long Dead GPS -> Adult floater)
      gamma[14,5,i,t] <- 0                                        # Pr(Long Dead GPS -> Adult breeder)
      gamma[14,6,i,t] <- 0                                        # Pr(Long Dead GPS -> Recently dead)
      gamma[14,7,i,t] <- 1                                        # Pr(Long Dead GPS -> Long dead)
      gamma[14,8,i,t] <-  0                                       # Pr(Long dead GPS -> Y1 GPS)
      gamma[14,9,i,t] <-  0                                       # Pr(Long dead GPS -> Y2 GPS)
      gamma[14,10,i,t] <- 0                                       # Pr(Long dead GPS -> Y3 GPS)
      gamma[14,11,i,t] <- 0                                       # Pr(Long dead GPS -> Adult floater GPS)
      gamma[14,12,i,t] <- 0                                       # Pr(Long dead GPS -> Adult breeder GPS)
      gamma[14,13,i,t] <- 0                                       # Pr(Long dead GPS -> Recently dead GPS)
      gamma[14,14,i,t] <- 0                                       # Pr(Long dead/not recovered  GPS -> Long dead GPS)
    } #t
  } #i
  
  # observation matrix: probabilities of y(t) given z(t)
  for (t in 1:(nyears_mr-1)){
    logit(p2[t]) <- alpha_mr[1] + zeta_f[t]
    logit(p3[t]) <- alpha_mr[2] + zeta_f[t]    
    logit(pAf[t]) <- alpha_mr[3] + zeta_f[t]
    logit(pAb[t]) <- alpha_mr[4] + zeta_b[t]
    omega[1,1,t] <- 0          # Pr(alive Y1 t --> non-detected t)
    omega[1,2,t] <- 1          # Pr(alive Y1 t --> detected Y1 t) # banding occasion for fledglings
    omega[1,3,t] <- 0          # Pr(alive Y1 t --> detected Y2 t)
    omega[1,4,t] <- 0          # Pr(alive Y1 t --> detected Y3 t)
    omega[1,5,t] <- 0          # Pr(alive Y1 t --> detected Adult floater t)
    omega[1,6,t] <- 0          # Pr(alive Y1 t --> detected Adult breeder t)
    omega[1,7,t] <- 0          # Pr(alive Y1 t --> detected dead t)
    omega[1,8,t] <- 0          # Pr(alive Y1 t --> Y1 GPS t)
    omega[1,9,t] <- 0          # Pr(alive Y1 t --> Y2 GPS t)
    omega[1,10,t] <- 0         # Pr(alive Y1 t --> Y3 GPS t)
    omega[1,11,t] <- 0         # Pr(alive Y1 t --> Adult floater GPS t)
    omega[1,12,t] <- 0         # Pr(alive Y1 t --> Adult breeder GPS t)
    omega[1,13,t] <- 0         # Pr(alive Y1 t --> detected dead GPS t)
    
    omega[2,1,t] <- 1 - p2[t]  # Pr(alive Y2 t --> non-detected t)
    omega[2,2,t] <- 0          # Pr(alive Y2 t --> detected Y1 t)
    omega[2,3,t] <- p2[t]      # Pr(alive Y2 t --> detected Y2 t)
    omega[2,4,t] <- 0          # Pr(alive Y2 t --> detected Y3 t)
    omega[2,5,t] <- 0          # Pr(alive Y2 t --> detected Adult floater t)
    omega[2,6,t] <- 0          # Pr(alive Y2 t --> detected Adult breeder t)
    omega[2,7,t] <- 0          # Pr(alive Y2 t --> detected dead t)
    omega[2,8,t] <- 0          # Pr(alive Y2 t --> detected Y1 GPS t)
    omega[2,9,t] <- 0          # Pr(alive Y2 t --> detected Y2 GPS t)
    omega[2,10,t] <- 0         # Pr(alive Y2 t --> detected Y3 GPS t)
    omega[2,11,t] <- 0         # Pr(alive Y2 t --> detected Adult floater GPS t)
    omega[2,12,t] <- 0         # Pr(alive Y2 t --> detected Adult breeder GPS t)
    omega[2,13,t] <- 0         # Pr(alive Y2 t --> detected detected dead GPS t)
    
    omega[3,1,t] <- 1 - p3[t]  # Pr(alive Y3  t -> non-detected t)
    omega[3,2,t] <- 0          # Pr(alive Y3  t -> detected Y1 t)
    omega[3,3,t] <- 0          # Pr(alive Y3  t -> detected Y2 t)
    omega[3,4,t] <- p3[t]      # Pr(alive Y3  t -> detected Y3 t)
    omega[3,5,t] <- 0          # Pr(alive Y3  t -> detected Adult floater t)
    omega[3,6,t] <- 0          # Pr(alive Y3  t -> detected Adult breeder t)
    omega[3,7,t] <- 0          # Pr(alive Y3  t -> detected dead t)
    omega[3,8,t] <- 0          # Pr(alive Y3 t -> detected Y1 GPS t)
    omega[3,9,t] <- 0          # Pr(alive Y3 t -> detected Y2 GPS t)
    omega[3,10,t] <- 0         # Pr(alive Y3 t -> detected Y3 GPS t)
    omega[3,11,t] <- 0         # Pr(alive Y3 t -> detected Adult floater GPS t)
    omega[3,12,t] <- 0         # Pr(alive Y3 t -> detected Adult breeder GPS t)
    omega[3,13,t] <- 0         # Pr(alive Y3 t -> detected detected dead GPS t)
    
    omega[4,1,t] <- 1 - pAf[t] # Pr(alive A floater t -> non-detected t)
    omega[4,2,t] <- 0          # Pr(alive A floater t -> detected Y1 t)
    omega[4,3,t] <- 0          # Pr(alive A floater t -> detected Y2 t)
    omega[4,4,t] <- 0          # Pr(alive A floater t -> detected Y3 t)
    omega[4,5,t] <- pAf[t]     # Pr(alive A floater t -> detected Adult floater t)
    omega[4,6,t] <- 0          # Pr(alive A floater t -> detected Adult breeder t)
    omega[4,7,t] <- 0          # Pr(alive A floater t -> detected dead t)
    omega[4,8,t] <- 0          # Pr(alive A floater t -> detected Y1 GPS t)
    omega[4,9,t] <- 0          # Pr(alive A floater t -> detected Y2 GPS t)
    omega[4,10,t] <- 0         # Pr(alive A floater t -> detected Y3 GPS t)
    omega[4,11,t] <- 0         # Pr(alive A floater t -> detected Adult floater GPS t)
    omega[4,12,t] <- 0         # Pr(alive A floater t -> detected Adult breeder GPS t)
    omega[4,13,t] <- 0         # Pr(alive A floater t -> detected detected dead GPS t)
    
    omega[5,1,t] <- 1 - pAb[t] # Pr(alive A breeder t -> non-detected t)
    omega[5,2,t] <- 0          # Pr(alive A breeder t -> detected Y1 t)
    omega[5,3,t] <- 0          # Pr(alive A breeder t -> detected Y2 t)
    omega[5,4,t] <- 0          # Pr(alive A breeder t -> detected Y3 t)
    omega[5,5,t] <- 0          # Pr(alive A breeder t -> detected Adult floater t)
    omega[5,6,t] <- pAb[t]     # Pr(alive A breeder t -> detected Adult breeder t)
    omega[5,7,t] <- 0          # Pr(alive A breeder t -> detected dead t)
    omega[5,8,t] <- 0          # Pr(alive A breeder t -> detected Y1 GPS t)
    omega[5,9,t] <- 0          # Pr(alive A breeder t -> detected Y2 GPS t)
    omega[5,10,t] <- 0         # Pr(alive A breeder t -> detected Y3 GPS t)
    omega[5,11,t] <- 0         # Pr(alive A breeder t -> detected Adult floater GPS t)
    omega[5,12,t] <- 0         # Pr(alive A breeder t -> detected Adult breeder GPS t)
    omega[5,13,t] <- 0         # Pr(alive A breeder t -> detected detected dead GPS t)
    
    omega[6,1,t] <- 0          # Pr(recently dead t -> non-detected t)
    omega[6,2,t] <- 0          # Pr(recently dead t -> detected Y1  t)
    omega[6,3,t] <- 0          # Pr(recently dead t -> detected Y2  t)
    omega[6,4,t] <- 0          # Pr(recently dead t -> detected Y3 floater  t)
    omega[6,5,t] <- 0          # Pr(recently dead t -> detected Adult floater  t)
    omega[6,6,t] <- 0          # Pr(recently dead t -> detected breeder t)
    omega[6,7,t] <- 1          # Pr(recently dead t -> detected dead t)
    omega[6,8,t] <- 0          # Pr( recently dead t -> detected Y1 GPS t)
    omega[6,9,t] <- 0          # Pr( recently dead t -> detected Y2 GPS t)
    omega[6,10,t] <- 0         # Pr( recently dead t -> detected Y3 GPS t)
    omega[6,11,t] <- 0         # Pr( recently dead t -> detected Adult floater GPS t)
    omega[6,12,t] <- 0         # Pr( recently dead t -> detected Adult breeder GPS t)
    omega[6,13,t] <- 0         # Pr( recently dead t -> detected detected dead GPS t)
    
    omega[7,1,t] <- 1          # Pr(dead not recovered t -> non-detected t)
    omega[7,2,t] <- 0          # Pr(dead not recovered t -> detected Y1  t)
    omega[7,3,t] <- 0          # Pr(dead not recovered t -> detected Y2  t)
    omega[7,4,t] <- 0          # Pr(dead not recovered t -> detected Y3 floater  t)
    omega[7,5,t] <- 0          # Pr(dead not recovered t -> detected Adult floater  t)
    omega[7,6,t] <- 0          # Pr(dead not recovered t -> detected Adult breeder t)
    omega[7,7,t] <- 0          # Pr(dead not recovered t -> detected dead t)
    omega[7,8,t] <- 0          # Pr(dead not recovered t -> detected Y1 GPS t)
    omega[7,9,t] <- 0          # Pr(dead not recovered t -> detected Y2 GPS t)
    omega[7,10,t] <- 0         # Pr(dead not recovered t -> detected Y3 GPS t)
    omega[7,11,t] <- 0         # Pr(dead not recovered t -> detected Adult floater GPS t)
    omega[7,12,t] <- 0         # Pr(dead not recovered t -> detected Adult breeder GPS t)
    omega[7,13,t] <- 0         # Pr(dead not recovered t -> detected detected dead GPS t)
    
    omega[8,1,t] <- 0          # Pr(alive Juv GPS t -> non-detected t)
    omega[8,2,t] <- 0          # Pr(alive Juv GPS t -> detected Y1 t) # banding occasion for fledglings
    omega[8,3,t] <- 0          # Pr(alive Juv GPS t -> detected Y2 t)
    omega[8,4,t] <- 0          # Pr(alive Juv GPS t -> detected Y3 t)
    omega[8,5,t] <- 0          # Pr(alive Juv GPS t -> detected Adult floater t)
    omega[8,6,t] <- 0          # Pr(alive Juv GPS t -> detected Adult breeder t)
    omega[8,7,t] <- 0          # Pr(alive Juv GPS t -> detected dead t)
    omega[8,8,t] <- 1          # Pr(alive Juv GPS t -> Y1 GPS t)
    omega[8,9,t] <- 0          # Pr(alive Juv GPS t -> Y2 GPS t)
    omega[8,10,t] <- 0         # Pr(alive Juv GPS t -> Y3 GPS t)
    omega[8,11,t] <- 0         # Pr(alive Juv GPS t -> Adult floater GPS t)
    omega[8,12,t] <- 0         # Pr(alive Juv GPS t -> Adult breeder GPS t)
    omega[8,13,t] <- 0         # Pr(alive Juv GPS t -> detected dead GPS t)
    
    omega[9,1,t] <- 0          # Pr(alive Y2 GPS t -> non-detected t)
    omega[9,2,t] <- 0          # Pr(alive Y2 GPS t -> detected Y1 t) # banding occasion for fledglings
    omega[9,3,t] <- 0          # Pr(alive Y2 GPS t -> detected Y2 t)
    omega[9,4,t] <- 0          # Pr(alive Y2 GPS t -> detected Y3 t)
    omega[9,5,t] <- 0          # Pr(alive Y2 GPS t -> detected Adult floater t)
    omega[9,6,t] <- 0          # Pr(alive Y2 GPS t -> detected Adult breeder t)
    omega[9,7,t] <- 0          # Pr(alive Y2 GPS t -> detected dead t)
    omega[9,8,t] <- 0          # Pr(alive Y2 GPS t -> Y1 GPS t)
    omega[9,9,t] <- 1          # Pr(alive Y2 GPS t -> Y2 GPS t)
    omega[9,10,t] <- 0         # Pr(alive Y2 GPS t -> Y3 GPS t)
    omega[9,11,t] <- 0         # Pr(alive Y2 GPS t -> Adult floater GPS t)
    omega[9,12,t] <- 0         # Pr(alive Y2 GPS t -> Adult breeder GPS t)
    omega[9,13,t] <- 0         # Pr(alive Y2 GPS t -> detected dead GPS t)
    
    omega[10,1,t] <- 0         # Pr(alive Y3 GPS t -> non-detected t)
    omega[10,2,t] <- 0         # Pr(alive Y3 GPS t -> detected Y1 t) # banding occasion for fledglings
    omega[10,3,t] <- 0         # Pr(alive Y3 GPS t -> detected Y2 t)
    omega[10,4,t] <- 0         # Pr(alive Y3 GPS t -> detected Y3 t)
    omega[10,5,t] <- 0         # Pr(alive Y3 GPS t -> detected Adult floater t)
    omega[10,6,t] <- 0         # Pr(alive Y3 GPS t -> detected Adult breeder t)
    omega[10,7,t] <- 0         # Pr(alive Y3 GPS t -> detected dead t)
    omega[10,8,t] <- 0         # Pr(alive Y3 GPS t -> Y1 GPS t)
    omega[10,9,t] <- 0         # Pr(alive Y3 GPS t -> Y2 GPS t)
    omega[10,10,t] <- 1        # Pr(alive Y3 GPS t -> Y3 GPS t)
    omega[10,11,t] <- 0        # Pr(alive Y3 GPS t -> Adult floater GPS t)
    omega[10,12,t] <- 0        # Pr(alive Y3 GPS t -> Adult breeder GPS t)
    omega[10,13,t] <- 0        # Pr(alive Y3 GPS t -> detected dead GPS t)
    
    omega[11,1,t] <- 0         # Pr(alive A floater GPS t -> non-detected t)
    omega[11,2,t] <- 0         # Pr(alive A floater GPS t -> detected Y1 t) # banding occasion for fledglings
    omega[11,3,t] <- 0         # Pr(alive A floater GPS t -> detected Y2 t)
    omega[11,4,t] <- 0         # Pr(alive A floater GPS t -> detected Y3 t)
    omega[11,5,t] <- 0         # Pr(alive A floater GPS t -> detected Adult floater t)
    omega[11,6,t] <- 0         # Pr(alive A floater GPS t -> detected Adult breeder t)
    omega[11,7,t] <- 0         # Pr(alive A floater GPS t -> detected dead t)
    omega[11,8,t] <- 0         # Pr(alive A floater GPS t -> Y1 GPS t)
    omega[11,9,t] <- 0         # Pr(alive A floater GPS t -> Y2 GPS t)
    omega[11,10,t] <- 0        # Pr(alive A floater GPS t -> Y3 GPS t)
    omega[11,11,t] <- 1        # Pr(alive A floater GPS t -> Adult floater GPS t)
    omega[11,12,t] <- 0        # Pr(alive A floater GPS t -> Adult breeder GPS t)
    omega[11,13,t] <- 0        # Pr(alive A floater GPS t -> detected dead GPS t)
    
    omega[12,1,t] <- 0         # Pr(alive A breeder GPS t -> non-detected t)
    omega[12,2,t] <- 0         # Pr(alive A breeder GPS t -> detected Y1 t) # banding occasion for fledglings
    omega[12,3,t] <- 0         # Pr(alive A breeder GPS t -> detected Y2 t)
    omega[12,4,t] <- 0         # Pr(alive A breeder GPS t -> detected Y3 t)
    omega[12,5,t] <- 0         # Pr(alive A breeder GPS t -> detected Adult floater t)
    omega[12,6,t] <- 0         # Pr(alive A breeder GPS t -> detected Adult breeder t)
    omega[12,7,t] <- 0         # Pr(alive A breeder GPS t -> detected dead t)
    omega[12,8,t] <- 0         # Pr(alive A breeder GPS t -> Y1 GPS t)
    omega[12,9,t] <- 0         # Pr(alive A breeder GPS t -> Y2 GPS t)
    omega[12,10,t] <- 0        # Pr(alive A breeder GPS t -> Y3 GPS t)
    omega[12,11,t] <- 0        # Pr(alive A breeder GPS t -> Adult floater GPS t)
    omega[12,12,t] <- 1        # Pr(alive A breeder GPS t -> Adult breeder GPS t)
    omega[12,13,t] <- 0        # Pr(alive A breeder GPS t -> detected dead GPS t)
    
    omega[13,1,t] <- 0         # Pr(recently dead GPS t -> non-detected t)
    omega[13,2,t] <- 0         # Pr(recently dead GPS t -> detected Y1 t) # banding occasion for fledglings
    omega[13,3,t] <- 0         # Pr(recently dead GPS t -> detected Y2 t)
    omega[13,4,t] <- 0         # Pr(recently dead GPS t -> detected Y3 t)
    omega[13,5,t] <- 0         # Pr(recently dead GPS t -> detected Adult floater t)
    omega[13,6,t] <- 0         # Pr(recently dead GPS t -> detected Adult breeder t)
    omega[13,7,t] <- 0         # Pr(recently dead GPS t -> detected dead t)
    omega[13,8,t] <- 0         # Pr(recently dead GPS t -> Y1 GPS t)
    omega[13,9,t] <- 0         # Pr(recently dead GPS t -> Y2 GPS t)
    omega[13,10,t] <- 0        # Pr(recently dead GPS t -> Y3 GPS t)
    omega[13,11,t] <- 0        # Pr(recently dead GPS t -> Adult floater GPS t)
    omega[13,12,t] <- 0        # Pr(recently dead GPS t -> Adult breeder GPS t)
    omega[13,13,t] <- 1        # Pr(recently dead GPS t -> detected dead GPS t)
    
    omega[14,1,t] <- 1         # Pr(long dead GPS t -> non-detected t)
    omega[14,2,t] <- 0         # Pr(long dead GPS t -> detected Y1 t) # banding occasion for fledglings
    omega[14,3,t] <- 0         # Pr(long dead GPS t -> detected Y2 t)
    omega[14,4,t] <- 0         # Pr(long dead GPS t -> detected Y3 t)
    omega[14,5,t] <- 0         # Pr(long dead GPS t -> detected Adult floater t)
    omega[14,6,t] <- 0         # Pr(long dead GPS t -> detected Adult breeder t)
    omega[14,7,t] <- 0         # Pr(long dead GPS t -> detected dead t)
    omega[14,8,t] <- 0         # Pr(long dead GPS t -> Y1 GPS t)
    omega[14,9,t] <- 0         # Pr(long dead GPS t -> Y2 GPS t)
    omega[14,10,t] <- 0        # Pr(long dead GPS t -> Y3 GPS t)
    omega[14,11,t] <- 0        # Pr(long dead GPS t -> Adult floater GPS t)
    omega[14,12,t] <- 0        # Pr(long dead GPS t -> Adult breeder GPS t)
    omega[14,13,t] <- 0        # Pr(long dead GPS t -> detected dead GPS t)
  }
  
  # likelihood 
  for (i in 1:nind){
    # latent state at first capture
    z[i,first[i]]  <- first_state[i]
    for (t in (first[i]+1):nyears_mr){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:14, i, t-1])
      # y(t) given z(t)
      OBS_mr[i,t] ~ dcat(omega[z[i,t], 1:13, t-1])
      
    }
  }
  
  
  
  #### * State-space model * ########################################################
  
  # Priors on process model
  N1start ~ dunif(1,20) # number of Y1
  N2start ~ dunif(1,15) # number of Y2
  N3start ~ dunif(1,10) # number of Y3
  N4start ~ dunif(1,10) # number of AY3 floaters: previously Y3
  N5start ~ dunif(1,20) # number of AY3 floaters: previously AY3 floaters
  N6start ~ dunif(0,5)  # number of AY3 floaters: previously AY3 breeders
  N7start ~ dunif(0,5)  # number of AY3 breeders: previously Y3
  N8start ~ dunif(0,5)  # number of AY3 breeders: previously AY3 floaters
  N9start ~ dunif(1,30) # number of AY3 breeders: previously AY3 breeders
  
  N[1,1,1] <- round(N1start)
  N[2,1,1] <- round(N2start)
  N[3,1,1] <- round(N3start)
  N[4,1,1] <- round(N4start)
  N[5,1,1] <- round(N5start)
  N[6,1,1] <- round(N6start)
  N[7,1,1] <- round(N7start)
  N[8,1,1] <- round(N8start)
  N[9,1,1] <- round(N9start)
  
  
  # Process model
  for (t in 1:(nyears-1+nyears_proj)){
    N[1,t+1,1] ~ dpois((N[7,t+1,1]+N[8,t+1,1]+N[9,t+1,1]) * (mean_fec_t[t+1] / sex_ratio)) # number of juveniles (flg to Y2)
    N[2,t+1,1] ~ dbin(phi1_f_t[t+6], N[1,t,1])                                # number of Y2 (last year's chicks)
    N[3,t+1,1] ~ dbin(phi2_f_t[t+6], N[2,t,1])                                # number of Y3 (can't breed yet)
    N[4,t+1,1] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,1])                      # number of floaters: previously Y3
    N[5,t+1,1] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,1]+N[5,t,1]+N[6,t,1])) # number of floaters: previously AY3 floaters
    N[6,t+1,1] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,1]+N[8,t,1]+N[9,t,1]))     # number of floaters: previously AY3 breeders
    N[7,t+1,1] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,1])                          # number of breeders: previously Y3
    N[8,t+1,1] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,1]+N[5,t,1]+N[6,t,1]))     # number of breeders: previously AY3 floaters
    N[9,t+1,1] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,1]+N[8,t,1]+N[9,t,1])) # number of breeders: previously AY3 breeders
  }
  
  ### Management scenarios
  
  ### Scenario 2: no nestwatchers
  
  # past
  for (t in 1:(nyears)){
    N[1,t,2] <- N[1,t,1]
    N[2,t,2] <- N[2,t,1]
    N[3,t,2] <- N[3,t,1] 
    N[4,t,2] <- N[4,t,1]
    N[5,t,2] <- N[5,t,1]
    N[6,t,2] <- N[6,t,1]
    N[7,t,2] <- N[7,t,1]
    N[8,t,2] <- N[8,t,1]
    N[9,t,2] <- N[9,t,1]
  }
  
  # future 
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,2] ~ dpois((N[7,t+1,2]+N[8,t+1,2]+N[9,t+1,2]) * (mean_fec_cl_t[t+1] / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,2] ~ dbin(phi1_f_t[t+6], N[1,t,2])                                # number of Y2 (last year's chicks)
    N[3,t+1,2] ~ dbin(phi2_f_t[t+6], N[2,t,2])                                # number of Y3 (can't breed yet)
    N[4,t+1,2] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,2])                      # number of floaters: previously Y3
    N[5,t+1,2] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,2]+N[5,t,2]+N[6,t,2])) # number of floaters: previously floaters
    N[6,t+1,2] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,2]+N[8,t,2]+N[9,t,2]))     # number of floaters: previously breeders
    N[7,t+1,2] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,2])                          # number of breeders: previously Y3
    N[8,t+1,2] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,2]+N[5,t,2]+N[6,t,2]))     # number of breeders: previously floaters
    N[9,t+1,2] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,2]+N[8,t,2]+N[9,t,2])) # number of breeders: previously breeders
  } #t
  
  
  ## Scenario 3: no closures
  
  # past
  for (t in 1:(nyears)){
    N[1,t,3] <- N[1,t,1]
    N[2,t,3] <- N[2,t,1]
    N[3,t,3] <- N[3,t,1]
    N[4,t,3] <- N[4,t,1]
    N[5,t,3] <- N[5,t,1]
    N[6,t,3] <- N[6,t,1]
    N[7,t,3] <- N[7,t,1]
    N[8,t,3] <- N[8,t,1]
    N[9,t,3] <- N[9,t,1]
  }
  # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,3] ~ dpois((N[7,t+1,3]+N[8,t+1,3]+N[9,t+1,3]) * (mean_fec_nw_t[t+1] / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,3] ~ dbin(phi1_f_t[t+6], N[1,t,3])                                # number of Y2 (last year's chicks)
    N[3,t+1,3] ~ dbin(phi2_f_t[t+6], N[2,t,3])                                # number of Y3 (can't breed yet)
    N[4,t+1,3] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,3])                      # number of floaters: previously Y3
    N[5,t+1,3] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,3]+N[5,t,3]+N[6,t,3])) # number of floaters: previously floaters
    N[6,t+1,3] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,3]+N[8,t,3]+N[9,t,3]))     # number of floaters: previously breeders
    N[7,t+1,3] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,3])                          # number of breeders: previously Y3
    N[8,t+1,3] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,3]+N[5,t,3]+N[6,t,3]))     # number of breeders: previously floaters
    N[9,t+1,3] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,3]+N[8,t,3]+N[9,t,3])) # number of breeders: previously breeders
  } #t
  
  ## Scenario 4: no managment 
  
  # past
  for (t in 1:(nyears)){
    N[1,t,4] <- N[1,t,1]
    N[2,t,4] <- N[2,t,1]
    N[3,t,4] <- N[3,t,1]
    N[4,t,4] <- N[4,t,1]
    N[5,t,4] <- N[5,t,1]
    N[6,t,4] <- N[6,t,1]
    N[7,t,4] <- N[7,t,1]
    N[8,t,4] <- N[8,t,1]
    N[9,t,4] <- N[9,t,1]
  }
  # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,4] ~ dpois((N[7,t+1,4]+N[8,t+1,4]+N[9,t+1,4]) * (mean_fec_nomgmt_t[t+1] / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,4] ~ dbin(phi1_f_t[t+6], N[1,t,4])                                # number of Y2 (last year's chicks)
    N[3,t+1,4] ~ dbin(phi2_f_t[t+6], N[2,t,4])                                # number of Y3 (can't breed yet)
    N[4,t+1,4] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,4])                      # number of floaters: previously Y3
    N[5,t+1,4] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,4]+N[5,t,4]+N[6,t,4])) # number of floaters: previously floaters
    N[6,t+1,4] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,4]+N[8,t,4]+N[9,t,4]))     # number of floaters: previously breeders
    N[7,t+1,4] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,4])                          # number of breeders: previously Y3
    N[8,t+1,4] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,4]+N[5,t,4]+N[6,t,4]))     # number of breeders: previously floaters
    N[9,t+1,4] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,4]+N[8,t,4]+N[9,t,4])) # number of breeders: previously breeders
  } #t
  
  ### Fec change scenarios (not linked to a particular mgmt action)
  
  ## Scenario 5: reduce fec 7%
  
  # past
  for (t in 1:(nyears)){
    N[1,t,5] <- N[1,t,1]
    N[2,t,5] <- N[2,t,1]
    N[3,t,5] <- N[3,t,1]
    N[4,t,5] <- N[4,t,1]
    N[5,t,5] <- N[5,t,1]
    N[6,t,5] <- N[6,t,1]
    N[7,t,5] <- N[7,t,1]
    N[8,t,5] <- N[8,t,1]
    N[9,t,5] <- N[9,t,1]
  }
  # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,5] ~ dpois((N[7,t+1,5]+N[8,t+1,5]+N[9,t+1,5]) * (mean_fec_t[t+1] * 0.93 / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,5] ~ dbin(phi1_f_t[t+6], N[1,t,5])                                # number of Y2 (last year's chicks)
    N[3,t+1,5] ~ dbin(phi2_f_t[t+6], N[2,t,5])                                # number of Y3 (can't breed yet)
    N[4,t+1,5] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,5])                      # number of floaters: previously Y3
    N[5,t+1,5] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,5]+N[5,t,5]+N[6,t,5])) # number of floaters: previously floaters
    N[6,t+1,5] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,5]+N[8,t,5]+N[9,t,5]))     # number of floaters: previously breeders
    N[7,t+1,5] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,5])                          # number of breeders: previously Y3
    N[8,t+1,5] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,5]+N[5,t,5]+N[6,t,5]))     # number of breeders: previously floaters
    N[9,t+1,5] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,5]+N[8,t,5]+N[9,t,5])) # number of breeders: previously breeders
  } #t
  
  ## Scenario 6: reduce fec 8%
  
  # past
  for (t in 1:(nyears)){
    N[1,t,6] <- N[1,t,1]
    N[2,t,6] <- N[2,t,1]
    N[3,t,6] <- N[3,t,1]
    N[4,t,6] <- N[4,t,1]
    N[5,t,6] <- N[5,t,1]
    N[6,t,6] <- N[6,t,1]
    N[7,t,6] <- N[7,t,1]
    N[8,t,6] <- N[8,t,1]
    N[9,t,6] <- N[9,t,1]
  }
  # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,6] ~ dpois((N[7,t+1,6]+N[8,t+1,6]+N[9,t+1,6]) * (mean_fec_t[t+1] * 0.92 / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,6] ~ dbin(phi1_f_t[t+6], N[1,t,6])                                # number of Y2 (last year's chicks)
    N[3,t+1,6] ~ dbin(phi2_f_t[t+6], N[2,t,6])                                # number of Y3 (can't breed yet)
    N[4,t+1,6] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,6])                      # number of floaters: previously Y3
    N[5,t+1,6] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,6]+N[5,t,6]+N[6,t,6])) # number of floaters: previously floaters
    N[6,t+1,6] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,6]+N[8,t,6]+N[9,t,6]))     # number of floaters: previously breeders
    N[7,t+1,6] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,6])                          # number of breeders: previously Y3
    N[8,t+1,6] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,6]+N[5,t,6]+N[6,t,6]))     # number of breeders: previously floaters
    N[9,t+1,6] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,6]+N[8,t,6]+N[9,t,6])) # number of breeders: previously breeders
  } #t
  
  ## Scenario 7: reduce fec 9%
  
  # past
  for (t in 1:(nyears)){
    N[1,t,7] <- N[1,t,1]
    N[2,t,7] <- N[2,t,1]
    N[3,t,7] <- N[3,t,1]
    N[4,t,7] <- N[4,t,1]
    N[5,t,7] <- N[5,t,1]
    N[6,t,7] <- N[6,t,1]
    N[7,t,7] <- N[7,t,1]
    N[8,t,7] <- N[8,t,1]
    N[9,t,7] <- N[9,t,1]
  }
  # # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,7] ~ dpois((N[7,t+1,7]+N[8,t+1,7]+N[9,t+1,7]) * (mean_fec_t[t+1] * 0.91 / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,7] ~ dbin(phi1_f_t[t+6], N[1,t,7])                                # number of Y2 (last year's chicks)
    N[3,t+1,7] ~ dbin(phi2_f_t[t+6], N[2,t,7])                                # number of Y3 (can't breed yet)
    N[4,t+1,7] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,7])                      # number of floaters: previously Y3
    N[5,t+1,7] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,7]+N[5,t,7]+N[6,t,7])) # number of floaters: previously floaters
    N[6,t+1,7] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,7]+N[8,t,7]+N[9,t,7]))     # number of floaters: previously breeders
    N[7,t+1,7] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,7])                          # number of breeders: previously Y3
    N[8,t+1,7] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,7]+N[5,t,7]+N[6,t,7]))     # number of breeders: previously floaters
    N[9,t+1,7] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,7]+N[8,t,7]+N[9,t,7])) # number of breeders: previously breeders
  } #t
  
  ## Scenario 8: reduce fec 10%
  
  # past
  for (t in 1:(nyears)){
    N[1,t,8] <- N[1,t,1]
    N[2,t,8] <- N[2,t,1]
    N[3,t,8] <- N[3,t,1]
    N[4,t,8] <- N[4,t,1]
    N[5,t,8] <- N[5,t,1]
    N[6,t,8] <- N[6,t,1]
    N[7,t,8] <- N[7,t,1]
    N[8,t,8] <- N[8,t,1]
    N[9,t,8] <- N[9,t,1]
  }
  # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,8] ~ dpois((N[7,t+1,8]+N[8,t+1,8]+N[9,t+1,8]) * (mean_fec_t[t+1] * 0.90 / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,8] ~ dbin(phi1_f_t[t+6], N[1,t,8])                                # number of Y2 (last year's chicks)
    N[3,t+1,8] ~ dbin(phi2_f_t[t+6], N[2,t,8])                                # number of Y3 (can't breed yet)
    N[4,t+1,8] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,8])                      # number of floaters: previously Y3
    N[5,t+1,8] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,8]+N[5,t,8]+N[6,t,8])) # number of floaters: previously floaters
    N[6,t+1,8] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,8]+N[8,t,8]+N[9,t,8]))     # number of floaters: previously breeders
    N[7,t+1,8] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,8])                          # number of breeders: previously Y3
    N[8,t+1,8] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,8]+N[5,t,8]+N[6,t,8]))     # number of breeders: previously floaters
    N[9,t+1,8] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,8]+N[8,t,8]+N[9,t,8])) # number of breeders: previously breeders
  } #t
  
  ## Scenario 9: reduce fec 11%
  
  # past
  for (t in 1:(nyears)){
    N[1,t,9] <- N[1,t,1]
    N[2,t,9] <- N[2,t,1]
    N[3,t,9] <- N[3,t,1]
    N[4,t,9] <- N[4,t,1]
    N[5,t,9] <- N[5,t,1]
    N[6,t,9] <- N[6,t,1]
    N[7,t,9] <- N[7,t,1]
    N[8,t,9] <- N[8,t,1]
    N[9,t,9] <- N[9,t,1]
  }
  # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,9] ~ dpois((N[7,t+1,9]+N[8,t+1,9]+N[9,t+1,9]) * (mean_fec_t[t+1] * 0.89 / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,9] ~ dbin(phi1_f_t[t+6], N[1,t,9])                                # number of Y2 (last year's chicks)
    N[3,t+1,9] ~ dbin(phi2_f_t[t+6], N[2,t,9])                                # number of Y3 (can't breed yet)
    N[4,t+1,9] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,9])                      # number of floaters: previously Y3
    N[5,t+1,9] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,9]+N[5,t,9]+N[6,t,9])) # number of floaters: previously floaters
    N[6,t+1,9] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,9]+N[8,t,9]+N[9,t,9]))     # number of floaters: previously breeders
    N[7,t+1,9] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,9])                          # number of breeders: previously Y3
    N[8,t+1,9] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,9]+N[5,t,9]+N[6,t,9]))     # number of breeders: previously floaters
    N[9,t+1,9] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,9]+N[8,t,9]+N[9,t,9])) # number of breeders: previously breeders
  } #t
  
  
  ## Scenario 10: reduce fec 12%
  
  # past
  for (t in 1:(nyears)){
    N[1,t,10] <- N[1,t,1]
    N[2,t,10] <- N[2,t,1]
    N[3,t,10] <- N[3,t,1]
    N[4,t,10] <- N[4,t,1]
    N[5,t,10] <- N[5,t,1]
    N[6,t,10] <- N[6,t,1]
    N[7,t,10] <- N[7,t,1]
    N[8,t,10] <- N[8,t,1]
    N[9,t,10] <- N[9,t,1]
  }
  # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,10] ~ dpois((N[7,t+1,10]+N[8,t+1,10]+N[9,t+1,10]) * (mean_fec_t[t+1] * 0.88 / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,10] ~ dbin(phi1_f_t[t+6], N[1,t,10])                                  # number of Y2 (last year's chicks)
    N[3,t+1,10] ~ dbin(phi2_f_t[t+6], N[2,t,10])                                  # number of Y3 (can't breed yet)
    N[4,t+1,10] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,10])                        # number of floaters: previously Y3
    N[5,t+1,10] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,10]+N[5,t,10]+N[6,t,10])) # number of floaters: previously floaters
    N[6,t+1,10] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,10]+N[8,t,10]+N[9,t,10]))     # number of floaters: previously breeders
    N[7,t+1,10] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,10])                            # number of breeders: previously Y3
    N[8,t+1,10] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,10]+N[5,t,10]+N[6,t,10]))     # number of breeders: previously floaters
    N[9,t+1,10] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,10]+N[8,t,10]+N[9,t,10])) # number of breeders: previously breeders
  } #t
  
  ## Scenario 11: reduce fec 13%
  
  # past
  for (t in 1:(nyears)){
    N[1,t,11] <- N[1,t,1]
    N[2,t,11] <- N[2,t,1]
    N[3,t,11] <- N[3,t,1]
    N[4,t,11] <- N[4,t,1]
    N[5,t,11] <- N[5,t,1]
    N[6,t,11] <- N[6,t,1]
    N[7,t,11] <- N[7,t,1]
    N[8,t,11] <- N[8,t,1]
    N[9,t,11] <- N[9,t,1]
  }
  # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,11] ~ dpois((N[7,t+1,11]+N[8,t+1,11]+N[9,t+1,11]) * (mean_fec_t[t+1] * 0.87 / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,11] ~ dbin(phi1_f_t[t+6], N[1,t,11])                                  # number of Y2 (last year's chicks)
    N[3,t+1,11] ~ dbin(phi2_f_t[t+6], N[2,t,11])                                  # number of Y3 (can't breed yet)
    N[4,t+1,11] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,11])                        # number of floaters: previously Y3
    N[5,t+1,11] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,11]+N[5,t,11]+N[6,t,11])) # number of floaters: previously floaters
    N[6,t+1,11] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,11]+N[8,t,11]+N[9,t,11]))     # number of floaters: previously breeders
    N[7,t+1,11] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,11])                            # number of breeders: previously Y3
    N[8,t+1,11] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,11]+N[5,t,11]+N[6,t,11]))     # number of breeders: previously floaters
    N[9,t+1,11] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,11]+N[8,t,11]+N[9,t,11])) # number of breeders: previously breeders
  } #t
  
  ## Scenario 12: reduce fec 14%
  
  # past
  for (t in 1:(nyears)){
    N[1,t,12] <- N[1,t,1]
    N[2,t,12] <- N[2,t,1]
    N[3,t,12] <- N[3,t,1]
    N[4,t,12] <- N[4,t,1]
    N[5,t,12] <- N[5,t,1]
    N[6,t,12] <- N[6,t,1]
    N[7,t,12] <- N[7,t,1]
    N[8,t,12] <- N[8,t,1]
    N[9,t,12] <- N[9,t,1]
  }
  # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,12] ~ dpois((N[7,t+1,12]+N[8,t+1,12]+N[9,t+1,12]) * (mean_fec_t[t+1] * 0.86 / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,12] ~ dbin(phi1_f_t[t+6], N[1,t,12])                                  # number of Y2 (last year's chicks)
    N[3,t+1,12] ~ dbin(phi2_f_t[t+6], N[2,t,12])                                  # number of Y3 (can't breed yet)
    N[4,t+1,12] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,12])                        # number of floaters: previously Y3
    N[5,t+1,12] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,12]+N[5,t,12]+N[6,t,12])) # number of floaters: previously floaters
    N[6,t+1,12] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,12]+N[8,t,12]+N[9,t,12]))     # number of floaters: previously breeders
    N[7,t+1,12] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,12])                            # number of breeders: previously Y3
    N[8,t+1,12] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,12]+N[5,t,12]+N[6,t,12]))     # number of breeders: previously floaters
    N[9,t+1,12] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,12]+N[8,t,12]+N[9,t,12])) # number of breeders: previously breeders
  } #t
  
  ## Scenario 13: reduce fec 15%
  
  # past
  for (t in 1:(nyears)){
    N[1,t,13] <- N[1,t,1]
    N[2,t,13] <- N[2,t,1]
    N[3,t,13] <- N[3,t,1]
    N[4,t,13] <- N[4,t,1]
    N[5,t,13] <- N[5,t,1]
    N[6,t,13] <- N[6,t,1]
    N[7,t,13] <- N[7,t,1]
    N[8,t,13] <- N[8,t,1]
    N[9,t,13] <- N[9,t,1]
  }
  # future
  for (t in nyears:(nyears-1+nyears_proj)){
    N[1,t+1,13] ~ dpois((N[7,t+1,13]+N[8,t+1,13]+N[9,t+1,13]) * (mean_fec_t[t+1] * 0.85 / sex_ratio))  # number of juveniles (flg to Y2)
    N[2,t+1,13] ~ dbin(phi1_f_t[t+6], N[1,t,13])                                  # number of Y2 (last year's chicks)
    N[3,t+1,13] ~ dbin(phi2_f_t[t+6], N[2,t,13])                                  # number of Y3 (can't breed yet)
    N[4,t+1,13] ~ dbin(phi3_f_t[t+6]*(1-psiFB), N[3,t,13])                        # number of floaters: previously Y3
    N[5,t+1,13] ~ dbin(phiAf_f_t[t+6]*(1-psiFB), (N[4,t,13]+N[5,t,13]+N[6,t,13])) # number of floaters: previously floaters
    N[6,t+1,13] ~ dbin(phiAb_f_t[t+6]*psiBF, (N[7,t,13]+N[8,t,13]+N[9,t,13]))     # number of floaters: previously breeders
    N[7,t+1,13] ~ dbin(phi3_f_t[t+6]*psiFB, N[3,t,13])                            # number of breeders: previously Y3
    N[8,t+1,13] ~ dbin(phiAf_f_t[t+6]*psiFB, (N[4,t,13]+N[5,t,13]+N[6,t,13]))     # number of breeders: previously floaters
    N[9,t+1,13] ~ dbin(phiAb_f_t[t+6]*(1-psiBF), (N[7,t,13]+N[8,t,13]+N[9,t,13])) # number of breeders: previously breeders
  } #t
  
  # Observation model
  
  # Prior
  sd_occ ~ dunif(0,5) 
  
  for (t in 1:nyears){  
    num_occ_OBS[t] ~ T(dnorm(N[7,t,1] + N[8,t,1] + N[9,t,1], sd = sd_occ), 0, )
  }
  
  
  #### * Nest productivity * ###################################################
  
  # Priors 
  rho[1] ~ dnorm(mean = 0, sd = 3) # intercept
  rho[2] ~ dnorm(mean = 0, sd = 3) # nestwatcher coefficient
  rho[3] ~ dnorm(mean = 0, sd = 3) # closure coefficient
  sdeps[1] ~ dunif(0,10)
  sdeps[2] ~ dunif(0,10)
  
  # Random effect for timestep, residual variation not explained by the covariates
  for (t in 1:nyears_prod+nyears_proj) {
    YReps[t] ~ dnorm(0, sd = sdeps[1]) 
  }
  
  # Random effect for breeding area (i.e. site)
  for (i in 1:nsites_prod) {
    SITEeps[i] ~ dnorm(0, sd = sdeps[2]) 
  }
  
  # Ecological model
  for (i in 1:nbrood){
    log(fec[i]) <- rho[1] + rho[2] * NW_prod[i] + rho[3] * closures[i] + 
      YReps[yrID_prod[i]] + SITEeps[siteID_prod[i]]
  } #i
  
  
  # Likelihood
  for (i in 1:nbrood){
    OBS_nestlings[i] ~ dpois(fec[i]) 
  }
  
  # Retro yearly means
  for (t in 1:nyears_prod){
    mean_fec_t[t] <- exp(rho[1] + rho[2] * propNW_t[t] + rho[3] * propCL_t[t] + YReps[t]) 
    mean_fec_nomgmt_t[t] <- exp(rho[1] + YReps[t]) 
    mean_fec_nw_t[t] <- exp(rho[1] + (rho[2] * propNW_t[t]) + YReps[t]) 
    mean_fec_cl_t[t] <- exp(rho[1] + (rho[3] * propCL_t[t]) + YReps[t]) 
  }
  
  # Future yearly means
  for (t in (nyears_prod+1):(nyears_prod+nyears_proj)){
    mean_fec_t[t] <- exp(rho[1] + rho[2] * propNW + rho[3] * propClos + YReps[t]) 
    mean_fec_nomgmt_t[t] <- exp(rho[1] + YReps[t]) 
    mean_fec_nw_t[t] <- exp(rho[1] + (rho[2] * propNW) + YReps[t]) 
    mean_fec_cl_t[t] <- exp(rho[1] + (rho[3] * propClos) + YReps[t]) 
  }
  
  # Overall mean (under current practices)
  mean_fec <-  exp(rho[1] + rho[2] * propNW + rho[3] * propClos) 
  # expected mean without closures
  mean_fec_nw <- exp(rho[1] + rho[2] * propNW)
  # expected mean without nestwatchers
  mean_fec_cl <- exp(rho[1] + rho[3] * propClos)
  # expected mean without management
  mean_fec_nomgmt <-  exp(rho[1])
  

  #### * Derived quantities * ##################################################
  
  # Number of breeders
  for (j in 1:scenarios){
    for (t in 1:(nyears+nyears_proj)){
      n_breeders[j,t] <- (N[7,t,j] + N[8,t,j]+ N[9,t,j])
    } #t
  } #j
  
  # Number of floaters
  for (j in 1:scenarios){
    for (t in 1:(nyears+nyears_proj)){
      n_floaters[t,j] <- (N[4,t,j] + N[5,t,j]+ N[6,t,j])
    } #t
  } #j
  
  # Number of effective breeders (i.e. sexually mature individuals, floaters and breeders)
  for (j in 1:scenarios){
    for (t in 1:nyears+nyears_proj){
      n_eff_breeders[t,j] <- N[4,t,j] + N[5,t,j] + N[6,t,j] + N[7,t,j] + N[8,t,j]+ N[9,t,j]
    } #t
  } #j
  
  # Number total
  for (t in 1:nyears+nyears_proj){
    for (j in 1:scenarios){
      n_total[t,j] <- (N[1,t,j] + N[2,t,j] + N[3,t,j] + N[4,t,j] + N[5,t,j] +
                         N[6,t,j] + N[7,t,j] + N[8,t,j] + N[9,t,j])
    } #j
  } #t
  
  # population growth rate
  for (t in 1:(nyears+nyears_proj-1)){
    for (j in 1:scenarios){
      lambda[t,j] <- n_total[t+1,j] / (n_total[t,j] + 0.000001)
    } #j
  } #t
  
  for (j in 1:scenarios){
    geomean_lambda[j] <- exp(mean(log(lambda[1:(nyears-1), j])))
    geomean_lambda_proj[j] <- exp(mean(log(lambda[nyears:(nyears+nyears_proj-1), j])))
  }
  
  
}) #end

