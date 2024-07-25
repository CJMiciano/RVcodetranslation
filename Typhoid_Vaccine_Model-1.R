######################## Question 6 ################################

######################## Typhoid Vaccine model  ########################

typhoid_vaccine_model <- function(times,yinit,pars){
  
  with(as.list(c(yinit,pars)), {
    # Initialize the model with the equilibrium state variables (at the end of the burn-in period)
    # Run the model for 10 years to get the "no vaccination" number of people in each state
    # Run the model for another 10 years with different vaccine coverage levels for different age classes
    
    # get total population
    N = sum(yinit)
    
    #set up vaccination vector
    if(times<10*12){
      kr = 0
      kc = 0
    } else {
      if(times<11*12) { #vaccination campaign is implemented over the first year
        kr = cov_routine
        kc = cov_campaign
      }
      else { #routine vaccination continues until end of the simulation
        kr = cov_routine
        kc = 0
      }
    }
    
    # Define each model state as a vector of states for each age group
    S = c(S1, S2, S3)
    I_S = c(I_S1, I_S2, I_S3)
    R = c(R1, R2, R3)
    C = c(C1, C2, C3)
    S_R = c(S_R1, S_R2, S_R3)
    I_A = c(I_A1, I_A2, I_A3)
    V_1= c(V_11, V_12, V_13)
    V_2 = c(V_21, V_22, V_23)
    
    # put demographic parameters into vectors for easier calculation
    births <- c(mu_0*N, 0, 0)
    age_in <- c(0, mu_1, mu_2) #demographic in-flow vector for non-susceptible states (multiplies by age group below)
    age_out = c(mu_1 ,mu_2, 0) #demographic outflow vector for all states
    die <- mu_3
    
    # Calculate force of infection vector
    lambda<- beta_0*sum(I_S, r_a*I_A, r_c*C)
    
    #differential equations
    # Susceptibles
    dS <- births - lambda*S  + w_v*V_1 + age_in*c(1,(1-kr*ve),1)*c(0,S[1],S[2]) - c(0,kc*ve,0)*S - age_out*S - die*S 
    # Infecteds
    dI_S <- lambda*S - delta*I_S + age_in*c(0,I_S[1],I_S[2]) - age_out*I_S - die*I_S
    # Recovereds
    dR <- delta*(1-theta)*I_S + delta*I_A - w*R + age_in*c(1,(1-kr*ve),1)*c(0,R[1],R[2]) - c(0,kc*ve,0)*R - age_out*R - die*R
    # Carriers
    dC <- delta*theta*I_S  + age_in*c(0,C[1],C[2]) - age_out*C - die*C
    # Susceptible -- Previously Infected
    dS_R <- w*R - lambda*S_R + w_v*V_2 + age_in*c(1,(1-kr*ve),1)*c(0,S_R[1],S_R[2]) - c(0,kc*ve,0)*S_R - age_out*S_R - die*S_R
    # Infected -- Previously Infected (Asymptomatic)
    dI_A <- lambda*S_R - delta*I_A + age_in*c(0,I_A[1],I_A[2]) - age_out*I_A - die*I_A
    # Vaccinated -- Never Infected
    dV_1 <- age_in*c(0,kr*ve,0)*c(0,S[1],S[2]) + c(0,kc*ve,0)*S - w_v*V_1 + age_in*c(0,V_1[1],V_1[2]) - age_out*V_1 - die*V_1
    # Vaccinated -- Previously Infected
    dV_2 <- age_in*c(0,kr*ve,0)*(c(0,S_R[1],S_R[2]) + c(0,R[1],R[2])) + c(0,kc*ve,0)*(S_R+R) - w_v*V_2 + age_in*c(0,V_2[1],V_2[2]) - age_out*V_2 - die*V_2
    
    return(list(c(dS, dI_S, dR, dC, dS_R, dI_A, dV_1, dV_2)))})	
  
}
