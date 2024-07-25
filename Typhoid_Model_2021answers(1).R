library(deSolve)

######################## Typhoid model (no age structure, no vaccination)  ########################

typhoid_model <- function(times,yinit,pars){
  
  with(as.list(c(yinit,pars)), {
    # We need to run the model for a "burn-in period" of 200 years 
    # to get the "equilibrium" number of people in each state

    # get total population
    N = sum(yinit)
    
    # define demographic parameters 
    births <- mu_0*N
    die <- mu_0
    
    # Calculate force of infection vector
    lambda<- beta_0*sum(I_S, r_a*I_A, r_c*C) # =(I_S + r_a*I_A + r_C*I_C) 
    
    #differential equations
    # Susceptibles
    dS <- births - lambda*S - die*S 
    # Infecteds
    dI_S <- lambda*S - delta*I_S - die*I_S
    # Recovereds
    dR <- delta*(1-theta)*I_S + delta*I_A - w*R - die*R
    # Carriers
    dC <- delta*theta*I_S - die*C
    # Susceptible -- Previously Infected
    dS_R <- w*R - lambda*S_R - die*S_R
    # Infected -- Previously Infected (Asymptomatic)
    dI_A <- lambda*S_R - delta*I_A - die*I_A

    return(list(c(dS, dI_S, dR, dC, dS_R, dI_A)))})	
  
}
	
################## Setting Parameters ##################

#Demographic parameters (times in months)
pop0<- 1.35e9 # population size
mu_0 = 0.0015  # birth rate = death rate

#Transmission Parameter
beta_0 = .75e-9 # baseline transmission rate

#Infection parameters
delta = 1 #1/duration of infectiousness
theta = .03 # proportion of infected who become carriers
r_c = .25 # relative infectiousness of chronic carriers
r_a = .33 # relative infectiousness of asymptomatic (subsequent) infections
w = 1/12 # 1/duration of immunity from natural infection (rate at which immunity wanes after infection)
rep = 0.1 # probability of symptomatic disease (first infections)

# Save the parameter values
pars <- cbind(mu_0, beta_0, delta, theta, r_c, r_a, w)

############### Setting Time Frame ####################

start_burnin = 0 # start date
end_burnin = 200*12 #NOTE: time is in months. Set to 100 for burn-in, and extend to run longer
times_burnin <- seq(start_burnin, end_burnin, by = 1) # gives a sequence from start to end in increments of 1

############ Setting Initial Conditions ###############

init_prev <- .0000001   #initial prevalence (assume epidemic begins with very low prevalence)

# Create Matrix of Initial State Variables: 6 Disease Compartments
yinit_in <- c(S =(1-init_prev)*pop0, # assume everyone except those infected are susceptible to begin with
              I_S = init_prev*pop0, 
              R=0, 
              C=0, 
              S_R=0, 
              I_A=0)

############# Running The Model ####################

results <- as.data.frame(ode(y=yinit_in, times=times_burnin, func=typhoid_model, parms=pars, method='rk4'))

## calculate incident cases

# Calculate the incidence rate of new infections (aka the "force of infection")
incidence_rate_lambda = beta_0*(results[,'I_S'] + r_a*results[,'I_A'] + r_c*results[,'C'])

# Calculate the number of new cases at each time point 
incidence_total = rep*results[,'S']*incidence_rate_lambda

############# Plotting Results ######################

## State Variable Plots for Each Age Group ##
plot(times_burnin, results$I_S, type= "l", col = 'black',ylab="Individuals",xlab="Time (Months)", main = "Typhoid States", ylim=c(0,1.01*max(results$I_S,results$R,results$C,results$I_A)))
lines(times_burnin, results$R, col = 'red')
lines(times_burnin, results$C, col = 'blue' )
lines(times_burnin, results$I_A, col = 'green')
#lines(times_burnin, results$S, col = 'orange')

legend("topleft",legend=c("Infected (Symptomatic)","Recovered","Carriers","Infected (Asymptomatic)"), 
       col=c("black","red","blue","green"),
       lty=1, ncol=1, cex =.3, pt.cex = 1)

## Incidence Plot ##
plot(times_burnin, incidence_total,type= "l", col='black', ylab="Cases",xlab="Time (Months)",main="Incident Cases, Total")

############ Question 3 ##############

## Incidence rate calculation
## rate per 100,000 is calculated as (cases/population)*100,000

# Get annual incidence rate (per 100K people) at end of burn-in period 
incidence =12*sum(incidence_total[end_burnin])/pop0*1e5
print(incidence)

######################## Typhoid age model (age structure, no vaccination)  ########################

typhoid_age_model <- function(times,yinit,pars){
  
  with(as.list(c(yinit,pars)), {
    # We need to run the model for a "burn-in period" of 100 years 
    # to get the "equilibrium" number of people in each state
    
    # get total population
    N = sum(yinit)
    
    # Define each model state as a vector of states for each age group
    S = c(S1, S2, S3)
    I_S = c(I_S1, I_S2, I_S3)
    R = c(R1, R2, R3)
    C = c(C1, C2, C3)
    S_R = c(S_R1, S_R2, S_R3)
    I_A = c(I_A1, I_A2, I_A3)
    
    # put demographic parameters into vectors for easier calculation
    births <- c(mu_0*N, 0, 0)
    age_in <- c(0, mu_1, mu_2) #demographic in-flow vector for non-susceptible states (multiplies by age group below)
    age_out = c(mu_1 ,mu_2, 0) #demographic outflow vector for all states
    die <- mu_3
    
    # Calculate force of infection vector
    lambda<- beta_0*sum(I_S, r_a*I_A, r_c*C)
    
    #differential equations
    # Susceptibles
    dS <- births - lambda*S  + age_in*c(0,S[1],S[2]) - age_out*S - die*S 
    # Infecteds
    dI_S <- lambda*S - delta*I_S + age_in*c(0,I_S[1],I_S[2]) - age_out*I_S - die*I_S
    # Recovereds
    dR <- delta*(1-theta)*I_S + delta*I_A - w*R + age_in*c(0,R[1],R[2]) - age_out*R - die*R
    # Carriers
    dC <- delta*theta*I_S  + age_in*c(0,C[1],C[2]) - age_out*C - die*C
    # Susceptible -- Previously Infected
    dS_R <- w*R - lambda*S_R + age_in*c(0,S_R[1],S_R[2]) - age_out*S_R - die*S_R
    # Infected -- Previously Infected (Asymptomatic)
    dI_A <- lambda*S_R - delta*I_A + age_in*c(0,I_A[1],I_A[2]) - age_out*I_A - die*I_A
    
    return(list(c(dS, dI_S, dR, dC, dS_R, dI_A)))})	
  
}

## Additional demographic parameters (times in months)
pop_dist <- c(.1,.3,.6)  #initial population age distribution (this will also equilibrate during the burn-in period)

mu_0 = 0.0015  # Birth rate
mu_1 = 1/9   # Rate at which people age out of first age class/into second age class
mu_2 = 1/(14.25*12) # Rate at which people age out of second age class/into third age class
mu_3 = mu_0   # set death rate equal to birth rate to keep population constant

# Save the parameter values
pars_age <- cbind(mu_0, mu_1, mu_2, mu_3, beta_0, delta, theta, r_c, r_a, w)

## Initial Conditions 

# Create Matrix of Initial State Variables: 3 Age Groups, 6 Disease Compartments
yinit_in <- c(S =(1-init_prev)*pop0*pop_dist, # assume everyone except those infected are susceptible to begin with
              I_S = init_prev*pop0*pop_dist, 
              R=rep(0,3), 
              C=rep(0,3), 
              S_R=rep(0,3), 
              I_A=rep(0,3))

############# Running The Model ####################

results <- as.data.frame(ode(y=yinit_in, times=times_burnin, func=typhoid_age_model, parms=pars_age,
	method='rk4'))

## calculate total and age-specific incident cases

# Calculate the incidence rate of new infections (aka the "force of infection")
incidence_rate_lambda = beta_0*(rowSums(results[,c('I_S1','I_S2','I_S3')]) + 
                                  r_a*rowSums(results[,c('I_A1','I_A2','I_A3')])+
                                  r_c*rowSums(results[,c('C1','C2','C3')]))

# Calculate the number of new cases at each time point overall and for each age group
incidence_total = rep*rowSums(results[,c('S1','S2','S3')])*incidence_rate_lambda

incidence1 = rep*results$S1*incidence_rate_lambda
incidence2 = rep*results$S2*incidence_rate_lambda
incidence3 = rep*results$S3*incidence_rate_lambda


############# Plotting Results ######################

## State Variable Plots for Each Age Group ##
plot(times_burnin, results$I_S1, type= "l", col = 'black',ylab="Individuals",
     xlab="Time (Months)", main = "Typhoid in Age Class 1", 
     ylim=c(0,1.01*max(results$I_S1,results$R1,results$C1,results$I_A1)))
lines(times_burnin, results$R1, col = 'red')
lines(times_burnin, results$C1, col = 'blue' )
lines(times_burnin, results$I_A1, col = 'green')
legend("topleft",legend=c("Infected (Symptomatic)","Recovered","Carriers","Infected (Asymptomatic)"), 
       col=c("black","red","blue","green"),
       lty=1, ncol=1, cex =.3, pt.cex = 1)

plot(times_burnin, results$I_S2, type= "l", col = 'black',ylab="Individuals",
     xlab="Time (Months)", main = "Typhoid in Age Class 2", 
     ylim=c(0,1.01*max(results$I_S2,results$R2,results$C2,results$I_A2)))
lines(times_burnin, results$R2, col = 'red')
lines(times_burnin, results$C2, col = 'blue' )
lines(times_burnin, results$I_A2, col = 'green')
legend("topleft",legend=c("Infected (Symptomatic)","Recovered","Carriers","Infected (Asymptomatic)"), 
       col=c("black","red","blue","green"),
       lty=1, ncol=1, cex =.3, pt.cex = 1)

plot(times_burnin, results$I_S3, type= "l", col = 'black',ylab="Individuals",
     xlab="Time (Months)", main = "Typhoid in Age Class 3", 
     ylim=c(0,1.01*max(results$I_S3,results$R3,results$C3,results$I_A3)))
lines(times_burnin, results$R3, col = 'red')
lines(times_burnin, results$C3, col = 'blue' )
lines(times_burnin, results$I_A3, col = 'green')
legend("topleft",legend=c("Infected (Symptomatic)","Recovered","Carriers","Infected (Asymptomatic)"), 
       col=c("black","red","blue","green"),
       lty=1, ncol=1, cex =.3, pt.cex = 1)

## Incidence Plots ##
plot(times_burnin, incidence_total,type= "l", col='black', ylab="Cases",xlab="Time (Months)",
     main="Incident Cases, Total")

plot(times_burnin, incidence1, type= "l", col='red', ylab="Cases",xlab="Time (Months)",
     main="Incident Cases, By Age", ylim=c(0,1.01*max(incidence1,incidence2,incidence3)))
lines(times_burnin, incidence2, col='green')
lines(times_burnin, incidence3, col='blue')
legend("topleft",legend=c("0-9mo","9mo-15y","15y+"), 
       col=c("red","green","blue"),
       lty=1, ncol=1, cex =.3, pt.cex = 1)

############ Question 4 ##############

## Incidence rate calculation
## rate per 100,000 is calculated as (cases/population)*100,000

# Get incidence at end of burn-in period (annual incidence ~ 12 months/year * monthly incidence)
incidence = 12*incidence_total[end_burnin]/pop0*1e5

# Average age of infection: age-group midpoint age * weight, weight = incident cases in age group/total incident cases
Avg_AoI = (.375*incidence1[end_burnin] + 7.5*incidence2[end_burnin] + (68.5/2)*incidence3[end_burnin])/incidence_total[end_burnin]

print(incidence)
print(Avg_AoI)

######################## Question 5 ################################

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

#Vaccination parameters
cov_routine = .9   # routine vaccination coverage 
cov_campaign = .1   # catch-up coverage for second age class  
ve = .875 # Initial vaccine efficacy
w_v = 1/(15*12) # 1/duration of immunity from vaccine

# Save the parameter values
pars_vac <- cbind(mu_0, mu_1, mu_2, mu_3, beta_0, delta, theta, r_c, r_a, w, cov_routine, cov_campaign, ve, w_v)

## Set Time Frame 

start_time = 1 # start date
vac_time = 10*12 # date of vaccine introduction
end_time = 20*12 # end date
times <- seq(start_time, end_time, by = 1) # gives a sequence from start to end in increments of 1

## Initial Conditions 
# Create Matrix of Initial State Variables: 3 Age Groups, 8 Disease Compartments
# Initialize the model with the values of the state variables at the end of the burn-in period
yinit_vac <- c(S = c(results$S1[end_burnin], results$S2[end_burnin], results$S3[end_burnin]),
              I_S = c(results$I_S1[end_burnin], results$I_S2[end_burnin], results$I_S3[end_burnin]), 
              R=c(results$R1[end_burnin], results$R2[end_burnin], results$R3[end_burnin]), 
              C=c(results$C1[end_burnin], results$C2[end_burnin], results$C3[end_burnin]), 
              S_R=c(results$S_R1[end_burnin], results$S_R2[end_burnin], results$S_R3[end_burnin]), 
              I_A=c(results$I_A1[end_burnin], results$I_A2[end_burnin], results$I_A3[end_burnin]),
              V_1=rep(0,3),
              V_2=rep(0,3))

## Run The Model 

results_vacRC <- as.data.frame(ode(y=yinit_vac, times=times, func=typhoid_vaccine_model, parms=pars_vac,
                             method='rk4'))

## calculate total and age-specific incident cases

# Calculate the incidence rate of new infections (aka the "force of infection")
incidence_rate_lambda_vacRC = beta_0*(rowSums(results_vacRC[,c('I_S1','I_S2','I_S3')]) + 
                                  r_a*rowSums(results_vacRC[,c('I_A1','I_A2','I_A3')])+
                                  r_c*rowSums(results_vacRC[,c('C1','C2','C3')]))

# Calculate the number of new cases at each time point overall and for each age group
incidence_total_vacRC = rep*rowSums(results_vacRC[,c('S1','S2','S3')])*incidence_rate_lambda_vacRC
incidence1_vacRC = rep*results_vacRC$S1*incidence_rate_lambda_vacRC
incidence2_vacRC = rep*results_vacRC$S2*incidence_rate_lambda_vacRC
incidence3_vacRC = rep*results_vacRC$S3*incidence_rate_lambda_vacRC


## Plotting Results ##

## State Variable Plots for Each Age Group ##
plot(times, results_vacRC$I_S1, type= "l", col = 'black',ylab="Individuals",
     xlab="Time (Months)", main = "Typhoid in Age Class 1", 
     ylim=c(0,1.01*max(results_vacRC$I_S1,results_vacRC$R1,results_vacRC$C1,results_vacRC$I_A1,results_vacRC$V_11,results_vacRC$V_21)))
lines(times, results_vacRC$R1, col = 'red')
lines(times, results_vacRC$C1, col = 'blue' )
lines(times, results_vacRC$I_A1, col = 'green')
lines(times, results_vacRC$V_11, col = 'orange')
lines(times, results_vacRC$V_21, col = 'magenta' )
legend("topleft",legend=c("Infected (Symptomatic)","Recovered","Carriers","Infected (Asymptomatic)",
                          "Vaccinated", "Vaccinated, Previously Infected"), 
       col=c("black","red","blue","green","orange","magenta"),
       lty=1, ncol=1, cex =.5, pt.cex = 1)

plot(times, results_vacRC$I_S2, type= "l", col = 'black',ylab="Individuals",
     xlab="Time (Months)", main = "Typhoid in Age Class 2", 
     ylim=c(0,1.01*max(results_vacRC$I_S2,results_vacRC$R2,results_vacRC$C2,results_vacRC$I_A2,results_vacRC$V_12,results_vacRC$V_22)))
lines(times, results_vacRC$R2, col = 'red')
lines(times, results_vacRC$C2, col = 'blue' )
lines(times, results_vacRC$I_A2, col = 'green')
lines(times, results_vacRC$V_12, col = 'orange')
lines(times, results_vacRC$V_22, col = 'magenta' )
legend("topleft",legend=c("Infected (Symptomatic)","Recovered","Carriers","Infected (Asymptomatic)",
                          "Vaccinated", "Vaccinated, Previously Infected"), 
       col=c("black","red","blue","green","orange","magenta"),
       lty=1, ncol=1, cex =.5, pt.cex = 1)

plot(times, results_vacRC$I_S3, type= "l", col = 'black',ylab="Individuals",
     xlab="Time (Months)", main = "Typhoid in Age Class 3", 
     ylim=c(0,1.01*max(results_vacRC$I_S3,results_vacRC$R3,results_vacRC$C3,results_vacRC$I_A3,results_vacRC$V_13,results_vacRC$V_23)))
lines(times, results_vacRC$R3, col = 'red')
lines(times, results_vacRC$C3, col = 'blue' )
lines(times, results_vacRC$I_A3, col = 'green')
lines(times, results_vacRC$V_13, col = 'orange')
lines(times, results_vacRC$V_23, col = 'magenta' )
legend("topleft",legend=c("Infected (Symptomatic)","Recovered","Carriers","Infected (Asymptomatic)",
                          "Vaccinated", "Vaccinated, Previously Infected"), 
       col=c("black","red","blue","green","orange","magenta"),
       lty=1, ncol=1, cex =.5, pt.cex = 1)

## Incidence Plots ##
plot(times, incidence_total_RC,type= "l", col='black', ylab="Cases",xlab="Time (Months)",
     main="Incident Cases, Total",
     ylim=c(0,1.01*max(incidence_total_RC)))



############ Question 6 ##############

# Probabilities
pr_inpat = 0.04  # probability of inpatient (hospitalization)
pr_outpat = 0.66 # probability of outpatient treatment
pr_no_treat = 0.30
pr_death_inpat = 0.05
pr_death_outpat = 0.005
pr_death_notreat = 0.01

# Calculate total symptomatic cases over 10 year period for:
# 1) No vaccination
# 2) Routine vaccination at 9 months
# 3) Routine vaccination plus catch up

#### No vaccination ###

#Vaccination coverage
cov_routine = 0  # routine vaccination coverage 
cov_campaign = 0  # catch-up coverage for second age class  

# Save the parameter values
pars_vac <- cbind(mu_0, mu_1, mu_2, mu_3, beta_0, delta, theta, r_c, r_a, w, cov_routine, cov_campaign, ve, w_v)

## Set Time Frame 
start_time = 1 # start date
vac_time = 10*12 # date of vaccine introduction
end_time = 20*12 # end date
times <- seq(start_time, end_time, by = 1) # gives a sequence from start to end in increments of 1

## Run The Model 
results_NoVac <- as.data.frame(ode(y=yinit_vac, times=times, func=typhoid_vaccine_model, parms=pars_vac,
                                   method='rk4'))

## calculate total and age-specific incident cases

incidence_rate_lambda_NoVac = beta_0*(rowSums(results_NoVac[,c('I_S1','I_S2','I_S3')]) + 
                                     r_a*rowSums(results_NoVac[,c('I_A1','I_A2','I_A3')])+
                                     r_c*rowSums(results_NoVac[,c('C1','C2','C3')]))

incidence_total_NoVac = rep*rowSums(results_NoVac[,c('S1','S2','S3')])*incidence_rate_lambda_NoVac
incidence1_NoVac = rep*results_NoVac$S1*incidence_rate_lambda_NoVac
incidence2_NoVac = rep*results_NoVac$S2*incidence_rate_lambda_NoVac
incidence3_NoVac = rep*results_NoVac$S3*incidence_rate_lambda_NoVac

## Calculate the total number of cases, hospitalizations, and deaths
Cases_NoVac = sum(incidence_total_NoVac[vac_time:end_time])
Hospitalizations_NoVac = Cases_NoVac*pr_inpat
Deaths_NoVac = Cases_NoVac*(pr_inpat*pr_death_inpat + pr_outpat*pr_death_outpat + pr_no_treat*pr_death_notreat)
  
#### Routine vaccination ###

#Vaccination coverage
cov_routine = 0.9  # routine vaccination coverage 
cov_campaign = 0  # catch-up coverage for second age class  

# Save the parameter values
pars_vac <- cbind(mu_0, mu_1, mu_2, mu_3, beta_0, delta, theta, r_c, r_a, w, cov_routine, cov_campaign, ve, w_v)

## Set Time Frame 
start_time = 1 # start date
vac_time = 10*12 # date of vaccine introduction
end_time = 20*12 # end date
times <- seq(start_time, end_time, by = 1) # gives a sequence from start to end in increments of 1

## Run The Model 
results_vacR0 <- as.data.frame(ode(y=yinit_vac, times=times, func=typhoid_vaccine_model, parms=pars_vac,
                                   method='rk4'))

## calculate total and age-specific incident cases

incidence_rate_lambda_vacR0 = beta_0*(rowSums(results_vacR0[,c('I_S1','I_S2','I_S3')]) + 
                                        r_a*rowSums(results_vacR0[,c('I_A1','I_A2','I_A3')])+
                                        r_c*rowSums(results_vacR0[,c('C1','C2','C3')]))

incidence_total_vacR0 = rep*rowSums(results_vacR0[,c('S1','S2','S3')])*incidence_rate_lambda_vacR0
incidence1_vacR0 = rep*results_vacR0$S1*incidence_rate_lambda_vacR0
incidence2_vacR0 = rep*results_vacR0$S2*incidence_rate_lambda_vacR0
incidence3_vacR0 = rep*results_vacR0$S3*incidence_rate_lambda_vacR0

## Calculate the total number of cases, hospitalizations, and deaths
Cases_vacR0=sum(incidence_total_vacR0[vac_time:end_time])
Hospitalizations_vacR0 = Cases_vacR0*pr_inpat
Deaths_vacR0 = Cases_vacR0*(pr_inpat*pr_death_inpat + pr_outpat*pr_death_outpat + pr_no_treat*pr_death_notreat)

#### Routine vaccination + catch-up campaign ###

## (See Question 5)

## Calculate the total number of cases, hospitalizations, and deaths
Cases_vacRC=sum(incidence_total_vacRC[vac_time:end_time])
Hospitalizations_vacRC = Cases_vacRC*pr_inpat
Deaths_vacRC = Cases_vacRC*(pr_inpat*pr_death_inpat + pr_outpat*pr_death_outpat + pr_no_treat*pr_death_notreat)

### Routine vaccination results ###

# Cases averted:
CasesAverted_Routine = Cases_NoVac - Cases_vacR0
# Hospitalizations averted:
HospAverted_Routine = Hospitalizations_NoVac - Hospitalizations_vacR0
# Deaths Averted
DeathsAverted_Routine = Deaths_NoVac - Deaths_vacR0

print(CasesAverted_Routine)
print(HospAverted_Routine)
print(DeathsAverted_Routine)

### Routine vaccination + catch-up results ###

# Cases averted:
CasesAverted_Camp = Cases_NoVac - Cases_vacRC
# Hospitalizations averted:
HospAverted_Camp = Hospitalizations_NoVac - Hospitalizations_vacRC
# Deaths Averted
DeathsAverted_Camp = Deaths_NoVac - Deaths_vacRC

print(CasesAverted_Camp)
print(HospAverted_Camp)
print(DeathsAverted_Camp)


############ Question 7 ##############

# Calculate number of people vaccinated in each scenario
# use same runs as in question 4

# Routine: infants are vaccinated as they age from age class 1 to age class 2
# vacc doses = sum of (aging rate from age class 1 to class 2 * vax coverage * (S + S_R + R + C + I_A)),
# with the last two states representing "wasted" doses
cov_routine = .9
Doses_vacR0 = sum(mu_1*cov_routine*results_vacR0[vac_time:end_time,c('S1','S_R1','R1','C1','I_A1')])
print(Doses_vacR0)

# Campaign: For first year, people in the second age class in all states except I_S are
# vaccinated each month from  at rate cov_campaign
# sum of (aging rate from class 1 to class 2 * vax coverage * vax efficacy* age class 1 ( S + S_R + R)) over 10 years
# + sum of vaccination rate * age class 2 ( S + S_R + R)
cov_routine = .9
cov_campaign = .1
Doses_Routine = sum(mu_1*cov_routine*results_vacRC[vac_time:end_time,c('S1','S_R1','R1','C1','I_A1')]) 
Doses_Campaign = sum(cov_campaign*results_vacRC[vac_time:(vac_time+11),c('S2','S_R2','R2','C2','I_A2')])
Doses_vacRC = Doses_Routine + Doses_Campaign
print(Doses_vacRC)


############ Question 8 ##############

# Calculate ICERs for different strategy pairs

# Routine vaccination vs no vaccination
# Calculate Cost of Routine - Cost of No vaccination = net cost
# Calculate DALYs with no vaccination -DALYs with Routine vaccination = DALYs averted (net benefit)

# Costs
cost_vacc_rout = 1.50 + .22 + 1.50  # routine vaccine 
cost_vacc_camp = 1.50 + .22 + .40 # campaign vaccine
cost_in_trt = 100  # inpatient treatment
cost_out_trt = 2.50 # outpatient treatment

# DALY info 
daly_wt = 0.13
dur_inf = 16/365
life_exp = 68.5
life_exp_lost = c(life_exp-.375,life_exp-7.5,life_exp-(68.5/2)) #life expectancy at death for each age group


### Costs of each strategy ###

# Cost of No vacc = 0 vaccine costs + cost of inpatient treatment * number of hospitalizations + cost of outpatient treatment * number outpatient treated
Cost_NoVac = cost_in_trt*Hospitalizations_NoVac + cost_out_trt*pr_outpat*Cases_NoVac

# Cost of Routine = Cost of Routine Vaccination * number of vaccinations + 
    # cost of inpatient treatment * number of hospitalizations + cost of outpatient treatment * number outpatient treated
Cost_vacR0 = cost_vacc_rout*Doses_vacR0 + cost_in_trt*Hospitalizations_vacR0 + cost_out_trt*pr_outpat*Cases_vacR0

# Cost of Routine + Campaign = Cost of Routine Vaccination * number of routine vaccinations + Cost of campaign vacc* number of camp vacc
# + cost of inpatient treatment * number of hospitalizations + cost of outpatient treatment * number outpatient treated
Cost_vacRC = cost_vacc_rout*Doses_Routine + cost_vacc_camp*Doses_Campaign +
  cost_in_trt*Hospitalizations_vacRC + cost_out_trt*pr_outpat*Cases_vacRC

### DALYs for each strategy ###

# DALY = YLD + YLL
# YLD (years lived with disability ) = Number of cases x duration of infection x disability weight
# YLL (years of life lost) = number of deaths*life expectancy at age of death
# # means we need to calculate separately for each age group

# DALYs for No Vacc
DALYs_NoVac=Cases_NoVac*dur_inf*daly_wt + 
  (pr_inpat*pr_death_inpat + pr_outpat*pr_death_outpat+pr_no_treat*pr_death_notreat)*
  (sum(incidence1_NoVac[vac_time:end_time])*(life_exp-.375) +
      sum(incidence2_NoVac[vac_time:end_time])*(life_exp-7.5) +
      sum(incidence3_NoVac[vac_time:end_time])*(life_exp-(68.5/2)))

# DALYs for Routine
DALYs_vacR0=Cases_vacR0*dur_inf*daly_wt + 
  (pr_inpat*pr_death_inpat + pr_outpat*pr_death_outpat+pr_no_treat*pr_death_notreat)*
  (sum(incidence1_vacR0[vac_time:end_time])*(life_exp-.375) +
     sum(incidence2_vacR0[vac_time:end_time])*(life_exp-7.5) +
     sum(incidence3_vacR0[vac_time:end_time])*(life_exp-(68.5/2)))

# DALYs for Routine vaccination + Campaign
DALYs_vacRC=Cases_vacRC*dur_inf*daly_wt + 
  (pr_inpat*pr_death_inpat + pr_outpat*pr_death_outpat+pr_no_treat*pr_death_notreat)*
  (sum(incidence1_vacRC[vac_time:end_time])*(life_exp-.375) +
     sum(incidence2_vacRC[vac_time:end_time])*(life_exp-7.5) +
     sum(incidence3_vacRC[vac_time:end_time])*(life_exp-(68.5/2)))

## part a) Routine vs No Vaccination
Routine_NoVac_ICER = (Cost_vacR0-Cost_NoVac)/(DALYs_NoVac-DALYs_vacR0)

print(Routine_NoVac_ICER)

## part b) Campaign vs No Vaccination
Camp_NoVac_ICER = (Cost_vacRC-Cost_NoVac)/(DALYs_NoVac-DALYs_vacRC)

print(Camp_NoVac_ICER)

## part c) Campaign vs Routine
Camp_Routine_ICER = (Cost_vacRC-Cost_vacR0)/(DALYs_vacR0-DALYs_vacRC)

print(Camp_Routine_ICER)


