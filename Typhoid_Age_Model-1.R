library(deSolve)

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

## Additional demographic parameters (times in months)
pop_dist <- c(.1,.3,.6)  #initial population age distribution (this will also equilibrate during the burn-in period)

mu_0 = 0.0015  # Birth rate
mu_1 = 1/9   # Rate at which people age out of first age class/into second age class
mu_2 = 1/(14.25*12) # Rate at which people age out of second age class/into third age class
mu_3 = mu_0   # set death rate equal to birth rate to keep population constant

# Save the parameter values
pars_age <- cbind(mu_0, mu_1, mu_2, mu_3, beta_0, delta, theta, r_c, r_a, w)

## Set the Time Frame 

start_burnin = 0 # start date
end_burnin = 200*12 #NOTE: time is in months. Set to 100 for burn-in, and extend to run longer
times_burnin <- seq(start_burnin, end_burnin, by = 1) # gives a sequence from start to end in increments of 1

## Initial Conditions 
init_prev <- .0000001   #initial prevalence (assume epidemic begins with very low prevalence)

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
       lty=1, ncol=1, cex =.9, pt.cex = 1)

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

############ Question 5 ##############

## Incidence rate calculation
## rate per 100,000 is calculated as (cases/population)*100,000

# Get incidence at end of burn-in period (annual incidence ~ 12 months/year * monthly incidence)

# Average age of infection: age-group midpoint age * weight, weight = incident cases in age group/total incident cases

incidence = 12*sum(incidence_total[end_burnin])/pop0*1e5

Avg_AoI = (0.375*incidence1[end_burnin]+7.5*incidence2[end_burnin]+(68.5/2)*incidence3[end_burnin])/(incidence_total[end_burnin])

print(incidence)
print(Avg_AoI)
