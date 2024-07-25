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

############ Question 4 ##############

## Incidence rate calculation
## rate per 100,000 is calculated as (cases/population)*100,000

# Get annual incidence rate (per 100K people) at end of burn-in period 

