rm(list=ls())                   # Clear all variables and functions

#### Parameter Values
delta <- 0.1   #
R1max <- 30    # maximum resource 1
R2max <- 30    # maximum resource 2
M <-   0.1   #maximum ingestion rate
WC <- 0.0001 #average adult body weight; grams
H <- 3.0 #
q1 <- 0.5 # preference for resource 1
sigma <- 0.5 # conversion efficiency
maint <- 0.01 # maintenance
xi <- 1.0 #maintenance cost of infection; immunity?; should be > 1.0
chi <- 1.0 # reproduction cost of infection; should be < 1.0
mu <- 0.0015 #background mortality
mui <- 0 #disease induced mortality
beta <- 0.01 # transmission rate
wcq = WC^-0.25 #power law with average adult body weight
  
pop.patch <- c(R1 = 1, # Resource 1
              R2 = 1, # Resource 2
              CS = 0, # Susceptible consumer
              CI = 0 # Infected consumer
  )     

values <- c(G1 = with(as.list(pop.patch), delta * (R1max- R1)), #
            G2 = with(as.list(pop.patch), delta * (R2max - R2)),   
            eta1 = with(as.list(pop.patch), (M * wcq* q1 * R1)/ (H + q1 * R1 + (1- q1)* R2)), 
            eta2=  with(as.list(pop.patch), (M * wcq* (1-q1) * R2)/ (H + q1 * R1 + (1- q1)* R2)),
            nus = sigma * (with(as.list(pop.patch), (M * wcq* q1 * R1)/ (H + q1 * R1 + (1- q1)* R2)) + with(as.list(pop.patch), (M * wcq* (1-q1) * R2)/ (H + q1 * R1 + (1- q1)* R2))) - maint*wcq ,
            nui = chi* (sigma * (with(as.list(pop.patch), (M * wcq* q1 * R1)/ (H + q1 * R1 + (1- q1)* R2)) + with(as.list(pop.patch), (M * wcq* (1-q1) * R2)/ (H + q1 * R1 + (1- q1)* R2)))) - xi*maint*wcq,
            I = beta,
            d = mu * wcq,
            di = (mu + mui) * wcq
            )            

patch.onehost <- function(t,y,parms){
  with(c(as.list(y),parms),{
    dR1dt <- G1 - (eta1 * (CS + CI))
    dR2dt <- G2 - (eta2 * (CS + CI))
    dCSdt <- (nus * CS) + (nui * CI) - (I * CS * CI) - (d *CS)
    dCIdt <- (I * CS * CI) - (di *CI)
    list(c(dR1dt,dR2dt,dCSdt,dCIdt))
  })
}


time= 1000
patch.onehost(t=time,y=pop.patch,parms=values)
delta.t <- 0.1                    # Set a small value for delta.t (0.1 day)

library(deSolve)                # Load libary to be used for numerical integration

time.out <- seq(0,1000,by=.1)     # INCLUDE THE APPROPRIATE VALUE TO GET A SEQUENCE
                               # FROM 0 to 365 BY STEPS OF 0.1


lsoda(
  y = pop.patch,               # Initial conditions for population
  times = time.out,             # Timepoints for evaluation
  func = patch.onehost,                   # Function to evaluate
  parms = values                # Vector of parameters
  )

pop.patch.out <- data.frame(lsoda(
  y = pop.patch,               # Initial conditions for population
  times = time.out,             # Timepoints for evaluation
  func = patch.onehost,                   # Function to evaluate
  parms = values                # Vector of parameters
  ))


tail(pop.patch.out)

plot(pop.patch.out$R1)