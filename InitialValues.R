set.seed(1) 


### FUNCTIONS ####

# function to get year from timepoint
getYear <- function(t){
  return(1 + (t - 1) %/% 12)}

# function to get month of year from timepoint
getMonth <- function(t){
  return(1 + (t - 1) %% 12)}


# function that calculates the log-density of the dirichlet distribution
lddirichlet <- function(x, a){
  wh <- which(a != 0)
  return(lgamma(sum(a[wh])) + sum((a[wh] - 1) * log(x[wh]) - lgamma(a[wh])))
}





### DATA ####

# surveillance and ips data
load("StartingValues.RData")

# divide the number of travellers (represented as n) by 1000 
n <- n/1000 






### PRIOR PARAMETERS ####

lambda_a = 1
lambda_b = 1
mu_alpha = rep(1, 12)
p_alpha = rep(0.8/gtypes, gtypes)
theta_a = 0.5
theta_b = 0.5





### PRIOR PARAMETER: PHI ###

phi_alpha <- matrix(0, dest + 2, dest + 1)

for (i in 1:dest) {
  
  # set the diagonal elements of the matrix
  # set the probability that a traveller who went to destination i was correctly recorded as having gone to destination i.
  phi_alpha[i, i] <- 0.8/3
  
  # set elements in the last two rows of column i 
  # set probability that a traveller who went to destination i was recorded as unknown or domestic respectively
  phi_alpha[dest + 1, i] <- 0.8/3 
  phi_alpha[dest + 2, i] <- 0.8/3  
  
}

# set elements in the last two rows of the last column
# set the probability that a traveller whose destination is unknown or domestic was actually recorded as one of those categories
phi_alpha[(dest + 1):(dest + 2), dest + 1] <- 0.8/2






### INITIAL VALUES ####

lambda_init <-  matrix((1/100000)/(dest + 1), dest + 1, yrs)

mu_init <-  matrix(1/12, dest + 1, 12)

p_init <-  matrix(1/(gtypes), dest + 1, gtypes) 

phi_init <- cbind(diag(rep(0.5, dest)), rep(0, dest))

phi_init <- rbind(phi_init, c(rep(0.4, dest), 0.9))

phi_init <- rbind(phi_init, rep(0.1, dest+1))

theta_init <- 0.9817398















