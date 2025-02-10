# This script initialises key variables, including prior distributions and starting values for MCMC sampling.

# Set random seed for reproducibility
set.seed(1) 

### DEFINE FUNCTIONS ####
# Function to get the year from a given timepoint
getYear <- function(t) {
  return(1 + (t - 1) %/% 12)
}

# Function to get the month from a given timepoint
getMonth <- function(t) {
  return(1 + (t - 1) %% 12)
}

# Function to calculate the log-density of the Dirichlet distribution
lddirichlet <- function(x, a) {
  wh <- which(a != 0)  # Consider only non-zero elements
  return(lgamma(sum(a[wh])) + sum((a[wh] - 1) * log(x[wh]) - lgamma(a[wh])))
}

### LOAD DATA ####
# Load preprocessed surveillance and travel data
load("StartingValues.RData")

# Scale the number of travelers (n) by dividing by 1000 for numerical stability
n <- n / 1000  

### DEFINE PRIORS ####
# Priors for lambda (travel rate)
lambda_a <- 1  
lambda_b <- 1  

# Priors for mu (seasonal variation in travel)
mu_alpha <- rep(1, 12)  

# Prior for p (genetic type proportions)
p_alpha <- rep(0.8 / gtypes, gtypes)  

# Priors for theta (underreporting rate)
theta_a <- 0.5  
theta_b <- 0.5  

### PRIOR FOR PHI (MISCLASSIFICATION) ####
# Initialise prior matrix for misclassification probabilities
phi_alpha <- matrix(0, dest + 2, dest + 1)

# Assign probabilities to diagonal elements (correct reporting)
for (i in 1:dest) {
  phi_alpha[i, i] <- 0.8 / 3  # Probability of correct destination classification
  
  # Assign probabilities for unknown or domestic misclassification
  phi_alpha[dest + 1, i] <- 0.8 / 3  # Probability of being misclassified as "Unknown"
  phi_alpha[dest + 2, i] <- 0.8 / 3  # Probability of being misclassified as "Domestic"
}

# Assign probabilities for travelers originally classified as "Unknown" or "Domestic"
phi_alpha[(dest + 1):(dest + 2), dest + 1] <- 0.8 / 2  

### INITIAL VALUES ####
# Initialise lambda (travel rate per year)
lambda_init <- matrix((1 / 100000) / (dest + 1), dest + 1, yrs)

# Initialise mu (seasonal travel proportions, uniform over months)
mu_init <- matrix(1 / 12, dest + 1, 12)

# Initialise p (genetic type proportions, uniform distribution)
p_init <- matrix(1 / gtypes, dest + 1, gtypes)

# Initialise phi (misclassification matrix)
phi_init <- cbind(diag(rep(0.5, dest)), rep(0, dest))  # Correctly classified
phi_init <- rbind(phi_init, c(rep(0.4, dest), 0.9))  # "Unknown" misclassification
phi_init <- rbind(phi_init, rep(0.1, dest + 1))  # "Domestic" misclassification

# Initialise theta (underreporting probability)
theta_init <- 0.9817398  









