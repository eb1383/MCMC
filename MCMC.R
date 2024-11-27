### MCMC SIMULATIONS TO FIT A BAYESIAN SPATIOTEMPORAL MODEL TO ESTIMATE THE TRANSMISSION DYNAMICS OF DISEASE 


# load("StartingValues_230613.RData")

set.seed(1) 




### LIBRARIES ### 
library(MCMCpack)
library(mvtnorm)
library(ggplot2)



### FUNCTIONS TO GET TIME ###

# two helper functions that extract the year and month from a given timepoint
# timepoint is represented as a number (1 = January, 2 = February, ...)

# function to get year from timepoint
getYear <- function(t){
  return(1 + (t - 1) %/% 12)}

# function to get month of year from timepoint
getMonth <- function(t){
  return(1 + (t - 1) %% 12)}



### LOG DENSITY FUNCTION OF DIRCIHLET DISTRIBUTION ###
# function that calculates the log-density of the dirichlet distribution
lddirichlet <- function(x, a){
  wh <- which(a != 0)
  return(lgamma(sum(a[wh])) + sum((a[wh] - 1) * log(x[wh]) - lgamma(a[wh])))
}


# load shigella data
source("LoadDataContinent.R")

# load ips data
source("LoadIPSContinent.R")

# save it so dont need to run again
save.image("StartingValues.RData")

# divide the number of travellers (represented as n) by 1000 
# unit of analysis becomes trips per 1000 people
n <- n/1000 



# create n for number of years not months
group_size <- 12

# calculate the number of new columns
num_new_columns <- ncol(n) %/% group_size

# create a new matrix to store the results
N <- matrix(0, nrow = nrow(n), ncol = num_new_columns)

# iterate over each group of columns and calculate the sums
for (i in 1:num_new_columns) {
  start_column <- (i - 1) * group_size + 1
  end_column <- i * group_size
  N[, i] <- rowSums(n[, start_column:end_column])
}



### MCMC PARAMETERS ###

# set number of MCMC iterations
samples = 10


### STARTING VALUES ###

# define the initial values of the model parameters

# lambda: the average rate of infection per year split by known travel + domestic

# 1 /100000 prob of individual being ill
lambda_init <-  matrix((1/100000)/(dest + 1), dest + 1, yrs)

# mu: proportion who travelled in given month 

# start with equal prob each month
mu_init <-  matrix(1/12, dest + 1, 12)

# p: the probability of going to a country and getting disease type j

# start with equal prob of getting each strain
p_init <-  matrix(1/(gtypes), dest + 1, gtypes) 

# phi: the probability of being recorded as going to a certain country
# phi[i,i] - prob went to a country, got recorded as going to that country

# stops starting rate being recorded as egypt when went to madagascar
phi_init <- cbind(diag(rep(0.5, dest)), rep(0, dest))

phi_init <- rbind(phi_init, c(rep(0.4, dest), 0.9))

phi_init <- rbind(phi_init, rep(0.1, dest+1))

# theta: the probability of a SNP being provided given that the individual was typed
theta_init <- 0.9817398




### LIKELIHOOD FUNCTION ####

# the likelihood function of the model calculates the log-likelihood of the observed data 
# given the model parameters

# this code is used to estimate the transmission dynamics of the disease based on genetic data 
# and travel patterns of the infected individuals 


# inputs

# dat: a three-dimensional array of observed data, where the first dimension represents the destination (including unknown and domestic), 
# the second dimension represents the genetic type, and the third dimension represents the time point.

# n: a matrix of the population sizes at each destination and time point.
# lambda: a matrix of the transmission rates at each destination and year.
# mu: a matrix of proportion of travellers at each destination and month.
# p: a matrix of the probabilities of each genetic type at each destination.
# phi: a matrix of the probabilities of travel between destinations.
# theta: the probability that the genetic data is non-missing.

llikelihood <- function(dat, n, lambda, mu, p, phi, theta, tims=1:timepoints, dsts=1:(dest+1)){
  
  # start by initialising the log likelihood to zero
  ll <- 0 
  # iterate over each timepoint t, destinationi i, and genetic type j 
  # to calculate the log-likelihood for the observed data for each combination
  
  for (t in tims){
    
    for (i in dsts){
      
      for (j in 1 : gtypes){
        
        
        # log likelihood contribution for complete data (when the genetic type and destination are both known)
        
        # use dpois to calculate the probability density function of the poisson distribution
        # to calculate the probability of observing the data dat given the expected number of cases 
        # by setting log = TRUE, dpois returns the log of the probability density function which is more numerically stable
        
        
        # the expected number of cases is calculated as the product of the following terms:
        
        # phi[i, i]: the probability of observing the destination i, given the individual travelled and got recorded 
        # n[i, t]: the number of travellers to destination i at the timepoint t
        # lambda[i, getYear(t)]: the transmission rate at destination i and the year of the timepoint t
        # mu[i, getMonth(t)]: the proportion of travellers at destination i in month of timepoint t
        # p[i, j]: the probability of observing genetic type j at destination i 
        # theta: the probability the genetic type is non-missing
        
        # the log-likelihood contribution is then added to the overall log-likelihood ll
        
        
        ll <- ll + dpois(dat[i, j, t], phi[i,i] * n[i, t] * lambda[i, getYear(t)] * mu[i, getMonth(t)] * p[i, j] * theta, log = TRUE)
      }
      
      # log likelihood for incomplete data (when the genetic type is unknown)
      
      j = gtypes + 1
      
      # the probability of the genetic type being missing is set to 1
      # the probability of travel between the destinations is multiplied by the complement of the probability that the genetic data is missing.

      ll <- ll + dpois(dat[i, j, t], phi[i,i] * n[i, t] * lambda[i, getYear(t)] * mu[i, getMonth(t)] * (1 - theta), log = TRUE)
      
    }
    
    # log likelihood for incomplete data (when the destination is unknown)
    
    # iterate over possible genetic types
    # calculating the mean number of cases for each genetic type
    for (j in 1 : gtypes){ 
      
      pois_mean <- 0
      
      k = dest + 2 # all destinations including unknown and domestic
      
      # iterate over the possible destinations including the unknown
      for (i in 1 : (dest + 1)){ 
        
        # poisson mean is calculated as the sum of the means for each possible genetic type 
        # weighted by the probability that a case with that genetic type would be observed at that destination
        
        pois_mean <- pois_mean + phi[k,i] * n[i, t] * lambda[i, getYear(t)] * mu[i, getMonth(t)] * p[i, j] * theta
        
        # phi[k, i] represents the probability of an individual traveling to destination i and being recorded as "unknown" rather than being assigned a specific destination
        # n[i, t] represents the number of individuals present in destination i at time point t
        # lambda[i, getYear(t)]: the annual incidence rate of the disease in location i in the year of time point t
        # mu[i, getMonth(t)]: the monthly variation in incidence rate of the disease in location i at the month of time point t
        # p[i, j] represents the probability of an individual with genetic type j arriving in destination i
        # theta: the probability the genetic type is non-missing
        # all of these terms are multiplied together to give the Poisson mean for the incomplete data.
        
      }
      
      # the contribution to the log likelihood is then calculated
      # this is the log probability of observing the data given the Poisson mean
    
      ll <- ll + dpois(dat[k, j, t], pois_mean, log = TRUE)
    }
    
    # log likelihood for incomplete data (when the destination and genetic type are unknown)
    
    pois_mean <- 0
    
    j = gtypes + 1 # include unknown genetic types
    
    # when both genetic type and destination are unknown 
    # the poisson mean is calculated as the sum of means for each possible combination of genetic type and destination
    # weighted by the probability of observing a case with genetic type at that destination
    # (which is assumed to be uniform across all destination and genetic types)
    
    for (i in 1 : (dest + 1)){
      
      pois_mean <- pois_mean + phi[k,i]*n[i, t]*lambda[i, getYear(t)]*mu[i, getMonth(t)]*(1 - theta)
      
    }
    
    # log likelihood for each time point is updated by adding the log-density of the observed data (assumed to be poisson distributed)
    # given the calculated poisson mean for each combination of destination and genetic type
    
    ll <- ll + dpois(dat[k, j, t], pois_mean, log = TRUE)
  }
  
  return(ll)
  
}





### PHI ###

# tweak phi so cannot be recorded at France when went to Europe

# create a matrix and initalise to zero
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





### PRIORS ####

# initalise a list of priors to specify the prior distributions for the model parameters
priors <- list(lambda_a = 1, lambda_b = 1, mu_alpha = rep(1, 12), p_alpha = rep(0.8/gtypes, gtypes), phi_alpha = phi_alpha, theta_a = 0.5, theta_b = 0.5)

# lambda_a: the shape parameter for the prior distribution of lambda 
# set to 1 to correspond to a uniform prior distribution 

# lambda_b: the scale parameter for the prior distribution of lambda
# set to 1 to correspond to a uniform prior distribution

# mu_alpha: a vector of shape parameters for the prior distribution of mu
# set to 0.8/12 to correspond to a uniform prior distribution 

# p_alpha: a vector of shape parameters for the prior distribution of p
# set to 0.8/gtypes to correspond to a uniform prior distribution 

# phi_alpha: a matrix of shape parameters for the prior distributions of the phi 
# the values in this matrix are set based on the assumption that 80% of the travel histories are recorded accurately 
# and are distributed uniformly across the different categories of travel history.

# theta_a: the shape parameter for the prior distribution of the theta
# set to 0.5 which corresponds to a uniform prior distribution.

# theta_b: the shape parameter for the prior distribution of the theta
# set to 0.5 which corresponds to a uniform prior distribution.


# priors set to be non-informative, the prior beliefs do not favour any particular values for the unknown parameters
# to allow minimal impact on the posterior distribution allowing the data to drive the results of the analysis 


# the prior probabilities are then updated based on the observed data to obtain the posterior probabilities 









### METROPOLIS ALGORITHM ###

# this code defines a function which implements the metropolis hastings algorithm 
# the goal of this algorithm is to generate a sequence of random samples from a probability distribution 
# the algorithm generates a markov chain whose stationary distribution is the target distribution of interest 
# it is designed to ensure that the samples drawn from the chain converge to a stationary distribution

# function inputs
# samples: the number of posterior samples to be generated
# dat: a 3-dimensional array containing the data to be used in the model
# priors: a list of prior distributions for the model parameters
# lambda_init, mu_init, p_init, phi_init, theta_init: initial values for the model parameters

metropolis_algorithm <- function(samples, dat, priors, lambda_init, mu_init, p_init, phi_init, theta_init){ 
  
  # pre-compute useful things
  # create array to store cases by month, gtype and year for each destination
  cases_by_month <- matrix(NA, dest + 1, 12) 
  cases_by_type <- matrix(NA, dest + 1, gtypes) 
  cases_by_year <- matrix(NA, dest + 1, yrs) 
  
  
  # compute summary statistics from the input data
  
  # iterate through each destination and calculate the number of cases by month and by genetic type
  
  for (i in 1 : (dest + 1)) {
    
    cases_by_month[i, ] <- rowSums(matrix(t(dat[i, , ]), 12, (gtypes + 1) * yrs)) 
    
    # dat[i, ] picks out the data corresponding to the ith destination
    # matrix(t(dat[i, , ]), 12) transposes the data and organises it into a matrix with 12 rows (one for each month)
    # and gtypes + 1 * yrs columns (one for each genetic type in each year)
    # rowSums computes the total number of cases in each month for that destination
    
    # calculate total number of cases for each genetic type in that destination, removing the unknown
    cases_by_type[i, ] <- apply(dat[i, 1:gtypes, ], 1, sum) 
  }
  
  # iterate through each year and calculate the number of cases by destination for that year
  
  for (y in 1 : yrs) {
    
    cases_by_year[ , y] <- apply(dat[1:(dest + 1), , 1:12 + (y - 1)*12], 1, sum)
    
    # 1:12 + (y-1)*12 selects the months for that year
    # apply(dat[1:dest....]) calculates the total number of cases for each destination in that year
  }
  
  # priors 
  lambda_a <- priors$lambda_a
  lambda_b <- priors$lambda_b
  mu_alpha <- priors$mu_alpha
  p_alpha <- priors$p_alpha
  phi_alpha <- priors$phi_alpha
  theta_a <- priors$theta_a
  theta_b <- priors$theta_b
  
  
  
  # initialise counters for acceptance/rejection rates of the metropolis hasting algorithm 
  lambda_accept <- matrix(0, dest + 1, yrs)
  lambda_reject <- matrix(0, dest + 1, yrs)
  rw_accept <- rep(0,yrs) # additional MH RW update to improve mixing of lambdas for each year
  rw_reject <- rep(0,yrs)
  mu_accept <- rep(0, dest + 1)
  mu_reject <- rep(0, dest + 1)
  p_accept <- rep(0, dest + 1)  
  p_reject <- rep(0, dest + 1) 
  phi_accept <- rep(0, dest + 1)
  phi_reject <- rep(0, dest + 1)
  theta_accept <- 0
  theta_reject <- 0
  
  
  
  # initial values for the parameters
  lambda <- lambda_init
  mu <- mu_init
  p <- p_init
  phi <- phi_init
  theta <- theta_init
  
  
  # create arrays to store the posterior samples for each parameter
  # arrays have dimensions based on the number of destinations, years and genetic types
  lambda_stored <- array(NA, c(dest + 1, yrs, samples))
  
  mu_stored <- array(NA, c(dest + 1, 12, samples))
  
  p_stored <- array(NA, c(dest + 1, gtypes, samples))
  
  phi_stored <- array(NA, c(dest + 2, dest + 1, samples)) 
  
  theta_stored <- rep(NA, samples)
  
  
  
  ### TUNING PARAMETERS ###
  
  # lambda 
  
  
  iterations_no_diminish <- 5000
  iterations_adaptive <- 15000
  min_beta <- 1
  
  
  
  # initialise tuning parameter with each element set to 1
  lambda_beta <- matrix(1, dest + 1, yrs)
  
  # target acceptance rate
  # algorithm should accept proposed samples with a probability close to 0.44
  lambda_target_acceptance <- 0.23
  
  # define the step size 
  lambda_delta <- 1/lambda_target_acceptance/(1-lambda_target_acceptance) + log(1 + N)
  # 1/lambda_target/(1-lambda_target) determines the scale of the step size 
  # the larger the desired acceptance rate the smaller the step size 
  # log(1+n) adds a log term to increase the step size as n increases
  # this allows the algorithm to explore the posterior distribution more efficiently 
  lambda_mean<-matrix(NA,dest+1,yrs)
  lambda_sigma <- array(NA,c(dest+1,dest+1,yrs))
  for (t in 1:yrs) {
    lambda_sigma[,,t]<-diag(rep(0.1,dest+1))
  }
  lambda_scale <- rep(1,yrs)
  rw_delta<-1/lambda_target_acceptance/(1-lambda_target_acceptance)
  
  
  
  
  # mu
  
  mu_beta <- rep(1, dest + 1)  
  # target acceptance rate (vector so 0.234)
  mu_target_acceptance <- 0.234   
  # define the step size
  mu_delta <- 1/mu_target_acceptance/(1-mu_target_acceptance)
  
  
  
  
  # p
  
  p_beta <- rep(1, dest + 1)
  
  # target acceptance rate (vector so 0.234)
  p_target_acceptance <- 0.234 
  
  # define the step size
  p_delta <- 1/p_target_acceptance/(1-p_target_acceptance)
  
  
  
  
  # phi
  
  phi_beta <- rep(1, dest + 1)
  
  # target acceptance rate (vector so 0.234)
  phi_target_acceptance <- 0.234 
  
  # define the step size
  phi_delta <- 1/phi_target_acceptance/(1-phi_target_acceptance)
  
  
  
  # implement metropolis hastings algorithm to sample from the posterior distribution
  
  # update one element at a time in a nested loop over years and destinations
  
  for (it in 1 : samples) {
    
    # records the current system time to keep track of the algorithm run time
    
    st <- Sys.time() 
    
    
    #### UPDATE LAMBDA ####
    
    # implement metropolis hastings algorithm to sample from the posterior distribution of the parameter lambda
    

    # start a loop over years
    for (y in 1 : yrs) {
      if (it<=iterations_no_diminish) {
        # start a loop over destinations (including domestic)
        for (i in 1 : (dest + 1)) {
          # create a copy of the current lambda parameter for the proposal
          lambda_prop <- lambda
              
          # proposal
          # update the ith destination in year y of the proposal lambda_prop 
          # drawing a random value from a gamma distribution with shape parameter lambda_a + cases_by_year[i, y] + lambda[i, y]*lambda_beta[i, y] 
          # and scale parameter lambda_b + n[i,y] + lambda_beta[i, y]
          lambda_prop[i,y] <- rgamma(1, lambda_a + cases_by_year[i, y] + lambda[i, y]*lambda_beta[i, y], lambda_b + N[i,y] + lambda_beta[i, y])
          
          # where beta is an adaptive constant 
          # bigger beta -> smaller move away from lambda (a and Y less important)
  
          
          
          # log posterior density of the proposal 
          # first term: log prior density of gamma distribution with shape parameter lambda_a and scale lambda_b
          # second term: log-likelihood of the data given by the proposal
          # third term; log prior density of the proposal
          posterior_prop <- dgamma(lambda_prop[i, y], lambda_a, lambda_b, log = TRUE) + llikelihood(dat, n, lambda_prop, mu, p, phi, theta,tims=(y-1)*12 + 1:12,dsts=i) - 
            dgamma(lambda_prop[i, y], lambda_a + cases_by_year[i, y] + lambda[i, y]*lambda_beta[i, y], lambda_b + N[i,y] + lambda_beta[i, y], log = TRUE)        
          
          # log posterior density of the current lambda[i, y]
          # terms same as above but with proposal and current values swapped
          posterior_curr <- dgamma(lambda[i, y], lambda_a, lambda_b, log = TRUE) + llikelihood(dat, n, lambda, mu, p, phi, theta,tims=(y-1)*12 + 1:12,dsts=i) - 
            dgamma(lambda[i, y], lambda_a + cases_by_year[i, y] + lambda_prop[i, y]*lambda_beta[i, y], lambda_b + N[i,y] + lambda_beta[i, y], log = TRUE)       
          
          # generate a random value from a uniform distribution between 0 and 1 to be used in the acceptance probability calculation
          rand_unif <- runif(n = 1)        
          
          # check if the proposal is valid (i.e., not NA) and if the acceptance probability is greater than a random uniform value 
          # if so, the proposal is accepted.       
          if(!is.na(posterior_prop) && posterior_prop - posterior_curr > log(rand_unif)) {          
            
            # increment the number of acceptances
            lambda_accept[i, y] <- lambda_accept[i, y] + 1  
            
            # if the proposal is accepted this line updates the value of lambda with the proposed value stored in lambda prop
            lambda[i,y] <- lambda_prop[i,y]
            
          } else {  
            
            # increment the number of rejections
            lambda_reject[i,y] <- lambda_reject[i, y] + 1            
          }
          
          # calculate the acceptance probability of the proposal distribution
          # if the proposed value has a posterior probability of NA, the acceptance probability is set to 0 
          # otherwise, it is calculated using the ratio of posterior probabilities of the proposed and current values, bounded between 0 and 1
          ap <- ifelse(is.na(posterior_prop), 0, min(1, exp(posterior_prop - posterior_curr))) # exp to get rid of log scale for acceptance prob
          
          # updates the lambda_beta parameter, which is used to adjust the step size of the random walk proposal distribution. 
          # formula uses the difference between the target acceptance probability (0.44) and the actual acceptance probability ap
          # and scales it by a factor that depends on the step size lambda_delta, the number of iterations it, and the number of times the parameter i has been updated
          if (it<=iterations_adaptive) {
            lambda_beta[i, y]<-lambda_beta[i, y]*exp(lambda_delta[i,y]/max(1,it-iterations_no_diminish)*(lambda_target_acceptance-ap))
          }
          
        }
      }
        
      # Second update to improve moves over lambda
      lambda_prop <-lambda
      lambda_prop[,y]<-rmvnorm(1,lambda_prop[,y],lambda_scale[y]*2.38^2/(dest+1)*lambda_sigma[,,y])
      if (min(lambda_prop[,y])>0) { 
        posterior_prop <- sum(dgamma(lambda_prop[, y], lambda_a, lambda_b, log = TRUE)) + llikelihood(dat, n, lambda_prop, mu, p, phi, theta,tims=(y-1)*12 + 1:12)        
      
        # log posterior density of the current lambda[i, y]
        # terms same as above but with proposal and current values swapped
        posterior_curr <- sum(dgamma(lambda[, y], lambda_a, lambda_b, log = TRUE)) + llikelihood(dat, n, lambda, mu, p, phi, theta,tims=(y-1)*12 + 1:12)      
      
        # generate a random value from a uniform distribution between 0 and 1 to be used in the acceptance probability calculation
        rand_unif <- runif(n = 1)        
      
        # check if the proposal is valid (i.e., not NA) and if the acceptance probability is greater than a random uniform value 
        # if so, the proposal is accepted.       
        if(!is.na(posterior_prop) && posterior_prop - posterior_curr > log(rand_unif)) {          
        
          # increment the number of acceptances
          rw_accept[y] <- rw_accept[y] + 1         
          # if the proposal is accepted this line updates the value of lambda with the proposed value stored in lambda prop
          lambda[,y] <- lambda_prop[,y]          
        } else {            
          # increment the number of rejections
          rw_reject[y] <- rw_reject[y] + 1            
        }
        ap<-min(1,exp(posterior_prop - posterior_curr))
      } else { # direct rejection
        # increment the number of rejections
        rw_reject[y] <- rw_reject[y] + 1
        ap<-0 
      }
      if (it<=iterations_adaptive) {
        if (is.na(ap)) {ap<-0}
        lambda_scale[y]<-lambda_scale[y]*exp(rw_delta/max(1,it-iterations_no_diminish)*(ap-lambda_target_acceptance))
        if (it==1) {
          lambda_mean[,y]<-(lambda[,y]+lambda_init[,y])/2
          lambda_sigma[,,y]<-(tcrossprod(lambda_init[,y])+tcrossprod(lambda[,y])-2*tcrossprod(lambda_mean[,y])+(dest+3)*lambda_sigma[,,y])/(dest+5)
        } else {
          new_mean<-it/(it+1)*lambda_mean[,y]+lambda[,y]/(it+1)
          lambda_sigma[,,y]<-((it+dest+3)*lambda_sigma[,,y]+tcrossprod(lambda[,y])+it*tcrossprod(lambda_mean[,y])-(it+1)*tcrossprod(new_mean))/(it+dest+4)
          lambda_mean[,y]<-new_mean
        }
      }
    }
    # print message that shows the minimum acceptance rate for lambda and gives the index of the lambda parameter with the lowest number of acceptances   
    cat("  Lambda min acceptance rate", min(lambda_accept/(lambda_accept + lambda_reject)), which.min(lambda_accept), "\n")
    cat("  RW min acceptance rate", min(rw_accept/(rw_accept + rw_reject)), which.min(rw_accept), "\n")      
    
    #### UPDATE MU ####
    
    # implement metropolis hastings algorithm to sample from the posterior distribution of the parameter mu
    
    # start a loop over destinations
    for (i in 1 : (dest + 1)){ 
      
      # create a copy of the current mu parameter for the proposal
      mu_prop <- mu
      
      # proposal
      
      # propose a new value for the ith element of the mu vector by sampling from a dircihlet distribution with updated hyperparameters
      # hyperparameters are calculated as the sum of the prior hyperparameters mu_alpha, the observed counts cases_by_month 
      # and a tuning parameter mu_beta*mu which controls the acceptance rate of the proposed values
      mu_prop[i,] <- c(rdirichlet(1, mu_alpha + cases_by_month[i, ] + mu_beta[i] * mu[i,]))
      
      # log posterior probability of the proposed value mu_prop
      # first term: log dirichlet prior probability given the prior hyperparameters mu_alpha
      # second term: the log likelihood of the data given the current values of the parameters 
      # third term: the log dircihlet prior probability of mu_prop[i, ]
      posterior_prop <- lddirichlet(mu_prop[i, ], mu_alpha) + llikelihood(dat, n, lambda, mu_prop, p, phi, theta,dsts=i) - 
        lddirichlet(mu_prop[i, ], mu_alpha + cases_by_month[i, ] + mu_beta[i] * mu[i, ])
      
      # log posterior probability of the current value mu[i, ]
      # same method as above but with the current value of mu instead of the proposed value
      posterior_curr <- lddirichlet(mu[i,], mu_alpha) + llikelihood(dat, n, lambda, mu, p, phi, theta,dsts=i) - 
        lddirichlet(mu[i,], mu_alpha + cases_by_month[i,] + mu_beta[i] * mu_prop[i,])
      
      # generate a random value from a uniform distribution between 0 and 1 to be used in the acceptance probability calculation
      rand_unif <- runif(n = 1) 
      
      # check if the proposal is valid (i.e., not NA) and if the acceptance probability is greater than a random uniform value 
      # if so, the proposal is accepted.
      if(!is.na(posterior_prop) && posterior_prop - posterior_curr > log(rand_unif)) { 
        
        # increment the number of acceptances
        mu_accept[i] <- mu_accept[i] + 1   
        
        # if the proposal is accepted this line updates the value of lambda with the proposed value stored in lambda prop
        mu[i,] <- mu_prop[i,]
            
      } else { 
        
        # increment the number of rejections
        mu_reject[i] <- mu_reject[i] + 1
        
      
      }
      
      # calculate the acceptance probability of the proposal distribution
      # calculated using the ratio of posterior probabilities of the proposed and current values
      # bounded between 0 and 1
      if (it<=iterations_adaptive) {
        ap <- min(1, exp(posterior_prop - posterior_curr))
        
        if (is.na(ap)) {ap<-0}
        
        # reset number of iterations
        mu_beta[i] <- max(min_beta,mu_beta[i]*exp(mu_delta/max(1, it - iterations_no_diminish) * (mu_target_acceptance - ap)))
      }
      ### test code
      if (it%%300==0 && i==which.min(mu_accept)) {        
        cols=rainbow(12)
        M<-max(mu_stored[i,,1:(it-1)])
        plot(1:it,t="n",ylim=c(0,M),main=paste("mu",i,round(c(posterior_prop,posterior_curr),3)),xlab=paste("beta =",mu_beta[i]),
            ylab=paste("acceptance =",mu_accept[i]/(mu_accept[i]+mu_reject[i])))
        for (j in 1:12) {
          lines(1:(it-1),mu_stored[i,j,1:(it-1)],col=cols[j])
          points(rep(it,12),mu[i,],col=cols)
          points(rep(it+1,12),mu_prop[i,],col=cols,pch=2)
        }
      }
      #cat("    ",i,mu_beta[i],round(mu[i,],3),"\n")
    }
    
    # print message that shows the minimum acceptance rate for mu and gives the index of the mu parameter with the lowest number of acceptances
    cat("  Mu min acceptance rate", min(mu_accept/(mu_accept + mu_reject)), which.min(mu_accept), "\n")
    #if (it<=iterations_adaptive) {
    #  cat("    Mu max beta", max(mu_beta),which.max(mu_beta),"\n")
    #  cat("    Mu min beta", min(mu_beta),which.min(mu_beta),"\n")
    #}    
       
       
    #### UPDATE P ####
    
    # implement metropolis hastings algorithm to sample from the posterior distribution of the parameter p
    
    # start a loop over destinations
    for (i in 1:(dest + 1)) {
      
      # create a copy of the current p parameter for the proposal
      p_prop <- p
      
      # proposal

      # update p using an adaptive Dirichlet random walk proposal with mean proportional to the current value of p and variance that is scaled based on a tuning parameter p_beta
      # where tuning parameter p_beta is adjusted after each iteration to achieve a target acceptance rate (p_target_acceptance)
      p_prop[i,] <- c(rdirichlet(1, p_alpha + cases_by_type[i, ] + p_beta[i] * p[i, ]))
      
      # log posterior density of the proposal
      # first term: log of the dircihlet prior density for the proposed p_prop using p_alpha
      # second term: log of the poisson likelihood of the data
      # third term: log of the dirichlet prior density for the proposed p_prop
      posterior_prop <- lddirichlet(p_prop[i, ], p_alpha) + llikelihood(dat, n, lambda, mu, p_prop, phi, theta,dsts=i) -
        lddirichlet(p_prop[i, ], p_alpha + cases_by_type[i, ] + p_beta[i] * p[i, ])
      
      # log posterior density of the current p
      # terms same as above but with proposal and current values swapped
      posterior_curr <- lddirichlet(p[i, ], p_alpha) + llikelihood(dat, n, lambda, mu, p, phi, theta,dsts=i) - 
        lddirichlet(p[i, ], p_alpha + cases_by_type[i, ] + p_beta[i] * p_prop[i, ])
      
      # generate a random value from a uniform distribution between 0 and 1 to be used in the acceptance probability calculation
      rand_unif <- runif(n = 1) 
      
      # check if the proposal is valid (i.e., not NA) and if the acceptance probability is greater than a random uniform value 
      # if so, the proposal is accepted.
      if(!is.na(posterior_prop) && posterior_prop - posterior_curr > log(rand_unif)) { 
        
        # increment the number of acceptances
        p_accept[i] <- p_accept[i] + 1
        
        # if the proposal is accepted this line updates the value of p with the proposed value stored in p prop
        p[i,] <- p_prop[i, ]
        
        
      } else {   
        
        # increment the number of rejections
        p_reject[i] <- p_reject[i] + 1  
        
      }
      
      if (it<=iterations_adaptive) {
        # calculate the acceptance probability of the proposal distribution
        # calculated using the ratio of posterior probabilities of the proposed and current values
        # bounded between 0 and 1
        
        ap <- min(1, exp(posterior_prop - posterior_curr))
        if (is.na(ap)) {ap<-0}
        # do 200 iterations without diminishing step size
        # divides by 1 for first 201 iterations then divides by 2, 3...
        
        p_beta[i] <- max(min_beta,p_beta[i]*exp(p_delta/max(1, it - iterations_no_diminish) * (p_target_acceptance - ap)))
      }
      ### test code
      if (it%%300==100 && i==which.min(p_accept)) {
        cols=rainbow(gtypes)
        M<-max(p_stored[i,,1:(it-1)])
        plot(1:it,t="n",ylim=c(0,M),main=paste("p",i,round(c(posterior_prop,posterior_curr),3)),xlab=paste("beta =",p_beta[i]),
              ylab=paste("acceptance =",p_accept[i]/(p_accept[i]+p_reject[i])))
        for (j in 1:gtypes) {
          lines(1:(it-1),p_stored[i,j,1:(it-1)],col=cols[j])
          points(rep(it,gtypes),p[i,],col=cols)
          points(rep(it+1,gtypes),p_prop[i,],col=cols,pch=2)
        }
      }   
      
    }
    
    # print message that shows the minimum acceptance rate for p and gives the index of the p parameter with the lowest number of acceptances
    
    cat("  p min acceptance rate", min(p_accept/(p_accept + p_reject)), which.min(p_accept), "\n")
    
    
    
    
    
    #### UPDATE PHI #### 
    
    # start a loop over the destinations
    for (i in 1:(dest + 1)) {
      
      # create a copy of the current phi parameter for the proposal
      phi_prop <- phi
      
      # proposal
      
      # update phi using an adaptive Dirichlet random walk proposal with mean proportional to the current value of phi and variance that is scaled based on a tuning parameter phi_beta
      phi_prop[ ,i] <- c(rdirichlet(1, phi_alpha[ ,i] + phi_beta[i] * phi[ ,i]))
      
      
  
      # log posterior density of the proposal
      # first term: log of the dircihlet prior density for the proposed phi_prop using phi_alpha
      # second term: log of the poisson likelihood of the data
      # third term: log of the dirichlet prior density for the proposed phi_prop
      posterior_prop <- lddirichlet(phi_prop[ ,i], phi_alpha[ ,i]) + llikelihood(dat, n, lambda, mu, p, phi_prop, theta,dsts=i) - 
        lddirichlet(phi_prop[ ,i], phi_alpha[ ,i]  + phi_beta[i] * phi[ ,i])
      
      # log posterior density of the current phi
      # terms same as above but with proposal and current values swapped
      posterior_curr <- lddirichlet(phi[ ,i], phi_alpha[ ,i]) + llikelihood(dat, n, lambda, mu, p, phi, theta,dsts=i) -
        lddirichlet(phi[ ,i], phi_alpha[ ,i]  + phi_beta[i] * phi_prop[ ,i])
      
      # generate a random value from a uniform distribution between 0 and 1 to be used in the acceptance probability calculation
      rand_unif <- runif(n = 1) 
      
      # check if the proposal is valid (i.e., not NA) and if the acceptance probability is greater than a random uniform value 
      # if so, the proposal is accepted    
      if(!is.na(posterior_prop) && posterior_prop - posterior_curr > log(rand_unif)) { 
        
      
        # increment the number of acceptances 
        phi_accept[i] <- phi_accept[i] + 1
        
        phi[ ,i] <- phi_prop[ ,i]
        

      } else {   
        
        # increment the number of rejections
        phi_reject[i] <- phi_reject[i] + 1 
        
      }
      if (it<=iterations_adaptive) {
        # calculate the acceptance probability of the proposal distribution
        # calculated using the ratio of posterior probabilities of the proposed and current values
        # bounded between 0 and 1
        
        ap <- min(1, exp(posterior_prop - posterior_curr))
        if (is.na(ap)) {ap<-0}
        
        phi_beta[i] <- max(min_beta,phi_beta[i]*exp(phi_delta/max(1, it - iterations_no_diminish) * (phi_target_acceptance - ap)))
      }
      ### test code
      if (it%%300==200 && i==which.min(phi_accept)) {        
        cols=rainbow(dest+2)
        M<-max(phi_stored[,i,1:(it-1)])
        plot(1:it,t="n",ylim=c(0,M),main=paste("phi",i,round(c(posterior_prop,posterior_curr),3)),xlab=paste("beta =",phi_beta[i]),
             ylab=paste("acceptance =",phi_accept[i]/(phi_accept[i]+phi_reject[i])))
        for (j in 1:(dest+2)) {
          lines(1:(it-1),phi_stored[j,i,1:(it-1)],col=cols[j])
          points(rep(it,dest+2),phi[,i],col=cols)
          points(rep(it+1,dest+2),phi_prop[,i],col=cols,pch=2)
        }
      }
    }

    # print message that shows the minimum acceptance rate for phi and gives the index of the phi parameter with the lowest number of acceptances
    cat("  phi min acceptance rate", min(phi_accept/(phi_accept + phi_reject)), which.min(phi_accept), "\n")
    
    

  
     
    
    
    #### UPDATE THETA ####
    
    # update the scale parameter theta
    theta <- rbeta(1, theta_a + sum(dat[ , 1:gtypes, ]), theta_b + sum(dat[ , gtypes + 1, ]))  # 1, success, failure
    
    # check if the current state of the model has zero posterior probability
    # if so print a warning message indicating that the starting state of the model is invalid
    if(posterior_curr == -Inf){
      
      cat("Starting state has zero probability \n")
      
    }
    
    
    
    # store everything
    lambda_stored[ , ,it] <- lambda
    mu_stored[ , , it] <- mu
    p_stored[ , , it] <- p
    phi_stored[ , , it] <- phi
    theta_stored[it] <- theta
    
    # print current iteration number and the time taken to run that iteration
    cat("Iteration", it, difftime(Sys.time(), st, units = "secs"), "seconds.\n" ) 
    
  }
  
  betas <- list(lambda_beta = lambda_beta, mu_beta = mu_beta, phi_beta = phi_beta, p_beta = p_beta)
  
  acceptances <- list(lambda_accept = lambda_accept, mu_accept = mu_accept, p_accept = p_accept, phi_accept = phi_accept)
  
  rejections <- list(lambda_reject = lambda_reject, mu_reject = mu_reject, p_reject = p_reject, phi_reject = phi_reject)
  
  return(list(lambda_stored = lambda_stored, mu_stored = mu_stored, p_stored = p_stored, phi_stored = phi_stored, theta_stored = theta_stored, betas = betas, acceptances = acceptances, rejections = rejections, cases_by_month = cases_by_month, cases_by_year = cases_by_year))
  
}

posterior_list <- metropolis_algorithm(samples, dat, priors, lambda_init, mu_init, p_init, phi_init, theta_init)


# save it so dont need to run again
save.image("Output.RData")
