set.seed(1) 

### LIBRARIES ### 
library(MCMCpack)
library(mvtnorm)
library(ggplot2)

### FUNCTIONS TO GET TIME ###
# function to get year from timepoint
getYear <- function(t){
  return(1 + (t - 1) %/% 12)}
# function to get month of year from timepoint
getMonth <- function(t){
  return(1 + (t - 1) %% 12)}

### LOG DENSITY FUNCTION OF DIRCIHLET DISTRIBUTION ###
lddirichlet <- function(x, a){
  wh <- which(a != 0)
  return(lgamma(sum(a[wh])) + sum((a[wh] - 1) * log(x[wh]) - lgamma(a[wh])))
}

# load data
source("")
source("")
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
phi_init <- cbind(diag(rep(0.5, dest)), rep(0, dest))
phi_init <- rbind(phi_init, c(rep(0.4, dest), 0.9))
phi_init <- rbind(phi_init, rep(0.1, dest+1))
# theta: the probability of a SNP being provided given that the individual was typed
theta_init <- 0.9817398

### LIKELIHOOD FUNCTION ####
llikelihood <- function(dat, n, lambda, mu, p, phi, theta, tims=1:timepoints, dsts=1:(dest+1)){
  # start by initialising the log likelihood to zero
  ll <- 0 
  # iterate over each timepoint t, destinationi i, and genetic type j 
  for (t in tims){
    for (i in dsts){
      for (j in 1 : gtypes){
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
      ll <- ll + dpois(dat[i, j, t], phi[i,i] * n[i, t] * lambda[i, getYear(t)] * mu[i, getMonth(t)] * (1 - theta), log = TRUE)
    }
    # log likelihood for incomplete data (when the destination is unknown)
    # iterate over possible genetic types, calculating the mean number of cases for each genetic type
    for (j in 1 : gtypes){ 
      pois_mean <- 0
      k = dest + 2 # all destinations including unknown and domestic
      for (i in 1 : (dest + 1)){ 
        # poisson mean is calculated as the sum of the means for each possible genetic type 
        # weighted by the probability that a case with that genetic type would be observed at that destination
        pois_mean <- pois_mean + phi[k,i] * n[i, t] * lambda[i, getYear(t)] * mu[i, getMonth(t)] * p[i, j] * theta
      }
      ll <- ll + dpois(dat[k, j, t], pois_mean, log = TRUE)
    }
    # log likelihood for incomplete data (when the destination and genetic type are unknown)
    pois_mean <- 0
    j = gtypes + 1 # include unknown genetic types
    for (i in 1 : (dest + 1)){
      pois_mean <- pois_mean + phi[k,i]*n[i, t]*lambda[i, getYear(t)]*mu[i, getMonth(t)]*(1 - theta) 
    }
    ll <- ll + dpois(dat[k, j, t], pois_mean, log = TRUE)
  } 
  return(ll)
}


### PHI ###
# create a matrix and initalise to zero
phi_alpha <- matrix(0, dest + 2, dest + 1)
for (i in 1:dest) {
  # set the probability that a traveller who went to destination i was correctly recorded as having gone to destination i.
  phi_alpha[i, i] <- 0.8/3
  # set probability that a traveller who went to destination i was recorded as unknown or domestic respectively
  phi_alpha[dest + 1, i] <- 0.8/3 
  phi_alpha[dest + 2, i] <- 0.8/3  
}
# set the probability that a traveller whose destination is unknown or domestic was actually recorded as one of those categories
phi_alpha[(dest + 1):(dest + 2), dest + 1] <- 0.8/2

### PRIORS ####
# initalise a list of priors to specify the prior distributions for the model parameters
priors <- list(lambda_a = 1, lambda_b = 1, mu_alpha = rep(1, 12), p_alpha = rep(0.8/gtypes, gtypes), phi_alpha = phi_alpha, theta_a = 0.5, theta_b = 0.5)


### METROPOLIS ALGORITHM ###
metropolis_algorithm <- function(samples, dat, priors, lambda_init, mu_init, p_init, phi_init, theta_init){ 
  # create array to store cases by month, gtype and year for each destination
  cases_by_month <- matrix(NA, dest + 1, 12) 
  cases_by_type <- matrix(NA, dest + 1, gtypes) 
  cases_by_year <- matrix(NA, dest + 1, yrs) 
  # iterate through each destination and calculate the number of cases by month and by genetic type
  for (i in 1 : (dest + 1)) {
    cases_by_month[i, ] <- rowSums(matrix(t(dat[i, , ]), 12, (gtypes + 1) * yrs)) 
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
  lambda_stored <- array(NA, c(dest + 1, yrs, samples))
  mu_stored <- array(NA, c(dest + 1, 12, samples))
  p_stored <- array(NA, c(dest + 1, gtypes, samples))
  phi_stored <- array(NA, c(dest + 2, dest + 1, samples)) 
  theta_stored <- rep(NA, samples)

  ### TUNING PARAMETERS ###
  iterations_no_diminish <- 5000
  iterations_adaptive <- 15000
  min_beta <- 1
  # initialise tuning parameter with each element set to 1
  lambda_beta <- matrix(1, dest + 1, yrs)
  # target acceptance rate
  lambda_target_acceptance <- 0.23
  # define the step size 
  lambda_delta <- 1/lambda_target_acceptance/(1-lambda_target_acceptance) + log(1 + N)
  # 1/lambda_target/(1-lambda_target) determines the scale of the step size 
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
  mu_target_acceptance <- 0.234   
  mu_delta <- 1/mu_target_acceptance/(1-mu_target_acceptance)
 
  # p
  p_beta <- rep(1, dest + 1)
  p_target_acceptance <- 0.234 
  p_delta <- 1/p_target_acceptance/(1-p_target_acceptance)

  # phi
  phi_beta <- rep(1, dest + 1)
  phi_target_acceptance <- 0.234 
  phi_delta <- 1/phi_target_acceptance/(1-phi_target_acceptance)
  
  # implement metropolis hastings algorithm to sample from the posterior distribution
  for (it in 1 : samples) {
    # records the current system time to keep track of the algorithm run time
    st <- Sys.time() 
    
    #### UPDATE LAMBDA ####
    # start a loop over years
    for (y in 1 : yrs) {
      if (it<=iterations_no_diminish) {
        # start a loop over destinations (including domestic)
        for (i in 1 : (dest + 1)) {
          lambda_prop <- lambda
          # update the ith destination in year y of the proposal lambda_prop 
          # drawing a random value from a gamma distribution with shape parameter lambda_a + cases_by_year[i, y] + lambda[i, y]*lambda_beta[i, y] 
          lambda_prop[i,y] <- rgamma(1, lambda_a + cases_by_year[i, y] + lambda[i, y]*lambda_beta[i, y], lambda_b + N[i,y] + lambda_beta[i, y])
          # beta is an adaptive constant 
          # bigger beta -> smaller move away from lambda (a and Y less important)
          posterior_prop <- dgamma(lambda_prop[i, y], lambda_a, lambda_b, log = TRUE) + llikelihood(dat, n, lambda_prop, mu, p, phi, theta,tims=(y-1)*12 + 1:12,dsts=i) - 
            dgamma(lambda_prop[i, y], lambda_a + cases_by_year[i, y] + lambda[i, y]*lambda_beta[i, y], lambda_b + N[i,y] + lambda_beta[i, y], log = TRUE)        
          posterior_curr <- dgamma(lambda[i, y], lambda_a, lambda_b, log = TRUE) + llikelihood(dat, n, lambda, mu, p, phi, theta,tims=(y-1)*12 + 1:12,dsts=i) - 
            dgamma(lambda[i, y], lambda_a + cases_by_year[i, y] + lambda_prop[i, y]*lambda_beta[i, y], lambda_b + N[i,y] + lambda_beta[i, y], log = TRUE)       
          # generate a random value from a uniform distribution between 0 and 1 to be used in the acceptance probability calculation
          rand_unif <- runif(n = 1)        
          # check if the proposal is valid (i.e., not NA) and if the acceptance probability is greater than a random uniform value 
          # if so, the proposal is accepted.       
          if(!is.na(posterior_prop) && posterior_prop - posterior_curr > log(rand_unif)) {          
            # increment the number of acceptances
            lambda_accept[i, y] <- lambda_accept[i, y] + 1  
            lambda[i,y] <- lambda_prop[i,y]
          } else {  
            # increment the number of rejections
            lambda_reject[i,y] <- lambda_reject[i, y] + 1            
          }
          # calculate the acceptance probability of the proposal distribution
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
        posterior_curr <- sum(dgamma(lambda[, y], lambda_a, lambda_b, log = TRUE)) + llikelihood(dat, n, lambda, mu, p, phi, theta,tims=(y-1)*12 + 1:12)      
        rand_unif <- runif(n = 1)        
        if(!is.na(posterior_prop) && posterior_prop - posterior_curr > log(rand_unif)) {          
          rw_accept[y] <- rw_accept[y] + 1         
          lambda[,y] <- lambda_prop[,y]          
        } else {            
          rw_reject[y] <- rw_reject[y] + 1            
        }
        ap<-min(1,exp(posterior_prop - posterior_curr))
      } else { # direct rejection
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
    # start a loop over destinations
    for (i in 1 : (dest + 1)){ 
      mu_prop <- mu
      mu_prop[i,] <- c(rdirichlet(1, mu_alpha + cases_by_month[i, ] + mu_beta[i] * mu[i,]))
      posterior_prop <- lddirichlet(mu_prop[i, ], mu_alpha) + llikelihood(dat, n, lambda, mu_prop, p, phi, theta,dsts=i) - 
        lddirichlet(mu_prop[i, ], mu_alpha + cases_by_month[i, ] + mu_beta[i] * mu[i, ])
      posterior_curr <- lddirichlet(mu[i,], mu_alpha) + llikelihood(dat, n, lambda, mu, p, phi, theta,dsts=i) - 
        lddirichlet(mu[i,], mu_alpha + cases_by_month[i,] + mu_beta[i] * mu_prop[i,])
      rand_unif <- runif(n = 1) 
      if(!is.na(posterior_prop) && posterior_prop - posterior_curr > log(rand_unif)) { 
        mu_accept[i] <- mu_accept[i] + 1   
        mu[i,] <- mu_prop[i,]   
      } else {  
        mu_reject[i] <- mu_reject[i] + 1
      }  
      if (it<=iterations_adaptive) {
        ap <- min(1, exp(posterior_prop - posterior_curr))
        if (is.na(ap)) {ap<-0}
        mu_beta[i] <- max(min_beta,mu_beta[i]*exp(mu_delta/max(1, it - iterations_no_diminish) * (mu_target_acceptance - ap)))
      }
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
    }
    cat("  Mu min acceptance rate", min(mu_accept/(mu_accept + mu_reject)), which.min(mu_accept), "\n")
    
    #### UPDATE P ####
    # loop over destinations
    for (i in 1:(dest + 1)) {
      p_prop <- p
      # update p using an adaptive Dirichlet random walk proposal with mean proportional to the current value of p and variance that is scaled based on a tuning parameter p_beta
      # where tuning parameter p_beta is adjusted after each iteration to achieve a target acceptance rate (p_target_acceptance)
      p_prop[i,] <- c(rdirichlet(1, p_alpha + cases_by_type[i, ] + p_beta[i] * p[i, ]))
      posterior_prop <- lddirichlet(p_prop[i, ], p_alpha) + llikelihood(dat, n, lambda, mu, p_prop, phi, theta,dsts=i) -
        lddirichlet(p_prop[i, ], p_alpha + cases_by_type[i, ] + p_beta[i] * p[i, ])
      posterior_curr <- lddirichlet(p[i, ], p_alpha) + llikelihood(dat, n, lambda, mu, p, phi, theta,dsts=i) - 
        lddirichlet(p[i, ], p_alpha + cases_by_type[i, ] + p_beta[i] * p_prop[i, ])
      rand_unif <- runif(n = 1) 
      if(!is.na(posterior_prop) && posterior_prop - posterior_curr > log(rand_unif)) { 
        p_accept[i] <- p_accept[i] + 1
        p[i,] <- p_prop[i, ]
      } else {   
        p_reject[i] <- p_reject[i] + 1  
      }
      if (it<=iterations_adaptive) {
        ap <- min(1, exp(posterior_prop - posterior_curr))
        if (is.na(ap)) {ap<-0}
        # do 200 iterations without diminishing step size
        # divides by 1 for first 201 iterations then divides by 2, 3...
        p_beta[i] <- max(min_beta,p_beta[i]*exp(p_delta/max(1, it - iterations_no_diminish) * (p_target_acceptance - ap)))
      }
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
    cat("  p min acceptance rate", min(p_accept/(p_accept + p_reject)), which.min(p_accept), "\n")
   
    #### UPDATE PHI #### 
    # loop over the destinations
    for (i in 1:(dest + 1)) {
      phi_prop <- phi
      phi_prop[ ,i] <- c(rdirichlet(1, phi_alpha[ ,i] + phi_beta[i] * phi[ ,i]))
      posterior_prop <- lddirichlet(phi_prop[ ,i], phi_alpha[ ,i]) + llikelihood(dat, n, lambda, mu, p, phi_prop, theta,dsts=i) - 
        lddirichlet(phi_prop[ ,i], phi_alpha[ ,i]  + phi_beta[i] * phi[ ,i])
      posterior_curr <- lddirichlet(phi[ ,i], phi_alpha[ ,i]) + llikelihood(dat, n, lambda, mu, p, phi, theta,dsts=i) -
        lddirichlet(phi[ ,i], phi_alpha[ ,i]  + phi_beta[i] * phi_prop[ ,i])
      rand_unif <- runif(n = 1)  
      if(!is.na(posterior_prop) && posterior_prop - posterior_curr > log(rand_unif)) { 
        phi_accept[i] <- phi_accept[i] + 1
        phi[ ,i] <- phi_prop[ ,i]
      } else {   
        phi_reject[i] <- phi_reject[i] + 1 
      }
      if (it<=iterations_adaptive) {
        ap <- min(1, exp(posterior_prop - posterior_curr))
        if (is.na(ap)) {ap<-0}
        phi_beta[i] <- max(min_beta,phi_beta[i]*exp(phi_delta/max(1, it - iterations_no_diminish) * (phi_target_acceptance - ap)))
      }
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
    cat("  phi min acceptance rate", min(phi_accept/(phi_accept + phi_reject)), which.min(phi_accept), "\n")
    
    #### UPDATE THETA ####
    theta <- rbeta(1, theta_a + sum(dat[ , 1:gtypes, ]), theta_b + sum(dat[ , gtypes + 1, ]))  # 1, success, failure
    # check if the current state of the model has zero posterior probability, if so print a warning message indicating that the starting state of the model is invalid
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
