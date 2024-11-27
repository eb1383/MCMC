# Bayesian Hierarchical MCMC Model for Estimating Infection Rates
This R code implements a Bayesian hierarchical model to estimate infection rates based on travel and genetic data, using the Metropolis-Hastings algorithm for MCMC sampling. It models the spread of infection, considering multiple destinations, timepoints, and genetic types, while also accounting for missing or incomplete data (such as unknown travel destinations or unrecorded genetic types). The code is organised into several key sections: data loading, model setup, likelihood computation, MCMC sampling, and parameter updates. The initial data manipulating processes are outlined below:

# SNPSort.R
R script summarised cases by hierarchical SNP clades:
- SNP4: First 4 segments of SNP.
- SNP5: First 3 segments.
- SNP6: First 2 segments.
- SNP7: First segment (broadest clade).

# LoadDataContinent.R
R script processes surveillance data into an array suitable for modelling, aggregating cases by month, year, destinations, and SNP types.
- The SNP Type Classification calculate the probability of SNPs being typed, matching SNP values to a predefined list (as in SNPSort.R).
- The Travel Type Classification assigns cases to domestic, travel or unknown. Joins continent information to destinations and aggregates cases by continent for destinations with fewer than 20 cases.

# LoadIPSContinent.R
R script processes International Passenger Survey (IPS) data to calculate the total number of visits for each destination across months and years. The data is integrated with continent and population data to compute visits to domestic and international destinations, storing the final results in a structured array.

# InitialValues.R
This section loads initial data from StartingValues.RData. This includes surveillance and IPS (International Passenger Survey) data.
Scales the traveler count (n) by dividing it by 1000 for easier numerical handling. The section further contains custom utility functions used in the script:
- getYear(t): Computes the year from a given timepoint t
- getMonth(t): Computes the month of the year for a given timepoint t
- lddirichlet(x, a): Computes the log-density of the Dirichlet distribution for a given input vector x and parameter vector a. Only elements where a is non-zero are included in calculations.

Prior distributions for model parameters are defined:
- Lambda: Beta prior
- Mu: Uniform distribution across 12 months
- P: Uniform weighted distribution over genetic types
- Theta: Beta prior

The structure of the phi_alpha matrix, which represents probabilities related to traveler destination recording accuracy, is defined:
- Diagonal elements represent probabilities of correct recording for each destination.
- Specific rows and columns represent probabilities of misclassifications into unknown or domestic categories.

The starting values for model parameters are specified:
- lambda_init: Uniform distribution over destinations
- mu_init: Uniform distribution over months
- p_init: Uniform distribution over genetic types
- phi_init: Matrix for destination probabilities with manually assigned probabilities for misclassification
- theta_init: Predefined initial value for theta

# MCMC.R

You will need to load your own data for analysis. The following two scripts are expected to be sourced:

LoadDataContinent.R: Loads the infection data across different continents or regions.
LoadIPSContinent.R: Loads additional data, such as travel information or genetic types across regions.
Data should be provided in matrices:

n: Matrix representing the number of cases per destination over time, with columns representing months and rows representing destinations.
dat: Infection data, typically with dimensions representing the number of destinations, genetic types, and timepoints.

The model uses several key parameters:

lambda: Rate of infection for each destination and year.
mu: Proportion of individuals traveling to each destination per month.
p: Probability of getting each genetic type upon infection at each destination.
phi: Probability of observing each destination correctly based on the traveler’s recorded data.
theta: Probability of observing a valid genetic type in the data.

The model uses several key parameters:

lambda: Rate of infection for each destination and year.
mu: Proportion of individuals traveling to each destination per month.
p: Probability of getting each genetic type upon infection at each destination.
phi: Probability of observing each destination correctly based on the traveler’s recorded data.
theta: Probability of observing a valid genetic type in the data.

Initial Setup:
The model initialization sets the initial values for the parameters lambda_init, mu_init, p_init, phi_init, and theta_init. These are essential for starting the MCMC chain.

Likelihood Function:
The llikelihood function calculates the log-likelihood for the model parameters given the data. It is used during the MCMC sampling to evaluate the fit of each proposed parameter set.

MCMC Sampling:
The core of the model is the Metropolis-Hastings sampling algorithm implemented in the metropolis_algorithm function. You can run the sampling by calling the function with your data and priors:

This will run the MCMC sampler for 10,000 iterations, updating the parameters with each step.

Convergence Diagnostics:
After running the MCMC, you can assess convergence by examining the acceptance rates for the Metropolis-Hastings updates. The acceptance rates for the parameters lambda, mu, p, phi, and theta are printed at each iteration, providing insights into the mixing and convergence of the chain.

Posterior Samples:
The posterior samples for each parameter are stored in arrays (lambda_stored, mu_stored, p_stored, phi_stored, and theta_stored). You can plot and analyze these samples using standard plotting techniques in R, e.g., with ggplot2.

To improve the efficiency of the MCMC algorithm, several tuning parameters are used:

lambda_beta, mu_beta, p_beta, phi_beta: Step size parameters for the random walk proposals.
lambda_target_acceptance, mu_target_acceptance, p_target_acceptance, phi_target_acceptance: Target acceptance rates for the MCMC updates. These can be adjusted to improve the sampling efficiency.
Output
At the end of the sampling process, the following outputs are available:

lambda_stored: Array containing the posterior samples of the infection rates (lambda).
mu_stored: Array containing the posterior samples for the travel proportions (mu).
p_stored: Array containing the posterior samples for the genetic type probabilities (p).
phi_stored: Array containing the posterior samples for the recording probabilities (phi).
theta_stored: Array containing the posterior samples for the genetic type validity (theta).

Convergence Issues: If the algorithm does not seem to converge, consider adjusting the priors, increasing the number of iterations, or changing the step sizes (lambda_beta, mu_beta, etc.).
