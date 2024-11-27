# MCMC
Markov Chain Monte Carlo (MCMC) model to analyse data given large-scale missingness

# LoadDataContinent.R
R script processes surveillance data into an array suitable for modelling, aggregating cases by month, year, destinations, and SNP types.

The SNP Type Classification calculate the probability of SNPs being typed, matching SNP values to a predefined list (aggregating SNP with the same cluster)

The Travel Type Classification assigns cases to domestic, travel or unknown. Joins continent information to destinations and aggregates cases by continent for destinations with fewer than 20 cases.

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
