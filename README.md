# Epidemiological Modeling with MCMC 🧬

This repository contains R scripts designed for **epidemiological modeling** of infectious diseases using **Markov Chain Monte Carlo (MCMC)** techniques. The models estimate critical parameters like infection rates, travel patterns, and disease transmission probabilities using statistical methods such as the **Dirichlet distribution**. 📊🔬

## 📦 Dependencies

To run the scripts in this repository, you'll need the following R packages:

- **MCMCpack** 📉: For **Metropolis-Hastings MCMC** sampling.
- **mvtnorm** 📐: For working with **multivariate normal distributions**.
- **ggplot2** 📈: For creating **data visualisations** and plots.

You can install them with these commands:

```r
install.packages("MCMCpack")
install.packages("mvtnorm")
install.packages("ggplot2")
```

## 🗂️ Files in this Repository

### 1. `MCMC.R` 🧮
This script runs the primary epidemiological model using MCMC sampling. It estimates key parameters like:

- **λ (lambda)**: Travel rate per year (split by destination and domestic travel) 🌍✈️
- **μ (mu)**: Seasonal variation in travel 📅
- **p**: Proportion of genetic types (susceptibility to disease) 🧬
- **θ (theta)**: Probability of underreporting travel and disease data 🚨
- **φ (phi)**: Probability of misclassifying a destination or travel data ❌

#### Functions:

- `getYear(t)`: Extracts the year from a given timepoint 🗓️
- `getMonth(t)`: Extracts the month from a given timepoint 📅
- `lddirichlet(x, a)`: Computes the log-density of the Dirichlet distribution 🔢
- `llikelihood()`: Computes the log-likelihood of observed data based on the model 📉
- `metropolis_algorithm()`: Executes Metropolis-Hastings MCMC for parameter sampling 🔄

### 2. `LoadDataContinent.R` 🌍
This script loads travel and surveillance data specific to each continent or destination. It's essential for the model as it provides the population and travel data used to simulate disease transmission. 🌏

### 3. `LoadIPSContinent.R` 🏛️
Loads data on the **International Political Systems (IPS)** affecting travel and migration patterns, helping to model how these systems influence the spread of diseases across borders. 🌐

### 4. `init_variables.R` 🔑
This script initializes the key variables for the model, including prior distributions and starting values for MCMC sampling. It sets up important parameters and ensures reproducibility by setting a random seed. 🔄

## 🚀 Instructions for Usage

### Running the Model

1. **Install Dependencies**: Ensure you have all required R packages installed. 💻
2. **Load Data**: Run the scripts `LoadDataContinent.R` and `LoadIPSContinent.R` to load necessary travel and political system data. 🌍🔍
3. **Initialise Variables**: Run `init_variables.R` to set initial values and prior distributions for the model. 🔑
4. **Run the Model**: Execute `main_model.R` to start the **MCMC sampling** and generate posterior parameter estimates. 🔄💻
5. **Visualise Results**: After running the model, use `ggplot2` to create plots or export the results for further analysis. 📊

### Customizing the Model
Modify the **prior distributions**, **MCMC parameters**, or **initial values** as needed. For example, you can adjust the **lambda**, **mu**, or **p** priors to reflect different assumptions about disease spread and travel patterns. 🛠️
