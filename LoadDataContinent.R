# This script processes Shigella case data by:
#  - Standardising SNP types
#  - Assigning travel categories
#  - Aggregating case counts by destination, month, and genetic type

# Load required libraries
library(readr)
library(tidyverse)

# Load Shigella case data
df <- read_csv("Shigella_16_19.csv")

### SORT SNP TYPES ####
# Calculate the probability of having a typed SNP
prob_snp_typed <- (length(df$SNP) - sum(is.na(df$SNP))) / length(df$SNP) 
print(prob_snp_typed)  # Probability: 0.9817398

# Replace missing SNP values with "Unknown"
df$SNP <- replace(df$SNP, is.na(df$SNP), "Unknown")

# Define a list of known SNP types
snp_list <- c(
  "1.3.11.394", "1.3.11.843", "1.3.431", "1.1.18", "1.1.1", "1.1.29",
  "1.3.3", "4.49.49", "1.3.20", "1.3.437", "1.60.74", "2.75", "12.32",
  "1.52", "3.152", "34.42", "35.43", "3.47", "33", "78", "Unknown"
)

# Standardize SNP values: Match those in snp_list, otherwise mark as "Other"
for (pattern in snp_list) {
  df$SNP <- ifelse(grepl(paste0("^", pattern), df$SNP), pattern, df$SNP)
}
df$SNP <- ifelse(df$SNP %in% snp_list, df$SNP, "Other")

### SORT TRAVEL CATEGORIES ####
# Label "Non-Travel" cases as "Domestic"
df$Destination[df$Travel == "Non-Travel"] <- "Domestic"

# Handle missing travel data
df$Destination[df$Travel == "Unknown"] <- "Unknown"
df$Destination[is.na(df$Destination)] <- "Unknown"

# Ensure all unknown SNPs are labeled correctly
df$SNP[df$SNP == "Unknown"] <- "Unknown"

# Replace long continent names with shorter versions for readability
df$Destination <- gsub("African continent", "Africa", df$Destination)
df$Destination <- gsub("Asian continent", "Asia", df$Destination)
df$Destination <- gsub("European continent", "Europe", df$Destination)
df$Destination <- gsub("South American continent", "South America", df$Destination)

# Load continent mapping data
continents <- read_csv("CountriesContinents.csv")
colnames(continents) <- c("Continent", "Destination")

# Merge continent information into dataset
df <- left_join(df, continents, by = "Destination")

# Compute destination case counts
df <- df %>%
  add_count(Destination, wt = Cases)

# If a destination has fewer than 20 cases, group it by continent
df <- df %>% mutate(
  Destination = case_when(
    n >= 20 ~ Destination,
    n < 20  ~ Continent
  )
)

# Display unique destinations
unique(df$Destination)

### AGGREGATE DATA ####
# Aggregate cases by Destination, Month, SNP, Year, and Travel category
df <- aggregate(Cases ~ Destination + Month + SNP + Year + Travel, data = df, sum, na.action = NULL)

# Define time parameters
years <- c(2016, 2017, 2018, 2019)
yrs <- length(years)
timepoints <- 12 * yrs  # Total months

# Define unique destinations
destinations <- unique(df$Destination)

# Reorder destinations so "Domestic" and "Unknown" appear last
destinations <- c(setdiff(destinations, c("Domestic", "Unknown")), "Domestic", "Unknown")

# Count real destinations (excluding Domestic/Unknown)
dest <- length(destinations) - 2

# Define unique genetic types
genetictypes <- unique(df$SNP)
gtypes <- length(genetictypes) - 1  # Excluding "Unknown"

# Initialize 3D array to store processed case data
dat <- array(0, c(dest + 2, gtypes + 1, timepoints))

### POPULATE CASES ARRAY ####
for (t in 1:timepoints) {
  for (i in 1:(dest + 2)) {
    for (j in 1:gtypes) {
      
      # Find matching cases for each (month, year, destination, SNP)
      if (i <= dest) {
        wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Destination == destinations[i] & df$SNP == genetictypes[j])
      } else if (i == dest + 1) {
        wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Travel == "Non-Travel" & df$SNP == genetictypes[j])
      } else {
        wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Travel == "Unknown" & df$SNP == genetictypes[j])
      }
      
      # Store case count in array
      if (length(wh) > 0) { 
        dat[i, j, t] <- sum(df$Cases[wh])
      }
      
      # Warning for multiple entries
      if (length(wh) > 1) {
        cat("Multiple entries for", destinations[i], genetictypes[j], t, "\n") 
      }
    }
    
    # Process "Unknown" SNPs separately
    j <- gtypes + 1
    if (i <= dest) {
      wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Destination == destinations[i] & df$SNP == "Unknown")
    } else if (i == dest + 1) {
      wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Travel == "Non-Travel" & df$SNP == "Unknown")
    } else {
      wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Travel == "Unknown" & df$SNP == "Unknown")
    }
    
    # Store case count for "Unknown" SNPs
    if (length(wh) > 0) {
      dat[i, j, t] <- sum(df$Cases[wh])
    }
    
    # Warning for multiple entries
    if (length(wh) > 1) {
      cat("Multiple entries for", destinations[i], genetictypes[j], t, "\n") 
    }
  }
}


