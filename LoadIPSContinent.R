# This script processes international travel data by:
#  - Summarising visits per country and month
#  - Mapping destinations to continents
#  - Adjusting data for domestic travel
#  - Storing processed data in an array format

# Load required libraries
library(tidyverse)

# Load International Passenger Survey (IPS) data
df <- read_csv("")

# Compute total number of visits for each country and month
df <- df %>%
  group_by(Year, Month, Country) %>%
  summarize(Visits = sum(Weight), .groups = "drop") %>% 
  rename(Destination = Country)  # Rename column for consistency

### MAP COUNTRIES TO CONTINENTS ####
# Load continent mapping data
continents <- read_csv("")
colnames(continents) <- c("Continent", "Destination")

# Merge continent data with IPS data
df <- left_join(df, continents, by = "Destination")

# If a destination is not in the predefined list, replace it with its continent
df$Destination <- ifelse(df$Destination %in% destinations, df$Destination, df$Continent)

# Recompute total visits, now aggregated by continent (or country if large enough)
df <- df %>%
  group_by(Year, Month, Destination) %>%
  summarize(Visits = sum(Visits), .groups = "drop")

### INCORPORATE POPULATION ####
# Load population data for the corresponding years
pop <- read_csv("Population_16_19.csv")

# Merge population data into the dataset based on the Year
df <- left_join(df, pop, by = "Year")

# Compute the total number of visits for each month
df_visits <- df %>%
  group_by(Year, Month) %>%
  summarize(Total_Visits = sum(Visits), .groups = "drop")

# Compute the number of people who stayed at home (domestic travelers)
df_domestic <- df_visits %>%
  left_join(pop, by = "Year") %>%
  mutate(Domestic = Population - Total_Visits,  # Calculate domestic stayers
         Destination = "Domestic", 
         Visits = Domestic) %>%
  dplyr::select(Year, Month, Destination, Visits)  # Select relevant columns

# Append domestic travel data to the dataset
df <- bind_rows(df, df_domestic)

### STORE DATA IN AN ARRAY ####
# Initialise a storage array for visits
n <- array(0, c(dest + 1, timepoints))

# Populate the array with visit data for each time point
for (t in 1:timepoints) {
  for (i in 1:(dest + 1)) {
    
    # Identify the corresponding data rows
    if (i <= dest) { 
      wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Destination == destinations[i])
    } else {
      wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Destination == "Domestic")
    }
    
    # Store visit counts in the array
    if (length(wh) == 1) {      
      n[i, t] <- df$Visits[wh]     
    } else if (length(wh) > 1) {
      print("Error: Multiple entries detected")
    }
  }
}
