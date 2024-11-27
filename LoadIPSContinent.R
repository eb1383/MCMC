# libraries
library(tidyverse)

# load data
df <- read_csv("IPS_16_19.csv")

# compute the total number of visits for each country and month
df <- df %>%
  group_by(Year, Month, Country) %>%
  summarize(Visits = sum(Weight)) %>% 
  rename(Destination = Country)

# load continent data 
continents <- read_csv("CountriesContinents.csv")
colnames(continents) <- c("Continent", "Destination")
# join continents to ips 
df <- left_join(df, continents, by = "Destination")
# replace destination with continent 
df$Destination <- ifelse(df$Destination %in% destinations, df$Destination, df$Continent)
# compute the new total number of visits for each destination and month
df <- df %>%
  group_by(Year, Month, Destination) %>%
  summarize(Visits = sum(Visits)) 


# load population data
pop <- read_csv("Population_16_19.csv")
df <- left_join(df, pop, by = "Year")
# compute the total number of visits for each month
df_visits <- df %>%
  group_by(Year, Month) %>%
  summarize(Total_Visits = sum(Visits))
# compute the number of people who stayed at home for each month
df_domestic <- df_visits %>%
  left_join(pop, by = "Year") %>%
  mutate(Domestic = Population - Total_Visits) %>%
  mutate(Destination = "Domestic", Visits = Domestic) %>%
  dplyr::select(Year, Month, Destination, Visits) 

df <- bind_rows(df, df_domestic)

# set up array to store data
n <- array(0, c(dest + 1, timepoints))

for(t in 1 : timepoints){
  for (i in 1 : (dest + 1)){
    if(i <= dest){ 
      wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Destination == destinations[i])
    }
    else {
      wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Destination == "Domestic")
    }    
    if(length(wh) == 1){      
      n[i,t] <- df$Visits[wh]     
    }    
    else if(length(wh) > 1){
      print('stop')
    }
  }
  
}


