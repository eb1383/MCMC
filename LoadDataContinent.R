# libraries
library(readr)
library(tidyverse)

# load data
df <- read_csv("Shigella_16_19.csv")

### SORT SNP TYPES ###
# probability of being typed
(length(df$SNP) - sum(is.na(df$SNP))) / length(df$SNP) # probability 0.9817398
# replace na snps
df$SNP <- replace(df$SNP, is.na(df$SNP), "Unknown")
# define snp list
snp_list <- c(
  "1.3.11.394",
  "1.3.11.843",
  "1.3.431",
  "1.1.18",
  "1.1.1",
  "1.1.29",
  "1.3.3",
  "4.49.49",
  "1.3.20",
  "1.3.437",
  "1.60.74",
  "2.75",
  "12.32",
  "1.52",
  "3.152",
  "34.42",
  "35.43",
  "3.47",
  "33",
  "78",
  "Unknown"
)
# marks SNPS in snp_list
for (pattern in snp_list) {
  df$SNP <- ifelse(grepl(paste0("^", pattern), df$SNP), pattern, df$SNP)
}
#  mark other for SNPs not in snp_list
df$SNP <- ifelse(df$SNP %in% snp_list, df$SNP, "Other")

### SORT TRAVEL TYPES ###
# which cases are non-travel
wh <- which(df$Travel == "Non-Travel")
# label non-travel cases as domestic
df$Destination[wh] = "Domestic"
# which cases are missing
wh <- which(df$Travel == "Unknown")
# label missing cases as unknown
df$Destination[wh] = "Unknown"
# which destinations are missing given they travelled
wh <- which(is.na(df$Destination))
# label missing cases as unknown
df$Destination[wh] = "Unknown"
# which snp are missing
wh <- which(df$SNP == "Unknown")
# label missing snp as unknown
df$SNP[wh] = "Unknown"

# replace continent names for readability 
df$Destination <- gsub("African continent", "Africa", df$Destination)
df$Destination <- gsub("Asian continent", "Asia", df$Destination)
df$Destination <- gsub("European continent", "Europe", df$Destination)
df$Destination <- gsub("South American continent", "South America", df$Destination)

# load continent data 
continents <- read_csv("CountriesContinents.csv")
colnames(continents) <- c("Continent", "Destination")
# join continents to cases 
df <- left_join(df, continents, by = "Destination")
# destination count
df <- df %>% 
  add_count(Destination, wt = Cases)
# replace destination with continent if under 10 cases
df <- df %>% mutate(
        Destination = case_when(
                n >= 20 ~ Destination,
                n < 20  ~ Continent
        )
)
# view destinations with most cases
unique(df$Destination)

### COMPILE DATA ###
# aggregate cases 
df <- aggregate(Cases ~ Destination + Month + SNP + Year + Travel, data = df, sum, na.action = NULL)
# years as list
years = c(2016, 2017, 2018, 2019)
# number of yrs
yrs = length(years)
# number of timepoints
timepoints = 12*yrs 
# number of destinations
destinations = unique(df$Destination)
# reorder destination
destinations <- c(setdiff(destinations, c("Domestic", "Unknown")), "Domestic", "Unknown")
# number of true destinations
dest = length(destinations) - 2
# genetictypes as list 
genetictypes = unique(df$SNP)
# number of genetic types
gtypes = length(genetictypes) - 1

# set up array to store data
dat <- array(0, c(dest + 2, gtypes + 1, timepoints))

for(t in 1 : timepoints){
  for (i in 1 : (dest + 2)){
    for (j in 1 : gtypes){
      # typed
      # picking out rows which match
      # all cases in same month as t
      if(i <= dest){
        wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Destination == destinations[i] & df$SNP == genetictypes[j])
      }
      else if(i == dest + 1){
        wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Travel == "Non-Travel" & df$SNP == genetictypes[j])
      }
      else {
        wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Travel == "Unknown" & df$SNP == genetictypes[j])
      }
      if(length(wh) > 0){ 
        dat[i,j,t] <- sum(df$Cases[wh])
      }
      if(length(wh) > 1){
        cat("Multiple entries for", destinations[i], genetictypes[j], t, "\n") 
      }
    }
    j = gtypes + 1
    # untyped
    if(i <= dest){
      wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Destination == destinations[i] & df$SNP == "Unknown")
    }
    else if(i == dest + 1){
      wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Travel == "Non-Travel" & df$SNP == "Unknown")
    }
    else {
      wh <- which(df$Month == getMonth(t) & df$Year == years[getYear(t)] & df$Travel == "Unknown" & df$SNP == "Unknown")
    }
    if(length(wh) > 0){
      dat[i,j,t] <- sum(df$Cases[wh])
    }
    if(length(wh) > 1){
      cat("Multiple entries for", destinations[i], genetictypes[j], t, "\n") 
    }
  }
}


