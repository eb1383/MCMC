# libraries
library(tidyverse)
library(reshape2)

### PREPARE THE DATA ####
# load data 
df <- read_csv("")
df$Destination <- ifelse(df$Travel == 'Non-Travel', 'Domestic', df$Destination)
# compute the total number of visits for each country and month
df <- df %>%
  group_by(Destination, SNP) %>%
  summarize(Cases = sum(Cases)) 
# fill nas with unknown
df$Destination[is.na(df$Destination)] <- "Unknown"
# frequency by destination 
df_freq_dest <- df %>% 
  group_by(Destination) %>%
  summarise(N = n()) %>% 
  arrange(desc(N))

# load continent data 
continents <- read_csv("")
colnames(continents) <- c("Continent", "Destination")
# join continents to data 
df <- left_join(df, continents, by = "Destination")
# replace missing continents with unknown
df$Continent[is.na(df$Continent)] <- "Unknown"
# add total count for each destination to dataframe 
df <- left_join(df, df_freq_dest, by = "Destination")
# filter destinations with below 20 values as contintents
df <- df %>% 
  mutate(Destination = ifelse(N < 20, Continent, Destination))
# select relevant columns
df <- df[ , c(1:3)]
# replace na snps
df$SNP <- replace(df$SNP, is.na(df$SNP), "Unknown")

### SNP CLADE 4 ####
df4 <- df
# look at SNP clade 
df4 <- df4 %>%
  mutate(SNP4 = word(SNP, 1, 4, sep = "\\."))
# aggregate data
df4 <- dcast(df4, SNP4 ~ Destination, sum, value.var = "Cases", fill = 0)
# add total
df4$Total <- rowSums(df4[ , -1])

### SNP CLADE 5 ####
df5 <- df
# look at SNP clade 
df5 <- df5 %>%
  mutate(SNP5 = word(SNP, 1, 3, sep = "\\."))
# aggregate data
df5 <- dcast(df5, SNP5 ~ Destination, sum, value.var = "Cases", fill = 0)
# add total
df5$Total <- rowSums(df5[ , -1])

### SNP CLADE 6 ####
df6 <- df
# look at SNP clade 
df6 <- df6 %>%
  mutate(SNP6 = word(SNP, 1, 2, sep = "\\."))
# aggregate data
df6 <- dcast(df6, SNP6 ~ Destination, sum, value.var = "Cases", fill = 0)
# add total
df6$Total <- rowSums(df6[ , -1])

### SNP CLADE 7 ####
df7 <- df
# look at SNP clade 
df7 <- df7 %>%
  mutate(SNP7 = word(SNP, 1, 1, sep = "\\."))
# aggregate data
df7 <- dcast(df7, SNP7 ~ Destination, sum, value.var = "Cases", fill = 0)
# add total
df7$Total <- rowSums(df7[ , -1])

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
  "78"
)

for (pattern in snp_list) {
  df$SNP <- ifelse(grepl(paste0("^", pattern), df$SNP), pattern, df$SNP)
}

df <- aggregate(Cases ~ SNP, df, sum)



