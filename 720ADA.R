# Initialise required packages.
required_packages <- c('tidyverse', 'readxl', 'lubridate')

# Install and load required packages.
for (p in required_packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

# Set working directory - please CHANGE ACCORDINGLY.
setwd("/Users/rafaelkoll/Desktop/4/Masters/c/Course/ADA/Practical/Breast-cancer-incidence-analysis")

##### POPULATION ESTIMATES - STANDARD POPULATION ####
# Read in the population estimates files for East Anglia for the 2004-2017 and 1971-2003
# periods. Then combine the two using bind_rows, and finally filter the
# result to include data from 1980 until 2017, and only female (sex=2).
pop_8017 <-
  read_csv("./Data/pop_ea_200417.csv") %>% bind_rows(readr::read_table(
    "./Data/population.txt",
    col_names = c("year", "sex", "age", "population"),
    skip = 5
  )) %>% filter(year >= 1980, sex == 2)


##### INCIDENCE RATE OF OUR POPULATION ####
# Load age-specific incidence rate of our population.
cases_0417 <- read_excel("./Data/brca_incidence_2004_17.xlsx",
                         skip=3,
                         n_max=19) %>% select(everything(), -1, -last_col())

# Add the numbers in the last row (90+ label) with the numbers in the previous to last
# since our age bands in the standard population are 1-18, with the 18th including everyone 85+.
age_90_over <- as.numeric(as.vector(cases_0417[18,]))
cases_0417[18,] <- cases_0417[18,] + age_90_over 
# After merging, drop the last column (over 90), since the data is now on the column for 85+.
cases_0417 <- cases_0417[1:18,]

# Convert the wide table to longer.
cases_0417 <- cases_0417 %>% pivot_longer(cols = 1:ncol(cases_0417), names_to = "year",
                                          values_to = "cases") %>% arrange(year)


# Read in the incidence rate for 1980 until 2003.
cases_8003 <- read.csv("./Data/tumour.tsv", sep = '\t')
cases_8003 <- cases_8003 %>% mutate(date_of_diagnosis = dmy(date_of_diagnosis))


# HAVE TO PROBABLY patient_number NUMBER WITH patient.csv FILE TO GET THEIR AGE.
# Create a new 'year' feature, keeping only the year from the date_of_diagnosis.
cases_8003 <- cases_8003 %>% mutate(year = as.numeric(format(date_of_diagnosis, '%Y')))
# Drop date_of_diagnosis
cases_8003 <- cases_8003 %>% select(everything(), -date_of_diagnosis)

cases_8003 %>% group_by(year) %>% summarise(number=n())

# TODO: now have all the count for year without considering age bands (basically all the age bands are in the same count).
# Have to seperate them, probably will be able to do that when I have the patient.csv file.
