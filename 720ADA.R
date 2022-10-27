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
setwd(
  "/Users/rafaelkoll/Desktop/4/Masters/c/Course/ADA/Practical/Breast-cancer-incidence-analysis"
)


##### POPULATION ESTIMATES AND EUROPEAN STANDARD POPULATION ####
# Read in the population estimates files for East Anglia for the 2004-2017 and 1971-2003
# periods. Then combine the two using bind_rows, and finally filter the
# result to include data from 1980 until 2017, and only female (sex=2).
pop_8017 <-
  read_csv("./Data/pop_ea_200417.csv") %>% bind_rows(readr::read_table(
    "./Data/population.txt",
    col_names = c("year", "sex", "age", "population"),
    skip = 5
  )) %>% filter(year >= 1980, sex == 2)

# Read in European standard population
european_standard <-
  read_excel("./Data/european-standard-pop.xlsx") %>% rename(age = ageband)

##### INCIDENCE RATE OF OUR POPULATION ####
# Load age-specific incidence rate of our population. Include everything except
# first and last columns.
cases_0417 <- read_excel("./Data/brca_incidence_2004_17.xlsx",
                         skip = 3,
                         n_max = 19) %>% select(everything(), -1, -last_col())

# Add the numbers in the last row (90+ label) with the numbers in the previous to last (85-89 label)
# since our age bands in the standard population are 1-18, with the 18th including everyone 85+.
age_90_over <- as.numeric(as.vector(cases_0417[18,]))
cases_0417[18,] <- cases_0417[18,] + age_90_over

# After merging, drop the last column (over 90) -> since the data is now on the column for 85+.
cases_0417 <- cases_0417[1:18,]

# Convert the wide table to longer.
cases_0417 <-
  cases_0417 %>% pivot_longer(
    cols = 1:ncol(cases_0417),
    names_to = "year",
    values_to = "cases"
  ) %>% arrange(year)

# Enter age band for each year. Creates an array from 1 to 18,
# repeats it n times -> n = length of unique years in the dataset.
cases_0417["age"] <-
  rep(c(1:18), times = length(unique(cases_0417$year)))

# Convert year from chr to int.
cases_0417 <- cases_0417 %>% mutate(year = as.integer(year))

# Read in the incidence rate for 1980 until 2003.
cases_8003 <-
  read.csv("./Data/tumour.tsv", sep = '\t') %>% mutate(date_of_diagnosis = dmy(date_of_diagnosis))

# Only keep the first date of diagnosis for each patient -> on later tumours they must be considered as prevalence and not incidence.
cases_8003 <-
  cases_8003 %>% group_by(patient_number) %>% summarise(date_of_diagnosis = min(date_of_diagnosis))

# Put the age band in the cases_8003 dataframe.
patient <-
  read_csv("./Data/patient.csv") %>% mutate(date_of_birth = dmy(date_of_birth))

# Merge the patient and cases_8003 dataframe. Get their age at diagnosis.
cases_8003 <-
  merge(patient, cases_8003, by = "patient_number") %>% mutate(age_diagnosed = as.integer(interval(date_of_birth, date_of_diagnosis) / years(1)))

# Get age group of patient.
cases_8003 <- cases_8003 %>%
  mutate(age = cut(
    age_diagnosed,
    breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, Inf),
    labels = FALSE,
    right = FALSE
  ))

# Get the year (e.g. 1980) from the date of diagnosis.
cases_8003 <-
  cases_8003 %>% mutate(year = as.numeric(format(date_of_diagnosis, '%Y')))

# Get cases number for each age band and each year by summarising.
cases_8003 <-
  cases_8003 %>% group_by(year, age) %>% summarise(cases = n())

# Combine cases from 80-17.
cases_combined <-
  base::rbind(cases_8003, cases_0417) %>% arrange(year, age)

#### YEARLY INCIDENCE RATE FROM 1980 TO 2017 ####
# Combine cases and population
combined_cases_population <-
  left_join(cases_combined, pop_8017, by = c("year", "age"))
combined_cases_population <-
  combined_cases_population %>% mutate(crude_rate = (cases / population) * 100000)

# Calculate age distribution proportions on standard population.
european_standard <-
  european_standard %>% mutate(age_distribution_proportions = europop / sum(europop))
# Combine the combined_cases_population dataframe with the european_standard.
combined_cases_population <-
  left_join(combined_cases_population, european_standard, by = "age")

# Calculate expected incidence, and then variance (v) and standard error (se).
combined_cases_population <-
  combined_cases_population %>% mutate(expected_incidence = crude_rate *
                                         age_distribution_proportions) %>% group_by(age) %>% mutate(v = ((europop / 100000) ^2) * (cases / (population) ^ 2)) %>% mutate(se = sqrt(sum(v)) * 100000)

# Calculate direct age standardised death rate.
direct_age_standardised_death_rate <-
  combined_cases_population %>% group_by(year) %>% summarise(DASDR = sum(expected_incidence))


