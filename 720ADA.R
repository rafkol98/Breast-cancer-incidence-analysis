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
                         n_max = 19) %>% select(everything(),-1,-last_col())

# Add the numbers in the last row (90+ label) with the numbers in the previous to last (85-89 label)
# since our age bands in the standard population are 1-18, with the 18th including everyone 85+.
age_90_over <- as.numeric(as.vector(cases_0417[18, ]))
cases_0417[18, ] <- cases_0417[18, ] + age_90_over

# After merging, drop the last column (over 90) -> since the data is now on the column for 85+.
cases_0417 <- cases_0417[1:18, ]

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

# Combine the cases and population file using full join.
# Replace the NA values introduced (e.g. when there are no recorded cases for an age group in a year) with 0.
combined_cases_population <-
  full_join(cases_combined, pop_8017, by = c("year", "age")) %>% replace(is.na(.), 0)

#### PART 1: YEARLY INCIDENCE RATE FROM 1980 TO 2017 ####
# Combine cases and population
combined_cases_population_part1 <-
  combined_cases_population %>% mutate(crude_rate = (cases / population) * 100000)

# Calculate age distribution proportions on standard population.
european_standard <-
  european_standard %>% mutate(age_distribution_proportions = europop / sum(europop))

# Combine the combined_cases_population_part1 dataframe with the european_standard.
combined_cases_population_part1 <-
  left_join(combined_cases_population_part1, european_standard, by = "age")

# Calculate expected incidence.
combined_cases_population_part1 <-
  combined_cases_population_part1 %>% mutate(expected_incidence = crude_rate *
                                               age_distribution_proportions)

# Calculate direct age standardised death rate (dasdr), standard error,
# and the lower and upper 95 confidence intervals.
dasdr_final_table <-
  combined_cases_population_part1 %>% group_by(year) %>% summarise(
    dasdr = sum(expected_incidence),
    standard_error = sqrt(sum(((
      europop / 100000
    ) ^ 2) * (
      cases / (population ^ 2)
    ))) * 100000,
    lower_95_CI = dasdr - 1.96 * (standard_error),
    upper_95_CI = dasdr + 1.96 * (standard_error)
  )

# Drop the standard error before exporting the table (the exercise only asks for
# incidence rate and the confidence intervals).
dasdr_final_table <-
  dasdr_final_table %>% select(everything(),-standard_error)

# Round the columns to two decimal places.
dasdr_final_table <-
  dasdr_final_table %>% mutate_if(is.numeric, round, digits = 2)

# Export the final dasdr table.
write.csv(
  dasdr_final_table,
  "./Generated data/female_breast_cancer_incidence_standardised_80_17.csv",
  row.names = FALSE
)


my_theme <- theme(
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text.x = element_text(size = 14, angle=45, hjust=1),
  axis.text.y = element_text(size = 14),
  plot.title = element_text(hjust = 0.5),
  panel.grid.minor.x = element_blank()
)

# Line plot with confidence intervals
plot_breast_cancer_incidence <-
  ggplot(dasdr_final_table, aes(year, dasdr)) +
  geom_line(col = 'red') +
  geom_point(
    colour = 'red',
    fill = 'red',
    size = 2,
    shape = 22
  ) +
  geom_ribbon(aes(ymin = lower_95_CI, ymax = upper_95_CI), alpha = 0.1) + labs(x = "Year",  y = "Age-standardised incidence rate \n (100,000 person years)") + scale_x_continuous(breaks = seq(1980, 2017, by = 5)) + scale_y_continuous(breaks = seq(0, 180, by = 10), limits =
                                                                                                                                                                                                                                           c(0, NA)) + ggtitle("Female breast cancer incidence in East Anglia (1980 to 2017)") + my_theme
plot_breast_cancer_incidence

# # Save plot.
# ggsave(
#   "./Generated data/female breast cancer incidence (1980-2017).png",
#   plot_breast_cancer_incidence
# )


# ALTERNATIVE 1: change in confidence intervals visualisation.
plot_breast_cancer_incidence_alt1 <-
  ggplot(dasdr_final_table, aes(year, dasdr)) +
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), colour = "grey65") +
  geom_line(colour = "red") +
  geom_point(
    colour = "red",
    fill = "red",
    size = 2,
    shape = 22
  ) + labs(x = "Year",  y = "Age-standardised incidence rate \n (100,000 person years)") + scale_x_continuous(breaks = seq(1980, 2017, by = 5)) + scale_y_continuous(breaks = seq(0, 180, by = 10), limits =
                                                                                                                                                                       c(0, NA)) + ggtitle("Female breast cancer incidence in East Anglia (1980 to 2017)") + my_theme
plot_breast_cancer_incidence_alt1                                                                                                                    

# ALTERNATIVE 2: Bar plots with condfidence intervals.
plot_breast_cancer_incidence_alt2 <-
  ggplot(dasdr_final_table, aes(year, dasdr)) +
  geom_bar(fill = "red",  stat = "identity", alpha = 0.7) +
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), colour = "grey30") +
  labs(x = "Year",  y = "Age-standardised incidence rate \n (100,000 person years)") + scale_x_continuous(breaks = seq(1980, 2017, by = 5)) + scale_y_continuous(breaks = seq(0, 180, by = 10), limits =  c(0, NA)) + ggtitle("Female breast cancer incidence in East Anglia (1980 to 2017)") + my_theme

plot_breast_cancer_incidence_alt2


#### PART 2 ####

# Create new age groups for the age groups 20-49, 50-69, 70+.
# If between the age of 20-49 that means they are currently in the 5-10 age groups. Assign them to new_age_gp = 1.
# If between the age of 50-69 that means they are currently in the 11-14 age groups. Assign them to new_age_gp = 2.
# If between the age of 70+ that means they are currently in the 15-18 age groups. Assign them to new_age_gp = 3.
combined_cases_population_part2 <-
  combined_cases_population %>% mutate(new_age_gp = case_when(
    between(age, 1, 4) ~ '<20',
    between(age, 5, 10) ~ '20-49',
    between(age, 11, 14) ~ '50-69',
    between(age, 15, 18) ~ '70+'
  )) %>% filter(new_age_gp != '<20') %>% select(everything(), -age)

# Group by year and new_age_gp. Summarise cases and population. Calculate age specific incidence rate for each year.
final_table_part2 <-
  combined_cases_population_part2  %>%  group_by(year, new_age_gp) %>% summarise(cases = sum(cases), population = sum(population)) %>% mutate(incidence = (cases / population) * 100000)

# Scatter plot for age specific incidence (part 2).
ggplot(final_table_part2, aes(x=year, y=incidence, color=new_age_gp, shape=new_age_gp)) + geom_point() +
  labs(x = "Year",  y = "Incidence") + scale_x_continuous(breaks = seq(1980, 2017, by = 5)) + scale_y_continuous(breaks = seq(0, 500, by = 50), limits =  c(0, NA)) + ggtitle("Age specific female breast cancer incidence (1980 to 2017)") + my_theme

