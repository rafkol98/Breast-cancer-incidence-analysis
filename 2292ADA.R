# Applied Data Analysis - Assignment.
# Author: 305542292
# Date: 10/11/2022

# Set working directory - please CHANGE ACCORDINGLY.
setwd(
  "/Users/rafaelkoll/Desktop/4/Masters/c/Course/ADA/Practical/Breast-cancer-incidence-analysis"
)

# Initialise required packages.
required_packages <- c("tidyverse", "readxl", "lubridate", "lemon")

# Install and load required packages.
for (p in required_packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

#### THEME - USED FOR PLOTS ####
# Create my own theme - to be used in all the plots I will generate.
my_theme <- theme(
  # Format x and y axis.
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
  axis.text.y = element_text(size = 14),
  plot.title = element_text(hjust = 0.5),
  # Remove border from plots.
  panel.border = element_blank(),
  # Eliminate minor grid lines.
  panel.grid.minor = element_blank(),
  # Reduce the opacity of grid lines.
  panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255))
)

##### POPULATION ESTIMATES AND EUROPEAN STANDARD POPULATION ####
# Read in the population estimates files for East Anglia for the 2004-2017 and 1971-2003
# periods. Then combine the two using bind_rows, and finally filter the
# result to include data from 1980 until 2017, and only female (sex=2).
pop_8017 <-
  read_csv("./Data/pop_ea_200417.csv") %>%
  bind_rows(readr::read_table(
    "./Data/population.txt",
    col_names = c("year", "sex", "age", "population"),
    skip = 5
  )) %>%
  filter(year >= 1980, sex == 2)

# Read in European standard population
european_standard <-
  read_excel("./Data/european-standard-pop.xlsx") %>% rename(age = ageband)

##### CASES IN EAST ANGLIA POPULATION ####
# Load age-specific breast cancer incidence of the East Anglia population. Include everything except
# first (Row Labels) and last (Grand Total) columns.
cases_0417 <- read_excel("./Data/brca_incidence_2004_17.xlsx",
  skip = 3,
  n_max = 19
) %>% select(everything(), -1, -last_col())

# Add the numbers in the last row (90+ label) with the numbers in the previous to last (85-89 label)
# since our age bands in the standard population are 1-18, with the 18th including everyone 85+.
cases_0417[18, ] <- cases_0417[18, ] + as.numeric(as.vector(cases_0417[19, ]))

# After adding the two rows, drop the last column (over 90) -> the data is now included on the row for 85+.
cases_0417 <- cases_0417[1:18, ]

# Convert the wide table to long format. Sort out the resulting data frame by year.
cases_0417 <-
  cases_0417 %>%
  pivot_longer(
    cols = 1:ncol(cases_0417),
    names_to = "year",
    values_to = "cases"
  ) %>%
  arrange(year)

# Enter age band for each year. Creates an array from 1 to 18,
# repeats it n times. Where n = length of unique years in the dataset.
cases_0417["age"] <-
  rep(c(1:18), times = length(unique(cases_0417$year)))

# Convert year from chr to int.
cases_0417 <- cases_0417 %>% mutate(year = as.integer(year))

# Read in the incidence rate for 1980 until 2003.
cases_8003 <-
  read.csv("./Data/tumour.tsv", sep = "\t") %>% mutate(date_of_diagnosis = dmy(date_of_diagnosis))

# Only keep the first date of diagnosis for each patient.
# On later tumours they must be considered as prevalence and not incidence.
cases_8003 <-
  cases_8003 %>%
  group_by(patient_number) %>%
  summarise(date_of_diagnosis = min(date_of_diagnosis))

# Put the age band in the cases_8003 dataframe.
patient <-
  read_csv("./Data/patient.csv") %>% mutate(date_of_birth = dmy(date_of_birth))

# Merge the patient and cases_8003 dataframe. Get their age at diagnosis.
cases_8003 <-
  merge(patient, cases_8003, by = "patient_number") %>%
  mutate(age_diagnosed = as.integer(interval(date_of_birth, date_of_diagnosis) / years(1)))

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
  cases_8003 %>% mutate(year = as.numeric(format(date_of_diagnosis, "%Y")))

# Get cases numbers for each age band and each year by grouping by year and age
# and then summarising.
cases_8003 <-
  cases_8003 %>%
  group_by(year, age) %>%
  summarise(cases = n())

# Combine cases from 80-17 by binding the rows of the two dataframes.
# Then sort the dataframe by year and age.
cases_combined <-
  base::rbind(cases_8003, cases_0417) %>% arrange(year, age)

# Join the cases and population data frames using full join.
# Replace the NA values introduced (e.g. when there are no recorded cases for an age group in a year) with 0.
combined_cases_population <-
  full_join(cases_combined, pop_8017, by = c("year", "age")) %>% replace(is.na(.), 0)





#### PART 1: YEARLY INCIDENCE RATE FROM 1980 TO 2017 ####

# Calculate crude rate (per 100000 person years) by dividng cases by population
# and then multiplying by 100,000.
combined_cases_population_part1 <-
  combined_cases_population %>% mutate(crude_rate = (cases / population) * 100000)

# Include a new column that showcases the population distribution proportions for
# each age group of the European standard population.
european_standard <-
  european_standard %>% mutate(age_distribution_proportions = europop / sum(europop))

# Combine the combined_cases_population_part1 dataframe with the european_standard.
combined_cases_population_part1 <-
  left_join(combined_cases_population_part1, european_standard, by = "age")

# Calculate the expected incidence by multiplying the crude rate with the age distribution
# proportions..
combined_cases_population_part1 <-
  combined_cases_population_part1 %>% mutate(expected_incidence = crude_rate *
    age_distribution_proportions)

# Calculate direct age standardised incidence rate, standard error,
# and the lower and upper 95 confidence intervals.
table_part1 <-
  combined_cases_population_part1 %>%
  group_by(year) %>%
  summarise(
    standardised_incidence = sum(expected_incidence),
    standard_error = sqrt(sum(((europop / 100000)^2) * (cases / (population^2)))) * 100000,
    lower_95_CI = standardised_incidence - 1.96 * (standard_error),
    upper_95_CI = standardised_incidence + 1.96 * (standard_error)
  )

# Drop the standard error before exporting the table (the exercise only asks for
# incidence rate and the confidence intervals).
table_part1 <-
  table_part1 %>% select(everything(), -standard_error)

# Round the numeric columns to two decimal places (if its integer then it ignores that).
final_table_part1 <-
  table_part1 %>% mutate_if(is.numeric, round, digits = 2)

# Export table to csv. Used to copy table in the report.
# write.csv(final_table_part1, "Generated tables/table_1.csv", row.names = FALSE)

# Figure for part 1. Line plot with confidence intervals shown with a light shadowed ribbon.
fig_part1 <-
  ggplot(final_table_part1, aes(year, standardised_incidence)) +
  # Include line and points for every year.
  geom_line(col = "red") +
  geom_point(
    colour = "red",
    fill = "red",
    size = 2,
    shape = 22
  ) +
  # Show the confidence intervals with a shadowed ribbon.
  geom_ribbon(aes(ymin = lower_95_CI, ymax = upper_95_CI), alpha = 0.03, linetype = 2, colour = "gray80") +
  labs(x = "Year", y = "Age-standardised incidence rate \n (100,000 person-yrs)") +
  # Intervals of 5 in the x-axis and 20 in the y-axis.
  scale_x_continuous(breaks = seq(1980, 2017, by = 5), minor_breaks = 0) +
  scale_y_continuous(
    breaks = seq(0, 200, by = 20), limits =
      c(0, NA), expand = c(0, 0)
  ) +
  theme_bw() +
  my_theme

# Show figure created for part 1.
fig_part1

# ALTERNATIVE 1: Bar plots with condfidence intervals.
fig_part1_alt1 <-
  ggplot(final_table_part1, aes(year, standardised_incidence)) +
  # Show the data using a bar plot.
  geom_bar(fill = "red", stat = "identity", alpha = 0.7) +
  # Show the confidence intervals with error bars.
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), colour = "grey30") +
  labs(x = "Year", y = "Age-standardised incidence rate \n (100,000 person years)") +
  # Intervals of 5 in the x-axis and 20 in the y-axis.
  scale_x_continuous(breaks = seq(1980, 2017, by = 5)) +
  scale_y_continuous(breaks = seq(0, 200, by = 20), limits = c(0, NA), expand = c(0, 10)) +
  theme_bw() +
  my_theme

# Show the first alternative created for part 1.
fig_part1_alt1


# ALTERNATIVE 2: change in confidence intervals visualisation - use errorbars instead of ribbon.
fig_part1_alt2 <-
  ggplot(final_table_part1, aes(year, standardised_incidence)) +
  # Show the confidence intervals with error bars.
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), colour = "grey65") +
  # Include line and points for every year. After error bars to be on top of them.
  geom_line(col = "red") +
  geom_point(
    colour = "red",
    fill = "red",
    size = 2,
    shape = 22
  ) +
  labs(x = "Year", y = "Age-standardised incidence rate \n (100,000 person-yrs)") +
  # Intervals of 5 in the x-axis and 20 in the y-axis.
  scale_x_continuous(breaks = seq(1980, 2017, by = 5), minor_breaks = 0) +
  scale_y_continuous(
    breaks = seq(0, 200, by = 20), limits =
      c(0, NA), expand = c(0, 0)
  ) +
  theme_bw() +
  my_theme

# Show the second alternative created for part 1.
fig_part1_alt2



#### PART 2 ####

# Create new age groups for the age groups 20-49, 50-69, 70+.
# If between the age of 20-49 that means they are currently in the 5-10 age groups. Assign them to new_age_gp = 1.
# If between the age of 50-69 that means they are currently in the 11-14 age groups. Assign them to new_age_gp = 2.
# If between the age of 70+ that means they are currently in the 15-18 age groups. Assign them to new_age_gp = 3.
combined_cases_population_part2 <-
  combined_cases_population %>%
  mutate(new_age_gp = case_when(
    between(age, 1, 4) ~ "<20",
    between(age, 5, 10) ~ "20-49",
    between(age, 11, 14) ~ "50-69",
    between(age, 15, 18) ~ "70+"
  )) %>%
  # Exclude people with age group less than 20.
  filter(new_age_gp != "<20") %>%
  # Drop the previous age column.
  select(everything(), -age)

# Group by year and new_age_gp. Summarise cases and population. Calculate age specific incidence rate for each year.
final_table_part2 <-
  combined_cases_population_part2 %>%
  group_by(year, new_age_gp) %>%
  summarise(cases = sum(cases), population = sum(population)) %>%
  mutate(incidence = (cases / population) * 100000)

# Round the columns to two decimal places.
final_table_part2 <-
  final_table_part2 %>% mutate_if(is.numeric, round, digits = 2)

# Convert final table to wide format. First keep only year, age group and incidence columns
# and then convert the table from long to wide format.
final_table_part2_wider <- final_table_part2 %>%
  select(year, new_age_gp, incidence) %>%
  pivot_wider(names_from = new_age_gp, values_from = incidence)

# Export table to csv.
# write.csv(final_table_part2_wider, "Generated tables/table_2_wide.csv", row.names = FALSE)

# Scatter plot combined with line plot for age specific incidence (part 2).
# Coloured considering colour vision defficiency. Used shapes to distinguish between categories.
fig_part2 <- ggplot(final_table_part2, aes(x = year, y = incidence, colour = new_age_gp, shape = new_age_gp)) +
  geom_line(lwd = 1) +
  geom_point(size = 2.2) +
  labs(x = "Year", y = "Crude Incidence \n (100,000 person-yrs)") +
  # Intervals of 5 in the x-axis and 50 in the y-axis.
  scale_x_continuous(breaks = seq(1980, 2017, by = 5)) +
  scale_y_continuous(breaks = seq(0, 470, by = 50), limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  theme_bw() +
  my_theme +
  # Enter the colours used for each category manually.
  scale_color_manual(name = "Age group", limits = c("70+", "50-69", "20-49"), values = c("#0072B2", "#E69F00", "#AA00AA")) +
  scale_shape_discrete(name = "Age group", limits = c("70+", "50-69", "20-49"))

# Show the chosen figure created for part 2.
fig_part2

# ALTERNATIVE 1: ggplot standard colours vs colours chosen to aid colour vision deficiency.
fig_part2_alt_1 <- ggplot(final_table_part2, aes(x = year, y = incidence, colour = new_age_gp, shape = new_age_gp)) +
  geom_line(lwd = 1) +
  geom_point(size = 2.2) +
  labs(x = "Year", y = "Crude Incidence \n (100,000 person-yrs)") +
  # Intervals of 5 in the x-axis and 50 in the y-axis.
  scale_x_continuous(breaks = seq(1980, 2017, by = 5)) +
  scale_y_continuous(breaks = seq(0, 470, by = 50), limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  theme_bw() +
  my_theme

# Show the first alternative created for part 2.
fig_part2_alt_1

# ALTENATIVE 2: Scatter plot for age specific incidence (part 2).
fig_part2_alt_2 <- ggplot(final_table_part2, aes(x = year, y = incidence, colour = new_age_gp)) +
  geom_point() +
  labs(x = "Year", y = "Crude Incidence \n (100,000 person-yrs)") +
  # Intervals of 5 in the x-axis and 50 in the y-axis.
  scale_x_continuous(breaks = seq(1980, 2017, by = 5)) +
  scale_y_continuous(breaks = seq(0, 470, by = 50), limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  theme_bw() +
  my_theme +
  # Enter the colours used for each category manually.
  scale_color_manual(name = "Age group", limits = c("70+", "50-69", "20-49"), values = c("#0072B2", "#E69F00", "#AA00AA")) +
  scale_shape_discrete(name = "Age group", limits = c("70+", "50-69", "20-49"))

# Show the second alternative created for part 2.
fig_part2_alt_2