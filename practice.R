# Load required libraries
library(tidyverse)
library(ggplot2)

# Load malaria data
malaria_data <- readRDS("malaria_facility_count_data.rds")

# Inspect the data
glimpse(malaria_data)
summary(malaria_data)
head(malaria_data)


# Reshape data to long format for age groups
malaria_long <- malaria_data %>%
	select(data_date, starts_with("malaria_rdt_")) %>%
	pivot_longer(
		cols = starts_with("malaria_rdt_"),
		names_to = "age_group",
		values_to = "prevalence"
	) %>%
	mutate(age_group = case_when(
		age_group == "malaria_rdt_0-4" ~ "0-4",
		age_group == "malaria_rdt_5-14" ~ "5-14",
		age_group == "malaria_rdt_15" ~ "15+"
	))

# Plot prevalence by age group over time
ggplot(malaria_long, aes(x = data_date, y = prevalence, color = age_group)) +
	geom_line(stat = "summary", fun = mean) +
	labs(title = "Malaria Prevalence by Age Group Over Time",
			 x = "Date",
			 y = "Mean Prevalence",
			 color = "Age Group") +
	theme_minimal()
