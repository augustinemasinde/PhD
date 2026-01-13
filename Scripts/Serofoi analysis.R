library(readxl)
library(dplyr)
library(here)
library(serofoi)
chikdata <- read_excel(here("data", "chikungunya_data_Uganda.xlsx"))
chikdata <- chikdata %>%  select(UniqueKey, Year, Age_Yrs,IgM_CHIK)
chikdata <- chikdata %>%
  mutate(
    IgM_CHIK = recode(IgM_CHIK, "Nengative" = "Negative", "NA" = NA_character_)
  ) %>%
  filter(!is.na(IgM_CHIK))
chikdata <- chikdata %>%
  mutate(
    age_cate = case_when(
      Age_Yrs >= 1 & Age_Yrs <= 11 ~ "1-11",
      Age_Yrs >= 12 & Age_Yrs <= 18 ~ "12-18",
      Age_Yrs >= 19 & Age_Yrs <= 49 ~ "19-49",
      Age_Yrs >= 50 & Age_Yrs <= 64 ~ "50-64",
      Age_Yrs >= 65 & Age_Yrs <= 97 ~ "65-97",
      TRUE ~ NA_character_  # catches anything outside range or NA
    )
  )

#chikungunya transmission in 2019
chikdata2022 <- data.frame(
  survey_year = c(2022,2022,2022,2022,2022,2022,2022,2022,2022,2022),
  n_sample = c(195,152,100,98,102,103,65,29,4,1),
  n_seropositive = c(13,13,23,29,36,50,37,21,3,1),
  age_min = c(0,10,20,30,40,50,60,70,80,100),
  age_max = c(9,19,29,39,49,59,69,79,89,110)
)

chikdata2022 <- chikdata2022 %>%
  mutate(n_sample = as.integer(n_sample))



# Implementation of the models
seromodel_constant <- fit_seromodel(
  serosurvey = chikdata2022,
  foi_sigma_rw = sf_normal(0.2,0.2),
  model_type = "constant",
  foi_prior = sf_normal(0.01,1),
  is_seroreversion = TRUE,
  seroreversion_prior = sf_normal(0.006,0.3),
  chains = 2,
  iter = 1000
)



foi_index <- get_foi_index(chikdata2022, group_size = 1, model_type = "time")
seromodel_time <- fit_seromodel(
  serosurvey = chikdata2022,
  model_type = "time",
  foi_prior = sf_normal(0.05, 1),
  foi_index = foi_index,
  iter = 2500,
  is_seroreversion = TRUE
)


foi_index <- get_foi_index(chikdata2022, group_size = 1, model_type = "time")
seromodel_log_time <- fit_seromodel(
  serosurvey = chikdata2022,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  is_log_foi = TRUE,
  foi_index = foi_index,
  iter = 2000,
  is_seroreversion = TRUE
)


# Visualisation of the results
plot_constant<-plot_seromodel(
  seromodel = seromodel_constant,
  serosurvey = chikdata2022,
  foi_max = 0.07,
  size_text = 6
)
plot_time<-plot_seromodel(
  seromodel = seromodel_time,
  serosurvey = chikdata2022,
  foi_max = 0.07,
  size_text = 6
)
plot_log_time <- plot_seromodel(
  seromodel = seromodel_log_time,
  serosurvey = chikdata2022,
  foi_max = 0.07,
  size_text = 6
)

cowplot::plot_grid(plot_constant, plot_time, plot_log_time, ncol = 3)







# Implementation of the models
seromodel_constant <- fit_seromodel(
  serosurvey = chikdata2019,
  model_type = "constant",
  iter = 1000
)



foi_index <- get_foi_index(chikdata2019, group_size = 5, model_type = "time")
seromodel_time <- fit_seromodel(
  serosurvey = chikdata2019,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  foi_index = foi_index,
  iter = 2500
)


foi_index <- get_foi_index(chikdata2019, group_size = 5, model_type = "time")
seromodel_log_time <- fit_seromodel(
  serosurvey = chikdata2019,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  is_log_foi = TRUE,
  foi_index = foi_index,
  iter = 2000
)


# Visualisation of the results
plot_constant <- plot_seromodel(
  seromodel = seromodel_constant,
  serosurvey = chikdata2019,
  foi_max = 0.07,
  size_text = 6
)
plot_time <- plot_seromodel(
  seromodel = seromodel_time,
  serosurvey = chikdata2019,
  foi_max = 0.07,
  size_text = 6
)
plot_log_time <- plot_seromodel(
  seromodel = seromodel_log_time,
  serosurvey = chikdata2019,
  foi_max = 0.07,
  size_text = 6
)

cowplot::plot_grid(plot_constant, plot_time, plot_log_time, ncol = 3)



#chikungunya transmission in 2020
chikdata2020 <- data.frame(
  survey_year = c(2020,2020,2020,2020,2020),
  n_sample = c(77, 31,292, 49,21),
  n_seropositive = c(17,9,57,13,2),
  age_min = c(1,12,19,50,65),
  age_max = c(11,18,49,64,97)
)

chikdata2020 <- chikdata2020 %>%
  mutate(n_sample = as.integer(n_sample))

# Implementation of the models
seromodel_constant <- fit_seromodel(
  serosurvey = chikdata2020,
  model_type = "constant",
  iter = 1000
)



foi_index <- get_foi_index(chikdata2020, group_size = 5, model_type = "time")
seromodel_time <- fit_seromodel(
  serosurvey = chikdata2020,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  foi_index = foi_index,
  iter = 2500
)


foi_index <- get_foi_index(chikdata2020, group_size = 5, model_type = "time")
seromodel_log_time <- fit_seromodel(
  serosurvey = chikdata2020,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  is_log_foi = TRUE,
  foi_index = foi_index,
  iter = 2000
)


# Visualisation of the results
plot_constant <- plot_seromodel(
  seromodel = seromodel_constant,
  serosurvey = chikdata2020,
  foi_max = 0.07,
  size_text = 6
)
plot_time <- plot_seromodel(
  seromodel = seromodel_time,
  serosurvey = chikdata2020,
  foi_max = 0.07,
  size_text = 6
)
plot_log_time <- plot_seromodel(
  seromodel = seromodel_log_time,
  serosurvey = chikdata2020,
  foi_max = 0.07,
  size_text = 6
)

cowplot::plot_grid(plot_constant, plot_time, plot_log_time, ncol = 3)


#chikungunya transmission in 2021
chikdata2021 <- data.frame(
  survey_year = c(2021,2021,2021,2021,2021),
  n_sample = c(51, 42,345, 48,22),
  n_seropositive = c(3,5,24,2,2),
  age_min = c(1,12,19,50,65),
  age_max = c(11,18,49,64,97)
)

chikdata2021 <- chikdata2021 %>%
  mutate(n_sample = as.integer(n_sample))

# Implementation of the models
seromodel_constant <- fit_seromodel(
  serosurvey = chikdata2021,
  model_type = "constant",
  iter = 1000
)



foi_index <- get_foi_index(chikdata2021, group_size = 5, model_type = "time")
seromodel_time <- fit_seromodel(
  serosurvey = chikdata2021,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  foi_index = foi_index,
  iter = 2500
)


foi_index <- get_foi_index(chikdata2021, group_size = 5, model_type = "time")
seromodel_log_time <- fit_seromodel(
  serosurvey = chikdata2021,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  is_log_foi = TRUE,
  foi_index = foi_index,
  iter = 2000
)


# Visualisation of the results
plot_constant <- plot_seromodel(
  seromodel = seromodel_constant,
  serosurvey = chikdata2021,
  foi_max = 0.07,
  size_text = 6
)
plot_time <- plot_seromodel(
  seromodel = seromodel_time,
  serosurvey = chikdata2021,
  foi_max = 0.07,
  size_text = 6
)
plot_log_time <- plot_seromodel(
  seromodel = seromodel_log_time,
  serosurvey = chikdata2021,
  foi_max = 0.07,
  size_text = 6
)

cowplot::plot_grid(plot_constant, plot_time, plot_log_time, ncol = 3)


#chikungunya transmission in 2022
chikdata2022 <- data.frame(
  survey_year = c(2022,2022,2022,2022,2022),
  n_sample = c(96, 195,306, 202,50),
  n_seropositive = c(6,18,68,102,32),
  age_min = c(0,5,15,45,65),
  age_max = c(4,14,44,64,101)
)

chikdata2022 <- chikdata2022 %>%
  mutate(n_sample = as.integer(n_sample))

# Implementation of the models
seromodel_constant <- fit_seromodel(
  serosurvey = chikdata2022,
  foi_sigma_rw = sf_normal(0.2,0.2),
  model_type = "constant",
  foi_prior = sf_normal(0.01,1),
  is_seroreversion = TRUE,
  seroreversion_prior = sf_normal(0.006,0.3),
  chains = 2,
  iter = 1000
)



foi_index <- get_foi_index(chikdata2022, group_size = 5, model_type = "time")
seromodel_time <- fit_seromodel(
  serosurvey = chikdata2022,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  foi_index = foi_index,
  iter = 2500,
  is_seroreversion = TRUE
)


foi_index <- get_foi_index(chikdata2022, group_size = 5, model_type = "time")
seromodel_log_time <- fit_seromodel(
  serosurvey = chikdata2022,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  is_log_foi = TRUE,
  foi_index = foi_index,
  iter = 2000,
  is_seroreversion = TRUE
)


# Visualisation of the results
plot_seromodel(
  seromodel = seromodel_constant,
  serosurvey = chikdata2022,
  foi_max = 0.07,
  size_text = 6
)
plot_time<-plot_seromodel(
  seromodel = seromodel_time,
  serosurvey = chikdata2022,
  foi_max = 0.07,
  size_text = 6
)
plot_log_time <- plot_seromodel(
  seromodel = seromodel_log_time,
  serosurvey = chikdata2022,
  foi_max = 0.07,
  size_text = 6
)

cowplot::plot_grid(plot_constant, plot_time, plot_log_time, ncol = 3)


