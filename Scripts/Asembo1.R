library(haven)
library(dplyr)
library(arsenal)
library(Rsero)
library(readxl)
library(rstan)
library(logistf)

# Robust data loader: update data_dir to where your .dta files live (or use here::here("..."))
data_dir <- "/Users/augustinemasinde/Desktop/PhD files/Global Health"  # <- update if needed

data_files <- list(
  Asembodata = file.path(data_dir, "Asembo PBIDS arbo.dta"),
  Kilifidata = file.path(data_dir, "KHDSS arbo.dta"),
  Kiberadata = file.path(data_dir, "Kibera PBIDS arbo.dta"),
  manyattadata = file.path(data_dir, "Manyatta HDSS arbo.dta"),
  NairobiUrbandata = file.path(data_dir, "NUHDSS arbo.dta")
)

# report missing files and stop early
missing_files <- names(data_files)[!vapply(data_files, file.exists, logical(1))]
if(length(missing_files) > 0){
  stop("Missing data files: ", paste(missing_files, collapse = ", "),
       ". Update data_dir or move files to the working directory.")
}

# read files into variables used in the script
Asembodata <- haven::read_dta(data_files$Asembodata)
Kilifidata <- haven::read_dta(data_files$Kilifidata)
Kiberadata <- haven::read_dta(data_files$Kiberadata)
manyattadata <- haven::read_dta(data_files$manyattadata)
NairobiUrbandata <- haven::read_dta(data_files$NairobiUrbandata)

#descriptive statistics
#unique participants in each study site
length(unique(Asembodata$sample_id))
length(unique(manyattadata$sample_id))
length(unique(Kilifidata$fk_study))
length(unique(Kiberadata$sample_id))
length(unique(NairobiUrbandata$fk_study))

#Distribution of sex in each study site
table(Asembodata$sex)
table(manyattadata$sex)
table(Kilifidata$sex)
table(Kiberadata$sex)
table(NairobiUrbandata$sex)


nairobi_data <-Nairobi_clean %>%
  mutate(agecat = case_when(
    ageyrs >= 0 & ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs <15 ~ "5-14",
    ageyrs >= 15 & ageyrs <45 ~ "15-44",
    ageyrs >= 45 & ageyrs <65 ~ "45-64",
    ageyrs >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))


Asembodata<- Asembodata %>% 
  mutate(agecat = case_when(
    ageyrs >= 0 & ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs <15 ~ "5-14",
    ageyrs >= 15 & ageyrs <45 ~ "15-44",
    ageyrs >= 45 & ageyrs <65 ~ "45-64",
    ageyrs >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))
Asembodata$agecat <- factor(
  Asembodata$agecat,
  levels = c("<5", "5-14", "15-44", "45-64", "65+"),
  ordered = TRUE
)

table(Asembodata$agecat)

manyattadata<- manyattadata %>% 
  mutate(agecat = case_when(
    ageyrs >= 0 & ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs <15 ~ "5-14",
    ageyrs >= 15 & ageyrs <45 ~ "15-44",
    ageyrs >= 45 & ageyrs <65 ~ "45-64",
    ageyrs >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))
manyattadata$agecat <- factor(
  manyattadata$agecat,
  levels = c("<5", "5-14", "15-44", "45-64", "65+"),
  ordered = TRUE
)

table(manyattadata$agecat)

Kilifidata<- Kilifidata %>% 
  mutate(agecat = case_when(
    age_y >= 0 & age_y < 5 ~ "<5",
    age_y >= 5 & age_y <15 ~ "5-14",
    age_y >= 15 & age_y <45 ~ "15-44",
    age_y >= 45 & age_y <65 ~ "45-64",
    age_y >= 65 ~ "65+",
    TRUE ~ NA_character_
    
  ))
Kilifidata$agecat <- factor(
  Kilifidata$agecat,
  levels = c("<5", "5-14", "15-44", "45-64", "65+"),
  ordered = TRUE
)
table(Kilifidata$agecat)


Kiberadata<- Kiberadata %>%
  mutate(agecat = case_when(
    ageyrs >= 0 & ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs <15 ~ "5-14",
    ageyrs >= 15 & ageyrs <45 ~ "15-44",
    ageyrs >= 45 & ageyrs <65 ~ "45-64",
    ageyrs >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))

Kiberadata$agecat <- factor(
  Kiberadata$agecat,
  levels = c("<5", "5-14", "15-44", "45-64", "65+"),
  ordered = TRUE
)
table(Kiberadata$agecat)

NairobiUrbandata<- NairobiUrbandata %>%
  mutate(agecat = case_when(
    age_y >= 0 & age_y < 5 ~ "<5",
    age_y >= 5 & age_y <15 ~ "5-14",
    age_y >= 15 & age_y <45 ~ "15-44",
    age_y >= 45 & age_y <65 ~ "45-64",
    age_y >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))
NairobiUrbandata$agecat <- factor(
  NairobiUrbandata$agecat,
  levels = c("<5", "5-14", "15-44", "45-64", "65+"),
  ordered = TRUE
)
table(NairobiUrbandata$agecat)


library(dplyr)

# Asembo
Asembodata_clean <- Asembodata %>%
  select(site, sex, ageyrs, agecat, CHKpos, DENpos, RVFpos) 

# Manyatta
manyattadata_clean <- manyattadata %>%
  select(site, sex, ageyrs, agecat, CHKpos, DENpos, RVFpos)

# Kilifi
Kilifidata_clean <- Kilifidata %>%
  rename(ageyrs = age_y,
         CHKpos = chke1_pos,
         DENpos = denv_pos,
         RVFpos = rvfgc_pos,
         site = location) %>%
  select(site, sex, ageyrs, agecat, CHKpos, DENpos, RVFpos) 

# Kibera
Kiberadata_clean <- Kiberadata %>%
  rename(CHKpos = chke1_pos,
         DENpos = denv_pos,
         RVFpos = rvfgc_pos) %>%
  select(site, sex, ageyrs, agecat, CHKpos, DENpos, RVFpos)

# Nairobi Urban
Nairobi_clean <- NairobiUrbandata %>%
  rename(ageyrs = age_y,
         CHKpos = chke1_pos,
         DENpos = denv_pos,
         RVFpos = rvfgc_pos,
         site = location) %>%
  select(site, sex, ageyrs, agecat, CHKpos, DENpos, RVFpos)


Nairobi_clean$agecat2 <- as.character(Nairobi_clean$agecat)
Nairobi_clean$agecat2[Nairobi_clean$agecat2 %in% c("<5", "5-14")] <- "<15"
Nairobi_clean$agecat2 <- factor(
  Nairobi_clean$agecat2,
  levels = c("<15", "15-44", "45-64", "65+")
)
Nairobi_clean$agegroup <- cut(
  Nairobi_clean$ageyrs,
  breaks = c(-Inf, 24, 44, Inf),
  labels = c("<25", "26-44", "45+"),
  right = TRUE
)

Kiberadata_clean$agecat2 <- as.character(Kiberadata_clean$agecat)
Kiberadata_clean$agecat2[Kiberadata_clean$agecat2 %in% c("<5", "5-14")] <- "<15"
Kiberadata_clean$agecat2 <- factor(
  Kiberadata_clean$agecat2,
  levels = c("<15", "15-44", "45-64", "65+")
)


Kiberadata_clean$agegroup <- cut(
  Kiberadata_clean$ageyrs,
  breaks = c(-Inf, 24, 44, Inf),
  labels = c("<25", "26-44", "45+"),
  right = TRUE
)
#kilifi
Kilifidata_clean$agegroup <- cut(
  Kilifidata_clean$ageyrs,
  breaks = c(-Inf, 24, 44, Inf),
  labels = c("<25", "26-44", "45+"),
  right = TRUE
)

#Manyatta
manyattadata_clean$agegroup <- cut(
  manyattadata_clean$ageyrs,
  breaks = c(-Inf, 24, 44, Inf),
  labels = c("<25", "26-44", "45+"),
  right = TRUE
)

#ASEMBO
Asembodata_clean$agegroup <- cut(
  Asembodata_clean$ageyrs,
  breaks = c(-Inf, 24, 44, Inf),
  labels = c("<25", "26-44", "45+"),
  right = TRUE
)

Asembodata_clean$agegroup2 <- cut(
  Asembodata_clean$ageyrs,
  breaks = c(-Inf, 44, Inf),
  labels = c("<45", "45+"),
  right = TRUE
)





Nairobi_clean$agegroup <- factor(Nairobi_clean$agegroup, ordered = FALSE)
Nairobi_clean$sex<-as.factor(Nairobi_clean$sex)
Nairobi_clean$CHKpos<-as.factor(Nairobi_clean$CHKpos)
Nairobi_clean$DENpos<-as.factor(Nairobi_clean$DENpos)
Nairobi_clean$RVFpos<-as.factor(Nairobi_clean$RVFpos)
Nairobi_clean$agecat2 <- factor(Nairobi_clean$agecat2, ordered = FALSE)
#Asembo
Asembodata_clean$sex<-as.factor(Asembodata_clean$sex)
Asembodata_clean$agegroup2<-as.factor(Asembodata_clean$agegroup2)
Asembodata_clean$agegroup<-as.factor(Asembodata_clean$agegroup)
Asembodata_clean$agecat<-as.factor(Asembodata_clean$agecat)
Asembodata_clean$CHKpos<-as.factor(Asembodata_clean$CHKpos)
Asembodata_clean$DENpos<-as.factor(Asembodata_clean$DENpos)
Asembodata_clean$RVFpos<-as.factor(Asembodata_clean$RVFpos)
Asembodata_clean$agecat <- factor(Asembodata_clean$agecat, ordered = FALSE)
#Manyatta
manyattadata_clean$sex<-as.factor(manyattadata_clean$sex)
manyattadata_clean$agegroup<-as.factor(manyattadata_clean$agegroup)
manyattadata_clean$agecat<-as.factor(manyattadata_clean$agecat)
manyattadata_clean$CHKpos<-as.factor(manyattadata_clean$CHKpos)
manyattadata_clean$DENpos<-as.factor(manyattadata_clean$DENpos)
manyattadata_clean$RVFpos<-as.factor(manyattadata_clean$RVFpos)
manyattadata_clean$agecat <- factor(manyattadata_clean$agecat, ordered = FALSE)
#Kilifi
Kilifidata_clean$sex<-as.factor(Kilifidata_clean$sex)
Kilifidata_clean$agegroup<-as.factor(Kilifidata_clean$agegroup)
Kilifidata_clean$CHKpos<-as.factor(Kilifidata_clean$CHKpos)
Kilifidata_clean$DENpos<-as.factor(Kilifidata_clean$DENpos)
Kilifidata_clean$RVFpos<-as.factor(Kilifidata_clean$RVFpos)
Kilifidata_clean$agecat <- factor(Kilifidata_clean$agecat, ordered = FALSE)
#Kibera
Kiberadata_clean$agegroup<- as.factor(Kiberadata_clean$agegroup)
Kiberadata_clean$sex<-as.factor(Kiberadata_clean$sex)
Kiberadata_clean$agecat<-as.factor(Kiberadata_clean$agecat)
Kiberadata_clean$CHKpos<-as.factor(Kiberadata_clean$CHKpos)
Kiberadata_clean$DENpos<-as.factor(Kiberadata_clean$DENpos)
Kiberadata_clean$RVFpos<-as.factor(Kiberadata_clean$RVFpos)
Kiberadata_clean$agecat <- factor(Kiberadata_clean$agecat, ordered = FALSE)

#reference category for sex
Asembodata_clean$sex<- relevel(Asembodata_clean$sex, ref = "M")
Nairobi_clean$sex <- relevel(Nairobi_clean$sex, ref = "Male")
manyattadata_clean$sex <- relevel(manyattadata_clean$sex, ref = "M")
Kilifidata_clean$sex <- relevel(Kilifidata_clean$sex, ref = "m")
Kiberadata_clean$sex <- relevel(Kiberadata_clean$sex, ref = "M")

#Bivariate analysis for Nairobi

model0<- glm(DENpos ~ 1, data = Asembodata_clean, family = binomial) #null model
model1<- glm(DENpos  ~ agegroup2, data = Asembodata_clean, family = binomial) # model with one variable
model2<- glm(DENpos  ~ sex + agegroup2, data = Asembodata_clean, family = binomial) #full model
summary(model2)


#effect size estimates
exp_coef <- exp(coef(model2))
exp_confint <- exp(confint(model2))
data.frame(Estimate = coef(model2), OR = exp_coef, CI_lower = exp_confint[,1], CI_upper = exp_confint[,2])



# checking significance of adding variable 1
anova(model0, model1, test="Chisq") #check if adding variable 2 is significant






df_all <- bind_rows(
  Asembodata_clean,
  manyattadata_clean,
  Kilifidata_clean,
  Kiberadata_clean,
  Nairobi_clean
)
str(df_all)

df_all$site <- as.factor(df_all$site)

df_all <- df_all %>%
  mutate(sex = case_when(
    sex %in% c("M", "m", "Male", "male") ~ "Male",
    sex %in% c("F", "f", "Female", "female") ~ "Female",
    TRUE ~ NA_character_
  )) %>%
  mutate(sex = factor(sex, levels = c("Male", "Female")))
str(df_all)

df_all$sex<- as.factor(df_all$sex)
df_all$age_cat<- as.factor(df_all$age_cat)
df_all$CHKpos<- as.factor(df_all$CHKpos)
df_all$DENpos<- as.factor(df_all$DENpos)
df_all$RVFpos<- as.factor(df_all$RVFpos)
df_all$age_cat <- factor(df_all$age_cat, ordered = FALSE)
table(df_all$age_cat)



#logistic regression for site and seropositivity
site_CHK<- glm(CHKpos ~ site, data = df_all, family = binomial)
summary(site_CHK)
exp(coef(site_CHK))
exp(confint(site_CHK))

site_DEN<- glm(DENpos ~ site, data = df_all, family = binomial)
summary(site_DEN)
exp(coef(site_DEN))
exp(confint(site_DEN))

site_RVF<- glm(RVFpos ~ site, data = df_all, family = binomial)
summary(site_RVF)
exp(coef(site_RVF))
exp(confint(site_RVF))


#logistic regression for age and seropositivity
age_CHK<- glm(CHKpos ~  age_cat, data = df_all, family = binomial)
summary(age_CHK)
exp(coef(age_CHK))
exp(confint(age_CHK))

age_DEN<- glm(DENpos~ age_cat, data = df_all, family = binomial)
summary(age_DEN)
exp(coef(age_DEN))
exp(confint(age_DEN))

age_RVF <- glm(RVFpos~ age_cat, data = df_all, family = binomial)
summary(age_RVF)
exp(coef(age_RVF))
exp(confint(age_RVF))
#logistic regression for sex and seropositivity

sex_CHK<- glm(CHKpos~sex, data = df_all, family = binomial)
exp(coef(sex_CHK))
exp(confint(sex_CHK))

sex_DEN<- glm(DENpos~sex, data = df_all, family = binomial)
exp(coef(sex_DEN))
exp(confint(sex_DEN))

sex_RVF<- glm(RVFpos~sex, data = df_all, family = binomial)
exp(coef(sex_RVF))
exp(confint(sex_RVF))




#Asembo
asembodata<- Asembodata %>%
  select(ageyrs,sex,agecat,CHKpos, DENpos, RVFpos)
  
asembodata$RVFpos <- ifelse(asembodata$RVFpos == 1, TRUE, FALSE)
asembodata$DENpos <- ifelse(asembodata$DENpos == 1, TRUE, FALSE)
asembodata$CHKpos <- ifelse(asembodata$CHKpos == 1, TRUE, FALSE)

asembodata$ageyrs=  as.integer(asembodata$ageyrs)
asembodata <- asembodata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
asembodata$sex<-as.character(asembodata$sex)

asembo_chik <-SeroData(age_at_sampling = asembodata$ageyrs,
                       Y= asembodata$CHKpos,
                       sampling_year = 2022)
seroprevalence(asembo_chik)
seroprevalence.plot(asembo_chik, age_class = 5)

piecewisemodel = FOImodel(type = 'outbreak', prioralpha1 = 0.013, prioralpha2 = 1, seroreversion = 1, priorRho1 = 0.2, priorRho2 = 1, K=1, priorT1 = 20, priorT2 = 10, sensitivity = 0.9, specificity = 0.99)
constantmodel = FOImodel(type = 'constant', priorC1 = 0.013, priorC2 = 1, seroreversion = 1, priorRho1 = 0.2, priorRho2 = 1)
print(piecewisemodel)
print(constantmodel)

piecewisemodel = fit(data = asembo_chik,  model = piecewisemodel, chains=1, iter=5000)
constantmodel = fit(data = asembo_chik,  model = constantmodel, chains=1, iter=5000)

seroprevalence.fit(piecewisemodel, YLIM = 1, age_class = 5)
seroprevalence.fit(constantmodel, YLIM = 1, age_class = 5)
#model comparison
m1 = compute_information_criteria(piecewisemodel)
print(m1) 
m2 = compute_information_criteria(constantmodel)
print(m2) 

plot(piecewisemodel, YLIM =1)
plot_posterior(piecewisemodel)
parameters_credible_intervals(piecewisemodel)
parameters_credible_intervals(constantmodel)
traceplot_Rsero(piecewisemodel)





#manyatta
manyattadata<- manyattadata_clean %>%
  select(ageyrs,sex,agecat,CHKpos, DENpos, RVFpos)
manyattadata$RVFpos <- ifelse(manyattadata$RVFpos == 1, TRUE, FALSE)
manyattadata$DENpos <- ifelse(manyattadata$DENpos == 1, TRUE, FALSE)
manyattadata$CHKpos <- ifelse(manyattadata$CHKpos == 1, TRUE, FALSE)

manyattadata$ageyrs=  as.integer(manyattadata$ageyrs)
manyattadata <- manyattadata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
manyattadata$sex<-as.character(manyattadata$sex)

manyatta_chik <-SeroData(age_at_sampling = manyattadata$ageyrs,
                       Y= manyattadata$CHKpos,
                       sampling_year = 2022)
seroprevalence(manyatta_chik)
seroprevalence.plot(manyatta_chik, age_class = 5)
#fit the model constant model
manyattapiecewisemodel = FOImodel(type = 'constant', priorC1 = 0.012, priorC2 = 1, seroreversion = 1, priorRho1 = 0.2, priorRho2 = 1)
print(manyattapiecewisemodel)

manyattapiecewisemodel = fit(data = manyatta_chik,  model = manyattapiecewisemodel, chains=1, iter=5000)

seroprevalence.fit(manyattapiecewisemodel, YLIM=1, age_class = 5)


plot_posterior(manyattapiecewisemodel)

#Posterior distribution of relevant model parameters

traceplot_Rsero(manyattapiecewisemodel)
plot(manyattapiecewisemodel, YLIM =1)


parameters_credible_intervals(manyattapiecewisemodel)


#kilifi
kilifidata<- Kilifidata_clean %>%
  select(ageyrs,sex,age_cat,CHKpos, DENpos, RVFpos)


kilifidata$RVFpos <- ifelse(kilifidata$RVFpos == 1, TRUE, FALSE)
kilifidata$DENpos <- ifelse(kilifidata$DENpos == 1, TRUE, FALSE)
kilifidata$CHKpos <- ifelse(kilifidata$CHKpos == 1, TRUE, FALSE)

kilifidata$ageyrs=  as.integer(kilifidata$ageyrs)
kilifidata <- kilifidata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
kilifidata$sex<-as.character(kilifidata$sex)

kilifidata_chik <-SeroData(age_at_sampling = kilifidata$ageyrs,
                         Y= kilifidata$CHKpos,
                         sampling_year = 2022)
seroprevalence(kilifidata_chik)
seroprevalence.plot(kilifidata_chik, age_class = 5)
#fit the model constant model
kilifipiecewisemodel = FOImodel(type = 'constant', priorC1 = 0.02,  priorC2 = 1, seroreversion = 1, priorRho1 = 0.2, priorRho2 = 0.5)
kilifipiecewisemodel = fit(data = kilifidata_chik,  model = kilifipiecewisemodel, chains=1, iter=5000)
seroprevalence.fit(kilifipiecewisemodel,age_class = 5)

#Posterior distribution of relevant model parameters
parameters_credible_intervals(kilifipiecewisemodel)

plot(kilifipiecewisemodel, YLIM =1)



#kibera
kiberadata<- Kiberadata_clean %>%
  select(ageyrs,sex,age_cat,CHKpos, DENpos, RVFpos)
kiberadata$RVFpos <- ifelse(kiberadata$RVFpos == 1, TRUE, FALSE)
kiberadata$DENpos <- ifelse(kiberadata$DENpos == 1, TRUE, FALSE)
kiberadata$CHKpos <- ifelse(kiberadata$CHKpos == 1, TRUE, FALSE)

kiberadata$ageyrs=  as.integer(kiberadata$ageyrs)
kiberadata <- kiberadata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
kiberadata$sex<-as.character(kiberadata$sex)

kiberadata_chik <-SeroData(age_at_sampling = kiberadata$ageyrs,
                           Y= kiberadata$RVFpos,
                           sampling_year = 2022)
seroprevalence(kiberadata_chik)
seroprevalence.plot(kiberadata_chik, age_class = 5)

kiberapiecewisemodel = FOImodel(type = 'piecewise', priorC1 = 0.3, priorC2 = 1, seroreversion = 1, priorRho1 = 0.06, priorRho2 = 0, K=2, priorT1 = c(5,50), priorT2 = c(5,5))
kiberapiecewisemodel = fit(data = kiberadata_chik,  model = kiberapiecewisemodel, chains=3, iter=5000)
seroprevalence.fit(kiberapiecewisemodel, YLIM= 1, age_class = 5)

parameters_credible_intervals(kiberapiecewisemodel)

plot(kiberapiecewisemodel, YLIM =1)

#Nairobi Urban

Nairobidata<- Nairobi_clean %>%
  select(ageyrs,sex,age_cat,CHKpos, DENpos, RVFpos)

Nairobidata$RVFpos <- ifelse(Nairobidata$RVFpos == 1, TRUE, FALSE)
Nairobidata$DENpos <- ifelse(Nairobidata$DENpos == 1, TRUE, FALSE)
Nairobidata$CHKpos <- ifelse(Nairobidata$CHKpos == 1, TRUE, FALSE)

Nairobidata$ageyrs=  as.integer(Nairobidata$ageyrs)
Nairobidata <- Nairobidata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
Nairobidata$sex<-as.character(Nairobidata$sex)

Nairobidata_chik <-SeroData(age_at_sampling = Nairobidata$ageyrs,
                           Y= Nairobidata$RVFpos,
                           sampling_year = 2022)
seroprevalence(Nairobidata_chik)
seroprevalence.plot(Nairobidata_chik, age_class = 5)

nairobipiecewisemodel = FOImodel(type = 'piecewise', priorC1 = 0.01, priorC2 = 1, seroreversion = 1, priorRho1 = 0.06, priorRho2 = 0, K=2, priorT1 = c(20,40), priorT2 = c(5,5))
nairobipiecewisemodel = fit(data = Nairobidata_chik,  model = nairobipiecewisemodel, chains=3, iter=5000)
seroprevalence.fit(nairobipiecewisemodel, YLIM=1, age_class = 5)


parameters_credible_intervals(nairobipiecewisemodel)
plot(nairobipiecewisemodel, YLIM =1)



#combined data
#estimating seroprevalence and visualizing aggregated data
df_all<- df_all %>% mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
df_all$CHKpos <- ifelse(df_all$CHKpos == 1, TRUE, FALSE)
df_all$DENpos <- ifelse(df_all$DENpos == 1, TRUE, FALSE)
df_all$RVFpos <- ifelse(df_all$RVFpos == 1, TRUE, FALSE)
df_all$ageyrs=  as.integer(df_all$ageyrs)

df_chik <-SeroData(age_at_sampling = round(as.numeric(df_all$ageyrs)),
                             Y= as.numeric(df_all$RVFpos),
                             category = as.character(df_all$site),
                            reference.category = "NAIROBI",
                             sampling_year = as.numeric(2022))
seroprevalence(df_chik)
seroprevalence.plot(df_chik, age_class = 5)

#model of force of infections
#constantmodel= FOImodel(type="constant",   priorC1 = 0.013, priorC2 = 1, seroreversion = 1, priorRho1 = 0.2, priorRho2 = 0)
revmodel= FOImodel(type="constant", priorC1 = 0.05, priorC2 = 1, cat_lambda = TRUE, seroreversion = 1, priorRho1 = 0.1, priorRho2 = 0.58)
#fit the model
#constantmodel= fit(data=df_chik, model=constantmodel, chains=3, iter=5000)
revmodel= fit(data=df_chik, model=revmodel, chains=3, iter=5000)


#visualize the best model fit
#seroprevalence.fit(constantmodel, YLIM=1, age_class = 5)
seroprevalence.fit(revmodel, YLIM=1, age_class = 5)




#Plotting the posterior distributions
plot_posterior(constantmodel)
plot_posterior(revmodel)
#Posterior distribution of relevant model parameters
#parameters_credible_intervals(constantmodel)
parameters_credible_intervals(revmodel)


#model comparison
#m1 = compute_information_criteria(constantmodel)
#print(m1)
m2 = compute_information_criteria(revmodel)



#plotting annual force of infection
plot(piecewise2, YLIM =1)
plot(piecewise, YLIM =1)
plot(piecewise1, YLIM =1)


traceplot_Rsero(piecewise2)


#adjusting age classes
SE= seroprevalence.fit(FOIfit.constant,
                       individual_samples = 0,
                       age_class = 5,
                       YLIM = 1,
                       fill.color = "lightblue",
                       line.color = "blue")
print(SE)

#Full model with different risks by site
df_all$site<- as.character(df_all$site)
df_all<- df_all %>%
  mutate(ageyrs = ifelse(ageyrs <=1, 1, ageyrs))
df_all$ageyrs=  as.integer(df_all$ageyrs)
df_all$CHKpos <- ifelse(df_all$CHKpos == 1, TRUE, FALSE)

chik_full <- SeroData(age_at_sampling = df_all$ageyrs,
                      Y= df_all$CHKpos,
                      #category= df_all$site,           
                      #reference.category  = "NAIROBI",
                      sampling_year = 2022)


seroprevalence(chik_full)
seroprevalence.plot(chik_full, age_class = 5)
#fit the model constant model
ConstantModel = FOImodel(type = 'constant', K=2, priorC1 = 0.15,
                         priorC2 = 1, seroreversion = 1,priorRho1 = 0.3, priorRho2 = .01, cat_lambda = TRUE)

FOIfit.constant = fit(data = chik_full,  model = ConstantModel, chains=1, iter=5000)
seroprevalence.fit(FOIfit.constant, YLIM=1, age_class = 5)


#Posterior distribution of relevant model parameters
parameters_credible_intervals(FOIfit.constant)














#Determine co-infections
Asembodata <- Asembodata %>%
  mutate(co_infection = case_when(
    CHKpos == 1 & RVFpos == 1 ~ "CHIKV and RVFV",
    CHKpos == 1 & DENpos == 1 ~ "CHIKV and DENV",
    RVFpos == 1 & DENpos == 1 ~ "RVFV and DENV",
    CHKpos == 1 & RVFpos == 1 & DENpos == 1 ~ "CHIKV, RVFV and DENV",
    TRUE ~ "No co-infection"
  ))
table(Asembodata$co_infection)
ass1<- Asembodata %>%
  filter(CHKpos == 1 & RVFpos == 1 & DENpos == 1) %>% summarise(count = n())

manyattadata <- manyattadata %>%
  mutate(co_infection = case_when(
    CHKpos == 1 & RVFpos == 1 ~ "CHIKV and RVFV",
    CHKpos == 1 & DENpos == 1 ~ "CHIKV and DENV",
    RVFpos == 1 & DENpos == 1 ~ "RVFV and DENV",
    CHKpos == 1 & RVFpos == 1 & DENpos == 1 ~ "CHIKV, RVFV and DENV",
    TRUE ~ "No co-infection"
  ))
table(manyattadata$co_infection)
ass2<- manyattadata %>%
  filter(CHKpos == 1 & RVFpos == 1 & DENpos == 1) %>% summarise(count = n())
Kilifidata <- Kilifidata %>%
  mutate(co_infection = case_when(
    chke1_pos == 1 & rvfgc_pos == 1 ~ "CHIKV and RVFV",
    chke1_pos == 1 & denv_pos == 1 ~ "CHIKV and DENV",
    rvfgc_pos == 1 & denv_pos == 1 ~ "RVFV and DENV",
    chke1_pos == 1 & rvfgc_pos == 1 & denv_pos == 1 ~ "CHIKV, RVFV and DENV",
    TRUE ~ "No co-infection"
  ))
table(Kilifidata$co_infection)
ass3<- Kilifidata %>%
  filter(chke1_pos == 1 & rvfgc_pos == 1 & denv_pos == 1) %>% summarise(count = n())
Kiberadata <- Kiberadata %>%
  mutate(co_infection = case_when(
    chke1_pos == 1 & rvfgc_pos == 1 ~ "CHIKV and RVFV",
    chke1_pos == 1 & denv_pos == 1 ~ "CHIKV and DENV",
    rvfgc_pos == 1 & denv_pos == 1 ~ "RVFV and DENV",
    chke1_pos == 1 & rvfgc_pos == 1 & denv_pos == 1 ~ "CHIKV, RVFV and DENV",
    TRUE ~ "No co-infection"
  ))
table(Kiberadata$co_infection)

NairobiUrbandata <- NairobiUrbandata %>%
  mutate(co_infection = case_when(
    chke1_pos == 1 & rvfgc_pos == 1 ~ "CHIKV and RVFV",
    chke1_pos == 1 & denv_pos == 1 ~ "CHIKV and DENV",
    rvfgc_pos == 1 & denv_pos == 1 ~ "RVFV and DENV",
    chke1_pos == 1 & rvfgc_pos == 1 & denv_pos == 1 ~ "CHIKV, RVFV and DENV",
    TRUE ~ "No co-infection"
  ))
table(NairobiUrbandata$co_infection)

#all the three viruses
Asembodata %>%
  filter(CHKpos == 1 & RVFpos == 1 & DENpos == 1) %>% summarise(count = n())
manyattadata %>% filter(CHKpos == 1 & RVFpos == 1 & DENpos == 1) %>%
  summarise(count = n())
Kilifidata %>% filter(chke1_pos == 1 & rvfgc_pos == 1 & denv_pos == 1) 
Kiberadata %>% filter(chke1_pos == 1 & rvfgc_pos == 1 & denv_pos == 1) %>%
  summarise(count = n())
NairobiUrbandata %>% filter(chke1_pos == 1 & rvfgc_pos == 1 & denv_pos == 1) %>%
  summarise(count = n())
#CHIKV and RVFV
ass<-Kilifidata %>%
  filter(chke1_pos == 1 & rvfgc_pos == 1 & denv_pos == 1)

View(df_all)
#Group co-infections by site and count the number of co-infections by virus
df_all %>%
  filter(CHKpos == 1 & RVFpos == 1) %>%
  group_by(site) %>%
  summarise(count = n())
df_all %>%
  filter(CHKpos == 1 & DENpos == 1) %>%
  group_by(site) %>%
  summarise(count = n())
df_all %>%
  filter(RVFpos == 1 & DENpos == 1) %>%
  group_by(site) %>%
  summarise(count = n())
df_all %>%
  filter(CHKpos == 1 & RVFpos == 1 & DENpos == 1) %>%
  group_by(site) %>%
  summarise(count = n())
#age-specific co-infections
df_all %>%
  filter(CHKpos == 1 & RVFpos == 1) %>%
  group_by(age_cat) %>%
  summarise(count = n())
df_all %>%
  filter(CHKpos == 1 & DENpos == 1) %>%
  group_by(age_cat) %>%
  summarise(count = n())
df_all %>%
  filter(RVFpos == 1 & DENpos == 1) %>%
  group_by(age_cat) %>%
  summarise(count = n())
df_all %>%
  filter(CHKpos == 1 & RVFpos == 1 & DENpos == 1) %>%
  group_by(age_cat) %>%
  summarise(count = n())


df_all <- df_all %>%
  mutate(sex = case_when(
    sex %in% c("M", "m", "Male", "male") ~ "Male",
    sex %in% c("F", "f", "Female", "female") ~ "Female",
    TRUE ~ NA_character_
  )) %>%
  mutate(sex = factor(sex, levels = c("Male", "Female")))
#sex-specific co-infections
df_all %>% 
  filter(CHKpos == 1 & RVFpos == 1) %>%
  group_by(sex) %>%
  summarise(count = n())
df_all %>%
  filter(CHKpos == 1 & DENpos == 1) %>%
  group_by(sex) %>%
  summarise(count = n())
df_all %>%
  filter(RVFpos == 1 & DENpos == 1) %>%
  group_by(sex) %>%
  summarise(count = n())
df_all %>%
  filter(CHKpos == 1 & RVFpos == 1 & DENpos == 1) %>%
         group_by(sex) %>%
         summarise(count = n())






#visualizations
library(ggplot2)
library(tidyverse)
summary_df <- df_all %>%
  group_by(age_cat) %>%
  summarise(
    n = n(),
    prop_positive = mean(RVFpos, na.rm = TRUE),        # TRUE treated as 1, FALSE as 0
    prop_negative = 1 - mean(RVFpos, na.rm = TRUE)
  )
summary_long <- summary_df %>%
  pivot_longer(
    cols = starts_with("prop_"),
    names_to = "status",
    values_to = "proportion"
  ) %>%
  mutate(status = factor(status,
                         levels = c("prop_negative", "prop_positive"),
                         labels = c("Negative", "Positive")))

ggplot(summary_long, aes(x = age_cat, y = proportion, fill = status)) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  geom_text(
    aes(x = age_cat, y = 1.03, label = paste0("n = ", n)),
    data = summary_df,
    size = 4,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c("grey90", "#6baed6"),
    name = "Serostatus"
  ) +
  scale_y_continuous(
    name = "Proportion seropositive",
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0)
  ) +
  labs(
    title = "RVFV Proportion Seropositive by Age Group",
    x = "Age (years)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

nairobi <- data.frame(
  site = "Nairobi Urban HDSS",
  agecat = c("<5","5–14","15–44","45–64","65+"),
  CHIKV = c(2.2, 6.0, 11.5, 11.4, 3.9),
  CHIKV_low = c(0.3, 3.1, 8.1, 7.4, 0.5),
  CHIKV_high = c(7.7, 10.2, 15.6, 16.7, 13.5),
  DENV = c(0.0, 0.0, 1.3, 2.5, 2.0),
  DENV_low = c(0.0, 0.0, 0.4, 0.8, 0.0),
  DENV_high = c(4.0, 1.8, 3.3, 5.7, 10.5),
  RVFV = c(1.1, 2.5, 3.0, 0.5, 0.0),
  RVFV_low = c(0.0, 0.8, 1.4, 0.0, 0.0),
  RVFV_high = c(6.0, 5.7, 5.5, 3.1, 7.0)
)

kibera <- data.frame(
  site = "Kibera PBIDS",
  agecat = c("<5","5–14","15–44","45–64","65+"),
  CHIKV = c(4.6, 10.0, 37.6, 51.6, 40.0),
  CHIKV_low = c(0.6, 6.5, 33.0, 41.1, 5.3),
  CHIKV_high = c(15.5, 14.6, 42.4, 62.0, 85.3),
  DENV = c(0.0, 2.6, 1.6, 1.1, 0.0),
  DENV_low = c(0.0, 0.9, 0.6, 0.0, 0.0),
  DENV_high = c(8.0, 5.6, 3.3, 5.7, 52.2),
  RVFV = c(6.8, 2.2, 3.0, 3.2, 20.0),
  RVFV_low = c(1.4, 0.7, 0.9, 0.7, 0.5),
  RVFV_high = c(18.7, 5.0, 9.2, 9.0, 71.6)
)

manyatta <- data.frame(
  site = "Manyatta HDSS",
  agecat = c("<5","5–14","15–44","45–64","65+"),
  CHIKV = c(4.3, 15.4, 34.4, 39.2, 45.5),
  CHIKV_low = c(0.5, 10.6, 30.1, 28.0, 18.6),
  CHIKV_high = c(14.5, 21.2, 38.9, 51.2, 76.6),
  DENV = c(0.0, 1.0, 3.2, 2.7, 0.0),
  DENV_low = c(0.0, 0.1, 1.8, 0.3, 0.0),
  DENV_high = c(7.6, 3.7, 5.2, 9.4, 28.5),
  RVFV = c(0.0, 1.5, 1.5, 2.7, 0.0),
  RVFV_low = c(0.0, 0.3, 0.6, 0.3, 0.0),
  RVFV_high = c(7.6, 4.4, 3.0, 9.4, 28.5)
)

kilifi <- data.frame(
  site = "Kilifi HDSS",
  agecat = c("<5","5–14","15–44","45–64","65+"),
  CHIKV = c(6.3, 9.5, 22.2, 50.5, 64.0),
  CHIKV_low = c(3.2, 6.5, 17.7, 43.4, 49.2),
  CHIKV_high = c(13.1, 14.2, 27.3, 57.6, 77.1),
  DENV = c(7.3, 14.9, 22.6, 70.8, 88.0),
  DENV_low = c(3.0, 10.2, 14.4, 64.0, 75.7),
  DENV_high = c(14.4, 20.7, 30.8, 83.9, 95.5),
  RVFV = c(1.0, 3.6, 3.3, 5.5, 0.0),
  RVFV_low = c(0.5, 1.5, 1.6, 2.8, 0.0),
  RVFV_high = c(5.7, 7.3, 5.9, 9.5, 7.1)
)

asembo <- data.frame(
  site = "Asembo PBIDS",
  agecat = c("<5","5–14","15–44","45–64","65+"),
  CHIKV = c(56.9, 74.1, 70.0, 73.5, 78.3),
  CHIKV_low = c(44.7, 68.5, 64.7, 64.3, 65.8),
  CHIKV_high = c(68.6, 79.1, 74.9, 81.3, 87.9),
  DENV = c(0.0, 0.0, 0.6, 0.9, 5.0),
  DENV_low = c(0.0, 0.0, 0.1, 0.0, 1.0),
  DENV_high = c(5.0, 1.3, 2.2, 4.8, 13.9),
  RVFV = c(0.0, 3.6, 6.7, 5.3, 0.0),
  RVFV_low = c(0.0, 1.7, 4.2, 2.0, 0.0),
  RVFV_high = c(5.0, 6.5, 9.9, 11.2, 6.0)
)

### --- Combine all sites --- ###

sero_all <- bind_rows(nairobi, kibera, manyatta, kilifi, asembo)

# Tidy long format
sero_long <- sero_all %>%
  pivot_longer(cols = c(CHIKV, DENV, RVFV),
               names_to = "virus", values_to = "seroprev") %>%
  left_join(
    sero_all %>%
      pivot_longer(cols = c(CHIKV_low, DENV_low, RVFV_low),
                   names_to = "virus_low", values_to = "low") %>%
      mutate(virus = gsub("_low", "", virus_low)) %>%
      select(site, agecat, virus, low),
    by = c("site", "agecat", "virus")
  ) %>%
  left_join(
    sero_all %>%
      pivot_longer(cols = c(CHIKV_high, DENV_high, RVFV_high),
                   names_to = "virus_high", values_to = "high") %>%
      mutate(virus = gsub("_high", "", virus_high)) %>%
      select(site, agecat, virus, high),
    by = c("site", "agecat", "virus")
  ) %>%
  mutate(
    virus = factor(virus, levels = c("CHIKV", "DENV", "RVFV")),
    agecat = factor(gsub("–", "-", agecat),
                    levels = c("<5","5-14","15-44","45-64","65+"),
                    ordered = TRUE)
  )

### --- Define reusable plotting function --- ###

plot_virus <- function(virus_name) {
  pd <- position_dodge(width = 0.8)
  
  ggplot(filter(sero_long, virus == virus_name),
         aes(x = agecat, y = seroprev, fill = site)) +
    geom_col(position = pd, width = 0.7) +
    geom_errorbar(aes(ymin = low, ymax = high), position = pd, width = 0.2) +
    geom_text(aes(label = sprintf("%.1f", seroprev)),
              position = position_dodge(width = 0.8),
              vjust = -0.6, size = 3) +
    labs(
      title = paste0("Age- and Site-Specific Crude Seroprevalence of ", virus_name),
      x = "Age Group",
      y = "Crude Seroprevalence (%)",
      fill = "Site"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
    )
}

### --- Generate individual virus plots --- ###

plot_virus("CHIKV")
plot_virus("DENV")
plot_virus("RVFV")



       







#Bayesian population weighting
#asembo population data
file_path <- "/Users/augustinemasinde/Desktop/PhD files/Global Health/PBIDS_sites_Ase_Man_Kib_midyr_popn_2022.xlsx"

pop_data <- read_excel(
  path = file_path,
  sheet = 1,
  col_types = "text"
)
pop_data<- pop_data %>%
  select(Ageband,Sex,N_mid2022)
pop_data <- pop_data %>% 
  filter(if_any(everything(), ~ !is.na(.) & . != "")) %>% 
  filter(!grepl("Total", .[[1]], ignore.case = TRUE))
#create age categories
library(dplyr)

asembo_pop_data <- pop_data %>%
  mutate(age_group = case_when(
    Ageband %in% c("0–4") ~ "<5",
    Ageband %in% c("5–9", "10–14") ~ "5-14",
    Ageband %in% c("15–24", "25–34", "35–44") ~ "15-44",
    Ageband %in% c("45–54", "55–64") ~ "45-64",
    Ageband == "65+" ~ "65+",
    TRUE ~ NA_character_
  ))
asembo_pop_data$site<- "ASEMBO"

asembo_pop_data <- asembo_pop_data %>%
  mutate(N_mid2022 = as.numeric(N_mid2022))

asembo_pop_data <- asembo_pop_data %>%
  group_by(Sex, age_group) %>%
  summarise(N_mid2022 = sum(N_mid2022, na.rm = TRUE), .groups = "drop")


#manyatta population data
file_path<-'/Users/augustinemasinde/Desktop/PhD files/Global Health/manyatta.xlsx'
manyatta_pop_data <- read_excel(
  path = file_path,
  sheet = 1,
  col_types = "text"
)
manyatta_pop_data<- manyatta_pop_data %>% select(Ageband, Sex, N_mid2022)
manyatta_pop_data <- manyatta_pop_data %>% 
  filter(if_any(everything(), ~ !is.na(.) & . != "")) %>% 
  filter(!grepl("Total", .[[1]], ignore.case = TRUE))
#create age categories
manyatta_pop_data <- manyatta_pop_data %>%
  mutate(age_group = case_when(
    Ageband %in% c("0–4") ~ "<5",
    Ageband %in% c("5–9", "10–14") ~ "5-14",
    Ageband %in% c("15–24", "25–34", "35–44") ~ "15-44",
    Ageband %in% c("45–54", "55–64") ~ "45-64",
    Ageband == "65+" ~ "65+",
    TRUE ~ NA_character_
  ))
manyatta_pop_data <- manyatta_pop_data %>%
  mutate(N_mid2022 = as.numeric(N_mid2022))
manyatta_pop_data<-manyatta_pop_data %>% select(Sex, age_group, N_mid2022)

manyatta_pop_data <- manyatta_pop_data %>%
  group_by(Sex, age_group) %>%
  summarise(N_mid2022 = sum(N_mid2022, na.rm = TRUE), .groups = "drop")

#kibera population data
file_path<-'/Users/augustinemasinde/Desktop/PhD files/Global Health/kibera.xlsx'
kibera_pop_data <- read_excel(
  path = file_path,
  sheet = 1,
  col_types = "text"
)
kibera_pop_data<- kibera_pop_data %>% select(Ageband, Sex, N_mid2022)
kibera_pop_data <- kibera_pop_data %>% 
  filter(if_any(everything(), ~ !is.na(.) & . != "")) %>% 
  filter(!grepl("Total", .[[1]], ignore.case = TRUE))
#create age categories
kibera_pop_data <- kibera_pop_data %>%
  mutate(age_group = case_when(
    Ageband %in% c("0–4") ~ "<5",
    Ageband %in% c("5–9", "10–14") ~ "5-14",
    Ageband %in% c("15–24", "25–34", "35–44") ~ "15-44",
    Ageband %in% c("45–54", "55–64") ~ "45-64",
    Ageband == "65+" ~ "65+",
    TRUE ~ NA_character_
  ))

kibera_pop_data <- kibera_pop_data %>%
  mutate(N_mid2022 = as.numeric(N_mid2022))
kibera_pop_data<-kibera_pop_data %>% select(Sex, age_group, N_mid2022)

kibera_pop_data <- kibera_pop_data %>%
  group_by(Sex, age_group) %>%
  summarise(N_mid2022 = sum(N_mid2022, na.rm = TRUE), .groups = "drop")


#nairobi population data
nairobi_pop_data <- read_dta('/Users/augustinemasinde/Desktop/PhD files/Global Health/stdized_pop_by_sex_loc_nrb.dta')
nairobi_pop_data<- nairobi_pop_data %>% rename(site = location)
#create ageband from age_y column
nairobi_pop_data <- nairobi_pop_data %>%
  mutate(age_group = case_when(
    age_y >= 0 & age_y < 5 ~ "<5",
    age_y >= 5 & age_y <15 ~ "5-14",
    age_y >= 15 & age_y <45 ~ "15-44",
    age_y >= 45 & age_y <65 ~ "45-64",
    age_y >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))
nairobi_pop_data<- nairobi_pop_data %>%
  group_by(sex, age_group) %>%          
  summarise(N_mid2022 = sum(pop), .groups = "drop") %>% arrange(age_group)
nairobi_pop_data$site<- "NAIROBI"
nairobi_pop_data<- nairobi_pop_data %>% rename(Sex = sex)
#Create a duplicate column age_group and rename to Ageband
nairobi_pop_data$Ageband <- nairobi_pop_data$age_group
nairobi_pop_data<- nairobi_pop_data %>% select(Ageband,site,Sex, age_group, N_mid2022)
#convert every thing to character
nairobi_pop_data <- nairobi_pop_data %>%
  mutate(across(everything(), as.character))

#kilifi population data
kilifi_pop_data<-read_dta('/Users/augustinemasinde/Desktop/PhD files/Global Health/stdized_pop_by_sex_loc_klf.dta')
kilifi_pop_data<- kilifi_pop_data %>% rename(site = location)
#Recode sex, 1=Female and 2=Male
kilifi_pop_data <- kilifi_pop_data %>%
  mutate(Sex = recode(as.numeric(sex),   # make sure it’s numeric first
                      `1` = "Female",
                      `2` = "Male"))
kilifi_pop_data <- kilifi_pop_data %>%
  mutate(pop = ifelse(is.na(pop), 0, pop))



#create ageband from age_y column
kilifi_pop_data <- kilifi_pop_data %>%
  mutate(age_group = case_when(
    age_y >= 0 & age_y < 5 ~ "<5",
    age_y >= 5 & age_y <15 ~ "5-14",
    age_y >= 15 & age_y <45 ~ "15-44",
    age_y >= 45 & age_y <65 ~ "45-64",
    age_y >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))

kilifi_pop_data<- kilifi_pop_data %>% group_by(Sex, age_group) %>%          
  summarise(N_mid2022 = sum(pop), .groups = "drop")




#Method : Using Stan
# Nairobi Urban
nrb_set <- nairobi_pop_data %>% select(Ageband,Sex, N_mid2022, age_group)
nrb_set <- nrb_set %>% rename(agecat = age_group, pop = N_mid2022,sex = Sex)
nrb_set$pop<- as.numeric(nrb_set$pop)
nrb_set<-nrb_set %>% mutate(pw = pop/sum(pop)) %>% select(agecat, sex, pop, pw)
nrb_set <- nrb_set %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)

nairobi_data <-Nairobi_clean %>%
  mutate(agecat = case_when(
    ageyrs >= 0 & ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs <15 ~ "5-14",
    ageyrs >= 15 & ageyrs <45 ~ "15-44",
    ageyrs >= 45 & ageyrs <65 ~ "45-64",
    ageyrs >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))
nairobi_data <- nairobi_data %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)

nrb_counts <- 
  nairobi_data %>%
  group_by(sex, agecat, .drop = FALSE) %>%
  summarise(y = sum(RVFpos), n = n()) %>%
  mutate(n=ifelse(is.na(y),1,n), # replacing values in the n strata with no data with 0
         y=ifelse(is.na(y),0,y)) %>% select(agecat, sex, y, n)

nairobi_data <- nrb_set %>% 
  left_join(nrb_counts, by = c("sex", "agecat"))


nrb_pw <- 
  nairobi_data %>%
  group_by(sex, agecat,.drop = FALSE) %>%
  summarise(pw = mean(pw))
y <- array(nrb_counts$y, dim = c(5, 2))
n <- array(nrb_counts$n, dim = c(5, 2))
pw <- array(nrb_pw$pw, dim = c(5, 2))
tot_pw_age <- rowSums(pw, dims = 1)
tot_pw_sex <- colSums(pw)

#Asembo
asembo_set <- asembo_pop_data %>% select(age_group,Sex, N_mid2022)
asembo_set <- asembo_set %>% rename(agecat = age_group,sex = Sex, pop = N_mid2022)
asembo_set$pop<- as.numeric(asembo_set$pop)
asembo_set <- asembo_set %>%
  mutate(sex = recode(sex,
                      "Female" = "F",
                      "Male"   = "M"))
asembo_set<-asembo_set %>% mutate(pw = pop/sum(pop)) %>% select(agecat, sex, pop, pw)

asembo_set <- asembo_set %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)

asembo_data <-Asembodata %>%
  mutate(agecat = case_when(
    ageyrs >= 0 & ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs <15 ~ "5-14",
    ageyrs >= 15 & ageyrs <45 ~ "15-44",
    ageyrs >= 45 & ageyrs <65 ~ "45-64",
    ageyrs >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))

asembo_data<- asembo_data %>% select(agecat, sex, CHKpos,DENpos,RVFpos)


asembo_counts <- 
  asembo_data %>%
  group_by(sex, agecat, .drop = FALSE) %>%
  summarise(y = sum(RVFpos), n = n()) %>%
  mutate(n=ifelse(is.na(y),1,n), # replacing values in the n strata with no data with 0
         y=ifelse(is.na(y),0,y))
asembo_data <- asembo_counts %>% 
  left_join(asembo_set, by = c("sex", "agecat"))


asembo_pw <- 
  asembo_data %>%
  group_by(sex, agecat,.drop = FALSE) %>%
  summarise(pw = mean(pw))
y <- array(asembo_counts$y, dim = c(5, 2))
n <- array(asembo_counts$n, dim = c(5, 2))
pw <- array(asembo_pw$pw, dim = c(5, 2))
tot_pw_age <- rowSums(pw, dims = 1)
tot_pw_sex <- colSums(pw)

#Kilifi 
kilifi_set <- kilifi_pop_data %>% select(age_group,Sex, N_mid2022)
kilifi_set <- kilifi_set %>% rename(agecat = age_group,sex = Sex, pop = N_mid2022)
kilifi_set$pop<- as.numeric(kilifi_set$pop)
kilifi_set <- kilifi_set %>%
  mutate(sex = recode(sex,
                      "Female" = "f",
                      "Male"   = "m"))
kilifi_set<-kilifi_set %>% mutate(pw = pop/sum(pop)) %>% select(agecat, sex, pop, pw)

kilifi_set <- kilifi_set %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)


kilifi_data <-Kilifidata %>%
  mutate(agecat = case_when(
    age_y >= 0 & age_y < 5 ~ "<5",
    age_y >= 5 & age_y <15 ~ "5-14",
    age_y >= 15 & age_y<45 ~ "15-44",
    age_y >= 45 & age_y <65 ~ "45-64",
    age_y >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))

kilifi_data <- kilifi_data %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)


kilifi_data<- kilifi_data %>% select(agecat, sex, chke1_pos,denv_pos,rvfgc_pos)

kilifi_counts <- 
  kilifi_data %>%
  group_by(sex, agecat, .drop = FALSE) %>%
  summarise(y = sum(rvfgc_pos), n = n()) %>%
  mutate(n=ifelse(is.na(y),1,n), # replacing values in the n strata with no data with 0
         y=ifelse(is.na(y),0,y))




kilifi_data <- kilifi_counts %>% 
  left_join(kilifi_set, by = c("sex", "agecat"))






kilifi_pw <- 
  kilifi_data %>%
  group_by(sex, agecat,.drop = FALSE) %>%
  summarise(pw = mean(pw))


y <- array(kilifi_counts$y, dim = c(5, 2))
n <- array(kilifi_counts$n, dim = c(5, 2))
pw <- array(kilifi_pw$pw, dim = c(5, 2))
tot_pw_age <- rowSums(pw, dims = 1)
tot_pw_sex <- colSums(pw)





#Manyatta
manyatta_set <- manyatta_pop_data %>% select(age_group,Sex, N_mid2022)
manyatta_set <- manyatta_set %>% rename(agecat = age_group,sex = Sex, pop = N_mid2022)
manyatta_set$pop<- as.numeric(manyatta_set$pop)
manyatta_set<-manyatta_set %>% mutate(pw = pop/sum(pop)) %>% select(agecat, sex, pop, pw)
manyatta_set <- manyatta_set %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)



manyatta_data <-manyattadata_clean %>%
  mutate(agecat = case_when(
    ageyrs >= 0 & ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs <15 ~ "5-14",
    ageyrs >= 15 & ageyrs <45 ~ "15-44",
    ageyrs >= 45 & ageyrs <65 ~ "45-64",
    ageyrs >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))
manyatta_data <- manyatta_data %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)

manyatta_counts <- 
  manyatta_data %>%
  group_by(sex, agecat, .drop = FALSE) %>%
  summarise(y = sum(RVFpos), n = n()) %>%
  mutate(n=ifelse(is.na(y),1,n), # replacing values in the n strata with no data with 0
         y=ifelse(is.na(y),0,y)) %>% select(agecat, sex, y, n)


manyatta_counts <- manyatta_counts %>%
  mutate(
    sex = ifelse(sex == "F", "Female",
                 ifelse(sex == "M", "Male", sex))
  )



manyatta_data <- manyatta_counts %>% 
  left_join(manyatta_set, by = c("sex", "agecat"))

manyatta_data$sex<-as.factor(manyatta_data$sex)
manyatta_data$agecat<-as.factor(manyatta_data$agecat)


manyatta_pw <- 
  manyatta_data %>%
  group_by(sex, agecat,.drop = FALSE) %>%
  summarise(pw = mean(pw))
y <- array(manyatta_counts$y, dim = c(5, 2))
n <- array(manyatta_counts$n, dim = c(5, 2))
pw <- array(manyatta_pw$pw, dim = c(5, 2))
tot_pw_age <- rowSums(pw, dims = 1)
tot_pw_sex <- colSums(pw)

#Kibera
kibera_set <- kibera_pop_data %>% select(age_group,Sex, N_mid2022)
kibera_set <- kibera_set %>% rename(agecat = age_group,sex = Sex, pop = N_mid2022)
kibera_set$pop<- as.numeric(kibera_set$pop)
kibera_set <- kibera_set %>%
  mutate(sex = recode(sex,
                      "Female" = "F",
                      "Male"   = "M"))
kibera_set<-kibera_set %>% mutate(pw = pop/sum(pop)) %>% select(agecat, sex, pop, pw)
kibera_set <- kibera_set %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)

kibera_data <-Kiberadata %>%
  mutate(agecat = case_when(
    ageyrs >= 0 & ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs <15 ~ "5-14",
    ageyrs >= 15 & ageyrs <45 ~ "15-44",
    ageyrs >= 45 & ageyrs <65 ~ "45-64",
    ageyrs >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))


kibera_data<- kibera_data %>% select(agecat, sex, chke1_pos,denv_pos,rvfgc_pos)
kibera_data <- kibera_data %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)

kibera_counts <- 
  kibera_data %>%
  group_by(sex, agecat, .drop = FALSE) %>%
  summarise(y = sum(rvfgc_pos), n = n()) %>%
  mutate(n=ifelse(is.na(y),1,n), 
         y=ifelse(is.na(y),0,y))


kibera_counts$sex<-as.factor(kibera_counts$sex)
kibera_counts$agecat<-as.factor(kibera_counts$agecat)

kibera_data <-kibera_counts %>% 
  left_join(kibera_set, by = c("sex", "agecat"))



kibera_pw <- 
  kibera_data %>%
  group_by(sex, agecat,.drop = FALSE) %>%
  summarise(pw = mean(pw))
y <- array(kibera_counts$y, dim = c(5, 2))
n <- array(kibera_counts$n, dim = c(5, 2))
pw <- array(kibera_pw$pw, dim = c(5, 2))
tot_pw_age <- rowSums(pw, dims = 1)
tot_pw_sex <- colSums(pw)
###########################################################################################################
#combined data
nrb_set <- nairobi_pop_data %>% select(Ageband,Sex, N_mid2022, age_group)
nrb_set <- nrb_set %>% rename(agecat = age_group, pop = N_mid2022,sex = Sex)
nrb_set$pop<- as.numeric(nrb_set$pop)
nrb_set$site<- "NAIROBI"
nrb_set<- nrb_set %>% select(agecat,sex, pop, site)

asembo_set <- asembo_pop_data %>% select(age_group,Sex, N_mid2022)
asembo_set <- asembo_set %>% rename(agecat = age_group,sex = Sex, pop = N_mid2022)
asembo_set$pop<- as.numeric(asembo_set$pop)
asembo_set$site<- "ASEMBO"

kilifi_set <- kilifi_pop_data %>% select(age_group,Sex, N_mid2022)
kilifi_set <- kilifi_set %>% rename(agecat = age_group,sex = Sex, pop = N_mid2022)
kilifi_set$pop<- as.numeric(kilifi_set$pop)
kilifi_set$site<- "KILIFI"

manyatta_set <- manyatta_pop_data %>% select(age_group,Sex, N_mid2022)
manyatta_set <- manyatta_set %>% rename(agecat = age_group,sex = Sex, pop = N_mid2022)
manyatta_set$pop<- as.numeric(manyatta_set$pop)
manyatta_set$site<- "MANYATTA"

kibera_set <- kibera_pop_data %>% select(age_group,Sex, N_mid2022)
kibera_set <- kibera_set %>% rename(agecat = age_group,sex = Sex, pop = N_mid2022)
kibera_set$pop<- as.numeric(kibera_set$pop)
kibera_set$site<- "KIBERA"

df_set <- rbind(asembo_set, manyatta_set, kilifi_set, kibera_set, nrb_set)


df_set<-df_set %>% mutate(pw = pop/sum(pop)) %>% select(site,agecat, sex, pop, pw)
df_set <- df_set %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)
df_set$sex<- as.factor(df_set$sex)
df_set<- df_set %>% select(site, agecat, sex, pop, pw)

dt_set <- df_set %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)
df_set$site<-as.factor(df_set$site)
df_set$sex<- as.factor(df_set$sex)


df_all<- df_all %>% mutate(site = case_when(
  site == "Asembo" ~ "ASEMBO",
  site == "Manyatta" ~ "MANYATTA",
  site == "Kilifi" ~ "KILIFI",
  site == "Kibera" ~ "KIBERA",
  site == "Nairobi Urban" ~ "NAIROBI",
  TRUE ~ site
))

df_data <- df_all %>%
  mutate(agecat = case_when(
    ageyrs >= 0 & ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs <15 ~ "5-14",
    ageyrs >= 15 & ageyrs <45 ~ "15-44",
    ageyrs >= 45 & ageyrs <65 ~ "45-64",
    ageyrs >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))
df_data$CHKpos<- as.numeric(df_data$CHKpos)
df_data <- df_data %>% mutate(CHKpos = ifelse(CHKpos == 2, 1, 0))

df_data$DENpos<- as.numeric(df_data$DENpos)
df_data <- df_data %>% mutate(DENpos = ifelse(DENpos == 2, 1, 0))
df_data$RVFpos<- as.numeric(df_data$RVFpos)
df_data <- df_data %>% mutate(RVFpos = ifelse(RVFpos == 2, 1, 0))
df_data<- df_data %>% select(agecat, sex, site, CHKpos,DENpos,RVFpos)


df_data <- df_data %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)


df_counts <- 
  df_data %>%
  group_by(agecat,sex,site, .drop = FALSE) %>%
  summarise(y = sum(DENpos), n = n()) %>%
  mutate(n=ifelse(is.na(y),1,n), 
         y=ifelse(is.na(y),0,y)) 

df_data <-df_counts %>% 
  left_join(df_set, by = c("sex", "agecat", "site"))


################################################
asembo_pw <- 
  asembo_data %>%
  group_by(agecat,sex, .drop = FALSE) %>%
  summarise(pw = mean(pw))
y <- array(asembo_counts$y, dim = c(5, 2))
n <- array(asembo_counts$n, dim = c(5, 2))
pw <- array(asembo_pw$pw, dim = c(5, 2))
tot_pw_age <- rowSums(pw, dims = 1)
tot_pw_sex <- colSums(pw)



mrp_adjusted <-"data {
  int N_se; // denominator sensitivity
  int N_sp; // denominator specificity
  int x; // numerator sensitivity
  int z; // numerator specificity
  int y[5, 2]; // no. seropositives
  int n[5, 2]; // no. samples
  real pw[5, 2]; // proportion of population in each demographic subgroup
  real tot_pw_age[5]; // proportion of population in each age group
  real tot_pw_sex[2]; // proportion female and male
}

parameters {
  real<lower=0,upper=1> se; 
  real<lower=0,upper=1> sp; 
  real bsex[2];  
  real bage[5]; 
  real<lower=0> sd_age;
}

transformed parameters {
  real<lower=0,upper=1> p[5, 2]; 
  real<lower=0,upper=1> p_obs[5, 2]; 
  
  for(a in 1:5){
    for(s in 1:2){
      p[a, s] = inv_logit(bage[a] + 
                            bsex[s]);
      
      p_obs[a, s] = se * p[a, s] + 
        (1 - sp) * (1 - p[a, s]); 
    }
  }
}

model {
  //priors
  se ~ beta(1,1); 
  sp ~ beta(1,1); 
  bsex ~ normal(0,10); 
  bage ~ normal(0, sd_age); 
  sd_age ~ normal(0, 0.5); 
  
  //likelihood
  
  for(a in 1:5){
    for(s in 1:2){ 
      y[a, s] ~ binomial(n[a, s], p_obs[a, s]); 
    }
  }
  x ~ binomial(N_se, se); 
  z ~ binomial(N_sp, sp); 
}

generated quantities {
  real p_site = 0;
  vector[5] p_age = rep_vector(0, 5);
  vector[2] p_sex = rep_vector(0, 2);
  
  
  for(a in 1:5){
    for(s in 1:2){
      p_site += p[a, s] * pw[a, s];
    }
  }
  
  for(a in 1:5){
    for(s in 1:2){
      p_age[a] += p[a, s] * pw[a, s]/tot_pw_age[a];
    }
  }
  
  
  for(s in 1:2){
    for(a in 1:5){
      p_sex[s] += p[a, s] * pw[a, s]/tot_pw_sex[s];
    }
  }
  
  
}"

mrp_model <- stan_model(model_code = mrp_adjusted)

fit <- sampling(
  object = mrp_model,     
  data = list(
    y = y,
    n = n,
    pw = pw,
    tot_pw_age = tot_pw_age,
    tot_pw_sex = tot_pw_sex,
    x = 2,
    z = 139,
    N_se = 2,
    N_sp = 145
  ),
  chains = 3,
  warmup = 1000,
  iter = 10000,
  cores = 3,
  seed = 111,
  pars = c("bage","bsex","sd_age","p_age","p_sex","p_site","se","sp"),
  refresh = 0
)


options(scipen=999)
out_adj <- signif(summary(fit)$summary, 3)
out_adj_1 <- as.data.frame(out_adj)
names <- (dimnames(out_adj)[[1]])[11:20]
out_adj <- out_adj[9:17, c(1, 4, 8)]




#################################################
#overall model
####################################################
df_pw <- df_data %>%
  group_by(agecat, sex, site, .drop = FALSE) %>%
  summarise(pw = mean(pw))


y <- array(df_counts$y, dim = c(5, 2, 5))
n <- array(df_counts$n, dim = c(5, 2, 5))
pw <- array(df_pw$pw, dim = c(5, 2, 5))

tot_pw_age <- apply(pw, 1, sum)
tot_pw_sex <- apply(pw, 2, sum)
tot_pw_region <- apply(pw, 3, sum)



mrp_model_region <- stan_model(model_code = "
data {
  int N_se; // denominator sensitivity
  int N_sp; // denominator specificity
  int x;    // numerator sensitivity
  int z;    // numerator specificity
  int y[5, 2, 5]; // no. seropositives
  int n[5, 2, 5]; // no. samples
  real pw[5, 2, 5]; // proportion of population in each subgroup
  real tot_pw_age[5];   // sum of weights by age
  real tot_pw_sex[2];   // sum of weights by sex
  real tot_pw_region[5];// sum of weights by region
}

parameters {
  real<lower=0,upper=1> se;
  real<lower=0,upper=1> sp;
  real bsex[2];
  real bage[5];
  real bregion[5];
  real<lower=0> sd_age;
  real<lower=0> sd_region;
}

transformed parameters {
  real<lower=0,upper=1> p[5,2,5];
  real<lower=0,upper=1> p_obs[5,2,5];
  for(a in 1:5)
    for(s in 1:2)
      for(r in 1:5){
        p[a,s,r] = inv_logit(bage[a] + bsex[s] + bregion[r]);
        p_obs[a,s,r] = se * p[a,s,r] + (1 - sp) * (1 - p[a,s,r]);
      }
}

model {
  // Priors
  se ~ beta(89.2,9.8);
  sp ~ beta(157.4,54);
  bsex ~ normal(0, 5);
  bage ~ normal(0, sd_age);
  bregion ~ normal(0, sd_region);
  sd_age ~ normal(0, 5);
  sd_region ~ normal(0,10);

  // Likelihood
  for(a in 1:5)
    for(s in 1:2)
      for(r in 1:5)
        y[a,s,r] ~ binomial(n[a,s,r], p_obs[a,s,r]);

  x ~ binomial(N_se, se);
  z ~ binomial(N_sp, sp);
}

generated quantities {
  real p_overall = 0;
  vector[5] p_age = rep_vector(0, 5);
  vector[2] p_sex = rep_vector(0, 2);
  vector[5] p_region = rep_vector(0, 5);

  // overall prevalence
  for(a in 1:5)
    for(s in 1:2)
      for(r in 1:5)
        p_overall += p[a,s,r] * pw[a,s,r];

  // Regional prevalence
  for(r in 1:5)
    for(a in 1:5)
      for(s in 1:2)
        p_region[r] += p[a,s,r] * pw[a,s,r] / tot_pw_region[r];

  // Age-specific prevalence
  for(a in 1:5)
    for(s in 1:2)
      for(r in 1:5)
        p_age[a] += p[a,s,r] * pw[a,s,r] / tot_pw_age[a];

  // Sex-specific prevalence
  for(s in 1:2)
    for(a in 1:5)
      for(r in 1:5)
        p_sex[s] += p[a,s,r] * pw[a,s,r] / tot_pw_sex[s];
}
")

fit_region <- sampling(
  object = mrp_model_region,
  data = list(
    y = y,
    n = n,
    pw = pw,
    tot_pw_age = tot_pw_age,
    tot_pw_sex = tot_pw_sex,
    tot_pw_region = tot_pw_region,
    x = 74,
    z = 67,  
    N_se = 82, 
    N_sp = 90 
  ),
  chains = 1,
  warmup = 1000,
  iter = 10000,
  cores = 3,
  seed = 111,
  pars = c(
    "p_overall", "se", "sp"
   
  ),
  refresh = 500
)



print(fit_region, digits=3)
 



#tableby()
df_all$site<- as.factor(df_all$site)
df_all$sex<- as.factor(df_all$sex)
df_all$age_group<- as.factor(df_all$age_group)
df_all$CHKpos<- as.factor(df_all$CHKpos)
df_all$DENpos<- as.factor(df_all$DENpos)
df_all$RVFpos<- as.factor(df_all$RVFpos)

df_all$age_group <- relevel(factor(df_all$age_group, ordered = FALSE), ref = "<5")
df_all$site <- relevel(factor(df_all$site), ref = "NAIROBI")

# Run logistic regression again
tab1 <- glm(RVFpos ~ age_group + sex + site, family = binomial, data = df_all)

# View summary
summary(tab1)

# Optional: get adjusted odds ratios and 95% CIs
exp(cbind(OR = coef(tab1), confint(tab1)))



tab<- tableby(data = df_all, CHKpos ~ site + sex + age_group)
summary(tab, text = TRUE)

tab1<-glm(RVFpos ~ age_group + sex + site, family = binomial, data = df_all)
summary(tab1)
exp(cbind(OR = coef(tab1), confint(tab1)))


# Crude  ORs
tab_crude_site <- glm(RVFpos ~ site, family = binomial, data = df_all)
exp(cbind(OR = coef(tab_crude_site), confint(tab_crude_site)))
summary(tab_crude_site)

# Crude OR for sex
tab_crude_sex <- glm(DENpos ~ sex, family = binomial, data = df_all)
exp(cbind(OR = coef(tab_crude_sex), confint(tab_crude_sex)))
summary(tab_crude_sex)
# Crude OR for age group
tab_crude_age <- glm(DENpos ~ age_group, family = binomial, data = df_all)
exp(cbind(OR = coef(tab_crude_age), confint(tab_crude_age)))
summary(tab_crude_age)




# Fit the full model for RVFpos
full_model_rvf <- glm(RVFpos ~ site + age_group + sex, family = binomial, data = df_all)
# Extract ORs and 95% CIs
OR <- exp(coef(full_model_rvf))
CI <- exp(confint(full_model_rvf))

forest_data <- data.frame(
  predictor = names(OR),
  OR = OR,
  lower = CI[,1],
  upper = CI[,2]
)

# Remove intercept and any rows with NA
forest_data <- forest_data %>%
  filter(predictor != "(Intercept)" & !is.na(OR) & !is.na(lower) & !is.na(upper))

# Add significance
forest_data$significant <- ifelse(summary(full_model_rvf)$coefficients[forest_data$predictor,4] < 0.05, "yes", "no")

# Create a grouping variable
forest_data <- forest_data %>%
  mutate(group = case_when(
    grepl("^site", predictor) ~ "Site",
    grepl("^age", predictor) ~ "Age group",
    grepl("^sex", predictor) ~ "Sex"
  ))

# Order predictors within each group
age_levels <- grep("^age", forest_data$predictor, value = TRUE)
age_levels <- age_levels[order(as.numeric(gsub(".*?(\\d+).*", "\\1", age_levels)))]
site_levels <- grep("^site", forest_data$predictor, value = TRUE)
sex_levels <- grep("^sex", forest_data$predictor, value = TRUE)

forest_data$predictor <- factor(forest_data$predictor,
                                levels = c(age_levels, site_levels, sex_levels))

# Nice x-axis breaks
max_val <- ceiling(max(forest_data$upper, na.rm = TRUE))
min_val <- floor(min(forest_data$lower, na.rm = TRUE))
x_breaks <- pretty(c(min_val, max_val), n = 6)

# Plot
ggplot(forest_data, aes(x = OR, y = predictor, color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = x_breaks) +
  xlab("Odds Ratio") +
  ylab("Predictor") +
  ggtitle("Forest Plot of RVFpos") +
  theme_minimal() +
  facet_grid(group ~ ., scales = "free_y", space = "free")




