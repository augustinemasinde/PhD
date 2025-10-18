library(haven)
library(dplyr)
library(arsenal)
Asembodata <- read_dta("Asembo PBIDS arbo.dta")
Kilifidata <- read_dta("KHDSS arbo.dta")
Kiberadata <- read_dta("Kibera PBIDS arbo.dta")
manyattadata <- read_dta("Manyatta HDSS arbo.dta")
NairobiUrbandata <- read_dta("NUHDSS arbo.dta")

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


Asembodata<- Asembodata %>% 
  mutate(age_cat = case_when(
    ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs < 18 ~ "5-17",
    ageyrs >= 18 & ageyrs < 46 ~ "18-45",
    ageyrs >= 46 & ageyrs <65 ~ "46-64",
    ageyrs >= 65 ~ "65+"
  ))
Asembodata$age_cat <- factor(
  Asembodata$age_cat,
  levels = c("<5", "5-17", "18-45", "46-64", "65+"),
  ordered = TRUE
)

table(Asembodata$age_cat)

manyattadata<- manyattadata %>% 
  mutate(age_cat = case_when(
    ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs <18 ~ "5-17",
    ageyrs >= 18 & ageyrs < 46 ~ "18-45",
    ageyrs >= 46 & ageyrs <65 ~ "46-64",
    ageyrs >= 65 ~ "65+"
  ))
manyattadata$age_cat <- factor(
  manyattadata$age_cat,
  levels = c("<5", "5-17", "18-45", "46-64", "65+"),
  ordered = TRUE
)

table(manyattadata$age_cat)

Kilifidata<- Kilifidata %>% 
  mutate(age_cat = case_when(
    age_y<5 ~ "<5",
    age_y >= 5 & age_y < 18 ~ "5-17",
    age_y >= 18 & age_y < 46 ~ "18-45",
    age_y >= 46 & age_y < 65 ~ "46-64",
    age_y >= 65 ~ "65+"
    
  ))
Kilifidata$age_cat <- factor(
  Kilifidata$age_cat,
  levels = c("<5", "5-17", "18-45", "46-64", "65+"),
  ordered = TRUE
)
table(Kilifidata$age_cat)


Kiberadata<- Kiberadata %>%
  mutate(age_cat = case_when(
    ageyrs < 5 ~ "<5",
    ageyrs >= 5 & ageyrs < 18 ~ "5-17",
    ageyrs >= 18 & ageyrs < 46 ~ "18-45",
    ageyrs >= 46 & ageyrs <65 ~ "46-64",
    ageyrs >= 65 ~ "65+"
  ))

Kiberadata$age_cat <- factor(
  Kiberadata$age_cat,
  levels = c("<5", "5-17", "18-45", "46-64", "65+"),
  ordered = TRUE
)
table(Kiberadata$age_cat)

NairobiUrbandata<- NairobiUrbandata %>%
  mutate(age_cat = case_when(
    age_y < 5 ~ "<5",
    age_y >= 5 & age_y < 18 ~ "5-17",
    age_y >= 18 & age_y < 46 ~ "18-45",
    age_y >= 46 & age_y < 65 ~ "46-64",
    age_y >= 65 ~ "65+"
  ))
NairobiUrbandata$age_cat <- factor(
  NairobiUrbandata$age_cat,
  levels = c("<5", "5-17", "18-45", "46-64", "65+"),
  ordered = TRUE
)
table(NairobiUrbandata$age_cat)


#Distribution of age in each study site
table(Asembodata$age_cat)
table(manyattadata$age_cat)
table(Kilifidata$age_cat)
table(Kiberadata$age_cat)
table(NairobiUrbandata$age_cat)

#CHIKV seropositivity in all sites
table(Asembodata$CHKpos)
table(manyattadata$CHKpos)
table(Kilifidata$chke1_pos)
table(Kiberadata$chke1_pos)
table(NairobiUrbandata$chke1_pos)

#RVFV seropositivity in all sites
table(Asembodata$RVFpos)
table(manyattadata$RVFpos)
table(Kilifidata$rvfgc_pos)
table(Kiberadata$rvfgc_pos)
table(NairobiUrbandata$rvfgc_pos)

#DENV seropositivity in all sites
table(Asembodata$DENpos)
table(manyattadata$DENpos)
table(Kilifidata$denv_pos)
table(Kiberadata$denv_pos)
table(NairobiUrbandata$denv_pos)

#Age distribution seropositivity in Asembo
table(Asembodata$age_cat, Asembodata$CHKpos)
table(Asembodata$age_cat, Asembodata$RVFpos)
table(Asembodata$age_cat, Asembodata$DENpos)

#age distribution seropositivity in Manyatta
table(manyattadata$age_cat, manyattadata$CHKpos)
table(manyattadata$age_cat, manyattadata$RVFpos)
table(manyattadata$age_cat, manyattadata$DENpos)
#age distribution seropositivity in Kilifi
table(Kilifidata$age_cat, Kilifidata$chke1_pos)
table(Kilifidata$age_cat, Kilifidata$rvfgc_pos)
table(Kilifidata$age_cat, Kilifidata$denv_pos)
#age distribution seropositivity in Kibera
table(Kiberadata$age_cat, Kiberadata$chke1_pos)
table(Kiberadata$age_cat, Kiberadata$rvfgc_pos)
table(Kiberadata$age_cat, Kiberadata$denv_pos)
#age distribution seropositivity in Nairobi Urban
table(NairobiUrbandata$age_cat, NairobiUrbandata$chke1_pos)
table(NairobiUrbandata$age_cat, NairobiUrbandata$rvfgc_pos)
table(NairobiUrbandata$age_cat, NairobiUrbandata$denv_pos)
#sex distribution seropositivity in Asembo
table(Asembodata$sex, Asembodata$CHKpos)
table(Asembodata$sex, Asembodata$RVFpos)
table(Asembodata$sex, Asembodata$DENpos)
#sex distribution seropositivity in Manyatta
table(manyattadata$sex, manyattadata$CHKpos)
table(manyattadata$sex, manyattadata$RVFpos)
table(manyattadata$sex, manyattadata$DENpos)
#sex distribution seropositivity in Kilifi
table(Kilifidata$sex, Kilifidata$chke1_pos)
table(Kilifidata$sex, Kilifidata$rvfgc_pos)
table(Kilifidata$sex, Kilifidata$denv_pos)
#sex distribution seropositivity in Kibera
table(Kiberadata$sex, Kiberadata$chke1_pos)
table(Kiberadata$sex, Kiberadata$rvfgc_pos)
table(Kiberadata$sex, Kiberadata$denv_pos)
#sex distribution seropositivity in Nairobi Urban
table(NairobiUrbandata$sex, NairobiUrbandata$chke1_pos)
table(NairobiUrbandata$sex, NairobiUrbandata$rvfgc_pos)
table(NairobiUrbandata$sex, NairobiUrbandata$denv_pos)


#Descriptive statistics using arsenal package
Asembodata$agecat<- as.factor(Asembodata$agecat)
Asembodata$sex<- as.factor(Asembodata$sex)
asembobiv<-tableby(CHKpos ~ agecat + sex , data = Asembodata) %>%
  summary(text = TRUE, title = "Asembo Arbovirus Study chikv")
write2word(asembobiv, file = "Asembo biv CHIKV.docx")

Asembodata$RVFpos<- as.factor(Asembodata$RVFpos)
Asembodata$sex<- as.factor(Asembodata$sex)
Asembodata$agecat<- as.factor(Asembodata$agecat)
asembobiv2 <-tableby(RVFpos~ agecat + sex , data = Asembodata) %>%
  summary(text = TRUE, title = "Asembo Arbovirus Study rvfv")
write2word(asembobiv2, file = "Asembo biv RVFV.docx")

Asembodata$DENpos<- as.factor(Asembodata$DENpos)
Asembodata$sex<- as.factor(Asembodata$sex)
Asembodata$agecat<- as.factor(Asembodata$agecat)
asembobiv3<- tableby(DENpos ~ agecat + sex , data = Asembodata) %>%
  summary(text = TRUE, title = "Asembo Arbovirus Study denv")
write2word(asembobiv3, file = "Asembo biv DENV.docx")

manyattadata$agecat<- as.factor(manyattadata$agecat)
manyattadata$sex<- as.factor(manyattadata$sex)
manyattadata$CHKpos<- as.factor(manyattadata$CHKpos)
manyatabiv1 <- tableby(CHKpos~ agecat+ sex , data = manyattadata) %>%
  summary(text = TRUE, title = "Manyatta Arbovirus Study chikv")
write2word(manyatabiv1, file = "Manyatta biv CHIKV.docx")

manyattadata$RVFpos<- as.factor(manyattadata$RVFpos)
manyattadata$sex<- as.factor(manyattadata$sex)
manyattadata$agecat<- as.factor(manyattadata$agecat)
manyatabiv2 <- tableby(RVFpos ~ agecat + sex , data = manyattadata) %>%
  summary(text = TRUE, title = "Manyatta Arbovirus Study rvfv")
write2word(manyatabiv2, file = "Manyatta biv RVFV.docx")
manyattadata$DENpos<- as.factor(manyattadata$DENpos)
manyattadata$sex<- as.factor(manyattadata$sex)
manyattadata$agecat<- as.factor(manyattadata$agecat)
manyatabiv3<- tableby(DENpos ~ agecat + sex , data = manyattadata) %>%
  summary(text = TRUE, title = "Manyatta Arbovirus Study denv")
write2word(manyatabiv3, file = "Manyatta biv DENV.docx")


Kilifidata$agecat<- as.factor(Kilifidata$agecat)
Kilifidata$sex<- as.factor(Kilifidata$sex)
Kilifidata$chke1_pos<- as.factor(Kilifidata$chke1_pos)
kilifibiv1<- tableby(chke1_pos~agecat+sex , data = Kilifidata) %>%
  summary(text = TRUE, title = "Kilifi Arbovirus Study chikv")
write2word(kilifibiv1, file = "Kilifi biv CHIKV.docx")

Kilifidata$rvfgc_pos<- as.factor(Kilifidata$rvfgc_pos)
Kilifidata$sex<- as.factor(Kilifidata$sex)
Kilifidata$agecat<- as.factor(Kilifidata$agecat)
kilifibiv2 <- tableby(rvfgc_pos ~ agecat + sex , data = Kilifidata) %>%
  summary(text = TRUE, title = "Kilifi Arbovirus Study rvfv")
write2word(kilifibiv2, file = "Kilifi biv RVFV.docx")

Kilifidata$denv_pos<- as.factor(Kilifidata$denv_pos)
Kilifidata$sex<- as.factor(Kilifidata$sex)
Kilifidata$agecat<- as.factor(Kilifidata$agecat)
kilifibiv3 <- tableby(denv_pos ~ agecat+sex , data = Kilifidata) %>%
  summary(text = TRUE, title = "Kilifi Arbovirus Study denv")
write2word(kilifibiv3, file = "Kilifi biv DENV.docx")


Kiberadata$agecat<- as.factor(Kiberadata$agecat)
Kiberadata$sex<- as.factor(Kiberadata$sex)
Kiberadata$chke1_pos<- as.factor(Kiberadata$chke1_pos)
kiberabiv1 <- tableby(chke1_pos~agecat+sex , data = Kiberadata) %>%
  summary(text = TRUE, title = "Kibera Arbovirus Study chikv")
write2word(kiberabiv1, file = "Kibera biv CHIKV.docx")

Kiberadata$rvfgc_pos<- as.factor(Kiberadata$rvfgc_pos)
Kiberadata$sex<- as.factor(Kiberadata$sex)
Kiberadata$agecat<- as.factor(Kiberadata$agecat)
kiberabiv2 <- tableby(rvfgc_pos ~ agecat+sex , data = Kiberadata) %>%
  summary(text = TRUE, title = "Kibera Arbovirus Study rvfv")
write2word(kiberabiv2, file = "Kibera biv RVFV.docx")

Kiberadata$denv_pos<- as.factor(Kiberadata$denv_pos)
Kiberadata$sex<- as.factor(Kiberadata$sex)
Kiberadata$agecat<- as.factor(Kiberadata$agecat)
kiberabiv3 <- tableby(denv_pos ~agecat+sex , data = Kiberadata) %>%
  summary(text = TRUE, title = "Kibera Arbovirus Study denv")
write2word(kiberabiv3, file = "Kibera biv DENV.docx")


NairobiUrbandata$agecat<- as.factor(NairobiUrbandata$agecat)
NairobiUrbandata$sex<- as.factor(NairobiUrbandata$sex)
NairobiUrbandata$chke1_pos<- as.factor(NairobiUrbandata$chke1_pos)
nairobibiv1 <- tableby(chke1_pos ~ agecat+sex , data = NairobiUrbandata) %>%
  summary(text = TRUE, title = "Nairobi Urban Arbovirus Study chikv")
write2word(nairobibiv1, file = "Nairobi Urban biv CHIKV.docx")
NairobiUrbandata$rvfgc_pos<- as.factor(NairobiUrbandata$rvfgc_pos)
NairobiUrbandata$sex<- as.factor(NairobiUrbandata$sex)
NairobiUrbandata$agecat<- as.factor(NairobiUrbandata$agecat)
nairobibiv2 <- tableby(rvfgc_pos~agecat+sex , data = NairobiUrbandata) %>%
  summary(text = TRUE, title = "Nairobi Urban Arbovirus Study rvfv")
write2word(nairobibiv2, file = "Nairobi Urban biv RVFV.docx")
NairobiUrbandata$denv_pos<- as.factor(NairobiUrbandata$denv_pos)
NairobiUrbandata$sex<- as.factor(NairobiUrbandata$sex)
NairobiUrbandata$agecat<- as.factor(NairobiUrbandata$agecat)
nairobibiv3<- tableby(denv_pos ~ agecat+sex , data = NairobiUrbandata) %>%
  summary(text = TRUE, title = "Nairobi Urban Arbovirus Study denv")
write2word(nairobibiv3, file = "Nairobi Urban biv DENV.docx")


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

manyattadata <- manyattadata %>%
  mutate(co_infection = case_when(
    CHKpos == 1 & RVFpos == 1 ~ "CHIKV and RVFV",
    CHKpos == 1 & DENpos == 1 ~ "CHIKV and DENV",
    RVFpos == 1 & DENpos == 1 ~ "RVFV and DENV",
    CHKpos == 1 & RVFpos == 1 & DENpos == 1 ~ "CHIKV, RVFV and DENV",
    TRUE ~ "No co-infection"
  ))
table(manyattadata$co_infection)

Kilifidata <- Kilifidata %>%
  mutate(co_infection = case_when(
    chke1_pos == 1 & rvfgc_pos == 1 ~ "CHIKV and RVFV",
    chke1_pos == 1 & denv_pos == 1 ~ "CHIKV and DENV",
    rvfgc_pos == 1 & denv_pos == 1 ~ "RVFV and DENV",
    chke1_pos == 1 & rvfgc_pos == 1 & denv_pos == 1 ~ "CHIKV, RVFV and DENV",
    TRUE ~ "No co-infection"
  ))
table(Kilifidata$co_infection)

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

prop.test(1867, 4165, conf.level = 0.95, correct = F)


#data frame for site, sex, age distribution and seropositivity
df<- data.frame(
  site= c("Asembo", "Manyatta", "Kilifi", "Kibera", "Nairobi Urban"),
  total_samples= c(3056, 1400, 4165, 2193, 1545),
  chikv_positive= c(233, 98, 1867, 634, 343),
  rvfv_positive= c(79, 45, 379, 150, 85),
  denv_positive= c(39, 28, 264, 95, 56)
)



df <- data.frame(
  site = c("Asembo", "Manyatta", "Kilifi", "Kibera", "Nairobi Urban"),
  total_samples = c(3056, 1400, 4165, 2193, 1545),
  chikv_positive = c(233, 98, 1867, 634, 343),
  rvfv_positive = c(79, 45, 379, 150, 85),
  denv_positive = c(39, 28, 264, 95, 56)
)


names(Asembodata)
names(manyattadata)
names(Kilifidata)
names(Kiberadata)
names(NairobiUrbandata)




library(dplyr)

library(dplyr)

# Asembo
Asembodata_clean <- Asembodata %>%
  select(site, sex, ageyrs, age_cat, CHKpos, DENpos, RVFpos) %>%
  mutate(across(c(CHKpos, DENpos, RVFpos), ~ as.numeric(as.character(.))))

# Manyatta
manyattadata_clean <- manyattadata %>%
  select(site, sex, ageyrs, age_cat, CHKpos, DENpos, RVFpos) %>%
  mutate(across(c(CHKpos, DENpos, RVFpos), ~ as.numeric(as.character(.))))

# Kilifi
Kilifidata$site<- "Kilifi"
Kilifidata_clean <- Kilifidata %>%
  rename(ageyrs = age_y,
         CHKpos = chke1_pos,
         DENpos = denv_pos,
         RVFpos = rvfgc_pos) %>%
  select(site, sex, ageyrs, age_cat, CHKpos, DENpos, RVFpos) %>%
  mutate(across(c(CHKpos, DENpos, RVFpos), ~ as.numeric(as.character(.))))

# Kibera
Kiberadata_clean <- Kiberadata %>%
  rename(CHKpos = chke1_pos,
         DENpos = denv_pos,
         RVFpos = rvfgc_pos) %>%
  select(site, sex, ageyrs, age_cat, CHKpos, DENpos, RVFpos) %>%
  mutate(across(c(CHKpos, DENpos, RVFpos), ~ as.numeric(as.character(.))))

# Nairobi Urban
Nairobi_clean <- NairobiUrbandata %>%
  rename(ageyrs = age_y,
         CHKpos = chke1_pos,
         DENpos = denv_pos,
         RVFpos = rvfgc_pos) %>%
  select(site, sex, ageyrs, age_cat, CHKpos, DENpos, RVFpos) %>%
  mutate(across(c(CHKpos, DENpos, RVFpos), ~ as.numeric(as.character(.))))


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

#CHKpos to factor
df_all$CHKpos <- as.factor(df_all$CHKpos)
df_all$DENpos <- as.factor(df_all$DENpos)
df_all$RVFpos <- as.factor(df_all$RVFpos)

#set reference level for site, sex, CHKpos, DENpos and RVFpos
df_all$site <- relevel(df_all$site, ref = "NAIROBI")
df_all$sex <- relevel(df_all$sex, ref = "Male")
df_all$age_cat <- relevel(df_all$age_cat, ref = "<5")
df_all$CHKpos <- relevel(df_all$CHKpos, ref = "0")
df_all$DENpos <- relevel(df_all$DENpos, ref = "0")
df_all$RVFpos <- relevel(df_all$RVFpos, ref = "0")
str(df_all)

df_plot <- df_all %>%
  group_by(age_cat, CHKpos) %>%
  summarise(count = n(), .groups = "drop")

# Plot
ggplot(df_plot, aes(x = age_cat, y = count, fill = as.factor(CHKpos))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("0" = "lightblue", "1" = "purple"),
                    name = "CHIKV",
                    labels = c("Negative", "Positive")) +
  labs(x = "Age Group (years)", y = "Number of samples",
       title = "CHIKV Seroprevalence by Age Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(ggplot2)
library(dplyr)

# Aggregate by site, age group, and CHIKV status
df_plot <- df_all %>%
  group_by(site, age_cat, CHKpos) %>%
  summarise(count = n(), .groups = "drop")

# Plot with one facet per site
ggplot(df_plot, aes(x = age_cat, y = count, fill = as.factor(CHKpos))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("0" = "lightblue", "1" = "purple"),
                    name = "CHIKV",
                    labels = c("Negative", "Positive")) +
  labs(x = "Age Group (years)", y = "Number of samples",
       title = "CHIKV Seroprevalence by Age Group and Site") +
  facet_wrap(~ site) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



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


#adjusting for seropositivity of the three viruses using bayesian model
df_CHIK<- df_all %>%  select(CHKpos, sex, age_cat, site, ageyrs)
df_CHIK$CHKpos<- ifelse(df_CHIK$CHKpos==1, TRUE, FALSE)
df_CHIK$sex<- as.character(df_CHIK$sex)
df_CHIK$ageyrs<- as.integer(df_CHIK$ageyrs)
df_CHIK$ageyrs[df_CHIK$ageyrs == 0] <- 1
df_CHIK$site<- as.character(df_CHIK$site)
df_CHIK$age_cat<- as.character(df_CHIK$age_cat)

chik.aggregated = SeroData(age_at_sampling = df_CHIK$ageyrs,
                                Y = df_CHIK$CHKpos,
                                category = df_CHIK$site,
                                reference.category = "Nairobi Urban",
                                sampling_year = 2022)
#seroprevalence estimation
seroprevalence(chik.aggregated) # Value of the seroprevalence
seroprevalence.plot(chik.aggregated) # plots of the seroprevalence vs age

#fit the model constant model
ConstantModel = FOImodel(type = 'piecewise', K=2)
FOIfit.constant = fit(data = chik.aggregated,  model = ConstantModel, chains=1)
seroprevalence.fit(FOIfit.constant, YLIM=0.5)



#Bivariate analysis between sex and seropositivity
#chikungunya
chik_tab <- tableby(CHKpos~ sex + age_cat + site, data = df_all)
summary(chik_tab)
write2word(chik_tab, file = "chik biv.docx")
#dengue
den_tab <- tableby(DENpos~sex + age_cat + site, data = df_all)
summary(den_tab)
write2word(den_tab, file = "den biv.docx")
#rvfv
rvf_tab <- tableby(RVFpos~sex + age_cat + site, data = df_all)
summary(rvf_tab)
write2word(rvf_tab, file = "rvf biv.docx")


