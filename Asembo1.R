# ensure working directory is the folder with the .dta/.xlsx files
setwd("/Users/augustinemasinde/Desktop/PhD files/Global Health")

library(haven)
library(dplyr)
library(arsenal)
library(Rsero)
library(readxl)
library(rstan)
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
Kilifidata_clean <- Kilifidata %>%
  rename(ageyrs = age_y,
         CHKpos = chke1_pos,
         DENpos = denv_pos,
         RVFpos = rvfgc_pos,
         site = location) %>%
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
         RVFpos = rvfgc_pos,
         site = location) %>%
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






#fit the model age-dependent force of infection using rsero package
#estimating site specific force of infection using rsero
asembodata<- Asembodata %>%
  select(ageyrs,sex,age_cat,CHKpos, DENpos, RVFpos)
  
#recode the IgM_CHIK variable, Positive to TRUE , Negative to FALSE
asembodata$RVFpos <- ifelse(asembodata$RVFpos == 1, TRUE, FALSE)
asembodata$DENpos <- ifelse(asembodata$DENpos == 1, TRUE, FALSE)
asembodata$CHKpos <- ifelse(asembodata$CHKpos == 1, TRUE, FALSE)

asembodata$ageyrs=  as.integer(asembodata$ageyrs)
asembodata <- asembodata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
asembodata$sex<-as.character(asembodata$sex)

asembo_chik <-SeroData(age_at_sampling = asembodata$ageyrs,
                       Y= asembodata$CHKpos,
                       #category = asembodata$sex,
                       #reference.category = "M",
                       sampling_year = 2022)
seroprevalence(asembo_chik)
seroprevalence.plot(asembo_chik, age_class = 5)
#fit the model constant model
ConstantModel = FOImodel(type = 'constant', priorC1 = 0.6, priorC2 = 1, seroreversion = 1, priorRho1 = 0.20, priorRho2 = 0)

FOIfit.constant = fit(data = asembo_chik,  model = ConstantModel, chains=1, iter=5000)

#customized function to adjust age classes
seroprevalence.fit2<- function(FOIfit,                   
                               individual_samples = 0,
                               age_class = 5,
                               YLIM=1,
                               fill.color  = "#d7def3", 
                               line.color  = "#5e6b91", 
                               mid.age.plot = TRUE,
                               ...){
  
  plots  <- list()
  data = FOIfit$data
  
  chains <- rstan::extract(FOIfit$fit)
  se = FOIfit$model$se 
  sp = FOIfit$model$sp 
  A <- FOIfit$data$A
  latest_sampling_year <- max(FOIfit$data$sampling_year)
  years <- seq(1,A)
  
  index.plot = 0
  unique.categories = data$unique.categories
  sorted.year = sort.int(unique(FOIfit$data$sampling_year), index.return = TRUE)
  Y = 0
  
  for(sampling_year in sorted.year$x ){
    Y = Y + 1
    for(cat in unique.categories){ 
      
      if(length(unique.categories) == 1){
        title =  ""
      } else {
        title = paste0('Category: ',cat)
      }
      
      index.plot = index.plot + 1
      age_group = data$age_group[which(data$sampling_year ==  sampling_year)][1]
      w = which(data$sampling_year ==  sampling_year & data$category == cat, arr.ind = TRUE)[,1]
      subdat = subset(data, sub = w)
      
      # compute the proportion of seropositive
      P = chains$P[,,sorted.year$ix[Y], 1]
      d = data$categoryindex[w]
      p1 = proportions.index(d)
      
      M = dim(chains$P)[1] 
      Pinf = matrix(0, nrow = M, ncol = FOIfit$data$A)
      
      # infection probability weighted on the categories      
      for(i in 1:length(p1$index)){
        Pinf =  Pinf +  p1$prop[i] * ( se - (se + sp - 1) * chains$P[,,sorted.year$ix[Y], p1$index[i]] ) 
      }
      
      par_out <- apply(Pinf, 2, function(x) c(mean(x), quantile(x, probs = c(0.025, 0.975))))
      par_out[par_out > YLIM] = YLIM # set to the upper limits for plotting
      
      # X axis (model fit coordinates)
      years.plotted =  seq(latest_sampling_year - sampling_year + 1, dim(chains$P)[2])
      years.plotted.normal = years.plotted - min(years.plotted) + 1
      meanFit <- data.frame(x = years.plotted.normal, y = par_out[1, years.plotted ])
      
      # create the envelope
      xpoly <- c(years.plotted.normal, rev(years.plotted.normal))
      ypoly <- c(par_out[3, years.plotted ], rev(par_out[2, years.plotted ]))
      DataEnvelope = data.frame(x = xpoly, y = ypoly)
      
      # histogram of data using the requested age_class (for >=10 groups)
      histdata_full <- sero.age.groups(dat = subdat, age_class = age_class, YLIM = YLIM, mid.age.plot = mid.age.plot) 
      histdata_full$labels_text <- as.character(histdata_full$labels)
      
      # histogram of data with single-year bins for the 0-9 region
      histdata_1 <- sero.age.groups(dat = subdat, age_class = 1, YLIM = YLIM, mid.age.plot = mid.age.plot)
      histdata_1$labels_text <- as.character(histdata_1$labels)
      
      # trim to last complete row (as in original)
      last_row_index_full <- tail(which(complete.cases(histdata_full)), 1)
      if(length(last_row_index_full) > 0){
        histdata_full <- histdata_full[1:last_row_index_full,]
      }
      last_row_index_1 <- tail(which(complete.cases(histdata_1)), 1)
      if(length(last_row_index_1) > 0){
        histdata_1 <- histdata_1[1:last_row_index_1,]
      }
      
      # maximum age to restrict plotting of fit
      max.age = max(histdata_full$age)
      DataEnvelope = subset(DataEnvelope, x <= max.age)
      meanFit = subset(meanFit, x <= max.age )
      
      # Split into <10 and >=10
      hist_under10 <- subset(histdata_1, age < 10)        # single-year groups 0-9
      hist_over10  <- subset(histdata_full, age >= 10)    # original 5-year groups (age_class)
      
      env_under10  <- subset(DataEnvelope, x < 10)
      env_over10   <- subset(DataEnvelope, x >= 10)
      
      mean_under10 <- subset(meanFit, x < 10)
      mean_over10  <- subset(meanFit, x >= 10)
      
      ## X-axis breaks & labels:
      # under-10: use the single-year midpoints from hist_under10
      if (nrow(hist_under10) > 0) {
        breaks_under10 <- hist_under10$age
        labels_under10 <- hist_under10$labels_text
      } else {
        breaks_under10 <- numeric(0); labels_under10 <- character(0)
      }
      
      # >=10: keep the same grouping logic as original (place tick roughly at start of group as original code did)
      if (nrow(hist_over10) > 0) {
        breaks_over10 <- hist_over10$age - round(age_class/2)  # same pattern as original function
        labels_over10 <- hist_over10$labels_text
      } else {
        breaks_over10 <- numeric(0); labels_over10 <- character(0)
      }
      
      # combine & sort (ensure ascending)
      breaks_all <- c(breaks_under10, breaks_over10)
      labels_all <- c(labels_under10, labels_over10)
      if(length(breaks_all) > 0){
        ord <- order(breaks_all)
        breaks_all <- breaks_all[ord]
        labels_all <- labels_all[ord]
      }
      
      # Build the single combined plot (continuous fit, ribbon and points)
      p <- ggplot2::ggplot() + 
        ggplot2::geom_polygon(data = DataEnvelope, ggplot2::aes(x = x, y = y), fill = fill.color) + 
        ggplot2::geom_line(data = meanFit, ggplot2::aes(x = x, y = y), linewidth = 1, color = line.color)
      
      # add individual samples if requested
      if(individual_samples > 0){
        Index_samples <- sample(nrow(Pinf), individual_samples)
        for (i in Index_samples){
          ind_foi <- data.frame(x = years.plotted.normal, y = Pinf[i, years.plotted])
          p <- p + ggplot2::geom_line(data = ind_foi, ggplot2::aes(x = x, y = y), linewidth = 0.8, colour = "#bbbbbb", alpha = 0.6)
        }
      }
      
      # add data points & error bars for <10 and >=10
      if(nrow(hist_under10) > 0){
        p <- p + ggplot2::geom_point(data = hist_under10, ggplot2::aes(x = age, y = mean)) +
          ggplot2::geom_segment(data = hist_under10, ggplot2::aes(x = age, y = lower, xend = age, yend = upper))
      }
      if(nrow(hist_over10) > 0){
        p <- p + ggplot2::geom_point(data = hist_over10, ggplot2::aes(x = age, y = mean)) +
          ggplot2::geom_segment(data = hist_over10, ggplot2::aes(x = age, y = lower, xend = age, yend = upper))
      }
      
      # final formatting + combined x-axis
      p <- p +
        ggplot2::scale_x_continuous(breaks = breaks_all, labels = labels_all) +
        ggplot2::xlab("Age (years)") + 
        ggplot2::ylab("Seroprevalence") + 
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, size = 14, vjust = 0.5),
                       axis.text.y = ggplot2::element_text(size = 14),
                       text = ggplot2::element_text(size = 14)) +
        ggplot2::ylim(0, YLIM) +
        ggplot2::ggtitle(title)
      
      plots[[index.plot]] <- p 
      plots[[index.plot]]$category <- cat 
      plots[[index.plot]]$year <- sampling_year 
      
    } # end for cat
  } # end for sampling_year
  return(plots)
}

#' @export
proportions.index <- function(d){
  b.x = c()
  b.y = c()
  ii = 0
  for(i in unique(d)){
    ii = ii + 1
    b.x[ii]  = i  
    b.y[ii]  = sum(d == i) / length(d)
  }
  return(list(index = b.x, prop = b.y ))
}

seroprevalence.fit2(FOIfit.constant, YLIM = 1, age_class = 5)


#Plotting the posterior distributions
plot_posterior(FOIfit.constant)

#Posterior distribution of relevant model parameters
parameters_credible_intervals(FOIfit.constant)
traceplot_Rsero(FOIfit.constant)





#manyatta
#fit the model age-dependent force of infection using rsero package
#estimating site specific force of infection using rsero
manyattadata<- manyattadata_clean %>%
  select(ageyrs,sex,age_cat,CHKpos, DENpos, RVFpos)

#recode the IgM_CHIK variable, Positive to TRUE , Negative to FALSE
manyattadata$RVFpos <- ifelse(manyattadata$RVFpos == 1, TRUE, FALSE)
manyattadata$DENpos <- ifelse(manyattadata$DENpos == 1, TRUE, FALSE)
manyattadata$CHKpos <- ifelse(manyattadata$CHKpos == 1, TRUE, FALSE)

manyattadata$ageyrs=  as.integer(manyattadata$ageyrs)
manyattadata <- manyattadata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
manyattadata$sex<-as.character(manyattadata$sex)

manyatta_chik <-SeroData(age_at_sampling = manyattadata$ageyrs,
                       Y= manyattadata$CHKpos,
                       age_class = 5,
                       #category = manyattadata$sex,
                       #reference.category = "M",
                       sampling_year = 2022)
seroprevalence(manyatta_chik)
seroprevalence.plot(manyatta_chik, age_class = 5)
#fit the model constant model
ConstantModel = FOImodel(type = 'constant', priorC1 = 0.079,priorC2 = 1, seroreversion = 1, priorRho1 = 0.20, priorRho2 = 0)

FOIfit.constant = fit(data = manyatta_chik,  model = ConstantModel, chains=1, iter=5000)
seroprevalence.fit(FOIfit.constant, YLIM=1, age_class = 5)

#Plotting the posterior distributions
plot_posterior(FOIfit.constant)

#Posterior distribution of relevant model parameters
parameters_credible_intervals(FOIfit.constant)
traceplot_Rsero(FOIfit.constant)




#kilifi
#fit the model age-dependent force of infection using rsero package
#estimating site specific force of infection using rsero
kilifidata<- Kilifidata_clean %>%
  select(ageyrs,sex,age_cat,CHKpos, DENpos, RVFpos)

#recode the IgM_CHIK variable, Positive to TRUE , Negative to FALSE
kilifidata$RVFpos <- ifelse(kilifidata$RVFpos == 1, TRUE, FALSE)
kilifidata$DENpos <- ifelse(kilifidata$DENpos == 1, TRUE, FALSE)
kilifidata$CHKpos <- ifelse(kilifidata$CHKpos == 1, TRUE, FALSE)

kilifidata$ageyrs=  as.integer(kilifidata$ageyrs)
kilifidata <- kilifidata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
kilifidata$sex<-as.character(kilifidata$sex)

kilifidata_chik <-SeroData(age_at_sampling = kilifidata$ageyrs,
                         Y= kilifidata$CHKpos,
                         #category = kilifidata$sex,
                         #reference.category = "M",
                         sampling_year = 2022)
seroprevalence(kilifidata_chik)
seroprevalence.plot(kilifidata_chik, age_class = 5)
#fit the model constant model
ConstantModel = FOImodel(type = 'constant', priorC1 = 0.02, priorC2 = 1, seroreversion = 1, priorRho2 = .2, priorRho1 = 0)
FOIfit.constant = fit(data = kilifidata_chik,  model = ConstantModel, chains=1, iter=5000)
seroprevalence.fit(FOIfit.constant,age_class = 5)

#Posterior distribution of relevant model parameters
parameters_credible_intervals(FOIfit.constant)




#kibera
#fit the model age-dependent force of infection using rsero package
#estimating site specific force of infection using rsero
kiberadata<- Kiberadata_clean %>%
  select(ageyrs,sex,age_cat,CHKpos, DENpos, RVFpos)

#recode the IgM_CHIK variable, Positive to TRUE , Negative to FALSE
kiberadata$RVFpos <- ifelse(kiberadata$RVFpos == 1, TRUE, FALSE)
kiberadata$DENpos <- ifelse(kiberadata$DENpos == 1, TRUE, FALSE)
kiberadata$CHKpos <- ifelse(kiberadata$CHKpos == 1, TRUE, FALSE)

kiberadata$ageyrs=  as.integer(kiberadata$ageyrs)
kiberadata <- kiberadata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
kiberadata$sex<-as.character(kiberadata$sex)

kiberadata_chik <-SeroData(age_at_sampling = kiberadata$ageyrs,
                           Y= kiberadata$CHKpos,
                           #category = kiberadata$sex,
                           #reference.category = "M",
                           sampling_year = 2022)
seroprevalence(kiberadata_chik)
seroprevalence.plot(kiberadata_chik, age_class = 5)
#fit the model constant model
ConstantModel = FOImodel(type = 'constant', priorC1 = 0.01, priorC2 = 1, seroreversion = 1, priorRho2 = .2, priorRho1 = 0)
FOIfit.constant = fit(data = kiberadata_chik,  model = ConstantModel, chains=1, iter=5000)
seroprevalence.fit(FOIfit.constant, YLIM= 1, age_class = 5)

#Posterior distribution of relevant model parameters
parameters_credible_intervals(FOIfit.constant)

#Nairobi Urban
#kibera
#fit the model age-dependent force of infection using rsero package
#estimating site specific force of infection using rsero
Nairobidata<- Nairobi_clean %>%
  select(ageyrs,sex,age_cat,CHKpos, DENpos, RVFpos)

#recode the IgM_CHIK variable, Positive to TRUE , Negative to FALSE
Nairobidata$RVFpos <- ifelse(Nairobidata$RVFpos == 1, TRUE, FALSE)
Nairobidata$DENpos <- ifelse(Nairobidata$DENpos == 1, TRUE, FALSE)
Nairobidata$CHKpos <- ifelse(Nairobidata$CHKpos == 1, TRUE, FALSE)

Nairobidata$ageyrs=  as.integer(Nairobidata$ageyrs)
Nairobidata <- Nairobidata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
Nairobidata$sex<-as.character(Nairobidata$sex)

Nairobidata_chik <-SeroData(age_at_sampling = Nairobidata$ageyrs,
                           Y= Nairobidata$CHKpos,
                           #category = kiberadata$sex,
                           #reference.category = "M",
                           sampling_year = 2022)
seroprevalence(Nairobidata_chik)
seroprevalence.plot(Nairobidata_chik, age_class = 5)
#fit the model constant model
ConstantModel = FOImodel(type = 'piecewise', K=2, priorC1 = 0.2, priorC2 = 1)
FOIfit.constant = fit(data = Nairobidata_chik,  model = ConstantModel, chains=1, iter=5000)
seroprevalence.fit(FOIfit.constant, YLIM=1, age_class = 5)

#Posterior distribution of relevant model parameters
parameters_credible_intervals(FOIfit.constant)





#estimating seroprevalence and visualizing aggregated data
df_all<- df_all %>% mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
df_all$CHKpos <- ifelse(df_all$CHKpos == 1, TRUE, FALSE)
df_all$DENpos <- ifelse(df_all$DENpos == 1, TRUE, FALSE)
df_all$RVFpos <- ifelse(df_all$RVFpos == 1, TRUE, FALSE)
df_all$ageyrs=  as.integer(df_all$ageyrs)

df_chik <-SeroData(age_at_sampling = round(as.numeric(df_all$ageyrs)),
                            Y= as.numeric(df_all$CHKpos),
                            category = as.character(df_all$site),
                            reference.category = "NAIROBI",
                            sampling_year = as.numeric(2022))
seroprevalence(df_chik)
seroprevalence.plot(df_chik, age_class = 5)

#model of force of infections
piecewisemodel= FOImodel(type="constant", cat_lambda = TRUE, priorC1 = 0.2, priorC2 = 1, seroreversion = 1, priorRho1 = 0.20, priorRho2 = 1, se=0.93, sp=0.902)

#fit the model
FOIfit.piecewise= fit(data=df_chik, model=piecewisemodel, chains=1, iter=5000)

#visualize the best model fit
seroprevalence.fit(FOIfit.piecewise, YLIM=1, age_class = 5)



#Plotting the posterior distributions
plot_posterior(FOIfit.piecewise)

#Posterior distribution of relevant model parameters
parameters_credible_intervals(FOIfit.piecewise)



#plotting annual force of infection
plot(FOIfit.constant, YLIM =1)
plot(FOIfit.outbreak, YLIM =1)
plot(FOIfit.piecewise, YLIM =1)


traceplot_Rsero(FOIfit.constant)
traceplot_Rsero(FOIfit.outbreak)
traceplot_Rsero(FOIfit.piecewise)


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

nairobi_data <- nrb_counts %>% 
  left_join(nrb_set, by = c("sex", "agecat"))


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
  summarise(y = sum(chke1_pos), n = n()) %>%
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
  bsex ~ normal(0,5); 
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
  object = mrp_model,      # <- compiled model, not raw text
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

