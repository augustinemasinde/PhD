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

#Distribution of age in each study site
table(Asembodata$agecat)
table(manyattadata$agecat)
table(Kilifidata$agecat)
table(Kiberadata$agecat)
table(NairobiUrbandata$agecat)

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
table(Asembodata$agecat, Asembodata$CHKpos)
table(Asembodata$agecat, Asembodata$RVFpos)
table(Asembodata$agecat, Asembodata$DENpos)

#age distribution seropositivity in Manyatta
table(manyattadata$agecat, manyattadata$CHKpos)
table(manyattadata$agecat, manyattadata$RVFpos)
table(manyattadata$agecat, manyattadata$DENpos)
#age distribution seropositivity in Kilifi
table(Kilifidata$agecat, Kilifidata$chke1_pos)
table(Kilifidata$agecat, Kilifidata$rvfgc_pos)
table(Kilifidata$agecat, Kilifidata$denv_pos)
#age distribution seropositivity in Kibera
table(Kiberadata$agecat, Kiberadata$chke1_pos)
table(Kiberadata$agecat, Kiberadata$rvfgc_pos)
table(Kiberadata$agecat, Kiberadata$denv_pos)
#age distribution seropositivity in Nairobi Urban
table(NairobiUrbandata$agecat, NairobiUrbandata$chke1_pos)
table(NairobiUrbandata$agecat, NairobiUrbandata$rvfgc_pos)
table(NairobiUrbandata$agecat, NairobiUrbandata$denv_pos)
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
#merged dataset


