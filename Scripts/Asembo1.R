library(haven)
library(dplyr)
library(arsenal)
library(Rsero)
library(readxl)
library(rstan)
library(brms)
library(tidyverse)
library(Rsero)
library(mclust)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(lme4)
library(serojump)
library(devtools)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



longitudinal_data <- read.csv("Paired_R2_R3(Sheet1).csv", stringsAsFactors = FALSE)

n_total <- nrow(longitudinal_data)
n_total
#Time between visits
longitudinal_data <- longitudinal_data %>%
  mutate(
    visit_date_r2 = as.POSIXct(visit_date_r2, tz = "UTC"),
    visit_date_r3 = as.POSIXct(visit_date_r3, tz = "UTC"),
    days_between_visits = as.numeric(
      difftime(visit_date_r3, visit_date_r2, units = "days")
    )
  )

#Descriptive stats on serostatus at baseline and end
baseline_pos <- sum(longitudinal_data$CHIKE1_pos_r2 == 1, na.rm = TRUE)
baseline_neg <- sum(longitudinal_data$CHIKE1_pos_r2 == 0, na.rm = TRUE)

end_pos <- sum(longitudinal_data$CHKpos_r3 == 1, na.rm = TRUE)
end_neg <- sum(longitudinal_data$CHKpos_r3 == 0, na.rm = TRUE)

baseline_pos
baseline_neg
end_pos
end_neg

table(
  Baseline = longitudinal_data$CHIKE1_pos_r2,
  Endline  = longitudinal_data$CHKpos_r3
)

longitudinal_data <- longitudinal_data %>%
  mutate(
    seroreverted = ifelse(CHIKE1_pos_r2 == 1 & CHKpos_r3 == 0, 1, 0)
  )

at_risk <- longitudinal_data %>%
  filter(CHIKE1_pos_r2 == 1)

n_reversions <- sum(at_risk$seroreverted, na.rm = TRUE)

total_person_days <- sum(at_risk$days_between_visits, na.rm = TRUE)

rho_per_day <- n_reversions / total_person_days
rho_per_year <- rho_per_day * 365





analysis_chikv <- longitudinal_data %>%
  filter(
    CHIKE1_pos_r2 == 1,
    CHK_au_r2 > 0,
    CHK_au_r3 > 0,
    CHK_au_r3 < CHK_au_r2,
    days_between_visits > 0
  ) %>%
  mutate(
    delta_t = days_between_visits,
    delta_log = log(CHK_au_r3) - log(CHK_au_r2),
    log_dt = log(delta_t)
  ) %>%
  select(pair_id, delta_t, log_dt, delta_log)

ggplot(analysis_chikv, aes(x = delta_t, y = delta_log)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic()

fit_exp <- lm(delta_log ~ delta_t - 1, data = analysis_chikv)
summary(fit_exp)

k_day <- -coef(fit_exp)["delta_t"]
half_life_days <- log(2) / k_day

k_day
half_life_days

ggplot(analysis_chikv, aes(x = log_dt, y = delta_log)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic()

fit_power <- lm(delta_log ~ log_dt - 1, data = analysis_chikv)
summary(fit_power)

AIC(fit_exp, fit_power)

t_grid <- seq(
  min(analysis_chikv$delta_t),
  max(analysis_chikv$delta_t),
  length.out = 200
)

k_day <- -coef(fit_exp)["delta_t"]

pred_exp <- data.frame(
  delta_t = t_grid,
  delta_log = -k_day * t_grid,
  model = "Exponential"
)


beta <- -coef(fit_power)["log_dt"]

pred_power <- data.frame(
  delta_t = t_grid,
  delta_log = -beta * log(t_grid),
  model = "Power-law"
)


pred_all <- bind_rows(pred_exp, pred_power)


ggplot() +
  geom_point(
    data = analysis_chikv,
    aes(x = delta_t, y = delta_log),
    alpha = 0.6
  ) +
  geom_line(
    data = pred_all,
    aes(x = delta_t, y = delta_log, color = model),
    linewidth = 1.2
  ) +
  scale_color_manual(
    values = c("Exponential" = "#0072B2", "Power-law" = "#D55E00")
  ) +
  labs(
    x = "Days between visits",
    y = "Log antibody change",
    color = "Model",
    title = "CHIKV antibody waning: exponential vs power-law"
  ) +
  theme_classic(base_size = 13)







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




library(ggplot2)
library(dplyr)
library(Cairo)
library(grid)
#forest plot
adjusted_df <- read.csv(
  "adjusted_seroprevalence_by_site_pathogen.csv",
  stringsAsFactors = FALSE)

     
labels_combined <- unlist(
  lapply(unique(adjusted_df$site), function(s) {
    c(s, paste0("  ", adjusted_df$pathogen[adjusted_df$site == s]), "")
  })
)

labels_combined <- labels_combined[-length(labels_combined)]
y_positions <- rev(seq_along(labels_combined))

plot_labels_df <- data.frame(
  y_label = labels_combined,
  y = y_positions,
  stringsAsFactors = FALSE
)

get_y_pos <- function(site, pathogen) {
  site_idx <- which(labels_combined == site)
  pathogen_label <- paste0("  ", pathogen)
  pathogen_idx <- which(labels_combined == pathogen_label)
  
  # Filter pathogen indices that come after site index
  pathogen_idx_after_site <- pathogen_idx[pathogen_idx > site_idx]
  
  if (length(pathogen_idx_after_site) == 0) return(NA_real_)
  return(y_positions[pathogen_idx_after_site[1]])
}

plot_data <- adjusted_df %>%
  rowwise() %>%
  mutate(
    y = get_y_pos(site, pathogen)
  ) %>%
  ungroup()

ggplot() +
  geom_point(data = plot_data, aes(x = estimate, y = y, color = pathogen), size = 3) +
  geom_errorbarh(data = plot_data, aes(y = y, xmin = lower, xmax = upper, color = pathogen), height = 0.3, linewidth = 0.8) +
  scale_y_continuous(
    breaks = y_positions,
    labels = labels_combined,
    limits = range(y_positions) + c(-0.5, 0.5)
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    labels = scales::percent_format(scale = 1)
  ) +
  scale_color_manual(values = c("CHIKV" = "#D55E00", "DENV" = "#0072B2", "RVFV" = "#009E73")) +
  labs(x = "Adjusted Seroprevalence (%)", y = NULL, color = "Pathogen") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(family = "Times New Roman", size = 10, hjust = 0),
    axis.text.x = element_text(family = "Times New Roman", size = 10),
    axis.title.x = element_text(family = "Times New Roman", size = 10),
    legend.position = "bottom",
    legend.title = element_text(family = "Times New Roman", size = 10),
    legend.text = element_text(family = "Times New Roman", size = 9),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_blank()
  ) +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", color = "grey50", linewidth = 0.6)

library(Cairo)

dpi <- 300
width_mm <- 107
height_mm <- 140  # Adjust height for 5 sites with spacing

width_px <- round((width_mm / 25.4) * dpi)
height_px <- round((height_mm / 25.4) * dpi)

Cairo::Cairo(
  file = "forest_plot_lancet.tiff",
  type = "tiff",
  width = width_px,
  height = height_px,
  units = "px",
  dpi = dpi,
  bg = "white"
)

print(
  last_plot() + 
    theme(
      text = element_text(family = "Times New Roman"),
      axis.text.y = element_text(family = "Times New Roman", size = 10, hjust = 0),
      axis.text.x = element_text(family = "Times New Roman", size = 10),
      axis.title.x = element_text(family = "Times New Roman", size = 10),
      legend.position = "bottom",
      legend.title = element_text(family = "Times New Roman", size = 10),
      legend.text = element_text(family = "Times New Roman", size = 9)
    )
)
dev.off()


ggsave(
  filename = "forest_plot_lancet_single1.tiff",
  plot = p,
  device = "tiff",
  width = 4.21,
  height = 3.5, # adjust as needed
  dpi = 300,
  units = "in"
)

ggsave("forest_plot_highres.png", plot = p, width = 10, height = 6, dpi = 300, units = "in")
ggsave("forest_plot_highres1.tiff", plot = p, width = 10, height = 6, dpi = 300, units = "in", device = "tiff")


#By pathogen

library(dplyr)
library(ggplot2)
library(scales)
library(Cairo)

adjusted_df <- read.csv(
  "adjusted_seroprevalence_by_site_pathogen.csv",
  stringsAsFactors = FALSE
)

labels_combined <- unlist(
  lapply(unique(adjusted_df$pathogen), function(p) {
    c(
      p,
      paste0("  ", adjusted_df$site[adjusted_df$pathogen == p]),
      ""
    )
  })
)

labels_combined <- labels_combined[-length(labels_combined)]
y_positions <- rev(seq_along(labels_combined))

get_y_pos <- function(pathogen, site) {
  pathogen_idx <- which(labels_combined == pathogen)
  site_label <- paste0("  ", site)
  site_idx <- which(labels_combined == site_label)
  site_idx_after_pathogen <- site_idx[site_idx > pathogen_idx]
  if (length(site_idx_after_pathogen) == 0) return(NA_real_)
  y_positions[site_idx_after_pathogen[1]]
}

plot_data <- adjusted_df %>%
  rowwise() %>%
  mutate(
    y = get_y_pos(pathogen, site)
  ) %>%
  ungroup()

p <- ggplot() +
  geom_point(data = plot_data, aes(x = estimate, y = y, color = pathogen), size = 3) +
  geom_errorbarh(
    data = plot_data,
    aes(y = y, xmin = lower, xmax = upper, color = pathogen),
    height = 0.3,
    linewidth = 0.8
  ) +
  scale_y_continuous(
    breaks = y_positions,
    labels = labels_combined,
    limits = range(y_positions) + c(-0.5, 0.5)
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    labels = percent_format(scale = 1)
  ) +
  scale_color_manual(
    values = c(
      "CHIKV" = "#D55E00",
      "DENV"  = "#0072B2",
      "RVFV"  = "#009E73"
    )
  ) +
  labs(
    x = "Adjusted Seroprevalence (%)",
    y = NULL,
    color = "Pathogen"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.y = element_text(size = 10, hjust = 0),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_blank()
  ) +
  geom_vline(
    xintercept = c(0, 100),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.6
  )

dpi <- 300
width_mm <- 107
height_mm <- 140

width_px <- round((width_mm / 25.4) * dpi)
height_px <- round((height_mm / 25.4) * dpi)

Cairo(
  file = "forest_plot_lancet_grouped_by_pathogen.tiff",
  type = "tiff",
  width = width_px,
  height = height_px,
  units = "px",
  dpi = dpi,
  bg = "white"
)

print(p)
dev.off()

ggsave(
  filename = "forest_plot_lancet_grouped_by_pathogen.png",
  plot = p,
  width = 4.21,
  height = 5.5,
  dpi = 300,
  units = "in"
)




library(ggplot2)
library(dplyr)
library(forcats)

seroprev_df_ageadj <- read.csv("age_adjusted_seroprevalence.csv", stringsAsFactors = FALSE)

seroprev_df_ageadj <- seroprev_df_ageadj %>%
  filter(age_group != "Overall") %>%  
  filter(!is.na(age_group))           

seroprev_df_ageadj <- seroprev_df_ageadj %>%
  mutate(
    age_group = factor(age_group, levels = c("<5", "5-14", "15-44", "45-64", "65+")),
    pathogen = factor(pathogen, levels = c("CHIKV", "DENV", "RVFV")),
    site = factor(site, levels = c("Asembo", "Kilifi", "Manyatta", "Kibera", "Nairobi"))
  )

ggplot(seroprev_df_ageadj, aes(x = age_group, y = estimate, color = pathogen)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~ site, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  labs(
    title = "Age-stratified Adjusted Seroprevalence by Site and virus",
    x = "Age Group",
    y = "Adjusted Seroprevalence (%)",
    color = "Virus"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("age_stratified_seroprevalence.png", plot = last_plot(), width = 10, height = 6, dpi = 300, units = "in")

# --------------------------------------------------
# Data preparation (DO THIS ONCE)
# --------------------------------------------------

Asembodata <- Asembodata %>%
  dplyr::filter(DENV2_AIU >= 0) %>%
  dplyr::mutate(log_DENV_AIU = log10(DENV2_AIU + 1))
Manyattadata <- manyattadata %>%
  dplyr::filter(DENV2_AIU >= 0) %>%
  dplyr::mutate(log_DENV_AIU = log10(DENV2_AIU + 1))
Kiberadata <- Kiberadata %>%
  dplyr::filter(denv2ns1_au >= 0) %>%
  dplyr::mutate(log_DENV_AIU = log10(denv2ns1_au + 1))
Kilifidata <- Kilifidata %>%
  dplyr::filter(denv2ns1_au >= 0) %>%
  dplyr::mutate(log_DENV_AIU = log10(denv2ns1_au + 1))
NairobiUrbandata<- NairobiUrbandata %>%
  dplyr::filter(denv2ns1_au >= 0) %>%
  dplyr::mutate(log_DENV_AIU = log10(denv2ns1_au + 1))
Asembodata_children <- Asembodata %>% filter(ageyrs <= 5) %>% mutate(log_CHKV_AIU = log10(CHKV_AIU + 1))
# --------------------------------------------------
# Violin + boxplot (NO WARNINGS)
# --------------------------------------------------
p<-ggplot(
  Asembodata_children,
  aes(
    x = factor(CHKpos),
    y = log_CHKV_AIU
  )
) +
  geom_violin(
    trim = FALSE,
    width = 0.45,
    fill = "grey80",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.6
  ) +
  scale_x_discrete(
    labels = c("0" = "Seronegative", "1" = "Seropositive")
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(ylim = c(0, 5)) +   
  labs(
    title= "Asembo CHIKV in Children",
    x = "Serostatus",
    y = "log10 Antibody concentration (AIU + 1)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.line = element_line(linewidth = 0.8)
  )+ theme_grey(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),   # no vertical grid
    panel.grid.minor = element_blank(),     # no minor grid
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )

p1<-ggplot(
  Manyattadata,
  aes(
    x = factor(DENpos),
    y = log_DENV_AIU
  )
) +
  geom_violin(
    trim = FALSE,
    width = 0.45,
    fill = "grey80",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.6
  ) +
  scale_x_discrete(
    labels = c("0" = "Seronegative", "1" = "Seropositive")
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(ylim = c(0, 5)) +   # zoom WITHOUT dropping data
  labs(
    title= "Manyatta DENV",
    x = "Serostatus",
    y = "log10 Antibody concentration (AIU + 1)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.line = element_line(linewidth = 0.8)
  )+ theme_grey(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),   # no vertical grid
    panel.grid.minor = element_blank(),     # no minor grid
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )



p2<-ggplot(
  Kiberadata,
  aes(
    x = factor(denv_pos),
    y = log_DENV_AIU
  )
) +
  geom_violin(
    trim = FALSE,
    width = 0.45,
    fill = "grey80",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.6
  ) +
  scale_x_discrete(
    labels = c("0" = "Seronegative", "1" = "Seropositive")
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(ylim = c(0, 5)) +   # zoom WITHOUT dropping data
  labs(
    title= "Kibera DENV",
    x = "Serostatus",
    y = "log10 Antibody concentration (AIU + 1)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.line = element_line(linewidth = 0.8)
  )+ theme_grey(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),   # no vertical grid
    panel.grid.minor = element_blank(),     # no minor grid
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )



p3<-ggplot(
  Kilifidata,
  aes(
    x = factor(denv_pos),
    y = log_DENV_AIU
  )
) +
  geom_violin(
    trim = FALSE,
    width = 0.45,
    fill = "grey80",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.6
  ) +
  scale_x_discrete(
    labels = c("0" = "Seronegative", "1" = "Seropositive")
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(ylim = c(0, 5)) +   # zoom WITHOUT dropping data
  labs(
    title= "Kilifi DENV",
    x = "Serostatus",
    y = "log10 Antibody concentration (AIU + 1)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.line = element_line(linewidth = 0.8)
  )+ theme_grey(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),   # no vertical grid
    panel.grid.minor = element_blank(),     # no minor grid
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )



p4<-ggplot(
  NairobiUrbandata,
  aes(
    x = factor(denv_pos),
    y = log_DENV_AIU
  )
) +
  geom_violin(
    trim = FALSE,
    width = 0.45,
    fill = "grey80",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.6
  ) +
  scale_x_discrete(
    labels = c("0" = "Seronegative", "1" = "Seropositive")
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(ylim = c(0, 5)) +   # zoom WITHOUT dropping data
  labs(
    title= "NairobiUrban DENV",
    x = "Serostatus",
    y = "log10 Antibody concentration (AIU + 1)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.line = element_line(linewidth = 0.8)
  )+ theme_grey(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),   # no vertical grid
    panel.grid.minor = element_blank(),     # no minor grid
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )
gridExtra::grid.arrange(p, p1, p2, p3, p4, nrow=2)

#save
ggsave("DENV_violin_boxplot.png", 
       plot = gridExtra::grid.arrange(p, p1, p2, p3, p4, nrow=2),
       width = 12, height = 8, dpi = 300)

ggsave("Asembo_CHKV_children_violin_boxplot.png", 
       width = 6, height = 6, dpi = 300)



#Finite mixture model fitting using MCLUST package
Asembodata<- Asembodata %>% mutate(DENV2_AIU= ifelse(DENV2_AIU <0, 0, DENV2_AIU))
manyattadata<- manyattadata %>% mutate(DENV2_AIU= ifelse(DENV2_AIU <0, 0, DENV2_AIU))
Kilifidata<- Kilifidata %>% mutate(denv2ns1_au= ifelse(denv2ns1_au <0, 0, denv2ns1_au))
Kiberadata<- Kiberadata %>% mutate(denv2ns1_au= ifelse(denv2ns1_au <0, 0, denv2ns1_au))
NairobiUrbandata<- NairobiUrbandata %>% mutate(denv2ns1_au= ifelse(denv2ns1_au <0, 0, denv2ns1_au))
Asembodata_children<- Asembodata_children %>% mutate(CHKV_AIU= ifelse(CHKV_AIU <0, 0, CHKV_AIU))
Asembodata$DENpos<- as.factor(Asembodata$DENpos)
manyattadata$DENpos<- as.factor(manyattadata$DENpos)
Kilifidata$denv_pos<- as.factor(Kilifidata$denv_pos)
Kiberadata$denv_pos<- as.factor(Kiberadata$denv_pos)
NairobiUrbandata$denv_pos<- as.factor(NairobiUrbandata$denv_pos)
Asembodata_children$CHKpos<- as.factor(Asembodata_children$CHKpos)

#Histogram plot
p <- ggplot(
  Asembodata,
  aes(x = log(DENV2_AIU + 1), fill = DENpos)
) +
  geom_histogram(bins = 20, color = "black") +
  scale_fill_discrete(labels = c("0" = "Seronegative", "1" = "Seropositive")) +
  labs(
    title = "Asembo DENV",
    x = "Log(AIU values + 1)",
    y = "Number of samples",
    fill = "DENV serostatus"
  ) + theme_grey(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),   
    panel.grid.minor = element_blank(),    
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )

p1 <- ggplot(
  manyattadata,
  aes(x = log(DENV2_AIU + 1), fill = DENpos)
) +
  geom_histogram(bins = 30, color = "black") +
  scale_fill_discrete(labels = c("0" = "Seronegative", "1" = "Seropositive")) +
  labs(
    title = "Manyatta DENV",
    x = "Log(AIU values + 1)",
    y = "Number of samples",
    fill = "DENV serostatus"
  ) + theme_grey(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),   # no vertical grid
    panel.grid.minor = element_blank(),     # no minor grid
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )


p2 <- ggplot(
  Kilifidata,
  aes(x = log(denv2ns1_au + 1), fill = denv_pos)
) +
  geom_histogram(bins = 30, color = "black") +
  scale_fill_discrete(labels = c("0" = "Seronegative", "1" = "Seropositive")) +
  labs(
    title = "Kilifi DENV",
    x = "Log(AIU values + 1)",
    y = "Number of samples",
    fill = "DENV serostatus"
  ) + theme_grey(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),   # no vertical grid
    panel.grid.minor = element_blank(),     # no minor grid
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )


p3 <- ggplot(
  Kiberadata,
  aes(x = log(denv2ns1_au + 1), fill = denv_pos)
) +
  geom_histogram(bins = 30, color = "black") +
  scale_fill_discrete(labels = c("0" = "Seronegative", "1" = "Seropositive")) +
  labs(
    title = "Kibera DENV",
    x = "Log(AIU values + 1)",
    y = "Number of samples",
    fill = "DENV serostatus"
  ) + theme_grey(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),   # no vertical grid
    panel.grid.minor = element_blank(),     # no minor grid
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )


p4 <- ggplot(
  NairobiUrbandata,
  aes(x = log(denv2ns1_au + 1), fill = denv_pos)
) +
  geom_histogram(bins = 30, color = "black") +
  scale_fill_discrete(labels = c("0" = "Seronegative", "1" = "Seropositive")) +
  labs(
    title = "NairobiUrban DENV",
    x = "Log(AIU values + 1)",
    y = "Number of samples",
    fill = "DENV serostatus"
  ) + theme_grey(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),   # no vertical grid
    panel.grid.minor = element_blank(),     # no minor grid
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )
gridExtra::grid.arrange(p, p1, p2, p3, p4, nrow=2)

#save
#save
ggsave("DENV_HISTOGRAM.png", 
       plot = gridExtra::grid.arrange(p, p1, p2, p3, p4, nrow=2),
       width = 12, height = 8, dpi = 300)

ggsave("CHKPLOT_Children_histogram.png", 
       width = 12, height = 8, dpi = 300)


#Finite mixture model fitting using MCLUST package
library(mclust)
#Asembo
asembo_AIU_Values <- log(Asembodata$CHKV_AIU + 1)
model_asembo <- Mclust(asembo_AIU_Values, G = 2)
summary(model_asembo)
plot(model_asembo, what = "BIC")
# Extract parameters
means <- model_asembo$parameters$mean
sds <- sqrt(model_asembo$parameters$variance$sigmasq)
props <- model_asembo$parameters$pro

# Create density curves for each component
x_vals <- seq(min(asembo_AIU_Values), max(asembo_AIU_Values), length.out = 1000)

dens_df <- data.frame(
  x = rep(x_vals, 2),
  density = c(
    dnorm(x_vals, mean = means[1], sd = sds[1]) * props[1],
    dnorm(x_vals, mean = means[2], sd = sds[2]) * props[2]
  ),
  component = factor(rep(1:2, each = length(x_vals)))
)

#cross over point
dens_diff <- function(x) {
  props[1] * dnorm(x, means[1], sds[1]) -
    props[2] * dnorm(x, means[2], sds[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(means),
  upper = max(means)
)$root


# Plot histogram + density curves
p<-ggplot() +
  geom_histogram(aes(x = asembo_AIU_Values, y = after_stat(density)), 
                 bins = 50, fill = "gray80", color = "white") +
  geom_line(data = dens_df, aes(x = x, y = density, color = component), size = 1.2) +
  labs(title = "Asembo CHKV",
       x = "log(AIU + 1)", y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() + geom_vline(
    xintercept = crossover,
    linetype = "dashed",
    linewidth = 1
  )


#Manyatta
manyatta_AIU_Values <- log(manyattadata$CHKV_AIU + 1)
model_manyatta <- Mclust(manyatta_AIU_Values, G = 2, modelNames = "V")
summary(model_manyatta)
# Extract parameters
means <- model_manyatta$parameters$mean
sds <- sqrt(model_manyatta$parameters$variance$sigmasq)
props <- model_manyatta$parameters$pro
# Create density curves for each component
x_vals <- seq(min(manyatta_AIU_Values), max(manyatta_AIU_Values), length.out = 1000)
dens_df <- data.frame(
  x = rep(x_vals, 2),
  density = c(
    dnorm(x_vals, mean = means[1], sd = sds[1]) * props[1],
    dnorm(x_vals, mean = means[2], sd = sds[2]) * props[2]
  ),
  component = factor(rep(1:2, each = length(x_vals)))
)


#cross over point
dens_diff <- function(x) {
  props[1] * dnorm(x, means[1], sds[1]) -
    props[2] * dnorm(x, means[2], sds[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(means),
  upper = max(means)
)$root

# Plot histogram + density curves
p1<-ggplot() +
  geom_histogram(aes(x = manyatta_AIU_Values, y = after_stat(density)), 
                 bins = 50, fill = "gray80", color = "white") +
  geom_line(data = dens_df, aes(x = x, y = density, color = component), size = 1.2) +
  labs(title = "Manyatta CHKV",
       x = "log(AIU + 1)", y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() + geom_vline(
    xintercept = crossover,
    linetype = "dashed",
    linewidth = 1
  )


#Kilifi
kilifi_AIU_Values <- log(Kilifidata$chkve1_au + 1)
model_kilifi <- Mclust(kilifi_AIU_Values, G = 2, modelNames = "V")
# Extract parameters
means <- model_kilifi$parameters$mean
sds <- sqrt(model_kilifi$parameters$variance$sigmasq)
props <- model_kilifi$parameters$pro
# Create density curves for each component
x_vals <- seq(min(kilifi_AIU_Values), max(kilifi_AIU_Values), length.out = 1000)
dens_df <- data.frame(
  x = rep(x_vals, 2),
  density = c(
    dnorm(x_vals, mean = means[1], sd = sds[1]) * props[1],
    dnorm(x_vals, mean = means[2], sd = sds[2]) * props[2]
  ),
  component = factor(rep(1:2, each = length(x_vals)))
)

#cross over point
dens_diff <- function(x) {
  props[1] * dnorm(x, means[1], sds[1]) -
    props[2] * dnorm(x, means[2], sds[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(means),
  upper = max(means)
)$root

# Plot histogram + density curves
p2<-ggplot() +
  geom_histogram(aes(x = kilifi_AIU_Values, y = after_stat(density)), 
                 bins = 50, fill = "gray80", color = "white") +
  geom_line(data = dens_df, aes(x = x, y = density, color = component), size = 1.2) +
  labs(title = "Kilifi CHKV",
       x = "log(AIU + 1)", y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() + geom_vline(
    xintercept = crossover,
    linetype = "dashed",
    linewidth = 1
  )

#Kibera
kibera_AIU_Values <- log(Kiberadata$chkve1_au + 1)
model_kibera <- Mclust(kibera_AIU_Values, G = 2, modelNames = "V")
# Extract parameters
means <- model_kibera$parameters$mean
sds <- sqrt(model_kibera$parameters$variance$sigmasq)
props <- model_kibera$parameters$pro
# Create density curves for each component
x_vals <- seq(min(kibera_AIU_Values), max(kibera_AIU_Values), length.out = 1000)
dens_df <- data.frame(
  x = rep(x_vals, 2),
  density = c(
    dnorm(x_vals, mean = means[1], sd = sds[1]) * props[1],
    dnorm(x_vals, mean = means[2], sd = sds[2]) * props[2]
  ),
  component = factor(rep(1:2, each = length(x_vals)))
)

#cross over point
dens_diff <- function(x) {
  props[1] * dnorm(x, means[1], sds[1]) -
    props[2] * dnorm(x, means[2], sds[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(means),
  upper = max(means)
)$root

# Plot histogram + density curves
p3<-ggplot() +
  geom_histogram(aes(x = kibera_AIU_Values, y = after_stat(density)), 
                 bins = 50, fill = "gray80", color = "white") +
  geom_line(data = dens_df, aes(x = x, y = density, color = component), size = 1.2) +
  labs(title = "Kibera CHKV",
       x = "log(AIU + 1)", y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() + geom_vline(
    xintercept = crossover,
    linetype = "dashed",
    linewidth = 1
  )

#Nairobi Urban
nairobi_AIU_Values <- log(NairobiUrbandata$chkve1_au + 1)
model_nairobi <- Mclust(nairobi_AIU_Values, G = 2, modelNames = "V")
# Extract parameters
means <- model_nairobi$parameters$mean
sds <- sqrt(model_nairobi$parameters$variance$sigmasq)
props <- model_nairobi$parameters$pro
# Create density curves for each component
x_vals <- seq(min(nairobi_AIU_Values), max(nairobi_AIU_Values), length.out = 1000)
dens_df <- data.frame(
  x = rep(x_vals, 2),
  density = c(
    dnorm(x_vals, mean = means[1], sd = sds[1]) * props[1],
    dnorm(x_vals, mean = means[2], sd = sds[2]) * props[2]
  ),
  component = factor(rep(1:2, each = length(x_vals)))
)


#cross over point
dens_diff <- function(x) {
  props[1] * dnorm(x, means[1], sds[1]) -
    props[2] * dnorm(x, means[2], sds[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(means),
  upper = max(means)
)$root



# Plot histogram + density curves
p4<- ggplot() +
  geom_histogram(aes(x = nairobi_AIU_Values, y = after_stat(density)), 
                 bins = 50, fill = "gray80", color = "white") +
  geom_line(data = dens_df, aes(x = x, y = density, color = component), size = 1.2) +
  labs(title = "Nairobi Urban CHKV",
       x = "log(AIU + 1)", y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() + geom_vline(
    xintercept = crossover,
    linetype = "dashed",
    linewidth = 1
  )
gridExtra::grid.arrange(p, p1, p2, p3, p4, nrow=2)

#save
ggsave("CHKV_Mclust_density.png", 
       plot = gridExtra::grid.arrange(p, p1, p2, p3, p4, nrow=2),
       width = 12, height = 8, dpi = 300)



#MIXR FIT package
library(mixR)
#Asembo
mfi_values_asembo <- Asembodata$CHKV_AIU
mfi_values_positive <- mfi_values_asembo[mfi_values_asembo > 0]
log_asembo <- log(mfi_values_positive + 1)
x_data <- log_asembo

asembo_model <- mixfit(
  log_asembo,
  family = "gamma",   
  ncomp = 2)

print(asembo_model)

pi <- asembo_model$pi
mu <- asembo_model$mu
sd <- asembo_model$sd

#identifiability check
muI<- max(mu)
muS<- min(mu)
SigmaI<- sd[2]
SigmaS<- sd[1]

sep <- abs(muI - muS) / sqrt(SigS^2 + SigI^2)
sep





x_vals <- seq(min(x_data), max(x_data), length.out = 1000)

dens1 <- pi[1] * dnorm(x_vals, mu[1], sd[1])
dens2 <- pi[2] * dnorm(x_vals, mu[2], sd[2])
dens_mix <- dens1 + dens2



dens_diff <- function(x) {
  pi[1] * dnorm(x, mu[1], sd[1]) -
    pi[2] * dnorm(x, mu[2], sd[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(mu),
  upper = max(mu)
)$root

seropositive <- log_asembo > crossover
seroprevalence <- mean(seropositive, na.rm = TRUE)
cat("Estimated seroprevalence based on crossover cutoff:", round(seroprevalence * 100, 2), "%\n")



tiff(
  filename = "Asembo_DENV.normal.tiff",
  width = 4.21,      
  height = 3.8,
  units = "in",
  res = 300
)

par(
  family = "Times",
  cex = 0.83,           
  mar = c(5, 5, 4, 6)      
)

hist(
  x_data,
  breaks = 30,
  freq = FALSE,
  col = "white",
  border = "grey60",
  main = "Asembo DENV Normal",
  xlab = "log(AIU + 1)",
  ylab = "Density",
  font.main = 2           
)


polygon(
  c(x_vals, rev(x_vals)),
  c(dens1, rep(0, length(dens1))),
  col = rgb(1, 0, 0, 0.35),
  border = NA
)


polygon(
  c(x_vals, rev(x_vals)),
  c(dens2, rep(0, length(dens2))),
  col = rgb(0, 0.6, 0.6, 0.35),
  border = NA
)


lines(x_vals, dens_mix, lwd = 2, col = "black")


abline(v = crossover, lty = 2, lwd = 2, col = "black")


text(
  x = crossover,
  y = max(dens_mix) * 0.9,
  labels = paste0("Crossover = ", round(crossover, 2)),
  pos = 4,
  cex = 0.83
)


legend(
  "topright",
  inset = c(-0.35, 0),
  legend = c("Seronegative", "Seropositive", "Mixture", "Crossover line"),
  fill = c(
    rgb(1, 0, 0, 0.35),
    rgb(0, 0.6, 0.6, 0.35),
    NA, NA
  ),
  lty = c(NA, NA, 1, 2),
  lwd = c(NA, NA, 2, 2),
  col = c(NA, NA, "black", "black"),
  border = NA,
  bty = "n",
  xpd = TRUE,
  cex = 0.83
)

dev.off()

#Manyatta
mfi_values_manyatta<- manyattadata$CHKV_AIU
mfi_values_manyatta_positive<- mfi_values_manyatta[mfi_values_manyatta >0]
log_manyatta <- log(mfi_values_manyatta_positive + 1)
x_data <- log_manyatta
manyatta_model <- mixfit(log_manyatta, family = "gamma", ncomp = 2)
print(manyatta_model)


pi <- manyatta_model$pi
mu <- manyatta_model$mu
sd <- manyatta_model$sd

#identifiability check
muI<- max(mu)
muS<- min(mu)
SigmaI<- sd[2]
SigmaS<- sd[1]

sep <- abs(muI - muS) / sqrt(SigS^2 + SigI^2)
sep



x_vals <- seq(min(x_data), max(x_data), length.out = 1000)

dens1 <- pi[1] * dnorm(x_vals, mu[1], sd[1])
dens2 <- pi[2] * dnorm(x_vals, mu[2], sd[2])
dens_mix <- dens1 + dens2



dens_diff <- function(x) {
  pi[1] * dnorm(x, mu[1], sd[1]) -
    pi[2] * dnorm(x, mu[2], sd[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(mu),
  upper = max(mu)
)$root

seropositive <- log_manyatta > crossover
seroprevalence <- mean(seropositive, na.rm = TRUE)
cat("Estimated seroprevalence based on crossover cutoff:", round(seroprevalence * 100, 2), "%\n")

tiff(
  filename = "Manyatta_CHKV_normal.tiff",
  width = 4.21,      
  height = 3.8,
  units = "in",
  res = 300
)

par(
  family = "Times",
  cex = 0.83,           
  mar = c(5, 5, 4, 6)      
)

hist(
  x_data,
  breaks = 40,
  freq = FALSE,
  col = "white",
  border = "grey60",
  main = "Manyatta CHKV Normal",
  xlab = "log(AIU + 1)",
  ylab = "Density",
  font.main = 2           
)


polygon(
  c(x_vals, rev(x_vals)),
  c(dens1, rep(0, length(dens1))),
  col = rgb(1, 0, 0, 0.35),
  border = NA
)


polygon(
  c(x_vals, rev(x_vals)),
  c(dens2, rep(0, length(dens2))),
  col = rgb(0, 0.6, 0.6, 0.35),
  border = NA
)


lines(x_vals, dens_mix, lwd = 2, col = "black")


abline(v = crossover, lty = 2, lwd = 2, col = "black")


text(
  x = crossover,
  y = max(dens_mix) * 0.9,
  labels = paste0("Crossover = ", round(crossover, 2)),
  pos = 4,
  cex = 0.83
)


legend(
  "topright",
  inset = c(-0.35, 0),
  legend = c("Seronegative", "Seropositive", "Mixture", "Crossover line"),
  fill = c(
    rgb(1, 0, 0, 0.35),
    rgb(0, 0.6, 0.6, 0.35),
    NA, NA
  ),
  lty = c(NA, NA, 1, 2),
  lwd = c(NA, NA, 2, 2),
  col = c(NA, NA, "black", "black"),
  border = NA,
  bty = "n",
  xpd = TRUE,
  cex = 0.83
)

dev.off()

#Kilifi
mfi_values_kilifi<- Kilifidata$chkve1_au
mfi_values_kilifi_positive<- mfi_values_kilifi[mfi_values_kilifi >0]
log_kilifi <- log(mfi_values_kilifi_positive + 1)
x_data <- log_kilifi

kilifi_model <- mixfit(log_kilifi, family = "weibull", ncomp = 2)
print(kilifi_model)

pi <- kilifi_model$pi
mu <- kilifi_model$mu
sd <- kilifi_model$sd

#identifiability check
muI<- max(mu)
muS<- min(mu)
SigmaI<- sd[2]
SigmaS<- sd[1]

sep <- abs(muI - muS) / sqrt(SigS^2 + SigI^2)
sep



x_vals <- seq(min(x_data), max(x_data), length.out = 1000)

dens1 <- pi[1] * dnorm(x_vals, mu[1], sd[1])
dens2 <- pi[2] * dnorm(x_vals, mu[2], sd[2])
dens_mix <- dens1 + dens2



dens_diff <- function(x) {
  pi[1] * dnorm(x, mu[1], sd[1]) -
    pi[2] * dnorm(x, mu[2], sd[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(mu),
  upper = max(mu)
)$root

seropositive <- log_kilifi > crossover
seroprevalence <- mean(seropositive, na.rm = TRUE)
cat("Estimated seroprevalence based on crossover cutoff:", round(seroprevalence * 100, 2), "%\n")

tiff(
  filename = "Kilifi_DENV_gamma.tiff",
  width = 4.21,      
  height = 3.8,
  units = "in",
  res = 300
)

par(
  family = "Times",
  cex = 0.83,           
  mar = c(5, 5, 4, 6)      
)

hist(
  x_data,
  breaks = 40,
  freq = FALSE,
  col = "white",
  border = "grey60",
  main = "Kilifi DENV Gamma",
  xlab = "log(AIU + 1)",
  ylab = "Density",
  font.main = 2           
)


polygon(
  c(x_vals, rev(x_vals)),
  c(dens1, rep(0, length(dens1))),
  col = rgb(1, 0, 0, 0.35),
  border = NA
)


polygon(
  c(x_vals, rev(x_vals)),
  c(dens2, rep(0, length(dens2))),
  col = rgb(0, 0.6, 0.6, 0.35),
  border = NA
)


lines(x_vals, dens_mix, lwd = 2, col = "black")


abline(v = crossover, lty = 2, lwd = 2, col = "black")


text(
  x = crossover,
  y = max(dens_mix) * 0.9,
  labels = paste0("Crossover = ", round(crossover, 2)),
  pos = 4,
  cex = 0.83
)


legend(
  "topright",
  inset = c(-0.35, 0),
  legend = c("Seronegative", "Seropositive", "Mixture", "Crossover line"),
  fill = c(
    rgb(1, 0, 0, 0.35),
    rgb(0, 0.6, 0.6, 0.35),
    NA, NA
  ),
  lty = c(NA, NA, 1, 2),
  lwd = c(NA, NA, 2, 2),
  col = c(NA, NA, "black", "black"),
  border = NA,
  bty = "n",
  xpd = TRUE,
  cex = 0.83
)

dev.off()

#Kibera
mfi_values_kibera<- Kiberadata$chkve1_au
mfi_values_kibera_positive<- mfi_values_kibera[mfi_values_kibera >0]
log_kibera <- log(mfi_values_kibera_positive + 1)
x_data <- log_kibera

Kibera_model <- mixfit(log_kibera, family = "gamma", ncomp = 2)
print(Kibera_model)


pi <- Kibera_model$pi
mu <- Kibera_model$mu
sd <- Kibera_model$sd

#identifiability check
muI<- max(mu)
muS<- min(mu)
SigmaI<- sd[2]
SigmaS<- sd[1]

sep <- abs(muI - muS) / sqrt(SigS^2 + SigI^2)
sep



x_vals <- seq(min(x_data), max(x_data), length.out = 1000)

dens1 <- pi[1] * dnorm(x_vals, mu[1], sd[1])
dens2 <- pi[2] * dnorm(x_vals, mu[2], sd[2])
dens_mix <- dens1 + dens2



dens_diff <- function(x) {
  pi[1] * dnorm(x, mu[1], sd[1]) -
    pi[2] * dnorm(x, mu[2], sd[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(mu),
  upper = max(mu)
)$root

seropositive <- log_kibera > crossover
seroprevalence <- mean(seropositive, na.rm = TRUE)
cat("Estimated seroprevalence based on crossover cutoff:", round(seroprevalence * 100, 2), "%\n")

tiff(
  filename = "Kibera_CHKV_Normal.tiff",
  width = 4.21,      
  height = 3.8,
  units = "in",
  res = 300
)

par(
  family = "Times",
  cex = 0.83,           
  mar = c(5, 5, 4, 6)      
)

hist(
  x_data,
  breaks = 40,
  freq = FALSE,
  col = "white",
  border = "grey60",
  main = "Kibera CHKV Normal",
  xlab = "log(AIU + 1)",
  ylab = "Density",
  font.main = 2           
)


polygon(
  c(x_vals, rev(x_vals)),
  c(dens1, rep(0, length(dens1))),
  col = rgb(1, 0, 0, 0.35),
  border = NA
)


polygon(
  c(x_vals, rev(x_vals)),
  c(dens2, rep(0, length(dens2))),
  col = rgb(0, 0.6, 0.6, 0.35),
  border = NA
)


lines(x_vals, dens_mix, lwd = 2, col = "black")


abline(v = crossover, lty = 2, lwd = 2, col = "black")


text(
  x = crossover,
  y = max(dens_mix) * 0.9,
  labels = paste0("Crossover = ", round(crossover, 2)),
  pos = 4,
  cex = 0.83
)


legend(
  "topright",
  inset = c(-0.35, 0),
  legend = c("Seronegative", "Seropositive", "Mixture", "Crossover line"),
  fill = c(
    rgb(1, 0, 0, 0.35),
    rgb(0, 0.6, 0.6, 0.35),
    NA, NA
  ),
  lty = c(NA, NA, 1, 2),
  lwd = c(NA, NA, 2, 2),
  col = c(NA, NA, "black", "black"),
  border = NA,
  bty = "n",
  xpd = TRUE,
  cex = 0.83
)

dev.off()


#Nairobi Urban
mfi_values_nairobi <- NairobiUrbandata$chkve1_au
mfi_values_nairobi_positive <- mfi_values_nairobi[mfi_values_nairobi > 0]
x_data <- log(mfi_values_nairobi_positive + 1)


Nairobi_model <- mixfit(x_data, family = "gamma", ncomp = 2)
print(Nairobi_model)

pi <- Nairobi_model$pi
mu <- Nairobi_model$mu
sd <- Nairobi_model$sd

#identifiability check
muI<- max(mu)
muS<- min(mu)
SigmaI<- sd[2]
SigmaS<- sd[1]

sep <- abs(muI - muS) / sqrt(SigS^2 + SigI^2)
sep



x_vals <- seq(min(x_data), max(x_data), length.out = 1000)

dens1 <- pi[1] * dnorm(x_vals, mu[1], sd[1])
dens2 <- pi[2] * dnorm(x_vals, mu[2], sd[2])
dens_mix <- dens1 + dens2



dens_diff <- function(x) {
  pi[1] * dnorm(x, mu[1], sd[1]) -
    pi[2] * dnorm(x, mu[2], sd[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(mu),
  upper = max(mu)
)$root

seropositive <- log_nairobi > crossover
seroprevalence <- mean(seropositive, na.rm = TRUE)
cat("Estimated seroprevalence based on crossover cutoff:", round(seroprevalence * 100, 2), "%\n")

tiff(
  filename = "Nairobi_CHKV_Normal.tiff",
  width = 4.21,      
  height = 3.8,
  units = "in",
  res = 300
)

par(
  family = "Times",
  cex = 0.83,           
  mar = c(5, 5, 4, 6)      
)

hist(
  x_data,
  breaks = 40,
  freq = FALSE,
  col = "white",
  border = "grey60",
  main = "Nairobi CHKV Normal",
  xlab = "log(AIU + 1)",
  ylab = "Density",
  font.main = 2           
)


polygon(
  c(x_vals, rev(x_vals)),
  c(dens1, rep(0, length(dens1))),
  col = rgb(1, 0, 0, 0.35),
  border = NA
)


polygon(
  c(x_vals, rev(x_vals)),
  c(dens2, rep(0, length(dens2))),
  col = rgb(0, 0.6, 0.6, 0.35),
  border = NA
)


lines(x_vals, dens_mix, lwd = 2, col = "black")


abline(v = crossover, lty = 2, lwd = 2, col = "black")


text(
  x = crossover,
  y = max(dens_mix) * 0.9,
  labels = paste0("Crossover = ", round(crossover, 2)),
  pos = 4,
  cex = 0.83
)


legend(
  "topright",
  inset = c(-0.35, 0),
  legend = c("Seronegative", "Seropositive", "Mixture", "Crossover line"),
  fill = c(
    rgb(1, 0, 0, 0.35),
    rgb(0, 0.6, 0.6, 0.35),
    NA, NA
  ),
  lty = c(NA, NA, 1, 2),
  lwd = c(NA, NA, 2, 2),
  col = c(NA, NA, "black", "black"),
  border = NA,
  bty = "n",
  xpd = TRUE,
  cex = 0.83
)

dev.off()

den_asembo<- Asembodata$DENV2_AIU
den_asembo_positive<- den_asembo[den_asembo >0]
log_asembo <- log(den_asembo_positive + 1)
asembo_model <- mixfit(
  log_asembo,
  family = "normal",   
  ncomp = 2
)
print(asembo_model)
bic= select(den_asembo_positive,ncomp=2:6, family="weibull")
plot(bic)

b1=bs.test(den_asembo_positive, ncomp=c(2,4), family="weibull")

b1$p.value





#MixSMSN package
library("mixsmsn")
#Asembo CHKV
y_asembo <- log(Asembodata$CHKV_AIU + 1)

model_asembo <- smsn.mix(
  y_asembo,
  nu = 3,
  g = 2,
  family = "Skew.normal",
  get.init = TRUE,
  criteria = TRUE,
  group = TRUE
)





# Grid
x_vals <- seq(min(y_log), max(y_log), length.out = 1000)

# Component densities
dens1 <- pi[1] * sn::dsn(x_vals, mu[1], sd[1], lam[1])
dens2 <- pi[2] * sn::dsn(x_vals, mu[2], sd[2], lam[2])
dens_mix <- dens1 + dens2

# Crossover
dens_diff <- function(x) {
  pi[1] * sn::dsn(x, mu[1], sd[1], lam[1]) -
    pi[2] * sn::dsn(x, mu[2], sd[2], lam[2])
}

crossover <- uniroot(
  dens_diff,
  lower = min(mu),
  upper = max(mu)
)$root









library(smsn)
library(sn)
#Asembo CHKV
y_asembo <- log(Asembodata$CHKV_AIU + 1)


model_asembo <- smsn.mix(
  y_asembo,
  nu = 3,
  g = 2,
  family = "Skew.normal",
  get.init = TRUE,
  criteria = TRUE,
  group = TRUE
)

mix.print(model_asembo)


pi  <- model_asembo$pii
mu  <- model_asembo$mu
sd  <- model_asembo$sigma
lam <- model_asembo$shape


x_vals <- seq(min(y_asembo), max(y_asembo), length.out = 2000)

dens1 <- pi[1] * sn::dsn(x_vals, mu[1], sd[1], lam[1])
dens2 <- pi[2] * sn::dsn(x_vals, mu[2], sd[2], lam[2])
dens_mix <- dens1 + dens2


dens_diff <- function(x) {
  pi[1] * sn::dsn(x, mu[1], sd[1], lam[1]) -
    pi[2] * sn::dsn(x, mu[2], sd[2], lam[2])
}

diff_vals <- dens_diff(x_vals)
idx <- which(diff_vals[-1] * diff_vals[-length(diff_vals)] < 0)

if (length(idx) == 0) {
  stop("No crossover found between mixture components.")
}

crossover <- uniroot(
  dens_diff,
  lower = x_vals[idx[1]],
  upper = x_vals[idx[1] + 1]
)$root


seropositive <- y_asembo > crossover
seroprevalence <- mean(seropositive, na.rm = TRUE)
cat("Estimated seroprevalence (skew-normal):", round(seroprevalence * 100, 2), "%\n")


tiff(
  filename = "Asembo_CHKV.skewnormal.tiff",
  width = 4.21,
  height = 3.8,
  units = "in",
  res = 300
)

par(
  family = "Times",
  cex = 0.83,
  mar = c(5, 5, 4, 6)
)

y_max <- max(
  dens_mix,
  hist(y_asembo, breaks = 40, plot = FALSE)$density
) * 1.1
# Histogram
hist(
  y_asembo,
  breaks = 40,
  freq = FALSE,
  col = "white",
  border = "grey60",
  main = "Asembo CHKV Skew-normal",
  xlab = "log(AIU + 1)",
  ylab = "Density",
  font.main = 2,
  ylim = c(0, y_max)
)

# Seronegative component
polygon(
  c(x_vals, rev(x_vals)),
  c(dens1, rep(0, length(dens1))),
  col = rgb(1, 0, 0, 0.35),
  border = NA
)

# Seropositive component
polygon(
  c(x_vals, rev(x_vals)),
  c(dens2, rep(0, length(dens2))),
  col = rgb(0, 0.6, 0.6, 0.35),
  border = NA
)

# Mixture density
lines(x_vals, dens_mix, lwd = 2, col = "black")

# Crossover line
abline(v = crossover, lty = 2, lwd = 2, col = "black")

# Crossover label
text(
  x = crossover,
  y = max(dens_mix) * 0.9,
  labels = paste0("Crossover = ", round(crossover, 2)),
  pos = 4,
  cex = 0.83
)

# Legend
legend(
  "topright",
  inset = c(-0.35, 0),
  legend = c(
    "Seronegative",
    "Seropositive",
    "Mixture",
    "Crossover line"
  ),
  fill = c(
    rgb(1, 0, 0, 0.35),
    rgb(0, 0.6, 0.6, 0.35),
    NA,
    NA
  ),
  lty = c(NA, NA, 1, 2),
  lwd = c(NA, NA, 2, 2),
  col = c(NA, NA, "black", "black"),
  border = NA,
  bty = "n",
  xpd = TRUE,
  cex = 0.83
)

dev.off()

#Manyatta CHKV
y_manyatta <- log(manyattadata$CHKV_AIU + 1)
model_manyatta <- smsn.mix(
  y_manyatta,
  nu = 3,
  g = 2,
  family = "Skew.normal",
  get.init = TRUE,
  criteria = TRUE,
  group = TRUE
)
mix.print(model_manyatta)

pi  <- model_manyatta$pii
mu  <- model_manyatta$mu
sd  <- model_manyatta$sigma
lam <- model_manyatta$shape
# Grid
x_vals <- seq(min(y_manyatta), max(y_manyatta), length.out = 1000)
# Component densities
dens1 <- pi[1] * sn::dsn(x_vals, mu[1], sd[1], lam[1])
dens2 <- pi[2] * sn::dsn(x_vals, mu[2], sd[2], lam[2])
dens_mix <- dens1 + dens2
# Crossover
dens_diff <- function(x) {
  pi[1] * sn::dsn(x, mu[1], sd[1], lam[1]) -
    pi[2] * sn::dsn(x, mu[2], sd[2], lam[2])
}
crossover <- uniroot(
  dens_diff,
  lower = min(mu),
  upper = max(mu)
)$root
seropositive <- y_manyatta > crossover
seroprevalence <- mean(seropositive, na.rm = TRUE)
cat("Estimated seroprevalence (skew-normal):", round(seroprevalence * 100, 2), "%\n")

tiff(
  filename = "Manyatta_CHKV.skewnormal.tiff",
  width = 4.21,
  height = 3.8,
  units = "in",
  res = 300
)

par(
  family = "Times",
  cex = 0.83,
  mar = c(5, 5, 4, 6)
)

y_max <- max(
  dens_mix,
  hist(y_manyatta, breaks = 40, plot = FALSE)$density
) * 1.1
# Histogram
hist(
  y_manyatta,
  breaks = 40,
  freq = FALSE,
  col = "white",
  border = "grey60",
  main = "Manyatta CHKV Skew-normal",
  xlab = "log(AIU + 1)",
  ylab = "Density",
  font.main = 2,
  ylim = c(0, y_max)
)

# Seronegative component
polygon(
  c(x_vals, rev(x_vals)),
  c(dens1, rep(0, length(dens1))),
  col = rgb(1, 0, 0, 0.35),
  border = NA
)

# Seropositive component
polygon(
  c(x_vals, rev(x_vals)),
  c(dens2, rep(0, length(dens2))),
  col = rgb(0, 0.6, 0.6, 0.35),
  border = NA
)

# Mixture density
lines(x_vals, dens_mix, lwd = 2, col = "black")

# Crossover line
abline(v = crossover, lty = 2, lwd = 2, col = "black")

# Crossover label
text(
  x = crossover,
  y = max(dens_mix) * 0.9,
  labels = paste0("Crossover = ", round(crossover, 2)),
  pos = 4,
  cex = 0.83
)

# Legend
legend(
  "topright",
  inset = c(-0.35, 0),
  legend = c(
    "Seronegative",
    "Seropositive",
    "Mixture",
    "Crossover line"
  ),
  fill = c(
    rgb(1, 0, 0, 0.35),
    rgb(0, 0.6, 0.6, 0.35),
    NA,
    NA
  ),
  lty = c(NA, NA, 1, 2),
  lwd = c(NA, NA, 2, 2),
  col = c(NA, NA, "black", "black"),
  border = NA,
  bty = "n",
  xpd = TRUE,
  cex = 0.83
)

dev.off()


# -------------------------------
# Kilifi CHKV skew-normal mixture
# -------------------------------

library(smsn)
library(sn)

# Log-transform data
y_kilifi <- log(Kilifidata$denv2ns1_au + 1)

# Fit skew-normal mixture
model_kilifi <- smsn.mix(
  y_kilifi,
  nu = 3,
  g = 2,
  family = "Skew.normal",
  get.init = TRUE,
  criteria = TRUE,
  group = TRUE
)

mix.print(model_kilifi)

# Extract parameters
pi  <- model_kilifi$pii
mu  <- model_kilifi$mu
sd  <- model_kilifi$sigma
lam <- model_kilifi$shape

# -------------------------------------------------
# FIX: Order components by mean (mu)
# -------------------------------------------------
ord <- order(mu)   # low mean first = seronegative

pi  <- pi[ord]
mu  <- mu[ord]
sd  <- sd[ord]
lam <- lam[ord]

# Grid
x_vals <- seq(min(y_kilifi), max(y_kilifi), length.out = 1000)

# Component densities
dens_neg <- pi[1] * sn::dsn(x_vals, mu[1], sd[1], lam[1])  # seronegative
dens_pos <- pi[2] * sn::dsn(x_vals, mu[2], sd[2], lam[2])  # seropositive
dens_mix <- dens_neg + dens_pos

# -------------------------------------------------
# Crossover calculation
# -------------------------------------------------
dens_diff <- function(x) {
  pi[1] * sn::dsn(x, mu[1], sd[1], lam[1]) -
    pi[2] * sn::dsn(x, mu[2], sd[2], lam[2])
}

crossover <- uniroot(
  dens_diff,
  lower = mu[1],
  upper = mu[2]
)$root

# Seroprevalence
seropositive <- y_kilifi > crossover
seroprevalence <- mean(seropositive, na.rm = TRUE)

cat("Estimated seroprevalence (skew-normal):",round(seroprevalence * 100, 2), "%\n")

# -------------------------------------------------
# Plot
# -------------------------------------------------

tiff(
  filename = "Kilifi_DENV.skewnormal.tiff",
  width = 4.21,
  height = 3.8,
  units = "in",
  res = 300
)

par(
  family = "Times",
  cex = 0.83,
  mar = c(5, 5, 4, 6),
  xpd = NA
)

y_max <- max(
  dens_mix,
  hist(y_kilifi, breaks = 40, plot = FALSE)$density
) * 1.1

# Histogram
hist(
  y_kilifi,
  breaks = 40,
  freq = FALSE,
  col = "white",
  border = "grey60",
  main = "Kilifi DENV Skew-normal",
  xlab = "log(AIU + 1)",
  ylab = "Density",
  font.main = 2,
  ylim = c(0, y_max)
)

# Seronegative component (LOW)
polygon(
  c(x_vals, rev(x_vals)),
  c(dens_neg, rep(0, length(dens_neg))),
  col = rgb(1, 0, 0, 0.35),
  border = NA
)

# Seropositive component (HIGH)
polygon(
  c(x_vals, rev(x_vals)),
  c(dens_pos, rep(0, length(dens_pos))),
  col = rgb(0, 0.6, 0.6, 0.35),
  border = NA
)

# Mixture density
lines(x_vals, dens_mix, lwd = 2, col = "black")

# Crossover
abline(v = crossover, lty = 2, lwd = 2, col = "black")

text(
  x = crossover,
  y = max(dens_mix) * 0.9,
  labels = paste0("Crossover = ", round(crossover, 2)),
  pos = 4,
  cex = 0.83
)

# Legend
legend(
  "topright",
  inset = c(-0.35, 0),
  legend = c(
    "Seronegative",
    "Seropositive",
    "Mixture",
    "Crossover line"
  ),
  fill = c(
    rgb(1, 0, 0, 0.35),
    rgb(0, 0.6, 0.6, 0.35),
    NA,
    NA
  ),
  lty = c(NA, NA, 1, 2),
  lwd = c(NA, NA, 2, 2),
  col = c(NA, NA, "black", "black"),
  border = NA,
  bty = "n",
  cex = 0.83
)

dev.off()


#Kibera CHKV
y_kibera <- log(Kiberadata$chkve1_au + 1)
model_kibera <- smsn.mix(
  y_kibera,
  nu = 3,
  g = 2,
  family = "Skew.normal",
  get.init = TRUE,
  criteria = TRUE,
  group = TRUE
)
mix.print(model_kibera)

pi  <- model_kibera$pii
mu  <- model_kibera$mu
sd  <- model_kibera$sigma
lam <- model_kibera$shape


# -------------------------------------------------
# FIX: Order components by mean (mu)
# -------------------------------------------------
ord <- order(mu)   # low mean first = seronegative

pi  <- pi[ord]
mu  <- mu[ord]
sd  <- sd[ord]
lam <- lam[ord]

# Grid
x_vals <- seq(min(y_kibera), max(y_kibera), length.out = 1000)
# Component densities
dens1 <- pi[1] * sn::dsn(x_vals, mu[1], sd[1], lam[1])
dens2 <- pi[2] * sn::dsn(x_vals, mu[2], sd[2], lam[2])
dens_mix <- dens1 + dens2
# Crossover
dens_diff <- function(x) {
  pi[1] * sn::dsn(x, mu[1], sd[1], lam[1]) -
    pi[2] * sn::dsn(x, mu[2], sd[2], lam[2])
}
crossover <- uniroot(
  dens_diff,
  lower = min(mu),
  upper = max(mu)
)$root
seropositive <- y_kibera > crossover
seroprevalence <- mean(seropositive, na.rm = TRUE)
cat("Estimated seroprevalence (skew-normal):", round(seroprevalence * 100, 2), "%\n")
tiff(
  filename = "Kibera_CHKV.skewnormal.tiff",
  width = 4.21,
  height = 3.8,
  units = "in",
  res = 300
)   
par(
  family = "Times",
  cex = 0.83,
  mar = c(5, 5, 4, 6)
)
y_max <- max(
  dens_mix,
  hist(y_kibera, breaks = 40, plot = FALSE)$density
) * 1.1
# Histogram
hist(
  y_kibera,
  breaks = 40,
  freq = FALSE,
  col = "white",
  border = "grey60",
  main = "Kibera CHKV Skew-normal",
  xlab = "log(AIU + 1)",
  ylab = "Density",
  font.main = 2,
  ylim = c(0, y_max)
)
# Seronegative component
polygon(
  c(x_vals, rev(x_vals)),
  c(dens1, rep(0, length(dens1))),
  col = rgb(1, 0, 0, 0.35),
  border = NA
)
# Seropositive component
polygon(
  c(x_vals, rev(x_vals)),
  c(dens2, rep(0, length(dens2))),
  col = rgb(0, 0.6, 0.6, 0.35),
  border = NA
)
# Mixture density
lines(x_vals, dens_mix, lwd = 2, col = "black")
# Crossover line
abline(v = crossover, lty = 2, lwd = 2, col = "black")
# Crossover label
text(
  x = crossover,
  y = max(dens_mix) * 0.9,
  labels = paste0("Crossover = ", round(crossover, 2)),
  pos = 4,
  cex = 0.83
)
# Legend
legend(
  "topright",
  inset = c(-0.35, 0),
  legend = c(
    "Seronegative",
    "Seropositive",
    "Mixture",
    "Crossover line"
  ),
  fill = c(
    rgb(1, 0, 0, 0.35),
    rgb(0, 0.6, 0.6, 0.35),
    NA,
    NA
  ),
  lty = c(NA, NA, 1, 2),
  lwd = c(NA, NA, 2, 2),
  col = c(NA, NA, "black", "black"),
  border = NA,
  bty = "n",
  xpd = TRUE,
  cex = 0.83
)
dev.off()

#Nairobi Urban CHKV
y_nairobi <- log(NairobiUrbandata$chkve1_au + 1)
model_nairobi <- smsn.mix(
  y_nairobi,
  nu = 3,
  g = 2,
  family = "Skew.normal",
  get.init = TRUE,
  criteria = TRUE,
  group = TRUE
)
mix.print(model_nairobi)
pi  <- model_nairobi$pii
mu  <- model_nairobi$mu
sd  <- model_nairobi$sigma
lam <- model_nairobi$shape

ord <- order(mu)   # low mean first = seronegative

pi  <- pi[ord]
mu  <- mu[ord]
sd  <- sd[ord]
lam <- lam[ord]
# Grid
x_vals <- seq(min(y_nairobi), max(y_nairobi), length.out = 1000)
# Component densities
dens1 <- pi[1] * sn::dsn(x_vals, mu[1], sd[1], lam[1])
dens2 <- pi[2] * sn::dsn(x_vals, mu[2], sd[2], lam[2])
dens_mix <- dens1 + dens2
# Crossover
dens_diff <- function(x) {
  pi[1] * sn::dsn(x, mu[1], sd[1], lam[1]) -
    pi[2] * sn::dsn(x, mu[2], sd[2], lam[2])
}
crossover <- uniroot(
  dens_diff,
  lower = min(mu),
  upper = max(mu)
)$root
seropositive <- y_nairobi > crossover
seroprevalence <- mean(seropositive, na.rm = TRUE)
cat("Estimated seroprevalence (skew-normal):", round(seroprevalence * 100, 2), "%\n")
tiff(
  filename = "Nairobi_CHKV.skewnormal.tiff",
  width = 4.21,
  height = 3.8,
  units = "in",
  res = 300
)
par(
  family = "Times",
  cex = 0.83,
  mar = c(5, 5, 4, 6)
)
y_max <- max(
  dens_mix,
  hist(y_nairobi, breaks = 40, plot = FALSE)$density
) * 1.1
# Histogram
hist(
  y_nairobi,
  breaks = 40,
  freq = FALSE,
  col = "white",
  border = "grey60",
  main = "Nairobi CHKV Skew-normal",
  xlab = "log(AIU + 1)",
  ylab = "Density",
  font.main = 2,
  ylim = c(0, y_max)
)
# Seronegative component
polygon(
  c(x_vals, rev(x_vals)),
  c(dens1, rep(0, length(dens1))),
  col = rgb(1, 0, 0, 0.35),
  border = NA
)
# Seropositive component
polygon(
  c(x_vals, rev(x_vals)),
  c(dens2, rep(0, length(dens2))),
  col = rgb(0, 0.6, 0.6, 0.35),
  border = NA
)
# Mixture density
lines(x_vals, dens_mix, lwd = 2, col = "black")
# Crossover line
abline(v = crossover, lty = 2, lwd = 2, col = "black")
# Crossover label     
text(
  x = crossover,
  y = max(dens_mix) * 0.9,
  labels = paste0("Crossover = ", round(crossover, 2)),
  pos = 4,
  cex = 0.83
)
# Legend
legend(
  "topright",
  inset = c(-0.35, 0),
  legend = c(
    "Seronegative",
    "Seropositive",
    "Mixture",
    "Crossover line"
  ),
  fill = c(
    rgb(1, 0, 0, 0.35),
    rgb(0, 0.6, 0.6, 0.35),
    NA,
    NA
  ),
  lty = c(NA, NA, 1, 2),
  lwd = c(NA, NA, 2, 2),
  col = c(NA, NA, "black", "black"),
  border = NA,
  bty = "n",
  xpd = TRUE,
  cex = 0.83
)
dev.off()







############################

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

kilifidata <- Kilifidata %>% 
  mutate(agecat = case_when(
    age_y >= 0 & age_y < 5 ~ "<5",
    age_y >= 5 & age_y <15 ~ "5-14",
    age_y >= 15 & age_y <45 ~ "15-44",
    age_y >= 45 & age_y <65 ~ "45-64",
    age_y >= 65 ~ "65+",
    TRUE ~ NA_character_
    
  ))
kilifidata$agecat <- factor(
  kilifidata$agecat,
  levels = c("<5", "5-14", "15-44", "45-64", "65+"),
  ordered = FALSE
)
table(kilifidata$agecat)


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
library(brms)

# Asembo
Asembodata_clean <- Asembodata %>%
  select(site, sex, ageyrs, agecat, CHKpos, DENpos, RVFpos) 

# Manyatta
manyattadata_clean <- manyattadata %>%
  select(site, sex, ageyrs, agecat, CHKpos, DENpos, RVFpos)

# Kilifi
kilifidata_clean <- kilifidata %>%
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

# Nairobi_clean (female labeled "Female")
Nairobi_clean$sex <- factor(Nairobi_clean$sex)
Nairobi_clean$sex <- relevel(Nairobi_clean$sex, ref = "Female")
Nairobi_clean$agecat <- factor(Nairobi_clean$agecat, levels = c("15-44", setdiff(levels(Nairobi_clean$agecat), "15-44")))

Nairobi_clean$CHKpos <- as.factor(Nairobi_clean$CHKpos)
Nairobi_clean$DENpos <- as.factor(Nairobi_clean$DENpos)
Nairobi_clean$RVFpos <- as.factor(Nairobi_clean$RVFpos)
Nairobi_clean$agecat <- factor(Nairobi_clean$agecat, ordered = FALSE)


# Asembodata_clean (female labeled "F")
Asembodata_clean$sex <- factor(Asembodata_clean$sex)
Asembodata_clean$sex <- relevel(Asembodata_clean$sex, ref = "F")
Asembodata_clean$agecat <- factor(Asembodata_clean$agecat, levels = c("15-44", setdiff(levels(Asembodata_clean$agecat), "15-44")))

Asembodata_clean$CHKpos <- as.factor(Asembodata_clean$CHKpos)
Asembodata_clean$DENpos <- as.factor(Asembodata_clean$DENpos)
Asembodata_clean$RVFpos <- as.factor(Asembodata_clean$RVFpos)
Asembodata_clean$agecat <- factor(Asembodata_clean$agecat, ordered = FALSE)


# manyattadata_clean (female labeled "F")
manyattadata_clean$sex <- factor(manyattadata_clean$sex)
manyattadata_clean$sex <- relevel(manyattadata_clean$sex, ref = "F")
manyattadata_clean$agecat <- factor(manyattadata_clean$agecat, levels = c("15-44", setdiff(levels(manyattadata_clean$agecat), "15-44")))

manyattadata_clean$CHKpos <- as.factor(manyattadata_clean$CHKpos)
manyattadata_clean$DENpos <- as.factor(manyattadata_clean$DENpos)
manyattadata_clean$RVFpos <- as.factor(manyattadata_clean$RVFpos)
manyattadata_clean$agecat <- factor(manyattadata_clean$agecat, ordered = FALSE)


# Kilifidata_clean (female labeled "f")
kilifidata_clean$sex <- factor(kilifidata_clean$sex)
kilifidata_clean$sex <- relevel(kilifidata_clean$sex, ref = "f")
kilifidata_clean$agecat <- factor(kilifidata_clean$agecat, levels = c("15-44", setdiff(levels(kilifidata_clean$agecat), "15-44")))

kilifidata_clean$CHKpos <- as.factor(kilifidata_clean$CHKpos)
kilifidata_clean$DENpos <- as.factor(kilifidata_clean$DENpos)
kilifidata_clean$RVFpos <- as.factor(kilifidata_clean$RVFpos)
kilifidata_clean$agecat <- factor(kilifidata_clean$agecat, ordered = FALSE)


# Kiberadata_clean (female labeled "F")
Kiberadata_clean$sex <- factor(Kiberadata_clean$sex)
Kiberadata_clean$sex <- relevel(Kiberadata_clean$sex, ref = "F")
Kiberadata_clean$agecat <- factor(Kiberadata_clean$agecat, levels = c("15-44", setdiff(levels(Kiberadata_clean$agecat), "15-44")))

Kiberadata_clean$CHKpos <- as.factor(Kiberadata_clean$CHKpos)
Kiberadata_clean$DENpos <- as.factor(Kiberadata_clean$DENpos)
Kiberadata_clean$RVFpos <- as.factor(Kiberadata_clean$RVFpos)
Kiberadata_clean$agecat <- factor(Kiberadata_clean$agecat, ordered = FALSE)


# Assuming df_all combines all datasets and has a 'site' column for clustering
df_all$sex <- factor(df_all$sex)
df_all$sex <- relevel(df_all$sex, ref = "F")  # or "Female" or "f" depending on df_all labels—adjust accordingly
df_all$agecat <- factor(df_all$agecat, levels = c("15-44", setdiff(levels(df_all$agecat), "15-44")))

df_all$CHKpos <- as.factor(df_all$CHKpos)
df_all$DENpos <- as.factor(df_all$DENpos)
df_all$RVFpos <- as.factor(df_all$RVFpos)

# Define priors
priors <- c(
  set_prior("normal(0, 1)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "b")
)

# Bayesian logistic regression with random intercept for site
library(brms)

bayesreg <- brm(
  formula = RVFpos ~ sex + agecat,
  data = kilifidata_clean,
  family = bernoulli(),
  prior = priors,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.95)
)

# Extract posterior summary and exponentiate to get ORs
post <- posterior_summary(bayesreg, probs = c(0.025, 0.975))
post_OR <- exp(post)
post_OR <- round(post_OR, 2)
colnames(post_OR) <- c("Estimate", "Est.Error", "2.5%", "97.5%")
post_OR_clean <- post_OR[grep("^b_", rownames(post_OR)), ]
rownames(post_OR_clean) <- gsub("^b_", "", rownames(post_OR_clean))
print(post_OR_clean)


Asembodata_clean$agecat   <- as.factor(as.character(Asembodata_clean$agecat))
manyattadata_clean$agecat <- as.factor(as.character(manyattadata_clean$agecat))
Kilifidata_clean$agecat   <- as.factor(as.character(Kilifidata_clean$agecat))
Kiberadata_clean$agecat   <- as.factor(as.character(Kiberadata_clean$agecat))
Nairobi_clean$agecat      <- as.factor(as.character(Nairobi_clean$agecat))


df_all <- bind_rows(
  Asembodata_clean,
  manyattadata_clean,
  Kilifidata_clean,
  Kiberadata_clean,
  Nairobi_clean
)

df_all<- df_all %>%  select(site, agecat, sex, CHKpos, DENpos, RVFpos)
df_all$site <- as.factor(df_all$site)

df_all <- df_all %>%
  mutate(sex = case_when(
    sex %in% c("M", "m", "Male", "male") ~ "Male",
    sex %in% c("F", "f", "Female", "female") ~ "Female",
    TRUE ~ NA_character_
  )) %>%
  mutate(sex = factor(sex, levels = c("Male", "Female")))
str(df_all)

df_all$sex <- relevel(df_all$sex, ref = "Female")
df_all$site <- relevel(df_all$site, ref = "NAIROBI")


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
#model of force of infections
constantmodel = FOImodel(type = 'constant', priorC1 = 0.2, priorC2 = 1, seroreversion = 1, priorRho1 = 1.38, priorRho2 = 1)

piecewisemodel = FOImodel(type = 'piecewise', priorC1 = 0.2, priorC2= 1, seroreversion = 1, priorRho1 = 1.38, priorRho2 = 1,
                          K=2, priorT1 = 30, priorT2 = 10) 



print(constantmodel)
print(piecewisemodel)

#fit the model
constantfit = fit(data = asembo_chik,  model = constantmodel, chains=1, iter=5000)
piecewisefit = fit(data = asembo_chik,  model = piecewisemodel, chains=1, iter=5000)


#visualize the best model fit
seroprevalence.fit(constantfit, YLIM = 1, age_class = 5)
seroprevalence.fit(piecewisefit, YLIM = 1, age_class = 5)


#model comparison
m1 = compute_information_criteria(constantfit)
print(m1)

m2 = compute_information_criteria(piecewisefit)
print(m2) 

#Plotting the posterior distributions
plot_posterior(piecewisefit)
plot_posterior(constantfit)
plot_posterior(outbreakfit)
#Posterior distribution of relevant model parameters
parameters_credible_intervals(constantfit)
parameters_credible_intervals(piecewisefit)

#convert annual probabilities to annual force of infection
prob_to_foi <- function(prob_percent, period_years = 1) {
  p <- prob_percent / 100
  lambda <- -log(1 - p) / period_years

  data.frame(
    Annual_Prob_Percent = prob_percent,
    FOI_per_year = lambda
  )
}
annual_probs <- c(1.166755e+00)
foi_results <- prob_to_foi(annual_probs)

print(foi_results)


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
#fit the model constant mo
constantmodel = FOImodel(type = 'constant', priorC1 = 0.01, priorC2 = 1, seroreversion = 1, priorRho1 = 0.75, priorRho2 = 1)
print(constantmodel)
piecewisemodel = FOImodel(type = 'piecewise', priorC1 = 0.2, priorC2= 1, seroreversion = 1, priorRho1 = 0.2, priorRho2 = 1,
                          K=2, priorT1 = 25, priorT2 = 0.5)


constantfit = fit(data = manyatta_chik,  model =constantmodel, chains=1, iter=5000)

piecewisefit = fit(data = manyatta_chik,  model = piecewisemodel, chains=1, iter=5000)




seroprevalence.fit(constantfit, YLIM = 1, age_class = 5)
seroprevalence.fit(piecewisefit, YLIM = 1, age_class = 5)


parameters_credible_intervals(constantfit)
parameters_credible_intervals(piecewisefit)

m1 = compute_information_criteria(constantfit)
print(m1)

m2 = compute_information_criteria(piecewisefit)
print(m2) 
#kilifi
kilifidata<- Kilifidata_clean %>%
  select(ageyrs,sex,agecat,CHKpos, DENpos, RVFpos)


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
constantmodel = FOImodel(type = 'constant', priorC1 = 0.09,  priorC2 = 1,seroreversion = 1, priorRho1 =0.0006 , priorRho2 = 1)
piecewisemodel = FOImodel(type = 'piecewise', priorC1 = 0.01, priorC2= 1, seroreversion = 1, priorRho1 = 0.06,
                                              priorRho2 = 1, K=2, priorT1 = 50, priorT2 = 40)
            
constantfit = fit(data = kilifidata_chik,  model = constantmodel, chains=1, iter=5000)
piecewisefit = fit(data = kilifidata_chik,  model = piecewisemodel,parallel_chains=1, iter=5000)

seroprevalence.fit(constantfit, YLIM = 1, age_class = 5)
seroprevalence.fit(piecewisefit,age_class = 5)
seroprevalence.fit(independentfit, YLIM = 1, age_class = 5)

m1 = compute_information_criteria(constantfit)
print(m1)

m2 = compute_information_criteria(piecewisefit)
print(m2) 



#Posterior distribution of relevant model parameters
parameters_credible_intervals(constantfit)
parameters_credible_intervals(piecewisefit)



#kibera
kiberadata<- Kiberadata_clean %>%
  select(ageyrs,sex,agecat,CHKpos, DENpos, RVFpos)
kiberadata$RVFpos <- ifelse(kiberadata$RVFpos == 1, TRUE, FALSE)
kiberadata$DENpos <- ifelse(kiberadata$DENpos == 1, TRUE, FALSE)
kiberadata$CHKpos <- ifelse(kiberadata$CHKpos == 1, TRUE, FALSE)

kiberadata$ageyrs=  as.integer(kiberadata$ageyrs)
kiberadata <- kiberadata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
kiberadata$sex<-as.character(kiberadata$sex)

kiberadata_chik <-SeroData(age_at_sampling = kiberadata$ageyrs,
                           Y= kiberadata$CHKpos,
                           sampling_year = 2022)
seroprevalence(kiberadata_chik)
seroprevalence.plot(kiberadata_chik, age_class = 5)

constantmodel = FOImodel(type = 'constant', priorC1 = 0.01, priorC2 = 1, seroreversion = 1, priorRho1 = 0.06, priorRho2 = 1)
piecewisemodel = FOImodel(type = 'piecewise', priorC1 = 0.01, priorC2 = 1, seroreversion = 1, priorRho1 = 0.06, priorRho2 = 1, K=2, priorT1 = c(17, 32), priorT2 = c(10,10))

constantfit = fit(data = kiberadata_chik,  model = constantmodel, chains=1, iter=5000)
piecewisefit = fit(data = kiberadata_chik,  model = piecewisemodel, chains=4,control= list(adapt_delta=0.95, max_treedepth =15), iter=2000)

seroprevalence.fit(constantfit, YLIM= 1, age_class = 5)
seroprevalence.fit(piecewisefit, YLIM= 1, age_class = 5)

parameters_credible_intervals(constantfit)
parameters_credible_intervals(piecewisefit)

m1 = compute_information_criteria(constantfit)
print(m1)

m2 = compute_information_criteria(piecewisefit)
print(m2) 

#Nairobi Urban

Nairobidata<- Nairobi_clean %>%
  select(ageyrs,sex,agecat,CHKpos, DENpos, RVFpos)

Nairobidata$RVFpos <- ifelse(Nairobidata$RVFpos == 1, TRUE, FALSE)
Nairobidata$DENpos <- ifelse(Nairobidata$DENpos == 1, TRUE, FALSE)
Nairobidata$CHKpos <- ifelse(Nairobidata$CHKpos == 1, TRUE, FALSE)

Nairobidata$ageyrs=  as.integer(Nairobidata$ageyrs)
Nairobidata <- Nairobidata %>%
  mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
Nairobidata$sex<-as.character(Nairobidata$sex)

Nairobidata_chik <-SeroData(age_at_sampling = Nairobidata$ageyrs,
                           Y= Nairobidata$CHKpos,
                           sampling_year = 2022)
seroprevalence(Nairobidata_chik)
seroprevalence.plot(Nairobidata_chik, age_class = 5)

constantmodel = FOImodel(type = 'constant', priorC1 = 0.01, priorC2 = 1, seroreversion = 1, priorRho1 = 0.2, priorRho2 = 1)
piecewisemodel = FOImodel(type = 'piecewise', priorC1 = 0.01, priorC2 = 1, seroreversion = 1, priorRho1 = 0.2, priorRho2 = 1,
                          K=2, priorT1 = 20, priorT2 = 0)
constantfit = fit(data = Nairobidata_chik,  model = constantmodel, chains=1, iter=5000)
piecewisefit = fit(data = Nairobidata_chik,  model = piecewisemodel, chains=1, iter=5000)

seroprevalence.fit(constantfit, YLIM=1, age_class = 5)
seroprevalence.fit(piecewisefit, YLIM=1, age_class = 5)

m1 = compute_information_criteria(constantfit)
print(m1)

m2 = compute_information_criteria(piecewisefit)
print(m2) 


parameters_credible_intervals(constantfit)
parameters_credible_intervals(piecewisefit)




#combined data
#estimating seroprevalence and visualizing aggregated data
df_all<- df_all %>% mutate(ageyrs = ifelse(ageyrs == 0, 1, ageyrs))
df_all$CHKpos <- ifelse(df_all$CHKpos == 1, TRUE, FALSE)
df_all$DENpos <- ifelse(df_all$DENpos == 1, TRUE, FALSE)
df_all$RVFpos <- ifelse(df_all$RVFpos == 1, TRUE, FALSE)
df_all$ageyrs=  as.integer(df_all$ageyrs)

df_chik <-SeroData(age_at_sampling = round(as.numeric(df_all$ageyrs)),
                             Y= as.numeric(df_all$DENpos),
                             sampling_year = as.numeric(2022))
seroprevalence(df_chik)
seroprevalence.plot(df_chik, age_class = 5)

#model of force of infections
#constantmodel= FOImodel(type="constant",   priorC1 = 0.013, priorC2 = 1, seroreversion = 1, priorRho1 = 0.2, priorRho2 = 0)
revmodel= FOImodel(type="constant", priorC1 = 0.2, priorC2 = 1, cat_lambda = TRUE, seroreversion = 1, priorRho1 = 0.3, priorRho2 = 0.58)
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
m1 = compute_information_criteria(revmodel)
print(m1)
m2 = compute_information_criteria(revmodel)
print(m2)


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
library(readxl)
library(tidyverse)
library(readstata13)

file_path <- '/Users/augustinemasinde/Desktop/PhD files/Global Health/Data/PBIDS_sites_Ase_Man_Kib_midyr_popn_2022.xlsx'

pop_asembo <- read_excel(file_path, sheet = "Asembo_mid_2022")
pop_asembo<- pop_asembo %>%
  select(Ageband,Sex,N_mid2022)
pop_asembo <- pop_asembo %>% 
  filter(if_any(everything(), ~ !is.na(.) & . != "")) %>% 
  filter(!grepl("Total", .[[1]], ignore.case = TRUE))

asembo_pop_data <- pop_asembo %>%
  mutate(age_group = case_when(
    Ageband %in% c("0–4") ~ "<5",
    Ageband %in% c("5–9", "10–14") ~ "5-14",
    Ageband %in% c("15–24", "25–34", "35–44") ~ "15-44",
    Ageband %in% c("45–54", "55–64") ~ "45-64",
    Ageband == "65+" ~ "65+",
    TRUE ~ NA_character_
  ))

asembo_pop_data <- asembo_pop_data %>%
  mutate(N_mid2022 = as.numeric(N_mid2022))
asembo_pop_data <- asembo_pop_data %>%
  group_by(Sex, age_group) %>%
  summarise(N_mid2022 = sum(N_mid2022, na.rm = TRUE), .groups = "drop")


#manyatta population data
pop_manyatta <- read_excel(file_path, sheet = "Manyatta_mid_2022")

pop_manyatta<- pop_manyatta %>%
  select(Ageband,Sex,N_mid2022)
pop_manyatta <- pop_manyatta %>% 
  filter(if_any(everything(), ~ !is.na(.) & . != "")) %>% 
  filter(!grepl("Total", .[[1]], ignore.case = TRUE))

#create age categories
manyatta_pop_data <- pop_manyatta %>%
  mutate(age_group = case_when(
    Ageband %in% c("0–4") ~ "<5",
    Ageband %in% c("5–9", "10–14") ~ "5-14",
    Ageband %in% c("15–24", "25–34", "35–44") ~ "15-44",
    Ageband %in% c("45–54", "55–64") ~ "45-64",
    Ageband == "65+" ~ "65+",
    TRUE ~ NA_character_
  ))

#convert N_mid2022 to numeric
manyatta_pop_data <- manyatta_pop_data %>%
  mutate(N_mid2022 = as.numeric(N_mid2022))
manyatta_pop_data <- manyatta_pop_data %>%
  group_by(Sex, age_group) %>%
  summarise(N_mid2022 = sum(N_mid2022, na.rm = TRUE), .groups = "drop")

#kibera population data
pop_kibera <- read_excel(file_path, sheet = "Kibera_mid_2022")
pop_kibera<- pop_kibera %>%
  select(Ageband,Sex,N_mid2022)
pop_kibera <- pop_kibera %>% 
  filter(if_any(everything(), ~ !is.na(.) & . != "")) %>% 
  filter(!grepl("Total", .[[1]], ignore.case = TRUE))


#create age categories
kibera_pop_data <- pop_kibera %>%
  mutate(age_group = case_when(
    Ageband %in% c("0–4") ~ "<5",
    Ageband %in% c("5–9", "10–14") ~ "5-14",
    Ageband %in% c("15–24", "25–34", "35–44") ~ "15-44",
    Ageband %in% c("45–54", "55–64") ~ "45-64",
    Ageband == "65+" ~ "65+",
    TRUE ~ NA_character_
  ))

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

nairobi_data <-NairobiUrbandata %>%
  mutate(agecat = case_when(
    age_y >= 0 & age_y < 5 ~ "<5",
    age_y >= 5 & age_y <15 ~ "5-14",
    age_y >= 15 & age_y <45 ~ "15-44",
    age_y >= 45 & age_y <65 ~ "45-64",
    age_y >= 65 ~ "65+",
    TRUE ~ NA_character_
  ))
nairobi_data <- nairobi_data %>%
  mutate(agecat = factor(agecat, levels = c("<5", "5-14", "15-44", "45-64", "65+"))) %>%
  arrange(agecat)

nrb_counts <- 
  nairobi_data %>%
  group_by(sex, agecat, .drop = FALSE) %>%
  summarise(
    y = sum(rvfgc_pos == "1"),
    n = n(),
    .groups = "drop"
  ) %>%
  select(agecat, sex, y, n)


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
asembo_data <- asembo_data %>%
  mutate(RVFpos = as.integer(as.character(RVFpos)))


asembo_counts <- 
  asembo_data %>%
  group_by(sex, agecat, .drop = FALSE) %>%
  summarise(y = sum(RVFpos), n = n()) %>%
  mutate(n=ifelse(is.na(y),1,n), 
         y=ifelse(is.na(y),0,y))
#arrange by agecat
library(dplyr)

asembo_counts <- asembo_counts %>%
  mutate(
    agecat = factor(
      agecat,
      levels = c("<5", "5-14", "15-44", "45-64", "65+")
    ),
    sex = factor(sex, levels = c("F", "M"))
  ) %>%
  arrange(sex, agecat)



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
  summarise(y = sum(denv_pos), n = n()) %>%
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
  summarise(
    y = sum(RVFpos == 1, na.rm = TRUE),  # count positives
    n = n(),                              # total in stratum
    .groups = "drop"
  ) %>%
  mutate(
    y = ifelse(is.na(y), 0, y),          # replace NA with 0
    n = ifelse(is.na(n), 1, n)           # replace NA with 1 if no samples
  ) %>%
  select(agecat, sex, y, n)

manyatta_counts <- manyatta_counts %>%
  mutate(
    sex = ifelse(sex == "F", "Female",
                 ifelse(sex == "M", "Male", sex))
  )



manyatta_data <- manyatta_counts %>% 
  left_join(manyatta_set, by = c("sex", "agecat"))




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

kibera_set <-kibera_set %>%
  group_by(agecat, sex) %>%
  summarise(
    pop = sum(pop, na.rm = TRUE),
    .groups = "drop"
  )

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
  int N_se;                // denominator for sensitivity
  int N_sp;                // denominator for specificity
  int x;                   // numerator for sensitivity
  int z;                   // numerator for specificity
  int y[5, 2];             // number of seropositives
  int n[5, 2];             // number of samples
  real pw[5, 2];           // proportion of population in each subgroup
  real tot_pw_age[5];      // proportion of population in each age group
  real tot_pw_sex[2];      // proportion female and male
}

parameters {
  real<lower=0,upper=1> se; 
  real<lower=0,upper=1> sp; 
  vector[2] bsex_raw;      // raw sex effects
  vector[5] bage_raw;      // raw age effects
  real<lower=0> sd_age;    // age effect scale
}

transformed parameters {
  vector[2] bsex;           // centered sex effects
  vector[5] bage;           // scaled age effects
  matrix[5, 2] p;           // true prevalence
  matrix[5, 2] p_obs;       // observed prevalence

  // non-centered age effects
  bage = bage_raw * sd_age;

  // sum-to-zero constraint for sex effects
  bsex = bsex_raw - mean(bsex_raw);

  // probability model
  for (a in 1:5){
    for (s in 1:2){
      p[a, s] = inv_logit(bage[a] + bsex[s]);
      p_obs[a, s] = se * p[a, s] + (1 - sp) * (1 - p[a, s]);
    }
  }
}

model {
  // Priors
  se ~ beta(8, 2);
  sp ~ beta(8, 2);
  bsex_raw ~ normal(0, 0.5);
  bage_raw ~ normal(0, 1);
  sd_age ~ normal(0, 1) T[0,];

  // Likelihood
  for (a in 1:5){
    for (s in 1:2){
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

  // Population-weighted prevalence
  for (a in 1:5){
    for (s in 1:2){
      p_site += p[a, s] * pw[a, s];
    }
  }

  for (a in 1:5){
    for (s in 1:2){
      p_age[a] += p[a, s] * pw[a, s] / tot_pw_age[a];
    }
  }

  for (s in 1:2){
    for (a in 1:5){
      p_sex[s] += p[a, s] * pw[a, s] / tot_pw_sex[s];
    }
  }
}

"

mrp_model <- stan_model(model_code = mrp_adjusted)

fit <- sampling(
  object = mrp_model,     
  data = list(
    y = y,
    n = n,
    pw = pw,
    tot_pw_age = tot_pw_age,
    tot_pw_sex = tot_pw_sex,
    x = 52,
    z = 65,
    N_se = 57,
    N_sp =76
    
  ),
  chains = 3,
  warmup = 2000,
  iter = 10000,
  cores = 4,
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


#Site forest plot
# Priors must already exist as `priors`
bayesreg_chk <- brm(
  formula = CHKpos ~ site,
  data = df_all,
  family = bernoulli(),
  prior = priors,
  cores = 4
)

bayesreg_den <- brm(
  formula = DENpos ~ site,
  data = df_all,
  family = bernoulli(),
  prior = priors,
  cores = 4
)

bayesreg_rvf <- brm(
  formula = RVFpos ~ site,
  data = df_all,
  family = bernoulli(),
  prior = priors,
  cores = 4
)


make_forest_df <- function(brms_model) {
  
  sm <- summary(brms_model)$fixed   # Extract fixed-effects table
  
  OR <- exp(sm[, "Estimate"])
  lower <- exp(sm[, "l-95% CI"])
  upper <- exp(sm[, "u-95% CI"])
  
  # Probability of direction (Bayesian significance)
  pd_vals <- p_direction(brms_model)
  pd_vals <- pd_vals$pd[match(rownames(sm), pd_vals$Parameter)]
  
  forest_data <- data.frame(
    predictor = rownames(sm),
    OR = OR,
    lower = lower,
    upper = upper,
    pd = pd_vals,
    significant = ifelse(pd_vals > 0.95, "yes", "no")
  )
  
  # Remove intercept
  forest_data <- forest_data %>% filter(predictor != "Intercept")
  
  # Grouping (only site effects)
  forest_data$group <- "Site"
  
  # Keep predictor order as is
  forest_data$predictor <- factor(forest_data$predictor,
                                  levels = forest_data$predictor)
  
  return(forest_data)
}


plot_forest <- function(forest_data, title_text) {
  
  max_val <- ceiling(max(forest_data$upper, na.rm = TRUE))
  min_val <- floor(min(forest_data$lower, na.rm = TRUE))
  x_breaks <- pretty(c(min_val, max_val), n = 6)
  
  ggplot(forest_data, aes(x = OR, y = predictor, color = significant)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = x_breaks) +
    xlab("Odds Ratio") +
    ylab("Site") +
    ggtitle(title_text) +
    facet_grid(group ~ ., scales = "free_y", space = "free") +
    theme_minimal()
}


forest_chk <- make_forest_df(bayesreg_chk)
forest_den <- make_forest_df(bayesreg_den)
forest_rvf <- make_forest_df(bayesreg_rvf)


plot_forest(forest_chk, "Forest Plot of CHKpos by Site")
plot_forest(forest_den, "Forest Plot of DENpos by Site")
plot_forest(forest_rvf, "Forest Plot of RVFpos by Site")
# End of Asembo1.R








library(dplyr)
library(ggplot2)
library(forcats)
library(posterior)
library(cmdstanr)
library(tidyr)

# --- 1. Prepare data (assuming Kilifidata_clean is loaded) ---
df_age <- Kilifidata_clean %>%
  mutate(age = as.integer(ageyrs)) %>%
  group_by(age) %>%
  summarise(
    positive = sum(CHKpos == 1, na.rm = TRUE),
    total = sum(!is.na(CHKpos))
  ) %>%
  arrange(age) %>%
  ungroup()

# --- 2. Exposure matrix ---
ages_vec <- df_age$age
age_max <- max(ages_vec)

n_obs <- nrow(df_age)
exposure_matrix <- matrix(0L, nrow = n_obs, ncol = age_max)
for (i in seq_len(n_obs)) {
  a <- ages_vec[i]
  if (a > 0) exposure_matrix[i, 1:a] <- 1L
}

# --- 3. Stan model code ---
stan_file <- "time_foi_seroreversion.stan"

stan_code <- '
data {
  int<lower=1> n_obs;
  array[n_obs] int<lower=0> n_pos;
  array[n_obs] int<lower=1> n_total;
  int<lower=1> age_max;
  array[n_obs] int<lower=0> ages;
  array[n_obs, age_max] int<lower=0, upper=1> exposure_matrix;
}
parameters {
  vector<lower=-12, upper=2>[age_max] log_foi;
  real<lower=0, upper=1> rho;
}
transformed parameters {
  vector<lower=0>[age_max] foi = exp(log_foi);
  vector[n_obs] prevalence;
  for (i in 1:n_obs) {
    real cum_hazard = 0;
    int a = ages[i];
    for (j in 1:a) {
      if (exposure_matrix[i, j] == 1) {
        cum_hazard += foi[j] * pow(1 - rho, a - j);
      }
    }
    prevalence[i] = 1 - exp(-cum_hazard);
  }
}
model {
  log_foi ~ normal(log(0.01), 1.5);
  rho ~ beta(2, 20);
  for (i in 1:n_obs) {
    real p = fmin(fmax(prevalence[i], 1e-6), 1 - 1e-6);
    n_pos[i] ~ binomial(n_total[i], p);
  }
}
'

writeLines(stan_code, stan_file)

# --- 4. Prepare Stan data ---
stan_data <- list(
  n_obs = n_obs,
  n_pos = as.integer(df_age$positive),
  n_total = as.integer(df_age$total),
  age_max = as.integer(age_max),
  ages = as.integer(ages_vec),
  exposure_matrix = exposure_matrix
)

# --- 5. Compile model ---
mod <- cmdstan_model(stan_file)

# --- 6. Initial values ---
init_fun <- function() {
  list(
    log_foi = rep(log(0.01), age_max),
    rho = 0.05
  )
}

# --- 7. Sampling ---
fit <- mod$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  refresh = 500,
  init = init_fun
)

# --- 8. Extract posterior draws ---
draws <- fit$draws(format = "df") # data.frame format for easier manipulation
draws_posterior <- fit$draws(format = "draws_df") # posterior::draws_df format

# Extract names of log_foi variables
log_foi_names <- grep("^log_foi\\[", names(draws_posterior), value = TRUE)

# --- 9. Compute posterior mean prevalence per age ---
n_draws <- nrow(draws_posterior)
p_hat_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)

for (d in seq_len(n_draws)) {
  log_foi_vec <- as.numeric(draws_posterior[d, log_foi_names])
  foi_vec <- exp(log_foi_vec)
  rho_d <- draws_posterior$rho[d]
  for (i in seq_len(n_obs)) {
    a <- ages_vec[i]
    cum <- 0
    if (a > 0) {
      for (j in 1:a) {
        if (exposure_matrix[i, j] == 1L) {
          cum <- cum + foi_vec[j] * (1 - rho_d)^(a - j)
        }
      }
    }
    p_hat_mat[d, i] <- 1 - exp(-cum)
  }
}

p_hat_mean <- apply(p_hat_mat, 2, mean)
p_hat_low <- apply(p_hat_mat, 2, quantile, probs = 0.025)
p_hat_high <- apply(p_hat_mat, 2, quantile, probs = 0.975)

# --- 10. Prepare data frame for plotting ---
plot_df <- df_age %>%
  mutate(
    observed_prev = 100 * positive / total,
    obs_lower = 100 * qbeta(0.025, 1 + positive, 1 + total - positive),
    obs_upper = 100 * qbeta(0.975, 1 + positive, 1 + total - positive),
    model_prev = 100 * p_hat_mean,
    model_low = 100 * p_hat_low,
    model_high = 100 * p_hat_high
  )

# --- 11. Create 5-year age bins for plotting ---
breaks <- seq(0, max(plot_df$age) + 5, by = 5)
labels <- paste(breaks[-length(breaks)], breaks[-1] - 1, sep = "-")

plot_df_binned <- plot_df %>%
  mutate(age_bin = cut(age,
                       breaks = breaks,
                       right = FALSE,
                       labels = labels)) %>%
  group_by(age_bin) %>%
  summarise(
    observed_prev = mean(observed_prev),
    obs_lower = mean(obs_lower),
    obs_upper = mean(obs_upper),
    model_prev = mean(model_prev),
    model_low = mean(model_low),
    model_high = mean(model_high),
    .groups = "drop"
  ) %>%
  mutate(age_bin = forcats::fct_inorder(age_bin))

# --- 12. Plot binned observed vs predicted seroprevalence ---
p1 <- ggplot(plot_df_binned, aes(x = age_bin)) +
  geom_line(aes(y = model_prev, group = 1), color = "#2C7FB8", size = 0.9) +
  geom_ribbon(aes(ymin = model_low, ymax = model_high, group = 1), fill = "#2C7FB8", alpha = 0.2) +
  geom_point(aes(y = observed_prev), size = 3, color = "black") +
  geom_errorbar(aes(ymin = obs_lower, ymax = obs_upper), width = 0.4) +
  labs(x = "Age group (years)", y = "Seroprevalence (%)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

# --- 13. Plot FOI reconstruction with uncertainty ---
foi_draws <- draws_posterior %>%
  select(all_of(log_foi_names)) %>%
  mutate_all(exp)

foi_mean <- apply(foi_draws, 2, mean)
foi_low <- apply(foi_draws, 2, quantile, probs = 0.025)
foi_high <- apply(foi_draws, 2, quantile, probs = 0.975)

foi_df <- tibble(
  year_index = 1:age_max,
  foi_mean = foi_mean,
  foi_low = foi_low,
  foi_high = foi_high
)

p2 <- ggplot(foi_df, aes(x = year_index)) +
  geom_line(aes(y = foi_mean), color = "#D95F02", size = 1) +
  geom_ribbon(aes(ymin = foi_low, ymax = foi_high), fill = "#D95F02", alpha = 0.3) +
  labs(x = "Years ago (1 = most recent)", y = "Force of Infection (FOI)") +
  theme_minimal(base_size = 14)

print(p2)

# Optionally save plots
ggsave("seroprev_fit_binned.png", p1, width = 7, height = 5)
ggsave("foi_reconstruction.png", p2, width = 7, height = 5)


# Extract log_foi variable names from posterior draws
log_foi_names <- grep("^log_foi\\[", names(draws_posterior), value = TRUE)

# Extract and exponentiate FOI draws
foi_draws <- draws_posterior %>%
  select(all_of(log_foi_names)) %>%
  mutate_all(exp)

# Calculate mean FOI over all draws and ages (flatten all values)
mean_foi <- mean(unlist(foi_draws))

cat("Overall mean FOI:", mean_foi, "\n")
# Calculate mean FOI for the most recent 5 years
recent_5yr_names <- paste0("log_foi[", 1:5, "]")
recent_5yr_foi <- foi_draws %>%
  select(all_of(recent_5yr_names)) %>%
  mutate(mean_recent_5yr = rowMeans(.))
mean_recent_5yr <- mean(recent_5yr_foi$mean_recent_5yr)
cat("Mean FOI for most recent 5 years:", mean_recent_5yr, "\n")
# Calculate mean FOI for the most recent 10 years
recent_10yr_names <- paste0("log_foi[", 1:10, "]")
recent_10yr_foi <- foi_draws %>%
  select(all_of(recent_10yr_names)) %>%
  mutate(mean_recent_10yr = rowMeans(.))
mean_recent_10yr <- mean(recent_10yr_foi$mean_recent_10yr)
cat("Mean FOI for most recent 10 years:", mean_recent_10yr, "\n")






#mixture models and FOI
library(dplyr)
library(mgcv)
library(mixR)

df <- Asembodata %>% dplyr::select(ageyrs, CHKV_AIU, DENV2_AIU)
df<- manyattadata %>% dplyr::select(ageyrs, CHKV_AIU, DENV2_AIU)
df<- Kilifidata %>% dplyr::select(age_y, chkve1_au,denv2ns1_au) %>% mutate(ageyrs= age_y, CHKV_AIU= chkve1_au, DENV2_AIU= denv2ns1_au) %>% dplyr::select(ageyrs, CHKV_AIU, DENV2_AIU)
df<- Kiberadata %>% dplyr::select(ageyrs, chkve1_au, denv2ns1_au) %>% mutate(ageyrs= ageyrs, CHKV_AIU= chkve1_au, DENV2_AIU= denv2ns1_au) %>% dplyr::select(ageyrs, CHKV_AIU, DENV2_AIU)

df <- df %>% filter(!is.na(CHKV_AIU), CHKV_AIU >= 0)
a <- df$ageyrs
z <- log(df$CHKV_AIU + 1)
mix_model <- mixfit(z, family = "gamma", ncomp = 2)
print(mix_model)
plot(mix_model)

ord <- order(mix_model$mu)
muS  <- mix_model$mu[ord[1]]
muI  <- mix_model$mu[ord[2]]


gam_fit <- gam(z ~ s(a, bs = "ps", k = 10), data = df, method = "REML")
gam.check(gam_fit)
plot(gam_fit, shade = TRUE)

age_grid <- seq(min(a), max(a), length.out = 200)
mu_a <- predict(gam_fit, newdata = data.frame(a = age_grid))

pi_a <- (mu_a - muS) / (muI - muS)
pi_a <- pmin(pmax(pi_a, 0), 1)

seroprev_df <- data.frame(age = age_grid, pi_a = pi_a)

dmu_da <- diff(mu_a) / diff(age_grid)
age_foi <- age_grid[-1]

FOI_age <- dmu_da / (muI - mu_a[-length(mu_a)])
FOI_age[FOI_age < 0] <- NA

FOI_df <- data.frame(age = age_foi, FOI = FOI_age)

set.seed(123)
B <- 1000

FOI_boot <- replicate(B, {
  idx <- sample(seq_along(z), replace = TRUE)
  boot_df <- data.frame(z = z[idx], a = a[idx])
  fit_b <- gam(z ~ s(a, bs = "cr", k = k_val), data = boot_df, method = "REML")
  mu_b <- predict(fit_b, newdata = data.frame(a = age_grid))
  dmu_b <- diff(mu_b) / diff(age_grid)
  foi_b <- dmu_b / (muI - mu_b[-length(mu_b)])
  foi_b[foi_b < 0] <- NA
  foi_b
})

FOI_median <- apply(FOI_boot, 1, median, na.rm = TRUE)
FOI_low    <- apply(FOI_boot, 1, quantile, 0.025, na.rm = TRUE)
FOI_upp    <- apply(FOI_boot, 1, quantile, 0.975, na.rm = TRUE)

FOI_ci_df <- data.frame(age = age_foi, FOI = FOI_median, FOI_low = FOI_low, FOI_upp = FOI_upp)

FOI_total <- data.frame(
  FOI = mean(FOI_ci_df$FOI, na.rm = TRUE),
  lower_ci = mean(FOI_ci_df$FOI_low, na.rm = TRUE),
  upper_ci = mean(FOI_ci_df$FOI_upp, na.rm = TRUE)
)
print(FOI_total)

pi_individual <- approx(x = age_grid, y = pi_a, xout = a, rule = 2)$y
total_seroprev <- mean(pi_individual, na.rm = TRUE)
cat("Overall seroprevalence:", total_seroprev, "\n")

ggplot() +
  geom_line(data = seroprev_df, aes(x = age, y = pi_a), color = "blue", size = 1) +
  labs(x = "Age (years)", y = "Seroprevalence", title = "Estimated Seroprevalence vs Age") +
  theme_minimal()








library(dplyr)
library(scam)
library(mixR)
library(mgcv)



df <- Asembodata %>% dplyr::select(ageyrs, CHKV_AIU, DENV2_AIU)
df<- manyattadata %>% dplyr::select(ageyrs, CHKV_AIU, DENV2_AIU)
df<- Kilifidata %>% dplyr::select(age_y, chkve1_au,denv2ns1_au) %>% mutate(ageyrs= age_y, CHKV_AIU= chkve1_au, DENV2_AIU= denv2ns1_au) %>% dplyr::select(ageyrs, CHKV_AIU, DENV2_AIU)
df<- Kiberadata %>% dplyr::select(ageyrs, chkve1_au, denv2ns1_au) %>% mutate(ageyrs= ageyrs, CHKV_AIU= chkve1_au, DENV2_AIU= denv2ns1_au) %>% dplyr::select(ageyrs, CHKV_AIU, DENV2_AIU)
df <- df %>% filter(!is.na(CHKV_AIU), CHKV_AIU >= 0)
a <- df$ageyrs
z <- log(df$CHKV_AIU + 1)


mix_model <- mixfit(z, family = "gamma", ncomp = 2)
print(mix_model)

ord <- order(mix_model$mu)
muS <- mix_model$mu[ord[1]]
muI <- mix_model$mu[ord[2]]

iso_fit <- scam(
  z ~ s(a, bs = "mpi", k = 30),
  data = df,
)

summary(iso_fit)
plot(iso_fit, shade = TRUE)

age_grid <- seq(min(a), max(a), length.out = 200)

mu_a <- predict(iso_fit, newdata = data.frame(a = age_grid))

pi_a <- (mu_a - muS) / (muI - muS)
pi_a <- pmin(pmax(pi_a, 0), 1)

dmu_da <- diff(mu_a) / diff(age_grid)
age_foi <- age_grid[-1]

FOI_age <- dmu_da / (muI - mu_a[-length(mu_a)])
FOI_age[FOI_age < 0] <- NA

B <- 1000

boot <- replicate(B, {
  
  idx <- sample(seq_along(z), replace = TRUE)
  
  boot_df <- data.frame(
    a = a[idx],
    z = z[idx]
  )
  
  fit_b <- scam(
    z ~ s(a, bs = "mpi", k = 30),
    data = boot_df
  )
  plot(iso_fit, shade = TRUE)
  
  mu_b <- predict(fit_b, newdata = data.frame(a = age_grid))
  
  pi_b <- (mu_b - muS) / (muI - muS)
  pi_b <- pmin(pmax(pi_b, 0), 1)
  
  dmu_b <- diff(mu_b) / diff(age_grid)
  
  foi_b <- dmu_b / (muI - mu_b[-length(mu_b)])
  foi_b[foi_b < 0] <- NA
  
  c(pi_b, foi_b)
})

n_age <- length(age_grid)
n_foi <- length(age_foi)

pi_boot <- boot[1:n_age, ]
foi_boot <- boot[(n_age + 1):(n_age + n_foi), ]

pi_med <- apply(pi_boot, 1, median, na.rm = TRUE)
pi_low <- apply(pi_boot, 1, quantile, 0.025, na.rm = TRUE)
pi_up  <- apply(pi_boot, 1, quantile, 0.975, na.rm = TRUE)

foi_med <- apply(foi_boot, 1, median, na.rm = TRUE)
foi_low <- apply(foi_boot, 1, quantile, 0.025, na.rm = TRUE)
foi_up  <- apply(foi_boot, 1, quantile, 0.975, na.rm = TRUE)

seroprev_df <- data.frame(
  age = age_grid,
  pi = pi_med,
  pi_low = pi_low,
  pi_up = pi_up
)

FOI_df <- data.frame(
  age = age_foi,
  FOI = foi_med,
  FOI_low = foi_low,
  FOI_up = foi_up
)

FOI_total <- data.frame(
  FOI = mean(FOI_df$FOI, na.rm = TRUE),
  lower_ci = mean(FOI_df$FOI_low, na.rm = TRUE),
  upper_ci = mean(FOI_df$FOI_up, na.rm = TRUE)
)

Seroprev_total <- data.frame(
  seroprev = mean(seroprev_df$pi, na.rm = TRUE),
  lower_ci = mean(seroprev_df$pi_low, na.rm = TRUE),
  upper_ci = mean(seroprev_df$pi_up, na.rm = TRUE)
)

FOI_total
Seroprev_total





library(mgcv)
library(scam)
library(mixdist)

a <- df$ageyrs
z <- log(df$CHKV_AIU + 1)

dat <- data.frame(a = a, z = z)

log_mfi <- log(df$CHKV_AIU[df$CHKV_AIU > 0] + 1)

mix_model <- mixfit(log_mfi, family = "gamma", ncomp = 2)

ord <- order(mix_model$mu)
muS <- mix_model$mu[ord[1]]
muI <- mix_model$mu[ord[2]]

iso_fit <- scam(
  z ~ s(a, bs = "mpi", k = 35),
  data = dat
)
summary(iso_fit)
plot(iso_fit, shade = TRUE)

age_grid <- seq(min(a), max(a), length.out = 200)

mu_hat <- predict(
  iso_fit,
  newdata = data.frame(a = age_grid)
)

pi_a <- (mu_hat - muS) / (muI - muS)
pi_a <- pmin(pmax(pi_a, 0), 1)

seroprev_df <- data.frame(
  age = age_grid,
  pi = pi_a
)

dmu_da <- diff(mu_hat) / diff(age_grid)
age_mid <- age_grid[-1]

FOI <- dmu_da / (muI - mu_hat[-length(mu_hat)])
FOI[FOI < 0] <- NA

FOI_df <- data.frame(
  age = age_mid,
  FOI = FOI
)

pi_individual <- approx(
  x = age_grid,
  y = pi_a,
  xout = a,
  rule = 2
)$y

total_seroprev <- mean(pi_individual, na.rm = TRUE)

plot(seroprev_df$age, seroprev_df$pi, type = "l")
plot(FOI_df$age, FOI_df$FOI, type = "l")

total_seroprev
total_foi <- mean(FOI_df$FOI, na.rm = TRUE)
total_foi






library(dplyr)
library(scam)
library(mixR)

set.seed(123)

df <- Asembodata %>% dplyr::select(ageyrs, CHKV_AIU)

a <- df$ageyrs
z <- log(df$CHKV_AIU + 1)

log_mfi <- log(df$CHKV_AIU[df$CHKV_AIU > 0] + 1)

mix_model <- mixfit(log_mfi, family = "gamma", ncomp = 2)

ord <- order(mix_model$mu)
muS <- mix_model$mu[ord[1]]
muI <- mix_model$mu[ord[2]]

iso_fit <- scam(
  z ~ s(a, bs = "mpi", k = 30),
  data = df
)

age_grid <- seq(min(a), max(a), length.out = 200)

mu_a <- predict(iso_fit, newdata = data.frame(a = age_grid))

pi_a <- (mu_a - muS) / (muI - muS)
pi_a <- pmin(pmax(pi_a, 0), 1)

dmu_da <- diff(mu_a) / diff(age_grid)
age_foi <- age_grid[-1]

FOI_age <- dmu_da / (muI - mu_a[-length(mu_a)])
FOI_age[FOI_age < 0] <- NA

B <- 1000

boot <- replicate(B, {
  
  idx <- sample(seq_along(z), replace = TRUE)
  
  boot_df <- data.frame(
    a = a[idx],
    z = z[idx]
  )
  
  fit_b <- scam(
    z ~ s(a, bs = "mpi", k = 30),
    data = boot_df
  )
  
  mu_b <- predict(fit_b, newdata = data.frame(a = age_grid))
  
  pi_b <- (mu_b - muS) / (muI - muS)
  pi_b <- pmin(pmax(pi_b, 0), 1)
  
  dmu_b <- diff(mu_b) / diff(age_grid)
  foi_b <- dmu_b / (muI - mu_b[-length(mu_b)])
  foi_b[foi_b < 0] <- NA
  
  c(pi_b, foi_b)
})

pi_boot <- boot[1:200, ]
foi_boot <- boot[201:399, ]

pi_med <- apply(pi_boot, 1, median, na.rm = TRUE)
pi_low <- apply(pi_boot, 1, quantile, 0.025, na.rm = TRUE)
pi_up  <- apply(pi_boot, 1, quantile, 0.975, na.rm = TRUE)

foi_med <- apply(foi_boot, 1, median, na.rm = TRUE)
foi_low <- apply(foi_boot, 1, quantile, 0.025, na.rm = TRUE)
foi_up  <- apply(foi_boot, 1, quantile, 0.975, na.rm = TRUE)

seroprev_df <- data.frame(
  age = age_grid,
  pi = pi_med,
  pi_low = pi_low,
  pi_up = pi_up
)

FOI_df <- data.frame(
  age = age_foi,
  FOI = foi_med,
  FOI_low = foi_low,
  FOI_up = foi_up
)

FOI_total <- data.frame(
  FOI = mean(FOI_df$FOI, na.rm = TRUE),
  lower_ci = mean(FOI_df$FOI_low, na.rm = TRUE),
  upper_ci = mean(FOI_df$FOI_up, na.rm = TRUE)
)

Seroprev_total <- data.frame(
  seroprev = mean(seroprev_df$pi, na.rm = TRUE),
  lower_ci = mean(seroprev_df$pi_low, na.rm = TRUE),
  upper_ci = mean(seroprev_df$pi_up, na.rm = TRUE)
)

print(FOI_total)
print(Seroprev_total)

