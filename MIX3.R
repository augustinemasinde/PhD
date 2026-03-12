library(dplyr)
library(rstan)
library(ggplot2)
library(bayesplot)
library(loo)
library(patchwork)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# -----------------------------
# CLEAN + COMBINE DATA
# -----------------------------

clean_data <- function(df, value_col, age_col, site_name){
  df %>%
    filter(
      !!rlang::sym(value_col) > 0,
      !is.na(!!rlang::sym(value_col)),
      !is.na(!!rlang::sym(age_col))
    ) %>%
    mutate(
      log_antibody = log(!!rlang::sym(value_col) + 1),
      age = pmin(floor(!!rlang::sym(age_col)), 80),
      site = site_name
    ) %>%
    select(log_antibody, age, site)
}

serodata_all <- bind_rows(
  clean_data(Asembodata, "CHKV_AIU", "ageyrs", "Asembo"),
  clean_data(manyattadata, "CHKV_AIU", "ageyrs", "Manyatta"),
  clean_data(Kilifidata, "chkve1_au", "age_y", "Kilifi"),
  clean_data(Kiberadata, "chkve1_au", "ageyrs", "Kibera"),
  clean_data(NairobiUrbandata, "chkve1_au", "age_y", "Nairobi")
) %>%
  mutate(
    site_id = as.integer(factor(site)),
    age_group = age + 1
  )

A <- length(unique(serodata_all$age_group))
S <- length(unique(serodata_all$site_id))

stan_data <- list(
  N = nrow(serodata_all),
  S = S,
  A = A,
  y = serodata_all$log_antibody,
  site = serodata_all$site_id,
  age_group = serodata_all$age_group
)

# -----------------------------
# STAN MODEL
# -----------------------------

stan_code_D <- "
data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> A;
  vector[N] y;
  int<lower=1,upper=S> site[N];
  int<lower=1,upper=A> age_group[N];
}

parameters {
  // ---- Seronegative ----
  real mu0_global;
  vector[S] mu0_site_raw;
  real<lower=0> sigma_mu0;
  vector<lower=0>[S] sigma0;

  // ---- Seropositive ----
  matrix[A,S] mu1_raw;
  real<lower=0> sigma1;

  // ---- Seroprevalence ----
  matrix[A,S] alpha;
  matrix[A,S] beta;
}

transformed parameters {
  vector[S] mu0_site;
  matrix[A,S] mu1;
  matrix[A,S] pi;

  for (s in 1:S)
    mu0_site[s] = mu0_global + mu0_site_raw[s] * sigma_mu0;

  for (a in 1:A)
    for (s in 1:S) {
      mu1[a,s] = mu0_site[s] + exp(mu1_raw[a,s]);
      pi[a,s]  = inv_logit(alpha[a,s] + beta[a,s]);
    }
}

model {
  // ---- Anchors ----
  mu0_global ~ normal(1.45, 0.15);
  mu0_site_raw ~ normal(0,1);
  sigma_mu0 ~ normal(0,0.2);

  sigma0 ~ normal(0.5,0.2);

  to_vector(mu1_raw) ~ normal(log(4),0.25);
  sigma1 ~ normal(0.8,0.25);

  to_vector(alpha) ~ normal(-2,1.5);
  to_vector(beta)  ~ normal(0,1);

  // ---- Likelihood ----
  for (n in 1:N) {
    int a = age_group[n];
    int s = site[n];

    target += log_sum_exp(
      log1m(pi[a,s]) +
        normal_lpdf(y[n] | mu0_site[s], sigma0[s]),
      log(pi[a,s]) +
        normal_lpdf(y[n] | mu1[a,s], sigma1)
    );
  }
}

generated quantities {
  vector[N] p_seropos;

  for (n in 1:N) {
    int a = age_group[n];
    int s = site[n];

    real lp_neg = log1m(pi[a,s]) +
                  normal_lpdf(y[n] | mu0_site[s], sigma0[s]);
    real lp_pos = log(pi[a,s]) +
                  normal_lpdf(y[n] | mu1[a,s], sigma1);

    p_seropos[n] = exp(lp_pos - log_sum_exp(lp_neg, lp_pos));
  }
}
"

fit_D <- stan(
  model_code = stan_code_D,
  data = stan_data,
  iter = 3000,
  warmup = 1500,
  chains = 4,
  seed = 123,
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

# -----------------------------
# POSTERIOR EXTRACTION
# -----------------------------

posterior_D <- rstan::extract(fit_D)

mu0_med    <- apply(posterior_D$mu0_site, 2, median)
sigma0_med <- apply(posterior_D$sigma0, 2, median)

mu1_med    <- apply(posterior_D$mu1, c(1,2), median)
sigma1_med <- median(posterior_D$sigma1)

pi_med     <- apply(posterior_D$pi, c(1,2), median)

# -----------------------------
# MIXTURE PLOTS
# -----------------------------

plot_site_mixture_D <- function(site_name, site_index) {
  
  df_site <- serodata_all %>% filter(site == site_name)
  
  age_weights <- table(df_site$age_group)
  age_weights <- age_weights / sum(age_weights)
  
  x_grid <- seq(min(df_site$log_antibody),
                max(df_site$log_antibody),
                length.out = 500)
  
  dens_neg <- numeric(length(x_grid))
  dens_pos <- numeric(length(x_grid))
  
  for (a in as.integer(names(age_weights))) {
    
    w <- age_weights[as.character(a)]
    
    dens_neg <- dens_neg +
      w * (1 - pi_med[a, site_index]) *
      dnorm(x_grid, mu0_med[site_index], sigma0_med[site_index])
    
    dens_pos <- dens_pos +
      w * pi_med[a, site_index] *
      dnorm(x_grid, mu1_med[a, site_index], sigma1_med)
  }
  
  dens_df <- data.frame(x=x_grid, neg=dens_neg, pos=dens_pos)
  
  ggplot(df_site, aes(x = log_antibody)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, fill="grey85", color="grey40") +
    geom_line(data=dens_df, aes(x,neg),
              color="#D55E00", linetype="dashed", linewidth=1) +
    geom_line(data=dens_df, aes(x,pos),
              color="#0072B2", linewidth=1) +
    labs(title=paste("Model D –", site_name),
         x="Log antibody titre", y="Density") +
    theme_classic(base_size = 13)
}

site_names <- levels(factor(serodata_all$site))

plots_D <- lapply(seq_along(site_names),
                  function(i) plot_site_mixture_D(site_names[i], i))

final_plot_D <- wrap_plots(plots_D, ncol = 2)

ggsave("Figure_CHIKV_ModelD_SiteMixtures.png",
       final_plot_D, width=10, height=12, dpi=300)
