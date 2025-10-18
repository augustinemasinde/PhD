#Chikungunya script to perform population-weighted seroprevalence estimation
if (!requireNamespace("cmdstanr", quietly = TRUE) &&
    !requireNamespace("rstan", quietly = TRUE)) {
  stop("Please install either 'cmdstanr' (recommended) or 'rstan' before running this script.")
}
library(dplyr)
library(tidyr)
library(ggplot2)


df_all <- df_all %>%
  mutate(
    CHKpos = as.numeric(as.character(CHKpos)),
    DENpos = as.numeric(as.character(DENpos)),
    RVFpos = as.numeric(as.character(RVFpos))
  )

# ---------------------------
# 1. Aggregate df_all into strata counts
# ---------------------------
strata <- df_all %>%
  group_by(site, sex, age_cat) %>%
  summarise(
    Ntest = n(),
    Npos  = sum(RVFpos, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(site, sex, age_cat)


# ---------------------------
# Here we create a pop_size for the 5 sites and then normalize to overall weights
set.seed(123)
pop_weights <- expand.grid(
  site = unique(strata$site),
  sex  = unique(strata$sex),
  age_cat = unique(strata$age_cat),
  stringsAsFactors = FALSE
) %>%
  arrange(site, sex, age_cat) %>%
  mutate(
    #pop sizes
    pop_size = case_when(
      age_cat == "0-11" ~ 146357,
      age_cat == "12-17" ~ 46852,
      age_cat == "18-29" ~ 111110,
      age_cat == "30-49" ~ 152159,
      age_cat == "50-64" ~ 52654,
      age_cat == "65+" ~ 17407,
      TRUE ~ 1000
    ),
    pop_size = pop_size * sample(80:120, n(), replace = TRUE)  # add small variation per cell
  )

# Merge weights into strata (left join) 
strata <- strata %>%
  left_join(pop_weights, by = c("site", "sex", "age_cat"))

# normalize to overall population proportion (sum to 1)
strata <- strata %>%
  mutate(pop_weight = pop_size / sum(pop_size))

# Quick checks
if (any(is.na(strata$pop_weight))) stop("Missing pop_weight for some strata — fill in population weights for all strata.")
sum(strata$pop_weight)  # should be 1

# ---------------------------
# 3. Prepare index variables for Stan
# ---------------------------
# Create consistent level ordering (Stan expects integer indices)
site_levels <- sort(unique(strata$site))
age_levels  <- sort(unique(strata$age_cat))
sex_levels  <- sort(unique(strata$sex))

strata <- strata %>%
  mutate(
    site_id = as.integer(factor(site, levels = site_levels)),
    age_id  = as.integer(factor(age_cat, levels = age_levels)),
    sex_id  = as.integer(factor(sex, levels = sex_levels))
  )

# Build stan list
stan_data <- list(
  I     = nrow(strata),
  Npos  = as.integer(strata$Npos),
  Ntest = as.integer(strata$Ntest),
  w     = as.numeric(strata$pop_weight),
  S     = length(site_levels),
  A     = length(age_levels),
  G     = length(sex_levels),
  site  = as.integer(strata$site_id),
  age   = as.integer(strata$age_id),
  sex   = as.integer(strata$sex_id)
)

# ---------------------------
# 4. Stan model code (hierarchical logistic with imperfect test)
# ---------------------------
stan_code <- "
data {
  int<lower=1> I;
  int<lower=0> Npos[I];
  int<lower=0> Ntest[I];
  vector[I] w;
  int<lower=1> S;
  int<lower=1> A;
  int<lower=1> G;
  int<lower=1,upper=S> site[I];
  int<lower=1,upper=A> age[I];
  int<lower=1,upper=G> sex[I];
}
parameters {
  real mu;
  vector[S] alpha_site;
  vector[A] alpha_age;
  vector[G] alpha_sex;
  real<lower=0,upper=1> sens;
  real<lower=0,upper=1> spec;
}
transformed parameters {
  vector[I] p_true;
  vector[I] p_obs;
  for (i in 1:I) {
    real lin = mu + alpha_site[site[i]] + alpha_age[age[i]] + alpha_sex[sex[i]];
    p_true[i] = inv_logit(lin);
    p_obs[i]  = p_true[i] * sens + (1 - p_true[i]) * (1 - spec);
  }
}
model {
  // Weakly informative priors
  mu ~ normal(0, 1.5);
  alpha_site ~ normal(0, 1);
  alpha_age  ~ normal(0, 1);
  alpha_sex  ~ normal(0, 1);

  // Priors for test characteristics: replace with your preferred Beta(a,b)
  spec ~ beta(2765.505, 194.95); #99.3 specificity
  sens ~ beta(1321.92,55.08); #96.0 sensitivity

  // Likelihood
  for (i in 1:I)
    Npos[i] ~ binomial(Ntest[i], p_obs[i]);
}
generated quantities {
  real pop_prev = 0;
  for (i in 1:I) pop_prev = pop_prev + w[i] * p_true[i];
  vector[I] p_stratum = p_true;
}
"

# Write Stan file to working directory
stan_file <- "sero_hier.stan"
writeLines(stan_code, con = stan_file)

# ---------------------------
# 5. Fit model: prefer cmdstanr if available, otherwise rstan
# ---------------------------

if (requireNamespace("cmdstanr", quietly = TRUE)) {
  library(cmdstanr)
  # ensure cmdstan installed (user may need to run cmdstanr::install_cmdstan() once)
  mod <- cmdstan_model(stan_file)
  fit <- mod$sample(data = stan_data,
                    seed = 123,
                    chains = 4, parallel_chains = 4,
                    iter_warmup = 1000, iter_sampling = 2000,
                    adapt_delta = 0.95)
  # extract draws as matrix
  draws <- fit$draws(variables = c("pop_prev", "p_stratum", "sens", "spec"), format = "draws_matrix")
  # use posterior package to convert if needed
  post_pop <- as.numeric(draws[, "pop_prev"])
  # p_stratum columns are p_stratum[1], p_stratum[2], ...
  p_cols <- grep("^p_stratum\\[", colnames(draws), value = TRUE)
  p_mat <- draws[, p_cols]   # draws x I matrix
} else {
  # fallback to rstan
  library(rstan)
  rstan_options(auto_write = TRUE)
  fit_rstan <- stan(file = stan_file, data = stan_data,
                    chains = 4, iter = 3000, warmup = 1000, cores = 4,
                    control = list(adapt_delta = 0.95), seed = 123)
  post <- rstan::extract(fit_rstan)
  post_pop <- post$pop_prev
  p_mat <- post$p_stratum   # array draws x I
}

# ---------------------------
# 6. Posterior summaries
# ---------------------------
# population-weighted prevalence
pop_mean <- mean(post_pop)
pop_ci   <- quantile(post_pop, probs = c(0.025, 0.975))

cat("Population-weighted prevalence (posterior mean):", round(pop_mean, 4), "\n")
cat("95% CrI:", round(pop_ci[1],4), "-", round(pop_ci[2],4), "\n\n")

# Stratum summaries (site × sex × age)
n_draws <- nrow(p_mat)
I <- ncol(p_mat)
stratum_post_mean <- colMeans(p_mat)
stratum_post_ci_l <- apply(p_mat, 2, quantile, 0.025)
stratum_post_ci_u <- apply(p_mat, 2, quantile, 0.975)

strata_results <- strata %>%
  mutate(post_mean = stratum_post_mean,
         post_l = stratum_post_ci_l,
         post_u = stratum_post_ci_u)

# Site-level aggregated (population-weighted within site)
site_summaries <- strata_results %>%
  group_by(site) %>%
  summarize(
    mean = sum(pop_weight * post_mean) / sum(pop_weight),   # weighted by pop_weight (same denominator per site)
    lower = NA_real_, upper = NA_real_,
    .groups = "drop"
  )

# derive site-level posterior distribution from draws
site_post_list <- lapply(unique(strata$site_id), function(sid) {
  idx <- which(strata$site_id == sid)
  
  # subset posterior draws for only those strata (columns)
  p_sub <- p_mat[, idx, drop = FALSE]  # n_draws × k
  
  # normalize weights for strata within this site
  wsub <- strata$pop_weight[idx]
  wsub <- wsub / sum(wsub)
  
  # weighted average prevalence per draw
  draws_site <- p_sub %*% wsub  # n_draws × 1
  
  c(
    mean  = mean(draws_site),
    lower = quantile(draws_site, 0.025),
    upper = quantile(draws_site, 0.975)
  )
})
site_post_df <- do.call(rbind, site_post_list) %>% as.data.frame()
site_post_df$site <- unique(strata$site)

print(site_post_df)

# ---------------------------
# 7. Plots
# ---------------------------
# Posterior density of overall prevalence
df_pop <- data.frame(pop_prev = post_pop)
ggplot(df_pop, aes(pop_prev)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = pop_mean, linetype = "dashed") +
  labs(title = "Posterior: population-weighted seroprevalence", x = "Prevalence", y = "Density")

# ---------------------------
# Site-level population-weighted posterior prevalence
# ---------------------------

# Build list of site-level posterior summaries
site_post_list <- lapply(unique(strata$site_id), function(sid) {
  idx   <- which(strata$site_id == sid)
  p_sub <- p_mat[, idx, drop = FALSE]          # subset draws for this site
  wsub  <- strata$pop_weight[idx]
  wsub  <- wsub / sum(wsub)                    # normalize weights within site
  
  draws_site <- p_sub %*% wsub                 # weighted prevalence per draw
  
  data.frame(
    site  = unique(strata$site[idx]),
    mean  = mean(draws_site),
    lower = quantile(draws_site, 0.025),
    upper = quantile(draws_site, 0.975)
  )
})

# Combine into single data frame
site_post_df <- dplyr::bind_rows(site_post_list)

# Forest plot of site-level posterior prevalence
library(ggplot2)
ggplot(site_post_df, aes(x = reorder(site, mean), y = mean)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), color = "darkorange") +
  coord_flip() +
  labs(
    title = "Site-level population-weighted posterior prevalence (95% CrI)",
    y = "Prevalence", x = ""
  ) +
  theme_minimal()

#Chikungunya script to perform population-weighted seroprevalence estimation
if (!requireNamespace("cmdstanr", quietly = TRUE) &&
    !requireNamespace("rstan", quietly = TRUE)) {
  stop("Please install either 'cmdstanr' (recommended) or 'rstan' before running this script.")
}
library(dplyr)
library(tidyr)
library(ggplot2)

# ---------------------------
# 1. Aggregate df_all into strata counts
# ---------------------------

