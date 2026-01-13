library(dplyr)
library(mgcv)
library(mixR)
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
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Robust data loader: update data_dir to where your .dta files live (or use here::here("..."))
data_dir <- "/Users/augustinemasinde/Desktop/PhD files/Global Health"  # <- update if needed

data_files <- list(
  Asembodata = file.path(data_dir, "Asembo PBIDS arbo.dta"),
  Kilifidata = file.path(data_dir, "KHDSS arbo.dta"),
  Kiberadata = file.path(data_dir, "Kibera PBIDS arbo.dta"),
  manyattadata = file.path(data_dir, "Manyatta HDSS arbo.dta"),
  NairobiUrbandata = file.path(data_dir, "NUHDSS arbo.dta")
)

# Load required libraries
library(haven)
library(dplyr)
library(mixR)
library(ggplot2)

# === 1. Load and preprocess data ===
Asembodata <- zap_labels(Asembodata)

df <- dplyr::select(Asembodata, CHKV_AIU, ageyrs, DENV2_AIU)
z_raw <- df$DENV2_AIU          # Raw titres (MFI or similar)
a <- df$ageyrs + 0.5          # Age centered (midpoint of age year)

# === 2. Log-transform titre data (critical!) ===
z_logged <- log(z_raw + 1)

# Prepare data frame for spline fitting
z_and_a <- data.frame(z = z_logged, a = a)

# === 3. Plot raw log titres vs age (to inspect data) ===
ggplot(z_and_a, aes(x = a, y = z)) +
  geom_point(alpha = 0.3) +
  labs(title = "Raw log(titre + 1) vs Age",
       x = "Age",
       y = "log(titre + 1)") +
  theme_minimal()

# === 4. Fit mixture model on logged titres > 0 ===
mfi_values_positive <- z_logged[z_logged > 0]

asembo_model <- mixfit(
  mfi_values_positive,
  family = "weibull",
  ncomp = 2
)

muS  <- asembo_model$mu[1]
muI  <- asembo_model$mu[2]
SigS <- asembo_model$sd[1]
SigI <- asembo_model$sd[2]

cat("Mixture model means (log titre): muS =", muS, ", muI =", muI, "\n")

# === 5. Define spline fitting function (mpspline.fit) ===
mpspline.fit <- function(response, x.var, ps.intervals = 40, degree = 3, order = 2, 
                         link = "identity", family = "gaussian", alpha = 0.01, kappa = 1e8) {
  y <- response
  x <- x.var
  wts <- rep(1, length(y))
  q <- degree
  d <- order
  ndx <- ps.intervals
  m.binomial <- rep(1, length(y))
  n <- length(y)
  xl <- min(x)
  xr <- max(x)
  xmax <- xr + 0.01 * (xr - xl)
  xmin <- xl - 0.01 * (xr - xl)
  dx <- (xmax - xmin)/ndx
  knots <- seq(xmin - q * dx, xmax + q * dx, by = dx)
  b <- splines::spline.des(knots, x, q + 1, 0 * x)$design
  n.col <- ncol(b)
  p <- sqrt(alpha)*diff(diff(diag(n.col)))
  mp <- sqrt(kappa)*diff(diag(n.col))
  nix <- rep(0, n.col - d)
  mnix <- rep(0, n.col - 1)
  b <- as.matrix(b)
  coef.est <- rep(1, ncol(b))
  ineqdiag <- diag(diff(coef.est) >= 1e-9)
  
  if(family == "gaussian") {
    mu <- rep(mean(y), length(y))
  }
  if(family == "binomial") {
    mu <- (y + 0.5 * m.binomial)/2
  }
  
  it <- 0
  repeat {
    if(it == 0) {
      if(link == "identity") eta <- mu
      if(link == "logit") eta <- log(mu/(m.binomial - mu))
      if(link == "probit") eta <- qnorm(mu/m.binomial)
      if(link == "cloglog") eta <- log(-log(1 - mu/m.binomial))
    }
    it <- it + 1
    if(it > 25) break
    if(link == "identity") {
      mu <- eta; h.prime <- 1
    }
    if(link == "logit") {
      mu <- m.binomial/(1 + exp(-eta))
      h.prime <- mu * (1 - mu/m.binomial)
    }
    if(link == "probit") {
      mu <- m.binomial * pnorm(eta)
      h.prime <- m.binomial * dnorm(eta)
    }
    if(link == "cloglog") {
      mu <- m.binomial * (1 - exp(-exp(eta)))
      h.prime <- m.binomial * exp(eta) * exp(-exp(eta))
    }
    if(family == "gaussian") w <- rep(1, length(y))
    if(family == "binomial") w <- h.prime^2/(mu * (1 - mu/m.binomial))
    u <- (y - mu)/h.prime + eta
    f <- lsfit(rbind(b, p, ineqdiag %*% mp), c(u, nix, mnix),
               wt = c(wts, nix + 1, mnix + 1) * c(w, (nix + 1), (mnix + 1)), intercept = FALSE)
    coef.old <- coef.est
    coef.est <- as.vector(f$coef)
    ineqdiag <- diag(diff(coef.est) < 0)
    d.coef <- max(abs((coef.est - coef.old)/coef.old))
    if(d.coef < 1e-20) break
    eta <- b %*% coef.old
  }
  
  w <- w * wts
  e <- 1e-9
  h <- hat(f$qr, intercept = FALSE)[1:n]
  trace <- sum(h) - 1
  if(family == "gaussian") {
    dev <- sum((y - eta)^2)
    dispersion.parm <- dev / (n - trace)
  }
  if(family == "binomial") {
    dev <- 2 * sum((y + e) * log((y + e)/(mu + e)) + (m.binomial - y + e) * log((m.binomial - y + e)/(m.binomial - mu + e)))
    dispersion.parm <- 1
  }
  aic <- dev + 2 * trace
  bic <- dev + log(n) * trace
  x.seq <- seq(xl, xr, by = 1)
  b.seq <- splines::spline.des(knots, x.seq, q + 1, 0 * x.seq)$design
  yhat <- b.seq %*% as.vector(coef.old)
  
  if(link == "logit") yhat <- 1/(1 + exp(-yhat))
  if(link == "probit") yhat <- apply(yhat, c(1, 2), pnorm)
  if(link == "cloglog") yhat <- 1 - exp(-exp(yhat))
  
  return(list(x = x.seq, yhat = yhat, aic = aic, bic = bic, dev = dev))
}

# === 6. Spline estimator with updated parameters ===
spline.estimator <- function(original) {
  fit <- mpspline.fit(original$z, original$a, ps.intervals = 40, degree = 3, order = 2, alpha = 0.01)
  
  fitted_means <- as.vector(fit$yhat)
  
  pred_at_data <- approx(fit$x, fitted_means, xout = original$a)$y
  residuals <- original$z - pred_at_data
  sd_resid <- sd(residuals)
  
  df_a <- data.frame(
    Age = fit$x,
    Mean = fitted_means,
    SD = sd_resid
  )
  
  derivs <- diff(fitted_means) / diff(fit$x)
  df_a_deriv <- data.frame(
    Age = fit$x[-1],
    Mean = derivs,
    SD = rep(NA, length(derivs))
  )
  
  list(df_a = df_a, df_a_deriv = df_a_deriv)
}

# === 7. Seroprevalence calculation function ===
seroprev_cis <- function(B, agemin, agemax, muS, muI, mu_a_data) {
  ages <- seq(agemin, agemax, length.out = nrow(mu_a_data))
  pi_a_samples <- matrix(NA, nrow = B, ncol = length(ages))
  
  for (b in 1:B) {
    for (i in 1:length(ages)) {
      pi_a_samples[b, i] <- (mu_a_data$Mean[i] - muS) / (muI - muS)
      pi_a_samples[b, i] <- max(min(pi_a_samples[b, i], 1), 0)
    }
  }
  
  pi_a_mean <- apply(pi_a_samples, 2, mean)
  pi_a_low <- apply(pi_a_samples, 2, quantile, probs = 0.025)
  pi_a_upp <- apply(pi_a_samples, 2, quantile, probs = 0.975)
  
  data.frame(
    Age = ages,
    pi_a = pi_a_mean,
    pi_a_low = pi_a_low,
    pi_a_upp = pi_a_upp
  )
}

# === 8. FOI calculation function ===
FOI_cis <- function(B, agemin, agemax, muI, muS, mu_a_data, mu_deriv_data) {
  ages <- seq(agemin, agemax, length.out = nrow(mu_deriv_data))
  foi_samples <- matrix(NA, nrow = B, ncol = length(ages))
  
  for (b in 1:B) {
    for (i in 1:length(ages)) {
      mu_a <- mu_a_data$Mean[i]
      mu_deriv <- mu_deriv_data$Mean[i]
      
      denom <- muI - mu_a
      if (denom <= 0) {
        foi_samples[b, i] <- NA
      } else {
        foi_samples[b, i] <- mu_deriv / denom
      }
      if (!is.na(foi_samples[b, i]) && foi_samples[b, i] < 0) {
        foi_samples[b, i] <- 0
      }
    }
  }
  
  foi_mean <- apply(foi_samples, 2, mean, na.rm = TRUE)
  foi_low <- apply(foi_samples, 2, quantile, probs = 0.025, na.rm = TRUE)
  foi_upp <- apply(foi_samples, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  data.frame(
    Age = ages,
    FOI = foi_mean,
    FOI_low = foi_low,
    FOI_upp = foi_upp
  )
}

# === 9. Main execution ===

agemin <- floor(min(a))
agemax <- ceiling(max(a))

set.seed(123)
spline_mu_output <- spline.estimator(z_and_a)

spline_mu <- spline_mu_output$df_a
spline_deriv <- spline_mu_output$df_a_deriv

B <- 500

seroprevalence <- seroprev_cis(
  B = B,
  agemin = agemin,
  agemax = agemax,
  muS = muS,
  muI = muI,
  mu_a_data = spline_mu
)

FOI <- FOI_cis(
  B = B,
  agemin = agemin,
  agemax = agemax,
  muI = muI,
  muS = muS,
  mu_a_data = spline_mu,
  mu_deriv_data = spline_deriv
)

# === 10. Output results ===
print("Seroprevalence estimates (head):")
print(head(seroprevalence))

print("FOI estimates (head):")
print(head(FOI))

cat("Max seroprevalence:", max(seroprevalence$pi_a), "\n")
cat("Mean FOI (excluding NA):", mean(FOI$FOI, na.rm = TRUE), "\n")

# === 11. Plot results ===

# Spline fit mean titres by age
ggplot(spline_mu, aes(x = Age, y = Mean)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, fill = "blue") +
  geom_hline(yintercept = muS, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = muI, linetype = "dashed", color = "green", size = 1) +
  labs(title = "Spline fit of Mean log(titre + 1) by Age",
       y = "Mean log(titre + 1)",
       x = "Age") +
  theme_minimal()

# Derivative of spline mean titres by age
ggplot(spline_deriv, aes(x = Age, y = Mean)) +
  geom_line(color = "purple", size = 1) +
  labs(title = "Derivative of Mean log(titre + 1) by Age",
       y = "Derivative",
       x = "Age") +
  theme_minimal()

# Seroprevalence by age
ggplot(seroprevalence, aes(x = Age, y = pi_a)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = pi_a_low, ymax = pi_a_upp), alpha = 0.2, fill = "blue") +
  labs(title = "Seroprevalence by Age",
       y = "Seroprevalence",
       x = "Age") +
  theme_minimal()

# FOI by age
ggplot(FOI, aes(x = Age, y = FOI)) +
  geom_line(color = "red") +
  geom_ribbon(aes(ymin = FOI_low, ymax = FOI_upp), alpha = 0.2, fill = "red") +
  labs(title = "Force of Infection (FOI) by Age",
       y = "FOI",
       x = "Age") +
  theme_minimal()

