# Load/install required packages with pacman
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dagitty, tidySEM, tidyverse, tidybayes, posterior, brms, bayestestR, gridExtra, patchwork, ggtext)

# Working directory: run this script with the "analysis" folder as working
# directory (RStudio: Session > Set Working Directory > To Source File Location).
# All file paths below are relative to that folder; the manuscript in the
# "results" folder reads the outputs via "../analysis/".
if (basename(getwd()) != "analysis")
  stop("Please set the working directory to the 'analysis' folder before running this script.")

set.seed(333)

# CHUNK MARKERS ---------------------------------------------------------------
# Lines of the form
#     ## ---- <label> ----
# are knitr chunk separators (knitr::read_chunk). 

# exemplary theory --------------------------------------------------------

#### plot ####

# Stressor–Strain–Outcome (SSO) model, originally proposed by Thomas Koeske and Gary Koeske (1993)
# conceptually consistent with the transactional stress model by Lazarus

# Define the DAG using dagitty syntax
sso_dag <- dagitty("
dag {
  Stressor -> Strain -> Outcome
}
")

# Arrange nodes in a left-to-right causal flow
lo <- get_layout(
  "Stressor", "", "Strain", "", "Outcome",
  rows = 1
)

# Plot basic model
sso_dgp <- graph_sem(sso_dag, layout = lo, 
                     edges = data.frame(
                       from = c("Stressor", "Strain"),
                       to = c("Strain", "Outcome"),
                       colour = "darkblue",
                       size = 1.2
                     ))
sso_dgp

ggsave("00_sso_dag.png", sso_dgp, width = 8, height = 4, dpi = 300)

# formal theory characteristics -------------------------------------------

n <- 500

#### simulation ####
sso_sim <- tibble(
  Stressor = runif(n, 0, 1),
  Resources = runif(n, 0, 1)
) %>%
  mutate(
    # stronger negative resource effect, less noise
    Appraisal = 0.6 + Stressor^2 - 1.8 * Resources + Resources^2 + rnorm(n, 0, 7/200),
    
    Strain = 0.0833 + 0.9 * Appraisal + rnorm(n, 0, 0.4 / 6),
    
    Outcome = 0.5 + 0.4 * Strain - 0.6 * Strain^2 + rnorm(n, 0, 0.4 / 6),
    
    Stressor2 = Stressor^2,
    Appraisal2 = Appraisal^2,
    Strain2 = Strain^2,
    Resources2 = Resources^2, 
    Intercept = 1
  )

summary(sso_sim)

#### data standardization (0–1 scaling) ####

scale_01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

sso_sim <- sso_sim %>%
  mutate(
    Stressor  = scale_01(Stressor),
    Resources = scale_01(Resources),
    Appraisal = scale_01(Appraisal),
    Strain    = scale_01(Strain),
    Outcome   = scale_01(Outcome)
  ) %>%
  mutate(
    # recompute derived terms AFTER scaling
    Stressor2  = Stressor^2,
    Resources2 = Resources^2,
    Appraisal2 = Appraisal^2,
    Strain2    = Strain^2,
    Intercept  = 1
  )

sso_sim
summary(sso_sim)

# causal test

# cf. mediation analysis

# DAG causal inference
impliedConditionalIndependencies(sso_dag)

# DAG causal inference
cis_dag <- impliedConditionalIndependencies(sso_dag)

# Bonferroni-adjusted alpha
alpha_bonf <- 0.05 / length(cis_dag)
conf_bonf <- 1 - alpha_bonf
conf_bonf

# run CI tests
ci_tests_dag <- localTests(
  sso_dag,
  sso_sim,
  type = "cis.loess",
  R = 5000,
  tests = cis_dag,
  conf.level = conf_bonf,
  abbreviate.names = FALSE
)

# helper: names of the two confidence-interval columns returned by localTests();
# they are labelled by their quantiles (e.g. "0.4%", "99.6%"), which depend on
# the confidence level, so they are selected by name rather than by position
# (fallback: the last two numeric columns)
ci_bound_cols <- function(x) {
  nm <- names(x)
  cols <- nm[str_detect(nm, "%$")]
  if (length(cols) != 2) cols <- tail(nm[sapply(x, is.numeric)], 2)
  stopifnot(length(cols) == 2)
  cols
}

# clean + interpretation
ci_cols_dag <- ci_bound_cols(ci_tests_dag)
ci_tests_dag %>%
  as.data.frame() %>%
  tibble::rownames_to_column("test") %>%
  mutate(
    ci_violated = .data[[ci_cols_dag[1]]] > 0 | .data[[ci_cols_dag[2]]] < 0,
    interpretation = ifelse(ci_violated, "CI violated", "CI supported")
  )

# Visualize data generation
p1 <- ggplot(sso_sim, aes(x = Stressor, y = Strain)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", colour = "black", linewidth = 1) +
  labs(x = "Stressor", y = "Strain",
       title = "Theory: Quadratic") +
  theme_minimal()

p2 <- ggplot(sso_sim, aes(x = Strain, y = Outcome)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", colour = "black", linewidth = 1) +
  labs(x = "Strain", y = "Outcome",
       title = "Theory: Yerkes-Dodson") +
  theme_minimal()

grid_ssodgp <- sso_dgp / (p1 + p2)
ggsave("01_sso_dgp.png", grid_ssodgp, width = 12, height = 6, dpi = 300)
grid_ssodgp

# simulation --------------------------------------------------------------

#### linear fit ####

## ---- code-sso-linear ----

library(tidyverse)
library(brms)

# Structural equations of the linear SSO formalization
eq1 <- bf(Strain ~ 0 + Intercept + Stressor)    # Stressor -> Strain
eq2 <- bf(Outcome ~ 0 + Intercept + Strain)     # Strain -> Outcome

# Parameter priors as part of the formal theory
priors <- c(
  prior(normal(0.2, 0.05), class = "b", coef = "Intercept", resp = "Strain"),
  prior(normal(0.5, 0.1), class = "b", coef = "Stressor",  resp = "Strain"),
  
  prior(normal(0.8,  0.05), class = "b", coef = "Intercept", resp = "Outcome"),
  prior(normal(-0.5, 0.1), class = "b", coef = "Strain",    resp = "Outcome"),
  
  prior(exponential(10), class = "sigma", resp = "Strain"),
  prior(exponential(10), class = "sigma", resp = "Outcome")
)

# Data simulation from the formal theory (sampling from the priors only)
bscm_linear <- brm(
  eq1 + eq2 + set_rescor(FALSE),
  data = sso_sim,
  family = gaussian(),
  prior = priors,
  sample_prior = "only",
  chains = 8, cores = 8, iter = 2500, warmup = 500,
  seed = 333
)

## ---- linear-figure-and-checks ----

summary(bscm_linear)

# Check default priors
get_prior(eq1 + eq2 + set_rescor(FALSE), data = sso_sim)

# visualize formal theory 
# Node positions
node_positions <- tibble(
  node = c("Stressor", "Strain", "Outcome"),
  x = c(0, 1, 2),
  y = c(0, 0, 0)
)

# Equations
edges_eq <- tibble(
  label = c(
    "Strain ~ N(α<sub>St</sub> + β<sub>St</sub>Stressor, σ<sub>St</sub>)",
    "Outcome ~ N(α<sub>O</sub> + β<sub>O</sub>Strain, σ<sub>O</sub>)"
  ),
  x_mid = c(0.5, 1.5),
  y_pos = c(0.0025, 0.0025)
)

# Exogenous distribution
edges_exo <- tibble(
  label = c(
    "Stressor, Strain, Outcome in [0,1]\n(min–max scaled using sample-specific bounds)",
    ""   # nothing for second edge
  ),
  x_mid = c(1, 2),
  y_pos = c(0.0015, 0.0015)   # between eq (0.02) and priors (0.01)
)

# Priors
edges_prior <- tibble(
  label = c(
    "α<sub>St</sub> ~ N(0.2, 0.05), β<sub>St</sub> ~ N(0.5, 0.1), σ<sub>St</sub> ~ Exp(10)",
    "α<sub>O</sub> ~ N(0.8, 0.05), β<sub>O</sub> ~ N(-0.5, 0.1), σ<sub>O</sub> ~ Exp(10)"
  ),
  x_mid = c(0.5, 1.5),
  y_pos = c(0.0002, 0.0002)
)

p0 <- ggplot() +
  geom_segment(
    data = tibble(x_mid = c(0.5, 1.5)),
    aes(
      x = x_mid - 0.35, y = 0,
      xend = x_mid + 0.35, yend = 0
    ),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    linewidth = 1.2,
    colour = "darkblue"
  ) +
  geom_text(
    data = node_positions,
    aes(x = x, y = y, label = node),
    size = 6,
    family = "sans"
  ) +
  geom_richtext(
    data = edges_eq,
    aes(x = x_mid, y = y_pos, label = label),
    family = "sans",
    size = 5.5,
    hjust = 0.5,
    vjust = 0,
    fill = NA,
    label.color = NA
  ) +
  geom_richtext(
    data = edges_exo,
    aes(x = x_mid, y = y_pos, label = label),
    family = "sans",
    size = 5.5,
    hjust = 0.5,
    vjust = 0,
    fill = NA,
    label.color = NA
  ) +
  geom_richtext(
    data = edges_prior,
    aes(x = x_mid, y = y_pos, label = label),
    family = "sans",
    size = 5,
    hjust = 0.5,
    vjust = 0,
    fill = NA,
    label.color = NA
  ) +
  coord_cartesian(
    xlim = c(-0.2, 2.2),   
    ylim = c(-0.001, 0.004), 
    expand = FALSE
  ) +
  theme_void() #+
#theme(
#   plot.margin = margin(-1000, 0, -1200, 0)
# )

p0

ggsave("00_sso_formaltheory_linear.png", p0, width = 10, height = 2, dpi = 300, scale = 1.1)

# Visualize prior predictive distributions
pp_check(bscm_linear, resp = "Strain", ndraws = 30)
pp_check(bscm_linear, resp = "Outcome", ndraws = 30)

Stressor_seq_plot <- expand_grid(
  Stressor = seq(min(sso_sim$Stressor), max(sso_sim$Stressor), length.out = 100)
) %>%
  mutate(
    Intercept = 1
  )

strain_seq_plot <- expand_grid(
  Strain = seq(min(sso_sim$Strain), max(sso_sim$Strain), length.out = 100)
) %>%
  mutate(
    Intercept = 1
  )

# Stressor -> Strain
epred_strain_prior <- posterior_epred(
  bscm_linear,
  newdata = Stressor_seq_plot,
  re_formula = NA,
  resp = "Strain"
)

p1 <- Stressor_seq_plot %>%
  mutate(
    estimate = apply(epred_strain_prior, 2, median, na.rm = TRUE),
    lower = apply(epred_strain_prior, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_strain_prior, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = Stressor, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") +
  geom_line(linewidth = 1) +
  geom_point(
    data = sso_sim,
    aes(x = Stressor, y = Strain),
    alpha = 0.3, size = 2, inherit.aes = FALSE
  ) +
  labs(
    x = "Stressor", y = "Strain",
    title = "Prior predictive: Stressor -> Strain (linear)"
  ) +
  theme_minimal()
p1
ggsave("02_sso_prior_linear1.png", p1, width = 8, height = 6, dpi = 300)

# Strain -> Outcome
epred_outcome_prior <- posterior_epred(
  bscm_linear,
  newdata = strain_seq_plot,
  re_formula = NA,
  resp = "Outcome"
)

p2 <- strain_seq_plot %>%
  mutate(
    estimate = apply(epred_outcome_prior, 2, median, na.rm = TRUE),
    lower = apply(epred_outcome_prior, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_outcome_prior, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = Strain, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") +
  geom_line(linewidth = 1) +
  geom_point(
    data = sso_sim,
    aes(x = Strain, y = Outcome),
    alpha = 0.3, size = 2, inherit.aes = FALSE
  ) +
  labs(
    x = "Strain", y = "Outcome",
    title = "Prior predictive: Strain -> Outcome (linear)"
  ) +
  theme_minimal()
p2
ggsave("02_sso_prior_linear2.png", p2, width = 8, height = 6, dpi = 300)

combined_plot <- (p1 + p2)
combined_plot
ggsave("02_sso_prior_linear_combined.png", combined_plot, width = 8,height = 3,dpi = 300, scale = 1.5
)

#### theory refinement - functional relationship ####

## ---- code-sso-nonlinear ----

# Refined structural equations with quadratic terms
eq1 <- bf(Strain ~ 0 + Intercept + Stressor2)
eq2 <- bf(Outcome ~ 0 + Intercept + Strain + Strain2)

priors_detail <- c(
  prior(normal(0.2, 0.05), class = "b", coef = "Intercept", resp = "Strain"),
  prior(normal(0.5,    0.1), class = "b", coef = "Stressor2",   resp = "Strain"),
  prior(exponential(10), class = "sigma", resp = "Strain"),
  
  prior(normal(0.5,   0.1), class = "b", coef = "Intercept", resp = "Outcome"),
  prior(normal(0.8,   0.1), class = "b", coef = "Strain",    resp = "Outcome"),
  prior(normal(-1.5,  0.1), class = "b", coef = "Strain2",   resp = "Outcome"),
  prior(exponential(10), class = "sigma", resp = "Outcome")
)

bscm_prior_detail <- brm(
  eq1 + eq2 + set_rescor(FALSE),
  data = sso_sim,
  family = gaussian(),
  prior = priors_detail,
  sample_prior = "only",
  chains = 8, cores = 8, iter = 2500, warmup = 500,
  seed = 333
)
## ---- nonlinear-figure-and-checks ----

summary(bscm_prior_detail)

# plotting theory
edges_eq <- tibble(
  label = c(
    "Strain ~ N(α<sub>St</sub> + β<sub>St</sub>Stressor<sup>2</sup>, σ<sub>St</sub>)",
    "Outcome ~ N(α<sub>O</sub> + β<sub>O1</sub>Strain + β<sub>O2</sub>Strain<sup>2</sup>, σ<sub>O</sub>)"
  ),
  x_mid = c(0.5, 1.5),
  y_pos = c(0.0025, 0.0025)
)

edges_prior <- tibble(
  label = c(
    "α<sub>St</sub> ~ N(0.2, 0.05), β<sub>St</sub> ~ N(0.5, 0.1),<br>σ<sub>St</sub> ~ Exp(10)",
    "α<sub>O</sub> ~ N(0.5, 0.1), β<sub>O1</sub> ~ N(0.8, 0.1),<br>β<sub>O2</sub> ~ N(-1.5, 0.1), σ<sub>O</sub> ~ Exp(10)"
  ),
  x_mid = c(0.5, 1.5),
  y_pos = c(0.0002, 0.0002)
)

p01 <- ggplot() +
  geom_segment(
    data = tibble(x_mid = c(0.5, 1.5)),
    aes(x = x_mid - 0.35, y = 0, xend = x_mid + 0.35, yend = 0),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    linewidth = 1.2,
    colour = "darkblue"
  ) +
  geom_text(
    data = node_positions,
    aes(x = x, y = y, label = node),
    size = 5.5, family = "sans"
  ) +
  geom_richtext(
    data = edges_eq,
    aes(x = x_mid, y = y_pos, label = label),
    family = "sans", size = 4.8,
    hjust = 0.5, vjust = 0, fill = NA, label.color = NA
  ) +
  geom_richtext(
    data = edges_exo,
    aes(x = x_mid, y = y_pos, label = label),
    family = "sans", size = 4.8,
    hjust = 0.5, vjust = 0, fill = NA, label.color = NA
  ) +
  geom_richtext(
    data = edges_prior,
    aes(x = x_mid, y = y_pos, label = label),
    family = "sans", size = 4.5,
    lineheight = 1.1,
    hjust = 0.5, vjust = 0, fill = NA, label.color = NA
  ) +
  coord_cartesian(
    xlim = c(-0.15, 2.15),
    ylim = c(-0.001, 0.004),
    expand = FALSE
  ) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2))

p01

ggsave("00_sso_formaltheory_nonlinear.png", p01,
       width = 10, height = 2.4, dpi = 300, scale = 1)
pp_check(bscm_prior_detail, resp = "Strain", ndraws = 30)
pp_check(bscm_prior_detail, resp = "Outcome", ndraws = 30)

# Create sequence for visualization

# grids for plotting
Stressor_seq_plot <- expand_grid(
  Stressor = seq(min(sso_sim$Stressor), max(sso_sim$Stressor), length.out = 100)
) %>%
  mutate(
    Intercept = 1,
    Stressor2 = Stressor^2
  )

strain_seq_plot <- expand_grid(
  Strain = seq(min(sso_sim$Strain), max(sso_sim$Strain), length.out = 100)
) %>%
  mutate(
    Intercept = 1,
    Strain2 = Strain^2
  )

# Stressor -> Strain
epred_strain_prior <- posterior_epred(
  bscm_prior_detail,
  newdata = Stressor_seq_plot,
  re_formula = NA,
  resp = "Strain"
)

p1 <- Stressor_seq_plot %>%
  mutate(
    estimate = apply(epred_strain_prior, 2, median, na.rm = TRUE),
    lower = apply(epred_strain_prior, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_strain_prior, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = Stressor, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") +
  geom_line(linewidth = 1) +
  geom_point(
    data = sso_sim,
    aes(x = Stressor, y = Strain),
    alpha = 0.3, size = 2, inherit.aes = FALSE
  ) +
  labs(
    x = "Stressor", y = "Strain",
    title = "Prior predictive: Stressor -> Strain (non-linear)"
  ) +
  theme_minimal()
p1
ggsave("02_sso_prior_nonlinear1.png", p1, width = 8, height = 6, dpi = 300)

# Strain -> Outcome
epred_outcome_prior <- posterior_epred(
  bscm_prior_detail,
  newdata = strain_seq_plot,
  re_formula = NA,
  resp = "Outcome"
)

p2 <- strain_seq_plot %>%
  mutate(
    estimate = apply(epred_outcome_prior, 2, median, na.rm = TRUE),
    lower = apply(epred_outcome_prior, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_outcome_prior, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = Strain, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") +
  geom_line(linewidth = 1) +
  geom_point(
    data = sso_sim,
    aes(x = Strain, y = Outcome),
    alpha = 0.3, size = 2, inherit.aes = FALSE
  ) +
  labs(
    x = "Strain", y = "Outcome",
    title = "Prior predictive: Strain -> Outcome (non-linear)"
  ) +
  theme_minimal()
p2
ggsave("02_sso_prior_nonlinear2.png", p2, width = 8, height = 6, dpi = 300)

combined_plot2 <- p01 / (p1 + p2)
combined_plot2
ggsave("02_sso_prior_nonlinear_combined.png", combined_plot2, width = 14,height = 6,dpi = 300, scale = 1
)

c <- combined_plot / combined_plot2
c
ggsave("02_sso_prior_lin_nonlin_combined.png", c, width = 8, height = 6,dpi = 300, scale = 1.7
)


# direct inference: fitting formal theory to data -------------------------

## ---- code-sso-update ----

# Fit the formal theory to the data (posterior updating);
# formula and priors as specified in the formal theory above
bscm_detail <- update(
  bscm_prior_detail,
  sample_prior = "yes",
  prior = priors_detail, 
  chains = 8, cores = 8, iter = 2500, warmup = 500,
  seed = 333
)
## ---- update-diagnostics ----

summary(bscm_detail)

# Visualize posterior predictive distributions
pp_check(bscm_detail, resp = "Strain", ndraws = 30)
pp_check(bscm_detail, resp = "Outcome", ndraws = 30)

# Create sequence for visualization
Stressor_seq_plot <- expand_grid(
  Stressor = seq(min(sso_sim$Stressor), max(sso_sim$Stressor), length.out = 100)
) %>%
  mutate(
    Intercept = 1,
    Stressor2 = Stressor^2
  )

strain_seq_plot <- expand_grid(
  Strain = seq(min(sso_sim$Strain), max(sso_sim$Strain), length.out = 100)
) %>%
  mutate(
    Intercept = 1,
    Strain2 = Strain^2
  )
# Stressor -> Strain
epred_strain_post <- posterior_epred(
  bscm_detail,
  newdata = Stressor_seq_plot,
  re_formula = NA,
  resp = "Strain"
)

p1 <- Stressor_seq_plot %>%
  mutate(
    estimate = apply(epred_strain_post, 2, median, na.rm = TRUE),
    lower = apply(epred_strain_post, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_strain_post, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = Stressor, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") +
  geom_line(linewidth = 1) +
  geom_point(
    data = sso_sim,
    aes(x = Stressor, y = Strain),
    alpha = 0.3, size = 2, inherit.aes = FALSE
  ) +
  labs(
    x = "Stressor", y = "Strain",
    title = "Posterior fitted: Stressor -> Strain"
  ) +
  theme_minimal()
p1
ggsave("03_sso_posterior_nonlinear1.png", p1, width = 8, height = 6, dpi = 300)

# Strain -> Outcome
epred_outcome_post <- posterior_epred(
  bscm_detail,
  newdata = strain_seq_plot,
  re_formula = NA,
  resp = "Outcome"
)

p2 <- strain_seq_plot %>%
  mutate(
    estimate = apply(epred_outcome_post, 2, median, na.rm = TRUE),
    lower = apply(epred_outcome_post, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_outcome_post, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = Strain, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") +
  geom_line(linewidth = 1) +
  geom_point(
    data = sso_sim,
    aes(x = Strain, y = Outcome),
    alpha = 0.3, size = 2, inherit.aes = FALSE
  ) +
  labs(
    x = "Strain", y = "Outcome",
    title = "Posterior fitted: Strain -> Outcome"
  ) +
  theme_minimal()
p2
ggsave("03_sso_posterior_nonlinear2.png", p2, width = 8, height = 6, dpi = 300)

c <- p1 + p2
c
ggsave("03_sso_posterior_nonlin_comb.png", c, width = 14,height = 6,dpi = 300, scale = 0.7
)


# indirect inference ------------------------------------------------------

#### comparing prior and posterior of the nonlinear model ####

# 1. SIDE-BY-SIDE PARAMETER COMPARISON: PRIOR VS POSTERIOR

draws_detail <- as_draws_df(bscm_detail) %>%
  as_tibble()

param_mapping <- c(
  "b_Strain_Intercept"  = "prior_b_Strain_Intercept",
  "b_Strain_Stressor2"    = "prior_b_Strain_Stressor2",
  "b_Outcome_Intercept" = "prior_b_Outcome_Intercept",
  "b_Outcome_Strain"    = "prior_b_Outcome_Strain",
  "b_Outcome_Strain2"   = "prior_b_Outcome_Strain2"
)

posterior_params <- draws_detail %>%
  select(all_of(names(param_mapping))) %>%
  pivot_longer(everything(), names_to = "name", values_to = "value") %>%
  mutate(type = "Posterior")

prior_params <- draws_detail %>%
  select(all_of(unname(param_mapping))) %>%
  pivot_longer(everything(), names_to = "name", values_to = "value") %>%
  mutate(
    type = "Prior",
    name = recode(name, !!!setNames(names(param_mapping), unname(param_mapping)))
  )

bind_rows(prior_params, posterior_params) %>%
  ggplot(aes(x = value, fill = type, colour = type)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  facet_wrap(~name, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Prior" = "blue", "Posterior" = "red")) +
  scale_colour_manual(values = c("Prior" = "darkblue", "Posterior" = "darkred")) +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(
    title = "Prior vs.Posterior Parameter Distributions",
    fill = "Distribution",
    colour = "Distribution",
    x = "Parameter value",
    y = "Density"
  )
ggsave("04_sso_indir1.png", width = 8, height = 6, dpi = 300)


# 2. PRIOR PREDICTIVE VS POSTERIOR PREDICTIVE CHECK

pp_check(bscm_prior_detail, resp = "Strain", ndraws = 100) +
  ggtitle("Prior predictive check: Strain")
ggsave("04_sso_indir_prior1.png", width = 8, height = 6, dpi = 300)

pp_check(bscm_detail, resp = "Strain", ndraws = 100) +
  ggtitle("Posterior predictive check: Strain")
ggsave("04_sso_indir_post1.png", width = 8, height = 6, dpi = 300)

pp_check(bscm_prior_detail, resp = "Outcome", ndraws = 100) +
  ggtitle("Prior predictive check: Outcome")
ggsave("04_sso_indir_prior2.png", width = 8, height = 6, dpi = 300)

pp_check(bscm_detail, resp = "Outcome", ndraws = 100) +
  ggtitle("Posterior predictive check: Outcome")
ggsave("04_sso_indir_post2.png", width = 8, height = 6, dpi = 300)

# 3. COMPARE PREDICTED CURVES: PRIOR VS POSTERIOR

Stressor_seq_plot <- expand_grid(
  Stressor = seq(min(sso_sim$Stressor), max(sso_sim$Stressor), length.out = 100)
) %>%
  mutate(
    Intercept = 1,
    Stressor2 = Stressor^2
  )

strain_seq_plot <- expand_grid(
  Strain = seq(min(sso_sim$Strain), max(sso_sim$Strain), length.out = 100)
) %>%
  mutate(
    Intercept = 1,
    Strain2 = Strain^2
  )

# Stressor -> Strain
epred_strain_prior <- posterior_epred(
  bscm_prior_detail,
  newdata = Stressor_seq_plot,
  re_formula = NA,
  resp = "Strain"
)

epred_strain_post <- posterior_epred(
  bscm_detail,
  newdata = Stressor_seq_plot,
  re_formula = NA,
  resp = "Strain"
)

prior_strain <- Stressor_seq_plot %>%
  mutate(
    estimate = apply(epred_strain_prior, 2, median, na.rm = TRUE),
    lower = apply(epred_strain_prior, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_strain_prior, 2, quantile, probs = 0.975, na.rm = TRUE),
    type = "Prior"
  )

posterior_strain <- Stressor_seq_plot %>%
  mutate(
    estimate = apply(epred_strain_post, 2, median, na.rm = TRUE),
    lower = apply(epred_strain_post, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_strain_post, 2, quantile, probs = 0.975, na.rm = TRUE),
    type = "Posterior"
  )

indir31 <- bind_rows(prior_strain, posterior_strain) %>%
  ggplot(aes(x = Stressor, y = estimate, fill = type, colour = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = sso_sim,
    aes(x = Stressor, y = Strain, fill = NULL, colour = NULL),
    alpha = 0.3, size = 2
  ) +
  scale_fill_manual(values = c("Prior" = "blue", "Posterior" = "red")) +
  scale_colour_manual(values = c("Prior" = "darkblue", "Posterior" = "darkred")) +
  labs(
    x = "Stressor", y = "Strain",
    title = "Prior vs.Posterior Expected Curve: Stressor -> Strain",
    fill = "Model", colour = "Model"
  ) +
  theme_minimal()
indir31
ggsave("04_sso_indir3_1.png", indir31, width = 8, height = 6, dpi = 300)

# Strain -> Outcome

epred_outcome_prior <- posterior_epred(
  bscm_prior_detail,
  newdata = strain_seq_plot,
  re_formula = NA,
  resp = "Outcome"
)

epred_outcome_post <- posterior_epred(
  bscm_detail,
  newdata = strain_seq_plot,
  re_formula = NA,
  resp = "Outcome"
)

prior_outcome <- strain_seq_plot %>%
  mutate(
    estimate = apply(epred_outcome_prior, 2, median, na.rm = TRUE),
    lower = apply(epred_outcome_prior, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_outcome_prior, 2, quantile, probs = 0.975, na.rm = TRUE),
    type = "Prior"
  )

posterior_outcome <- strain_seq_plot %>%
  mutate(
    estimate = apply(epred_outcome_post, 2, median, na.rm = TRUE),
    lower = apply(epred_outcome_post, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_outcome_post, 2, quantile, probs = 0.975, na.rm = TRUE),
    type = "Posterior"
  )

indir32 <- bind_rows(prior_outcome, posterior_outcome) %>%
  ggplot(aes(x = Strain, y = estimate, fill = type, colour = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = sso_sim,
    aes(x = Strain, y = Outcome, fill = NULL, colour = NULL),
    alpha = 0.3, size = 2
  ) +
  scale_fill_manual(values = c("Prior" = "blue", "Posterior" = "red")) +
  scale_colour_manual(values = c("Prior" = "darkblue", "Posterior" = "darkred")) +
  labs(
    x = "Strain", y = "Outcome",
    title = "Prior vs.Posterior Expected Curve: Strain -> Outcome",
    fill = "Model", colour = "Model"
  ) +
  theme_minimal()
indir32
ggsave("04_sso_indir3_2.png", indir32, width = 8, height = 6, dpi = 300)

# 4. PARAMETER-UPDATING SUMMARY TABLE
param_names <- names(param_mapping)
prior_names <- unname(param_mapping)

prior_summary <- draws_detail %>%
  select(all_of(prior_names)) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>%
  mutate(Parameter = recode(Parameter, !!!setNames(param_names, prior_names))) %>%
  group_by(Parameter) %>%
  summarise(
    Prior_Mean = mean(value, na.rm = TRUE),
    Prior_SD = sd(value, na.rm = TRUE),
    Prior_Lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
    Prior_Upper_95 = quantile(value, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

posterior_summary_tbl <- draws_detail %>%
  select(all_of(param_names)) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>%
  group_by(Parameter) %>%
  summarise(
    Posterior_Mean = mean(value, na.rm = TRUE),
    Posterior_SD = sd(value, na.rm = TRUE),
    Posterior_Lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
    Posterior_Upper_95 = quantile(value, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

coef_comparison <- prior_summary %>%
  left_join(posterior_summary_tbl, by = "Parameter") %>%
  mutate(
    Parameter = recode(
      Parameter,
      "b_Strain_Intercept" = "Strain_Intercept",
      "b_Strain_Stressor2" = "Strain_Stressor2",
      "b_Outcome_Intercept" = "Outcome_Intercept",
      "b_Outcome_Strain" = "Outcome_Strain",
      "b_Outcome_Strain2" = "Outcome_Strain2"
    )
  ) %>%
  select(
    Parameter,
    Prior_Mean, Prior_SD, 
    Posterior_Mean, Posterior_SD
  )

print(coef_comparison)


#### comparing linear and nonlinear ####

# PRIOR model (linear)
draws_prior <- as_draws_df(bscm_linear) %>%
  as_tibble()

# POSTERIOR model (quadratic)
draws_post <- as_draws_df(bscm_detail) %>%
  as_tibble()
# tentative comparison: the two models do not share the same parameters

# 2. PRIOR PREDICTIVE VS POSTERIOR PREDICTIVE CHECK

pp_check(bscm_linear, resp = "Strain", ndraws = 100) +
  ggtitle("Prior predictive check (linear model): Strain")

pp_check(bscm_detail, resp = "Strain", ndraws = 100) +
  ggtitle("Posterior predictive check (quadratic model): Strain")

pp_check(bscm_linear, resp = "Outcome", ndraws = 100) +
  ggtitle("Prior predictive check (linear model): Outcome")

pp_check(bscm_detail, resp = "Outcome", ndraws = 100) +
  ggtitle("Posterior predictive check (quadratic model): Outcome")


# 3. COMPARE PREDICTED CURVES: LINEAR PRIOR MODEL VS QUADRATIC POSTERIOR MODEL

Stressor_seq_plot <- expand_grid(
  Stressor = seq(min(sso_sim$Stressor), max(sso_sim$Stressor), length.out = 100)
) %>%
  mutate(
    Intercept = 1,
    Stressor2 = Stressor^2
  )

strain_seq_plot <- expand_grid(
  Strain = seq(min(sso_sim$Strain), max(sso_sim$Strain), length.out = 100)
) %>%
  mutate(
    Intercept = 1,
    Strain2 = Strain^2
  )

# Stressor -> Strain

epred_strain_prior <- posterior_epred(
  bscm_linear,
  newdata = Stressor_seq_plot,
  re_formula = NA,
  resp = "Strain"
)

epred_strain_post <- posterior_epred(
  bscm_detail,
  newdata = Stressor_seq_plot,
  re_formula = NA,
  resp = "Strain"
)

prior_strain <- Stressor_seq_plot %>%
  mutate(
    estimate = apply(epred_strain_prior, 2, median, na.rm = TRUE),
    lower = apply(epred_strain_prior, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_strain_prior, 2, quantile, probs = 0.975, na.rm = TRUE),
    type = "Prior model"
  )

posterior_strain <- Stressor_seq_plot %>%
  mutate(
    estimate = apply(epred_strain_post, 2, median, na.rm = TRUE),
    lower = apply(epred_strain_post, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_strain_post, 2, quantile, probs = 0.975, na.rm = TRUE),
    type = "Posterior model"
  )

indir41 <- bind_rows(prior_strain, posterior_strain) %>%
  ggplot(aes(x = Stressor, y = estimate, fill = type, colour = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = sso_sim,
    aes(x = Stressor, y = Strain, fill = NULL, colour = NULL),
    alpha = 0.3, size = 2
  ) +
  scale_fill_manual(values = c("Prior model" = "blue", "Posterior model" = "red")) +
  scale_colour_manual(values = c("Prior model" = "darkblue", "Posterior model" = "darkred")) +
  labs(
    x = "Stressor", y = "Strain",
    title = "Linear Prior Model vs Quadratic Posterior Model: Stressor -> Strain",
    fill = "Model", colour = "Model"
  ) +
  theme_minimal()
indir41
ggsave("04_sso_indir4_1.png", indir41, width = 8, height = 6, dpi = 300)

# Strain -> Outcome

epred_outcome_prior <- posterior_epred(
  bscm_linear,
  newdata = strain_seq_plot,
  re_formula = NA,
  resp = "Outcome"
)

epred_outcome_post <- posterior_epred(
  bscm_detail,
  newdata = strain_seq_plot,
  re_formula = NA,
  resp = "Outcome"
)

prior_outcome <- strain_seq_plot %>%
  mutate(
    estimate = apply(epred_outcome_prior, 2, median, na.rm = TRUE),
    lower = apply(epred_outcome_prior, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_outcome_prior, 2, quantile, probs = 0.975, na.rm = TRUE),
    type = "Prior model"
  )

posterior_outcome <- strain_seq_plot %>%
  mutate(
    estimate = apply(epred_outcome_post, 2, median, na.rm = TRUE),
    lower = apply(epred_outcome_post, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(epred_outcome_post, 2, quantile, probs = 0.975, na.rm = TRUE),
    type = "Posterior model"
  )

indir42 <- bind_rows(prior_outcome, posterior_outcome) %>%
  ggplot(aes(x = Strain, y = estimate, fill = type, colour = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = sso_sim,
    aes(x = Strain, y = Outcome, fill = NULL, colour = NULL),
    alpha = 0.3, size = 2
  ) +
  scale_fill_manual(values = c("Prior model" = "blue", "Posterior model" = "red")) +
  scale_colour_manual(values = c("Prior model" = "darkblue", "Posterior model" = "darkred")) +
  labs(
    x = "Strain", y = "Outcome",
    title = "Linear Prior Model vs Quadratic Posterior Model: Strain -> Outcome",
    fill = "Model", colour = "Model"
  ) +
  theme_minimal()
indir42
ggsave("04_sso_indir4_2.png", indir42, width = 8, height = 6, dpi = 300)
c <- (indir41 + indir42) / (indir31 + indir32)
c
ggsave("04_sso_indir_all_combined.png", c, width = 16,height = 12,dpi = 300, scale = 0.7
)

# 4. MODEL-SPECIFIC PARAMETER SUMMARY TABLE

prior_summary_tbl <- draws_prior %>%
  select(
    b_Strain_Intercept,
    b_Strain_Stressor,
    b_Outcome_Intercept,
    b_Outcome_Strain
  ) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>%
  group_by(Parameter) %>%
  summarise(
    Prior_Mean = mean(value, na.rm = TRUE),
    Prior_SD = sd(value, na.rm = TRUE),
    Prior_Lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
    Prior_Upper_95 = quantile(value, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

posterior_summary_tbl <- draws_post %>%
  select(
    b_Strain_Intercept,
    b_Strain_Stressor2,
    b_Outcome_Intercept,
    b_Outcome_Strain,
    b_Outcome_Strain2
  ) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>%
  group_by(Parameter) %>%
  summarise(
    Posterior_Mean = mean(value, na.rm = TRUE),
    Posterior_SD = sd(value, na.rm = TRUE),
    Posterior_Lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
    Posterior_Upper_95 = quantile(value, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

coef_comparison <- full_join(prior_summary_tbl, posterior_summary_tbl, by = "Parameter") %>%
  mutate(
    Parameter = recode(
      Parameter,
      "b_Strain_Intercept" = "Strain_Intercept",
      "b_Strain_Stressor" = "Strain_Stressor",
      "b_Strain_Stressor2" = "Strain_Stressor2",
      "b_Outcome_Intercept" = "Outcome_Intercept",
      "b_Outcome_Strain" = "Outcome_Strain",
      "b_Outcome_Strain2" = "Outcome_Strain2"
    )
  ) %>%
  select(
    Parameter,
    Prior_Mean, Prior_SD,
    Posterior_Mean, Posterior_SD
  )

print(coef_comparison)


# latent variable resp. missing path -------------------------------------------------------

## ---- code-sso-latent ----

# Latent strain as a fully missing variable
sso_latent <- sso_sim %>%
  mutate(
    L_Strain = NA_real_   # latent variable, unobserved for all observations
  )

priors_latent <- c(
  # Latent equation: narrow priors serve identification, not theoretical commitment
  prior(normal(0.2, 0.01),  class = "b", coef = "Intercept",     resp = "LStrain"),
  prior(normal(0.5, 0.01),  class = "b", coef = "Stressor2",     resp = "LStrain"),
  
  # Outcome module: theoretical priors as in the refined SSO formalization
  prior(normal(0.5,  0.1),  class = "b", coef = "Intercept",     resp = "Outcome"),
  prior(normal(0.8,  0.1),  class = "b", coef = "miL_Strain",    resp = "Outcome"),
  prior(normal(-1.5, 0.1),  class = "b", coef = "ImiL_StrainE2", resp = "Outcome"),
  prior(exponential(10),    class = "sigma",                     resp = "Outcome")
)

# Standard deviation of the latent variable fixed to a constant directly in
# the model formula (required for identification)
eq0 <- bf(L_Strain | mi() ~ 0 + Intercept + Stressor2,
          sigma = 0.1)

eq1 <- bf(Outcome ~ 0 + Intercept + mi(L_Strain) + I(mi(L_Strain)^2))

bscm_latent_detail <- brm(
  eq0 + eq1 + set_rescor(FALSE),
  data    = sso_latent,
  family  = gaussian(),
  prior   = priors_latent,
  chains  = 8, cores = 8, iter = 3000, warmup = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed    = 333
)
## ---- latent-diagnostics ----

summary(bscm_latent_detail)

pp_check(bscm_latent_detail, resp = "Outcome",  ndraws = 30)

# Posterior-inferred latent strain per observation: the imputed values
# Ymi_LStrain[i] are informed both by the stressor (via the structural equation
# for L_Strain) and by the outcome (via the likelihood of the Outcome equation).
# Note: posterior_epred(..., resp = "LStrain") would instead return
# E[L_Strain | Stressor], i.e., the regression line, which does not use the
# information carried by the outcome.
strain_posterior_median <- as_draws_df(bscm_latent_detail) %>%
  as_tibble() %>%
  select(starts_with("Ymi_LStrain")) %>%
  apply(2, median)

# Scatterplot
# observed
scatter_p1 <- tibble(
  Stressor = sso_sim$Stressor,
  Strain   = sso_sim$Strain        
)

# inferred, observed
scatter_p2 <- tibble(
  L_Strain = strain_posterior_median,  
  Outcome  = sso_sim$Outcome           
)

# sequences
Stressor_seq_plot <- expand_grid(
  Stressor = seq(min(sso_sim$Stressor), max(sso_sim$Stressor), length.out = 100)
) %>%
  mutate(
    Intercept = 1,
    Stressor2 = Stressor^2,
    L_Strain  = NA_real_
  )

L_Strain_seq_plot <- tibble(
  L_Strain  = seq(0, 1, length.out = 200),
  Intercept = 1,
  Stressor2 = NA_real_
)

# posteriors
epred_strain_post <- posterior_epred(
  bscm_latent_detail,
  newdata    = Stressor_seq_plot,
  re_formula = NA,
  resp       = "LStrain"
)

epred_outcome_post <- posterior_epred(
  bscm_latent_detail,
  newdata    = L_Strain_seq_plot,
  re_formula = NA,
  resp       = "Outcome"
)

p1 <- Stressor_seq_plot %>%
  mutate(
    estimate = apply(epred_strain_post, 2, median,   na.rm = TRUE),
    lower    = apply(epred_strain_post, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper    = apply(epred_strain_post, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = Stressor, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") +
  geom_line(linewidth = 1) +
  geom_point(
    data = scatter_p1,
    aes(x = Stressor, y = Strain),   # observed
    alpha = 0.3, size = 2, inherit.aes = FALSE
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    x       = "Stressor",
    y       = "L_Strain (latent, posterior-inferred)",
    title   = "Posterior: Stressor -> L_Strain"
    #,caption = "Circles: Observed Strain"
  ) +
  theme_minimal()

p2 <- L_Strain_seq_plot %>%
  mutate(
    estimate = apply(epred_outcome_post, 2, median,   na.rm = TRUE),
    lower    = apply(epred_outcome_post, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper    = apply(epred_outcome_post, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = L_Strain, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") +
  geom_line(linewidth = 1) +
  geom_point(
    data = scatter_p2,
    aes(x = L_Strain, y = Outcome),  # posterior-inferred L_Strain and observed Outcome
    alpha = 0.3, size = 2, inherit.aes = FALSE
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    x       = "L_Strain (latent, posterior-inferred)",
    y       = "Outcome",
    title   = "Posterior: L_Strain -> Outcome (non-linear)"
  ) +
  theme_minimal()

# save
p1
ggsave("sso_posterior_latent_strain1.png", p1, width = 8, height = 6, dpi = 300)

p2
ggsave("sso_posterior_latent_strain2.png", p2, width = 8, height = 6, dpi = 300)

c <- p1 + p2
c
ggsave("sso_posterior_latent_strain_comb.png", c, width = 14, height = 6, dpi = 300, scale = 0.7)



# modularity - SSO_LAZ --------------------------------------------------------------

# replace one theory node by another theory

# DAG of the augmented theory (used for the graphical representation below)
sso_dag_lazarus <- dagitty("
dag {
  Stressor -> Appraisal
  Resources -> Appraisal
  Appraisal -> Strain
  Strain -> Outcome
}
")


#### causality ####

## ---- code-sso-ci-tests ----

library(dagitty)

sso_dag_lazarus_rev <- dagitty("
dag {
  Stressor -> Appraisal
  Appraisal -> Resources
  Appraisal -> Strain
  Strain -> Outcome
}
")

# Conditional independencies implied by the DAG
cis_dag_laz <- impliedConditionalIndependencies(sso_dag_lazarus_rev)

# Bonferroni-adjusted confidence level
alpha_bonf <- 0.05 / length(cis_dag_laz)
conf_bonf <- 1 - alpha_bonf

# Test the implied conditional independencies in the data
# using nonparametric (LOESS-based) regression
ci_tests_dag_laz <- localTests(
  sso_dag_lazarus_rev,
  sso_sim,
  type = "cis.loess",
  R = 5000,
  tests = cis_dag_laz,
  conf.level = conf_bonf,
  abbreviate.names = FALSE
)

## ---- ci-table-formatting ----

cis_dag_laz
conf_bonf
ci_tests_dag_laz

ci_cols_laz <- ci_bound_cols(ci_tests_dag_laz)

# Format a dagitty test label "A _||_ B | C" as math for the manuscript table,
# matching the notation used in the text (renders in both docx and PDF and
# leaves no ASCII pipe that could break markdown tables)
fmt_ci_test <- function(x) {
  parts <- str_match(x, "^(\\S+) _\\|\\|_ (\\S+) \\| (.+)$")
  stopifnot(!anyNA(parts[, 1]))
  sprintf("$\\text{%s} \\perp \\text{%s} \\mid \\text{%s}$",
          parts[, 2], parts[, 3], parts[, 4])
}

ci_daglaz_clean <- ci_tests_dag_laz %>%
  as.data.frame() %>%
  rownames_to_column("test") %>%
  mutate(test = fmt_ci_test(test))
ci_daglaz_clean <- ci_daglaz_clean %>%
  mutate(
    ci_lower = .data[[ci_cols_laz[1]]],
    ci_upper = .data[[ci_cols_laz[2]]],
    ci_violated = ci_lower > 0 | ci_upper < 0,
    interpretation = ifelse(ci_violated, "CI violated", "CI supported"),
    
    `CI (Bonf. adj.)` = sprintf("[%.3f, %.3f]", ci_lower, ci_upper)
  ) %>%
  select(
    Test = test,
    Estimate = estimate,
    `CI (Bonf. adj.)`,
    Interpretation = interpretation
  ) %>%
  mutate(Estimate = round(Estimate, 3))
ci_daglaz_clean

p_table <- tableGrob(ci_daglaz_clean, rows = NULL, theme = ttheme_default())
p_table
ggsave("05_ssolaz_causal.png", p_table, width = 14, height = 6, dpi = 300)


#### DAG plot ####

# Arrange nodes in a left-to-right causal flow
lo <- get_layout(
  "Stressor", "Appraisal", "Strain", "Outcome",
  "",         "Resources", "",       "",
  rows = 2
)

# Plot refined model
ssolaz_dgp <- graph_sem(
  sso_dag_lazarus,
  layout = lo,
  edges = data.frame(
    from = c("Stressor", "Resources", "Appraisal", "Strain"),
    to   = c("Appraisal", "Appraisal", "Strain", "Outcome"),
    colour = "darkblue",
    size = 1.2
  )
)
ssolaz_dgp


# Visualize refined data generation

# Create individual plots
p1 <- ggplot(sso_sim, aes(x = Stressor, y = Appraisal, colour = Resources)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, colour = "black", linewidth = 1) +
  labs(
    x = "Stressor",
    y = "Appraisal",
    title = "Refined theory: Stressor -> Appraisal",
    colour = "Resources"
  ) +
  theme_minimal() +
  scale_colour_gradient(low = "cornflowerblue", high = "black")

p2 <- ggplot(sso_sim, aes(x = Resources, y = Appraisal, colour = Stressor)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, colour = "black", linewidth = 1) +
  labs(
    x = "Resources",
    y = "Appraisal",
    title = "Refined theory: Resources -> Appraisal",
    colour = "Stressor"
  ) +
  theme_minimal() +
  scale_colour_gradient(low = "cornflowerblue", high = "black")

p3 <- ggplot(sso_sim, aes(x = Appraisal, y = Strain, colour = Stressor)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", colour = "black", linewidth = 1) +
  labs(
    x = "Appraisal",
    y = "Strain",
    title = "Refined theory: Appraisal -> Strain",
    colour = "Stressor"
  ) +
  theme_minimal() +
  scale_colour_gradient(low = "cornflowerblue", high = "black")

p4 <- ggplot(sso_sim, aes(x = Stressor, y = Strain)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", colour = "black", linewidth = 1) +
  labs(
    x = "Stressor",
    y = "Strain",
    title = "Coarse theory still approximately true: Stressor -> Strain"
  ) +
  theme_minimal()

p5 <- ggplot(sso_sim, aes(x = Strain, y = Outcome, colour = Appraisal)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", colour = "black", linewidth = 1) +
  labs(
    x = "Strain",
    y = "Outcome",
    title = "Theory: Strain -> Outcome (Yerkes-Dodson)",
    colour = "Appraisal"
  ) +
  theme_minimal() +
  scale_colour_gradient(low = "cornflowerblue", high = "black")

# Combine into grid and save
grid_plot <- (ssolaz_dgp + p1) / (p2 + p3) / (p4 + p5)

ggsave("05_ssolaz_dgp2.png", grid_plot, width = 8, height = 6, dpi = 300, scale = 1.5)

# Optionally display
grid_plot

#### indirect inference ####

## ---- code-ssolaz-modularity ----

# Augmented structural equations: the direct Stressor -> Strain path is
# replaced by an appraisal module; the Strain -> Outcome module is unchanged
eq1 <- bf(Appraisal ~ 0 + Intercept + Stressor2 + Resources + Resources2)
eq2 <- bf(Strain ~ 0 + Intercept + Appraisal)
eq3 <- bf(Outcome ~ 0 + Intercept + Strain + Strain2)

# Priors reflecting rough theoretical knowledge about the new module
priors_ssolaz_detail <- c(
  # Appraisal ~ Intercept + Stressor2 + Resources + Resources2
  prior(normal(0.75, 0.1), class = "b", coef = "Intercept",  resp = "Appraisal"),
  prior(normal(1, 0.3), class = "b", coef = "Stressor2",    resp = "Appraisal"),
  prior(normal(-2, 0.3), class = "b", coef = "Resources", resp = "Appraisal"),
  prior(normal(1, 0.3), class = "b", coef = "Resources2", resp = "Appraisal"),
  prior(exponential(16), class = "sigma", resp = "Appraisal"),
  
  # Strain ~ Intercept + Appraisal
  prior(normal(0.10, 0.1), class = "b", coef = "Intercept",  resp = "Strain"),
  prior(normal(1, 0.3), class = "b", coef = "Appraisal",  resp = "Strain"),
  prior(exponential(16), class = "sigma", resp = "Strain"),
  
  # Outcome ~ Intercept + Strain + Strain2
  # (unchanged from the refined SSO formalization above)
  prior(normal(0.5,  0.1), class = "b", coef = "Intercept", resp = "Outcome"),
  prior(normal(0.8,  0.1), class = "b", coef = "Strain",    resp = "Outcome"),
  prior(normal(-1.5, 0.1), class = "b", coef = "Strain2",   resp = "Outcome"),
  prior(exponential(10), class = "sigma", resp = "Outcome")
)

# Single fit; sample_prior = "yes" retains prior draws next to the posterior,
# allowing the prior vs. posterior comparisons shown below
bscm_ssolaz_detail <- brm(
  eq1 + eq2 + eq3 + set_rescor(FALSE),
  data = sso_sim,
  family = gaussian(),
  prior = priors_ssolaz_detail,
  sample_prior = "yes",
  chains = 8, cores = 8, iter = 2500, warmup = 500,
  seed = 333
)
## ---- modularity-diagnostics ----

summary(bscm_ssolaz_detail)

pp_check(bscm_ssolaz_detail, resp = "Appraisal", ndraws = 30)
pp_check(bscm_ssolaz_detail, resp = "Strain", ndraws = 30)
pp_check(bscm_ssolaz_detail, resp = "Outcome", ndraws = 30)


# draw table with posterior and prior samples
draws_detail <- as_draws_df(bscm_ssolaz_detail) %>% 
  as_tibble()

# grids
appraisal_seq_plot <- expand_grid(
  Stressor = seq(min(sso_sim$Stressor), max(sso_sim$Stressor), length.out = 100),
  Resources = mean(sso_sim$Resources)
) %>%
  mutate(
    Intercept = 1,
    Stressor2 = Stressor^2,
    Resources2 = Resources^2
  )

resources_seq_plot <- expand_grid(
  Resources = seq(min(sso_sim$Resources), max(sso_sim$Resources), length.out = 100),
  Stressor = mean(sso_sim$Stressor)
) %>%
  mutate(
    Intercept = 1,
    Stressor2 = Stressor^2,
    Resources2 = Resources^2
  )

strain_seq_plot <- expand_grid(
  Appraisal = seq(min(sso_sim$Appraisal), max(sso_sim$Appraisal), length.out = 100)
) %>%
  mutate(
    Intercept = 1
  )

outcome_seq_plot <- expand_grid(
  Strain = seq(min(sso_sim$Strain), max(sso_sim$Strain), length.out = 100)
) %>%
  mutate(
    Intercept = 1,
    Strain2 = Strain^2
  )

# posterior epreds
epred_appraisal_post <- posterior_epred(
  bscm_ssolaz_detail,
  newdata = appraisal_seq_plot,
  re_formula = NA,
  resp = "Appraisal"
)

epred_resources_post <- posterior_epred(
  bscm_ssolaz_detail,
  newdata = resources_seq_plot,
  re_formula = NA,
  resp = "Appraisal"
)

epred_strain_post <- posterior_epred(
  bscm_ssolaz_detail,
  newdata = strain_seq_plot,
  re_formula = NA,
  resp = "Strain"
)

epred_outcome_post <- posterior_epred(
  bscm_ssolaz_detail,
  newdata = outcome_seq_plot,
  re_formula = NA,
  resp = "Outcome"
)

# prior predictions computed manually from prior draws
prior_appraisal_mat <- with(
  draws_detail,
  outer(prior_b_Appraisal_Intercept, appraisal_seq_plot$Intercept) +
    outer(prior_b_Appraisal_Stressor2, appraisal_seq_plot$Stressor2) +
    outer(prior_b_Appraisal_Resources, appraisal_seq_plot$Resources) +
    outer(prior_b_Appraisal_Resources2, appraisal_seq_plot$Resources2)
)

prior_resources_mat <- with(
  draws_detail,
  outer(prior_b_Appraisal_Intercept, resources_seq_plot$Intercept) +
    outer(prior_b_Appraisal_Stressor2, resources_seq_plot$Stressor2) +
    outer(prior_b_Appraisal_Resources, resources_seq_plot$Resources) +
    outer(prior_b_Appraisal_Resources2, resources_seq_plot$Resources2)
)

prior_strain_mat <- with(
  draws_detail,
  outer(prior_b_Strain_Intercept, strain_seq_plot$Intercept) +
    outer(prior_b_Strain_Appraisal, strain_seq_plot$Appraisal)
)

prior_outcome_mat <- with(
  draws_detail,
  outer(prior_b_Outcome_Intercept, outcome_seq_plot$Intercept) +
    outer(prior_b_Outcome_Strain, outcome_seq_plot$Strain) +
    outer(prior_b_Outcome_Strain2, outcome_seq_plot$Strain2)
)

# helper
summarise_curve <- function(grid, mat, type_label) {
  grid %>%
    mutate(
      estimate = apply(mat, 2, median, na.rm = TRUE),
      lower = apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      upper = apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE),
      type = type_label
    )
}

curve_appraisal_prior <- summarise_curve(appraisal_seq_plot, prior_appraisal_mat, "Prior")
curve_appraisal_post  <- summarise_curve(appraisal_seq_plot, epred_appraisal_post, "Posterior")

curve_resources_prior <- summarise_curve(resources_seq_plot, prior_resources_mat, "Prior")
curve_resources_post  <- summarise_curve(resources_seq_plot, epred_resources_post, "Posterior")

curve_strain_prior <- summarise_curve(strain_seq_plot, prior_strain_mat, "Prior")
curve_strain_post  <- summarise_curve(strain_seq_plot, epred_strain_post, "Posterior")

curve_outcome_prior <- summarise_curve(outcome_seq_plot, prior_outcome_mat, "Prior")
curve_outcome_post  <- summarise_curve(outcome_seq_plot, epred_outcome_post, "Posterior")

# Create legend data for prior/posterior
legend_data <- tibble(x = NA, y = NA, type = c("Prior", "Posterior"))

# 1 Stressor -> Appraisal
p1 <- ggplot() +
  geom_point(data = sso_sim, aes(x = Stressor, y = Appraisal, colour = Resources), alpha = 0.6, size = 2) +
  geom_ribbon(data = curve_appraisal_prior, aes(x = Stressor, ymin = lower, ymax = upper), alpha = 0.18, fill = "blue") +
  geom_line(data = curve_appraisal_prior, aes(x = Stressor, y = estimate, linetype = "Prior"), colour = "darkblue", linewidth = 1) +
  geom_ribbon(data = curve_appraisal_post, aes(x = Stressor, ymin = lower, ymax = upper), alpha = 0.18, fill = "red") +
  geom_line(data = curve_appraisal_post, aes(x = Stressor, y = estimate, linetype = "Posterior"), colour = "darkred", linewidth = 1) +
  labs(x = "Stressor", y = "Appraisal", title = "Prior vs.Posterior: Stressor -> Appraisal", colour = "Resources", linetype = NULL) +
  theme_minimal() + scale_colour_gradient(low = "cornflowerblue", high = "black") + scale_linetype_manual(values = c("Prior" = "solid", "Posterior" = "solid"))

# 2 Resources -> Appraisal
p2 <- ggplot() +
  geom_point(data = sso_sim, aes(x = Resources, y = Appraisal, colour = Stressor), alpha = 0.6, size = 2) +
  geom_ribbon(data = curve_resources_prior, aes(x = Resources, ymin = lower, ymax = upper), alpha = 0.18, fill = "blue") +
  geom_line(data = curve_resources_prior, aes(x = Resources, y = estimate, linetype = "Prior"), colour = "darkblue", linewidth = 1) +
  geom_ribbon(data = curve_resources_post, aes(x = Resources, ymin = lower, ymax = upper), alpha = 0.18, fill = "red") +
  geom_line(data = curve_resources_post, aes(x = Resources, y = estimate, linetype = "Posterior"), colour = "darkred", linewidth = 1) +
  labs(x = "Resources", y = "Appraisal", title = "Prior vs.Posterior: Resources -> Appraisal", colour = "Stressor", linetype = NULL) +
  theme_minimal() + scale_colour_gradient(low = "cornflowerblue", high = "black") + scale_linetype_manual(values = c("Prior" = "solid", "Posterior" = "solid"))

# 3 Appraisal -> Strain
p3 <- ggplot() +
  geom_point(data = sso_sim, aes(x = Appraisal, y = Strain, colour = Stressor), alpha = 0.6, size = 2) +
  geom_ribbon(data = curve_strain_prior, aes(x = Appraisal, ymin = lower, ymax = upper), alpha = 0.18, fill = "blue") +
  geom_line(data = curve_strain_prior, aes(x = Appraisal, y = estimate, linetype = "Prior"), colour = "darkblue", linewidth = 1) +
  geom_ribbon(data = curve_strain_post, aes(x = Appraisal, ymin = lower, ymax = upper), alpha = 0.18, fill = "red") +
  geom_line(data = curve_strain_post, aes(x = Appraisal, y = estimate, linetype = "Posterior"), colour = "darkred", linewidth = 1) +
  labs(x = "Appraisal", y = "Strain", title = "Prior vs.Posterior: Appraisal -> Strain", colour = "Stressor", linetype = NULL) +
  theme_minimal() + scale_colour_gradient(low = "cornflowerblue", high = "black") + scale_linetype_manual(values = c("Prior" = "solid", "Posterior" = "solid"))

# 4 Strain -> Outcome
p4 <- ggplot() +
  geom_point(data = sso_sim, aes(x = Strain, y = Outcome, colour = Appraisal), alpha = 0.6, size = 2) +
  geom_ribbon(data = curve_outcome_prior, aes(x = Strain, ymin = lower, ymax = upper), alpha = 0.18, fill = "blue") +
  geom_line(data = curve_outcome_prior, aes(x = Strain, y = estimate, linetype = "Prior"), colour = "darkblue", linewidth = 1) +
  geom_ribbon(data = curve_outcome_post, aes(x = Strain, ymin = lower, ymax = upper), alpha = 0.18, fill = "red") +
  geom_line(data = curve_outcome_post, aes(x = Strain, y = estimate, linetype = "Posterior"), colour = "darkred", linewidth = 1) +
  labs(x = "Strain", y = "Outcome", title = "Prior vs.Posterior: Strain -> Outcome", colour = "Appraisal", linetype = NULL) +
  theme_minimal() + scale_colour_gradient(low = "cornflowerblue", high = "black") + scale_linetype_manual(values = c("Prior" = "solid", "Posterior" = "solid"))

# Combine and save
grid_plot <- ssolaz_dgp / (p1 + p2) / (p3 + p4)
ggsave("06_ssolaz_indir.png", grid_plot, width = 8, height = 7, dpi = 300, 
       scale = 1.5)
grid_plot


#### parameter updating: how far did the data move the theory? ####

priors_ssolaz_detail

## ---- code-ssolaz-update-table ----
# Table labels of the structural parameters (order = order of rows in the table)
param_labels_ssolaz <- tribble(
  ~Parameter,                ~Label,
  "b_Appraisal_Intercept",   "$\\alpha_{A}$",
  "b_Appraisal_Stressor2",   "$\\beta_{A,St^{2}}$",
  "b_Appraisal_Resources",   "$\\beta_{A,R}$",
  "b_Appraisal_Resources2",  "$\\beta_{A,R^{2}}$",
  "b_Strain_Intercept",      "$\\alpha_{S}$",
  "b_Strain_Appraisal",      "$\\beta_{S,A}$",
  "b_Outcome_Intercept",     "$\\alpha_{O}$",
  "b_Outcome_Strain",        "$\\beta_{O,S}$",
  "b_Outcome_Strain2",       "$\\beta_{O,S^{2}}$"
)

# Analytic prior specification (exact, avoids MC error in the prior draws),
# read from the fitted model so that the table cannot drift from the priors
# that were actually used
prior_spec_ssolaz <- prior_summary(bscm_ssolaz_detail) %>%
  as_tibble() %>%
  filter(class == "b", coef != "") %>%
  mutate(
    Parameter  = paste0("b_", resp, "_", coef),
    Prior_Mean = as.numeric(str_match(prior, "normal\\(\\s*([-0-9.]+)\\s*,")[, 2]),
    Prior_SD   = as.numeric(str_match(prior, ",\\s*([-0-9.]+)\\s*\\)")[, 2])
  ) %>%
  select(Parameter, Prior_Mean, Prior_SD)

prior_spec_ssolaz <- param_labels_ssolaz %>%
  left_join(prior_spec_ssolaz, by = "Parameter")

# every labelled parameter must have a user-specified normal prior
stopifnot(!anyNA(prior_spec_ssolaz$Prior_Mean), !anyNA(prior_spec_ssolaz$Prior_SD))

posterior_spec_ssolaz <- as_draws_df(bscm_ssolaz_detail) %>%
  as_tibble() %>%
  select(all_of(prior_spec_ssolaz$Parameter)) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>%
  group_by(Parameter) %>%
  summarise(
    Posterior_Mean = mean(value, na.rm = TRUE),
    Posterior_SD   = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

update_table_ssolaz <- prior_spec_ssolaz %>%
  left_join(posterior_spec_ssolaz, by = "Parameter") %>%
  mutate(
    # shift of the posterior mean in units of prior SD (= units of theoretical commitment)
    Shift_in_prior_SD = (Posterior_Mean - Prior_Mean) / Prior_SD
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  select(Label, Prior_Mean, Prior_SD, Posterior_Mean, Posterior_SD, Shift_in_prior_SD) %>%
  set_names(c("Parameter",
              "$M_{prior}$", "$SD_{prior}$",
              "$M_{post}$",  "$SD_{post}$",
              "$z$"))

## ---- update-table-output ----
knitr::kable(
  update_table_ssolaz,
  format    = "latex",
  escape    = FALSE,
  booktabs  = TRUE,
  align     = c("l", rep("r", 5)),
  linesep   = ""
)

#### session info ####

# Record the computational environment (R and package versions) for
# reproducibility; the text file accompanies the analysis code on OSF
sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#### save image ####

save.image("bscms-runningsexample.RData")