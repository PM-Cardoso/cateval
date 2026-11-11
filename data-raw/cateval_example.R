# Simulated example data for cateval
#
# Builds `cateval_example`: 3000 patients who each received one of three
# treatments (A, B or C), with a continuous outcome where lower is better, and
# a prediction model's expected outcome on each treatment.
#
# The data are designed so that:
#   * treatment effects vary between patients (B suits patients with a high
#     `marker`; C suits younger patients), so there is real heterogeneity to
#     validate;
#   * treatment choice depends on patient characteristics, as in routine care,
#     so matching on covariates matters;
#   * the prediction model gets the direction of the effects right but
#     overstates their size (its predicted differences are 1.25 times the true
#     ones, plus noise), so the calibration slope is below 1.
#
# Run from the package root to regenerate data/cateval_example.rda:
#   source("data-raw/cateval_example.R")

set.seed(20261111)
n <- 3000

# ---- Patient characteristics ------------------------------------------------
age      <- round(pmin(pmax(rnorm(n, 60, 12), 25), 90))
sex      <- sample(c("Female", "Male"), n, replace = TRUE)
baseline <- round(rnorm(n, 60, 10), 1)
marker   <- round(pmin(pmax(rnorm(n, 5, 2), 0), 10), 2)

# ---- True expected outcome on each treatment --------------------------------
# A is the reference. tau_B and tau_C are the true differences from A
# (negative = better than A).
mu_A  <- 20 + 0.6 * baseline + 0.1 * (age - 60) + 1 * (sex == "Male")
tau_B <- 0.5 - 1.2 * (marker - 5)     # B is better when the marker is high
tau_C <- 0.25 + 0.15 * (age - 60)     # C is better in younger patients
mu <- cbind(A = mu_A, B = mu_A + tau_B, C = mu_A + tau_C)

# ---- Treatment received -----------------------------------------------------
# Clinicians lean towards the right choice, but only partly: B is more likely
# with a high marker, C more likely in younger patients.
lin <- cbind(
  A = 0,
  B = 0.30 * (marker - 5),
  C = -0.04 * (age - 60)
)
prob <- exp(lin) / rowSums(exp(lin))
treatment <- apply(prob, 1, function(p) sample(colnames(prob), 1, prob = p))

# ---- Observed outcome -------------------------------------------------------
outcome <- round(mu[cbind(seq_len(n), match(treatment, colnames(mu)))] + rnorm(n, 0, 5), 1)

# ---- Model predictions ------------------------------------------------------
# A shared prognostic error (the model's baseline is not perfect), and
# treatment differences that are 25% too large plus some noise.
prog_err <- rnorm(n, 0, 1.5)
pred_A <- mu_A + prog_err
pred_B <- mu_A + prog_err + 1.25 * tau_B + rnorm(n, 0, 1)
pred_C <- mu_A + prog_err + 1.25 * tau_C + rnorm(n, 0, 1)

cateval_example <- data.frame(
  id        = seq_len(n),
  treatment = treatment,
  age       = age,
  sex       = sex,
  baseline  = baseline,
  marker    = marker,
  outcome   = outcome,
  pred_A    = round(pred_A, 2),
  pred_B    = round(pred_B, 2),
  pred_C    = round(pred_C, 2),
  stringsAsFactors = FALSE
)

if (sys.nframe() == 0L || identical(Sys.getenv("CATEVAL_SAVE"), "true")) {
  usethis::use_data(cateval_example, overwrite = TRUE)
}
