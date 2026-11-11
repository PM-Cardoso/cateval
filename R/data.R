#' Simulated example data: three treatments and a continuous outcome
#'
#' A simulated cohort of 3000 patients, each treated with one of three
#' treatments (`A`, `B` or `C`). The outcome is a continuous score measured
#' after treatment, where lower is better. For every patient a prediction model
#' gives the expected outcome on each of the three treatments, which is what
#' cateval validates.
#'
#' The data are built to have a known answer:
#' * Treatment effects vary between patients. `B` works better than `A` when
#'   `marker` is high, and `C` works better than `A` in younger patients.
#' * Treatment choice depends on `age` and `marker`, as it would in routine
#'   care. Patients on different treatments therefore differ at the start, and
#'   a fair comparison needs matching or adjustment.
#' * The prediction model ranks the treatments well but overstates the size of
#'   the differences between them: its predicted differences are 1.25 times the
#'   true ones, plus noise. A well-run validation should therefore find a
#'   calibration slope below 1.
#'
#' Nothing in the data is specific to a disease or a drug. The same layout
#' (one row per patient, the treatment received, the observed outcome and one
#' prediction column per treatment) is all cateval needs.
#'
#' @format A data frame with 3000 rows and 10 columns:
#' \describe{
#'   \item{id}{Patient identifier.}
#'   \item{treatment}{Treatment received: `"A"`, `"B"` or `"C"`.}
#'   \item{age}{Age in years.}
#'   \item{sex}{`"Female"` or `"Male"`.}
#'   \item{baseline}{Score before treatment.}
#'   \item{marker}{A biomarker, from 0 to 10.}
#'   \item{outcome}{Score after treatment. Lower is better.}
#'   \item{pred_A, pred_B, pred_C}{The model's predicted outcome on treatment
#'     `A`, `B` and `C`.}
#' }
#'
#' @source Simulated. The code is in `data-raw/cateval_example.R` in the
#'   package's GitHub repository.
#'
#' @examples
#' head(cateval_example)
#' table(cateval_example$treatment)
"cateval_example"
