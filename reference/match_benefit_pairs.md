# Match Concordant and Discordant Patients for Overall-Benefit Calibration

Builds the matched patient pairs used to validate a treatment-selection
model's predicted *overall benefit*. Each patient is labelled concordant
(received the predicted-best drug) or discordant (received a different
drug), and concordant patients are matched to clinically similar
discordant patients by nearest-neighbour Mahalanobis distance, exactly
on the predicted-best drug and forced to differ on the drug actually
received.

## Usage

``` r
match_benefit_pairs(
  data,
  drug_var,
  outcome_var,
  pred_cols,
  matching_var,
  match.exact = NULL,
  match.antiexact = NULL,
  caliper = NULL,
  conc_tolerance = NULL,
  symmetric = TRUE,
  seed = NULL
)
```

## Arguments

- data:

  A data frame with one row per patient, containing the received drug,
  the observed outcome, the predicted outcome for each drug and the
  matching covariates.

- drug_var:

  Character. Name of the column holding the drug actually received.

- outcome_var:

  Character. Name of the observed outcome column (e.g. the 12-month
  HbA1c). Benefits are computed on this absolute scale.

- pred_cols:

  Character vector. Names of the predicted-outcome columns, one per
  drug, sharing a common prefix (e.g. `c("pred.SGLT2", "pred.GLP1")`).

- matching_var:

  Character vector. Covariates used for Mahalanobis distance matching.

- match.exact:

  Optional character vector. Additional variables to match on exactly.
  The predicted-best drug is always included.

- match.antiexact:

  Optional character vector. Additional variables on which a matched
  pair must differ. The received-drug variable is always included.

- caliper:

  Optional named numeric vector passed to
  [`MatchIt::matchit`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
  (e.g. `c(prehba1c = 0.4)`); interpreted in standard-deviation units.

- conc_tolerance:

  Optional numeric. If supplied, a patient is concordant when they
  received any drug within this tolerance of the best predicted outcome,
  and a discordant match must differ from every drug in that set.

- symmetric:

  Logical. If `TRUE` (default), match in both directions and return the
  pooled, de-duplicated set in `combined`. If `FALSE`, only the
  concordant-as-treated direction is run.

- seed:

  Optional integer. Seed set before each matching call so the result is
  reproducible.

## Value

An object of class `cateval_benefit_match`: a list with elements

- concordant:

  Matched pairs, concordant patients as treated.

- discordant:

  Matched pairs, discordant patients as treated (`NULL` if
  `symmetric = FALSE`).

- combined:

  The concordant and discordant sets pooled, with duplicate pairs
  removed.

- matchit:

  A list of the underlying `matchit` objects (`$concordant` and
  `$discordant`) for balance / love plots.

Each pairs data frame has one row per matched pair, including
`concordant_drug`, `discordant_drug`, `index_drug` (the drug received by
the index/treated patient of the pair), `pred.benefit`, `obs.benefit`
and the retained matching covariates.

## Details

Matching is run in both directions by default (`symmetric = TRUE`): once
with concordant patients as the treated group and once with the roles
reversed. This returns three matched populations that can each be
inspected or summarised separately:

- `concordant` - concordant patients matched to discordant controls;

- `discordant` - discordant patients matched to concordant controls;

- `combined` - the two sets pooled, with pairs that appear in both
  directions counted once.

The predicted and observed benefits are computed on the scale of
`outcome_var` (i.e. absolute outcomes): for each pair, `pred.benefit` is
the concordant patient's predicted outcome on the discordant drug minus
the predicted outcome on their own (recommended) drug, and `obs.benefit`
is the observed outcome of the discordant patient minus that of the
concordant patient. Both are positive when the model's recommended drug
is expected to, and does, give the better outcome.

The matched pairs frames are the input to
[`benefit_deciles`](https://PM-Cardoso.github.io/cateval/reference/benefit_deciles.md),
[`benefit_overall`](https://PM-Cardoso.github.io/cateval/reference/benefit_overall.md),
[`benefit_calibration`](https://PM-Cardoso.github.io/cateval/reference/benefit_calibration.md)
and
[`benefit_by_drug`](https://PM-Cardoso.github.io/cateval/reference/benefit_by_drug.md).

**Limitation (confidence intervals under matching with replacement).**
Matching is done with replacement, so a single control patient can
appear in several matched pairs, which makes pairs that share a control
not independent. The confidence intervals from the summary functions –
both the analytical interval and the optional pair-level bootstrap –
treat matched pairs as independent, and so they do not widen to reflect
this dependence. They are therefore likely to be somewhat too narrow
(anti-conservative) when control reuse is substantial. Valid variance
estimation for matching with replacement is a known hard problem; the
intervals here should be read as approximate, and this limitation stated
when they are reported. The extent of reuse can be inspected directly,
e.g.
`table(table(c(pairs$cateval_row_id_conc, pairs$cateval_row_id_disc)))`.

## See also

[`benefit_deciles`](https://PM-Cardoso.github.io/cateval/reference/benefit_deciles.md),
[`benefit_overall`](https://PM-Cardoso.github.io/cateval/reference/benefit_overall.md),
[`benefit_calibration`](https://PM-Cardoso.github.io/cateval/reference/benefit_calibration.md),
[`benefit_by_drug`](https://PM-Cardoso.github.io/cateval/reference/benefit_by_drug.md)

## Examples

``` r
if (FALSE) { # \dontrun{
matched <- match_benefit_pairs(
  data = analysis_dataset,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.", c("SGLT2", "GLP1", "DPP4", "SU", "TZD")),
  matching_var = c("agetx", "prebmi", "prehba1c", "sex"),
  match.exact = c("sex"),
  caliper = c(prehba1c = 0.4),
  seed = 19840503
)

# inspect balance
plot(summary(matched$matchit$concordant))
# summarise the pooled set
benefit_deciles(matched$combined, cal_groups = 10)
} # }
```
