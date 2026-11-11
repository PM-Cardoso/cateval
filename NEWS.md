# cateval 2.5.0

## Example data and a worked example

* New dataset `cateval_example`: 3000 simulated patients on three treatments
  (`A`, `B`, `C`) with a continuous outcome and a prediction model's expected
  outcome on each treatment. The model is built to rank treatments well but
  overstate their differences, so a validation has a known answer (calibration
  slope about 0.8). The simulation code is in `data-raw/`.
* New article, "Worked example", on the website. It runs every function on
  `cateval_example`, with short plotting code and notes on how to read each
  result.
* Every help page now has examples that run on `cateval_example`, so the
  reference pages show real output. None are wrapped in `\dontrun{}` any more.
* Help pages no longer use diabetes-specific wording or examples; the package
  applies to any continuous outcome.

## `plot_optimal_drugs()`

* New arguments `title_count`, `title_sets` and `subtitle`. The default titles
  are neutral; the previous ones referred to HbA1c and mmol/mol.
* Clearer layout: percentage labels on both charts, readable axis labels, and
  the second chart split into the groups given in `groups`. The inputs and the
  tables returned with `plot = FALSE` are unchanged.

## Bug fix

* Matching failed if the data had a column called `id`, `weights` or
  `distance`, because `MatchIt::get_matches()` adds columns with those names.
  cateval now gives MatchIt's columns a `cateval_` prefix when the name is
  taken. A data column called `subclass` gets a clear error, since cateval uses
  that name to link the patients of each matched set.

## Website

* New "Worked example" tab, a shorter home page with a calibration figure,
  reorganised reference index (with an "Example data" section), and a cleaner
  theme.

# cateval 2.4.0

## Pairwise matching is always anti-exact on the received drug (behaviour change)

* In `pairwise_calibration()` / `pairwise_calibration_pair()` with
  `matching = TRUE`, patients predicted to do better on the first drug are still
  matched to patients predicted to do better on the second, but the match is now
  always anti-exact on the received drug, so every matched set compares patients
  who took different drugs. This is the same requirement as in
  `match_benefit_pairs()`, and it is protected in the same way: it is re-applied
  after `match_args`, so `match_args = list(antiexact = NULL)` cannot remove it
  and an `antiexact` entry adds to it. **Matched pairwise results therefore
  differ from 2.3.0.** Unmatched results (`matching = FALSE`) are unchanged.
* Matching arguments (`matching_var`, `match.exact`, `match.antiexact`,
  `caliper`, `match_args`) can now refer to any column, including `drug_var`,
  by its original name. The working columns are added as copies instead of
  renaming the originals; the matched cohort returned with
  `return_matched = TRUE` keeps the original columns as well as
  `dataset_drug_var`, `dataset_outcome_var` and `dataset_benefit`.
* A drug pair for which no matches can be formed (e.g. a very tight caliper) now
  gives a warning and no rows for that pair, instead of stopping the whole
  `pairwise_calibration()` run.

## Prediction columns and drug names of any kind

The package makes no assumption about what the treatments are called (drug
classes, molecules, trial arms, numeric codes, ...). This release makes that
hold throughout.

* `match_benefit_pairs()` no longer derives drug names by stripping the longest
  common prefix of `pred_cols`. That prefix could run into the drug names: for
  `c("pred.SGLT2", "pred.SU")` it was `"pred.S"`, giving drugs `GLT2` and `U`, no
  concordant patients and an opaque MatchIt error. Drug names are now matched
  against the values observed in `drug_var`:
  * any prefix (`"pred.A"`, `"yhat_outcome_A"`) or suffix (`"A_pred"`); where
    several names fit, the longest is taken (`"pred.BA"` is `BA`, not `A`);
  * failing that, ignoring upper/lower case, keeping the spelling in `drug_var`
    (not used when two observed names differ by case alone);
  * a drug that has a prediction column but that nobody received is read from
    the naming pattern shared by the other columns.
  The last two routes print the mapping they used. Otherwise the new `drugs`
  argument, or a named `pred_cols` vector, sets the mapping explicitly, and the
  mapping used is always returned in the new `drugs` element of the result. For
  the usual `pred.<drug>` columns the matched pairs are identical to 2.3.0.
* Numeric drug codes now work in `pairwise_calibration()`. Previously the code
  was used as a position when looking up the prediction column (drug `10` looked
  for the 10th column), giving "subscript out of bounds" or the wrong column.
* `drugs` arguments accept character, factor or numeric values.
* Patients whose received drug has no prediction column (or is missing) are now
  excluded with a warning; previously they were matched and produced pairs with
  a missing `pred.benefit`.
* Drug names containing spaces, hyphens, dots, brackets or other
  regular-expression characters (e.g. `"GLP-1 RA"`) now work with
  `conc_tolerance`; set membership was previously tested with a word-boundary
  regular expression.
* `get_best_drugs()` and `get_ranked_or_tolerant_drugs()` gain a `drugs`
  argument (drug name per column, by position), and remove `final_var_name` /
  `prediction_column` from column names literally rather than as a regular
  expression (so the `.` in `"pred."` no longer matches any character).

## Small calibration groups

* `benefit_by_group()` and `benefit_by_drug()`: a group with fewer than two
  pairs now gets `NA` confidence bounds (analytical and bootstrap) with one
  clear warning per `cal_groups` value, instead of `NaN` bounds and a warning
  from `qt()`. Groups that receive no pairs are now reported in the same warning
  (they were silently missing). `benefit_by_drug()` prefixes these warnings with
  the drug class.
* `n` in `benefit_by_group()` now counts pairs with an observed benefit, which is
  the sample size the mean and its interval are based on (consistent with
  `benefit_overall()`). It only differs from before when `obs.benefit` has
  missing values.
* `benefit_overall()` (fewer than two pairs) and `benefit_citl_slope()` (fewer
  than three complete pairs, or constant `pred.benefit`) return `NA` intervals
  with a warning instead of `NaN`.
* `pairwise_calibration_pair()`: a group whose outcome model has no residual
  degrees of freedom keeps its point estimate but gets `NA` Wald bounds, with a
  warning naming the group.

## Documentation

* `match_benefit_pairs()`: `pred.benefit` is computed from the *index* (treated)
  patient's predictions — the concordant patient in the concordant-as-treated
  direction and the discordant patient in the discordant-as-treated direction.
  The documentation previously said it always used the concordant patient's
  predictions; the computation itself is unchanged.
* `match_benefit_pairs()`: new section on `ratio > 1`. Each treated patient is
  paired with each of its controls, so it contributes up to `ratio` rows, all
  weighted equally by the summary functions. The internal comment that each
  matched set held exactly one concordant and one discordant patient was wrong
  for `ratio > 1` and has been corrected.
* README quick start: `plot_optimal_drugs()` used a tolerance column from a
  `get_best_drugs(rank = 1)` call; the call now uses `tolerance = 3`.

## Input checks

* `drug_var` and `outcome_var` are rejected as `matching_var`, `match.exact` or
  `adjustment_var` (and `outcome_var` as `match.antiexact`) with a clear message.
* A variable that is both exact and anti-exact after all overrides is rejected
  before calling MatchIt.
* Non-pair matching methods (e.g. `"full"`, `"subclass"`, `"exact"`) fail with a
  message explaining that each treated patient must be linked to its own
  control(s).
* Missing values in `match.exact` / `match.antiexact` variables now exclude the
  patient before matching, as missing `matching_var` values already did.

# cateval 2.3.0

## Multiple `cal_groups` in `benefit_by_group()`

* `benefit_by_group()` now accepts a vector of `cal_groups` (e.g. `c(3, 5, 10)`),
  as `benefit_by_drug()` and `pairwise_calibration()` already did. Each value is
  summarised in turn on the same matched pairs and the results are stacked.
* The output gains an `n_groups` column (first column) recording which
  `cal_groups` value produced each row, matching `benefit_by_drug()`. This is the
  one visible change for existing code: a scalar `cal_groups` returns the same
  rows and values as before, plus this column.
* Internally the single-setting summary moved to an unexported worker
  (`.benefit_by_group_one()`); `benefit_by_drug()` is unchanged and continues
  to call `benefit_by_group()` once per setting and drug.

# cateval 2.2.0

All changes below are additive: every new argument defaults to the behaviour of
the previous release, so existing code produces the same results.

## Per-direction matching in `match_benefit_pairs()`

* `match_args_conc` and `match_args_disc` apply `MatchIt::matchit()` arguments to
  a single matching direction. `match_args` continues to apply to both, and the
  direction-specific lists are layered on top of it, so the order of precedence
  is: package defaults, then the named arguments (`method`, `distance`,
  `replace`, `ratio`, `caliper`, `std.caliper`), then `match_args`, then that
  direction's own list.
* Because `utils::modifyList()` treats `NULL` as removal, a direction list can
  also unset an argument for one direction only, e.g.
  `match_args_disc = list(caliper = NULL)`.
* `match_args_disc` is ignored with a warning when `symmetric = FALSE`, since
  the discordant-as-treated direction is not run.
* All three lists must be fully named; an unnamed element is now rejected with a
  clear message rather than failing inside `matchit()`.

## The required matching constraints can no longer be dropped

* Matching is always exact on the predicted-best drug and always anti-exact on
  the received-drug variable (and on the tolerance-set columns when
  `conc_tolerance` is used). These are re-applied after any override, so passing
  `exact` or `antiexact` through `match_args` now *adds* to them instead of
  replacing them. Previously such an override silently removed the constraint
  that makes a matched pair interpretable.

## Maintenance

* `pairwise_calibration()` no longer creates the intermediate `drug_1_pred` and
  `drug_2_pred` columns, which were used once to compute `benefit` and then
  carried unused through the matched cohort returned by `return_matched = TRUE`.
* Removed two vestigial `requireNamespace()` guards in
  `pairwise_calibration_pair()` for packages that are already hard dependencies.
* `globalVariables()` gained `grouping` and `n_groups`, which are used as bare
  column names but were undeclared, and lost `drug_1_pred` and `drug_2_pred`.
* Corrected the comment on `.build_match_formula()`, which claimed that all
  exact- and anti-exact-matched variables are dropped from the distance formula;
  categorical ones are in fact retained, which is harmless because they are
  constant within an exact stratum.

# cateval 2.1.1

Documentation-only fixes; no change to any function's behaviour or output.

* `pairwise_calibration()`: the `@return` documentation listed `cal_groups` and
  `grouping` columns that the function does not produce. It now documents the
  columns actually returned (`mean`, `coef`, `coef_low`, `coef_high`, `n_groups`,
  `drug1`, `n_drug1`, `drug2`, `n_drug2`). Note that the output has no per-group
  identifier column; rows are ordered by group within each `n_groups` block.
* `benefit_citl_slope()`: the `@return` documentation omitted the `metric` column
  that identifies each row (`calibration_in_the_large` / `calibration_slope`); it
  is now documented alongside `estimate`, `lower.ci` and `upper.ci`.

# cateval 2.1.0

All changes below are additive: every new argument defaults to the behaviour of
the previous release, so existing code produces the same results.

## Configurable matching

* `match_benefit_pairs()`, `pairwise_calibration()` and
  `pairwise_calibration_pair()` now expose the matching options that were
  previously fixed: `method`, `distance`, `replace` and `ratio`, together with a
  `match_args` list that is applied last and can supply any other
  `MatchIt::matchit()` argument.
* `caliper` and `std.caliper` are now available in the pairwise functions, which
  previously had no caliper support.
* `seed` is now available in the pairwise functions, so a matched cohort can be
  reproduced exactly.

## Matching is performed once per drug pair

* In the pairwise functions, matching previously happened inside the
  `cal_groups` loop and was therefore repeated for each requested number of
  groups. It is now performed once per drug pair and the same matched cohort is
  reused for every value of `cal_groups`, so the different groupings are directly
  comparable. Unseeded results may differ slightly from earlier versions.

## Inspecting the matched cohort

* `return_matched = TRUE` makes the pairwise functions return a list with
  `results`, `matched` and `matchit`, the latter two holding one entry per drug
  pair. This allows the matching to be checked (for example with a balance plot),
  as was already possible for the overall-benefit workflow.

## Bootstrap confidence intervals for the pairwise analysis

* `bootstrap`, `n_boot`, `boot_seed` and `conf` are new in the pairwise
  functions. With `bootstrap = TRUE` the intervals are percentile bootstrap
  intervals instead of the model's Wald intervals; the point estimates are
  unchanged. Calibration groups are fixed before resampling, and matched sets are
  resampled as a unit.
* `conf` also controls the level of the Wald intervals (default 0.95, as before).

# cateval 2.0.0

## Breaking: functions renamed for clarity and consistency

Several exported functions were renamed so the names describe what they do and
follow a consistent scheme (`benefit_*` for the overall-benefit family,
`pairwise_*` for the drug-pair family, `plot_*` for plotting). **Update any code
that calls the old names.**

| Old name | New name | Why |
|----------|----------|-----|
| `benefit_deciles()` | `benefit_by_group()` | It splits into `cal_groups` bins (any number), not only deciles; now parallels `benefit_by_drug()`. |
| `benefit_calibration()` | `benefit_citl_slope()` | Names the two metrics it returns (calibration-in-the-large and slope) and removes the clash with the calibration-curve data. |
| `unified_validation()` | `pairwise_calibration()` | Describes what it does — pairwise benefit calibration across all drug pairs. |
| `calibration_hte()` | `pairwise_calibration_pair()` | Reads as the single-pair engine behind `pairwise_calibration()`. |
| `optimal_drug_comparison_plot()` | `plot_optimal_drugs()` | Clearer, and groups plotting functions under a `plot_` prefix. |

No arguments, return values or behaviour changed — only the function names.

# cateval 1.6.0

## `benefit_by_drug()` supports multiple calibration groupings

* `benefit_by_drug()` now accepts a vector of `cal_groups` values (e.g.
  `c(3, 5, 10)`) as well as a single value. Each requested number of groups is
  computed for every drug class and the results are stacked into one data frame.
* The output gains an `n_groups` column recording which `cal_groups` setting
  produced each row, matching the convention already used by
  `calibration_hte()` and `unified_validation()`. The identifying columns
  (`index_drug`, `n_groups`) now lead the returned data frame.
* Existing single-value calls are unchanged apart from the added `n_groups`
  column.

# cateval 1.5.0

## Overall-benefit calibration reworked

* `match_benefit_pairs()` replaces `compute_overall_benefit()` as the entry point.
  It labels each patient concordant or discordant, runs symmetric 1:1 matching
  (concordant-as-treated and discordant-as-treated), and returns the
  `concordant`, `discordant` and de-duplicated `combined` matched populations
  together with the underlying `matchit` objects.
* New summary functions take a matched-pairs frame and return the calibration
  numbers: `benefit_deciles()` (mean observed vs predicted benefit by group),
  `benefit_by_drug()` (the same, split by the index/treated patient's drug),
  `benefit_overall()` (mean observed benefit of following the model) and
  `benefit_calibration()` (calibration-in-the-large and calibration slope).
* Each summary offers analytical confidence intervals by default or an optional
  pair-level percentile bootstrap (`bootstrap = TRUE`). Note that the bootstrap
  resamples matched pairs as independent, so it does not by itself account for
  the dependence introduced by matching with replacement.
* Predicted benefit is anchored on the index (treated) patient of each pair, and
  observed benefit is computed on the absolute outcome scale.

## Removed

* `compute_overall_benefit()`, `compute_overall_benefit_performance()` and
  `plot_benefit_calibration()` have been removed. Plotting is left to the caller;
  the package returns the underlying numbers only.

## Other

* Added a test suite and continuous integration (`R-CMD-check`).
