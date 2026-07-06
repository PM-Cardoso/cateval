# Package index

## Overall-benefit calibration

Match concordant and discordant patients, then summarise how well the
predicted benefit of following the model matches the observed benefit.

- [`match_benefit_pairs()`](https://PM-Cardoso.github.io/cateval/reference/match_benefit_pairs.md)
  : Match Concordant and Discordant Patients for Overall-Benefit
  Calibration
- [`benefit_deciles()`](https://PM-Cardoso.github.io/cateval/reference/benefit_deciles.md)
  : Group Matched Pairs by Predicted Benefit and Summarise Observed
  Benefit
- [`benefit_by_drug()`](https://PM-Cardoso.github.io/cateval/reference/benefit_by_drug.md)
  : Overall-Benefit Calibration Split by Drug Class
- [`benefit_overall()`](https://PM-Cardoso.github.io/cateval/reference/benefit_overall.md)
  : Mean Observed Benefit of Following the Model
- [`benefit_calibration()`](https://PM-Cardoso.github.io/cateval/reference/benefit_calibration.md)
  : Calibration-in-the-Large and Calibration Slope of Predicted Benefit

## Pairwise drug-class calibration

Compare predicted and observed treatment effects between pairs of drugs.

- [`unified_validation()`](https://PM-Cardoso.github.io/cateval/reference/unified_validation.md)
  : Perform Pairwise Heterogeneous Treatment Effect Calibration Across
  Multiple Drugs
- [`calibration_hte()`](https://PM-Cardoso.github.io/cateval/reference/calibration_hte.md)
  : Estimate Heterogeneous Treatment Effects Across Calibration Groups

## Predicted-best drug utilities

Identify the recommended drug (or the set of near-best drugs) for each
patient and visualise the resulting recommendations.

- [`get_best_drugs()`](https://PM-Cardoso.github.io/cateval/reference/get_best_drugs.md)
  : Assign Best Drugs Based on Prediction Rankings or Tolerance
- [`get_ranked_or_tolerant_drugs()`](https://PM-Cardoso.github.io/cateval/reference/get_ranked_or_tolerant_drugs.md)
  : Get Ranked or Tolerance-Based Drug Recommendation
- [`optimal_drug_comparison_plot()`](https://PM-Cardoso.github.io/cateval/reference/optimal_drug_comparison_plot.md)
  : Plotting optimal drug combinations
