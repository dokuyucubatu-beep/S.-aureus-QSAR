# Metric definitions and prespecified diagnostics

Let y be observed pMIC, p the prediction and n the evaluated sample size.

- **Predictive test R²** = 1 - sum((y-p)²)/sum((y-mean(y))²). It can be negative.
  Pearson correlation squared is reported separately and cannot replace this.
- **RMSE** = sqrt(mean((y-p)²)); **MAE** = mean(abs(y-p)). Units are pMIC.
- **5-fold Q²** uses all pooled out-of-fold training predictions and the overall
  training mean in the denominator. Every molecule is predicted by a model that
  did not fit its scaffold group. It is not the mean of fold R² scores.
- **10-fold R² mean ± SD** uses the ten separate validation-fold predictive R²
  values; SD has ddof=1. Fold sizes and target variances differ, so this differs
  from pooled Q². These are fixed-protocol CV estimates, not nested tuning results.
- **Bootstrap 95% CI:** sample the 37 held-out scaffold groups with replacement,
  retaining all compounds of each sampled group, 1,000 times. Recalculate predictive
  R² for the fixed test predictions. Report undefined resamples explicitly. This
  interval covers test-set sampling uncertainty conditionally on the fitted model;
  it does not include model-refitting or training-data uncertainty.

## Golbraikh–Tropsha diagnostics

The direction convention is explicit:

```
k  = sum(y*p)/sum(p*p)       # observed ~ k * predicted, through origin
k' = sum(y*p)/sum(y*y)       # predicted ~ k' * observed, through origin
r0²  = 1 - sum((y-k*p)²)  / sum((y-mean(y))²)
r0'² = 1 - sum((p-k'*y)²) / sum((p-mean(p))²)
r² = Pearson(y,p)²
```

The prespecified diagnostic requires positive r, r²>0.6, |r0²-r0'²|<0.3, and
either the k direction or k' direction satisfying 0.85<slope<1.15 and
abs((r²-r0²)/r²)<0.1 in the corresponding direction. The absolute difference
rule is a conservative explicitly stated implementation. Q²>0.5 is reported as
a separate conventional robustness flag. These historical criteria have known
statistical limitations and are **not OECD acceptance thresholds**. See Alexander
et al. (2015); predictive R², RMSE and MAE remain primary error measures.

## Y-scrambling

One hundred independent permutations shuffle y within the training set globally.
Each repeats the complete selector/imputer/scaler/RF fit and the same five scaffold
folds. No test labels are used. The Monte Carlo p-value for actual CV Q² is
`(1 + count(Q²_scrambled >= Q²_actual)) / (100 + 1)`.

Corrected `cR²p = R * sqrt(R² - mean(Rr)²)`, where R is the actual training
Pearson correlation and Rr are training correlations of the shuffled-label fits.
This is distinct from actual/scrambled Q²; both are reported. A negative radicand
is undefined, not forced to zero. cR²p>0.5 is a conventional diagnostic only;
flexible models may fit randomized labels and therefore are penalized by it.

## Applicability domain

On selected descriptors standardized using training means/SDs, augment the design
matrix with an intercept. For a query x, `h = x' pinv(X'X) x`, pseudoinverse
relative cutoff 1e-12; `h* = 3*(p+1)/n_train`. Also compute the maximum ECFP4
(radius 2, 2048 bits) Tanimoto similarity to any training molecule. Its cutoff is
the fifth percentile of training leave-self-out nearest-neighbor similarities.

A query is in the defined AD only if **both** thresholds pass. Neither threshold
uses test labels or residuals. For a nonlinear RF, leverage is a descriptor-space
extrapolation diagnostic, not a calibrated uncertainty interval. The Williams-style
plot uses residuals divided by the training out-of-fold RMSE; its ±3 lines are
diagnostic and do not define membership or discard test errors.

The molecule:descriptor ratio is reported for all 542 compounds and separately
for the actual fitting set. It is not a parameter-count guarantee for a forest.
