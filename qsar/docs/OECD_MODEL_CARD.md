# OECD-aligned model documentation and unresolved limitations

This model is a research baseline. **OECD approval, regulatory acceptance, or full
scientific validation is not claimed.** The OECD 2004/2007 principles and the 2024
QAF second edition guide this evidence mapping. Numerical diagnostics alone do
not establish conformity or fitness for a specific regulatory decision.

| OECD principle | Evidence supplied | Limitations/status |
|---|---|---|
| 1. Defined endpoint | pMIC, molar conversion, ATCC 29213, explicit 24 h microdilution, individual ChEMBL assays, exact-value/repeat policies | Partially addressed: full-paper medium, inoculum, temperature and readout harmonization not established; censored compounds excluded |
| 2. Unambiguous algorithm | Pinned dependencies, raw/curated snapshots, fixed scaffold IDs, descriptor list/selection logs, config, source code, serialized model and hashes | Implemented and reproducible in the recorded environment; neither GNN nor mechanism model |
| 3. Defined applicability domain | Training-only leverage + ECFP4 threshold; per-test h/similarity/AD flags; coverage and within-AD metrics | Defined empirically; accuracy inside AD is not guaranteed; prospective calibration absent |
| 4. Fit, robustness and predictivity | Test errors, grouped 5/10-fold CV, scaffold bootstrap CI, Y permutations, GT diagnostics, baseline comparisons | All values and failures retained; scaffold holdout is not independent-study/laboratory validation; no test-based optimization |
| 5. Mechanistic interpretation, if possible | Global SHAP, exact feature definitions, training background, held-out explanations and additivity check | Exploratory descriptor associations only; no causal antibacterial mechanism, target specificity, or proof of DHFR inhibition |

## Intended use and domain

Prioritization/research comparison for single-component organic compounds against
the stated strain under the description-filtered assay setting. Not a clinical
susceptibility test, treatment recommendation, MRSA/general-strain predictor,
toxicity/selectivity model, or validated regulatory replacement for experiments.

## Data and chemical identity

Canonical isomeric SMILES preserve recorded charge and stereochemistry. No salt
stripping, pH-neutralization or tautomer merging is applied. Thus measured forms
are not silently replaced, but cross-source tautomer/protonation inconsistencies
may remain. Scaffold grouping ignores stereochemistry to keep stereoisomers in
one split. All acyclic molecules share one ACYCLIC group. Exact duplicates are
aggregated before splitting; conflict exclusions are logged before modeling.

## Scientific release assessment

Read `RESULTS.md` and `results/reference/metrics.json` before interpreting the model.
Failure of a diagnostic remains visible and is not repaired by selecting a new
split or retuning on the held-out test. Even strong results cannot resolve the
endpoint-harmonization and prospective-validation gaps above.

The card supplies evidence relevant to a QMRF-style assessment. It is not a
completed regulatory submission or an official OECD QMRF/QPRF certification.
