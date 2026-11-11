# S. aureus ATCC 29213 antibacterial QSAR

Reproducible **phenotypic pMIC regression**, using DeepChem 2.8.0 RDKit descriptors
and its `SklearnModel` interface with a fixed Random Forest. This is separate from
the DHFR docking experiment: MIC measures whole-cell growth inhibition and is not
a DHFR binding endpoint. This release uses a descriptor model, not a graph neural
network or a pretrained antibacterial predictor.

## Data and endpoint

ChEMBL snapshot: literature MIC records for **S. aureus ATCC 29213**, with explicit
24-hour microdilution descriptions. Broad collection contained 989 records;
curation retained **542 compounds**, split before fitting into **439 training and
103 test molecules** (146 and 37 nonoverlapping Bemis–Murcko scaffold groups).

The response is `pMIC = -log10(MIC [mol/L])`; larger values mean lower MIC.
Exact `=` values only; no replacement of censored `>` or `<` values. Multi-component
structures/salts, unsupported elements and flagged records are excluded. Mass MIC
is converted using the molecular weight of the same single-component structure:
`MIC [M] = MIC [ug/mL] * 1e-3 / MW`. Canonical isomeric SMILES define duplicates;
replicates are aggregated by median pMIC only when their range is <=1 log unit.
Every rejected activity and its reason is retained.

Assay metadata were reviewed and bactericidal, MIC90, staining-based and
single-concentration descriptions excluded in addition to resistance variants,
co-treatment, serum, biofilm and special-condition experiments. Curation is based
on ChEMBL descriptions, not full-paper verification of all medium, inoculum and
readout details. Thus this is a **description-filtered, cross-publication endpoint**;
complete laboratory protocol harmonization is not established.

## Leakage controls and descriptor elimination

The primary RF settings were fixed in [config.json](config.json), committed with
the curated split in `74bd292` before training. There is no hyperparameter search
and no model selection based on the test set. Ridge(alpha=10) and training-mean
predictions are prespecified comparators.

Within **every training fold and every Y permutation**, the selector:

1. Removes descriptors with >5% missing/nonfinite values.
2. Learns median imputation from that fold's training rows.
3. Removes variance <=1e-12 and greedily removes Pearson-correlated descriptors
   with |r|>0.90, using stable descriptor order.
4. Selects up to 15 descriptors by univariate regression F score, also enforcing
   at least ten fitting molecules per selected descriptor.
5. Learns standardization from those training rows.

SHAP is a post-fit explanation, not a test-informed feature selector. The held-out
test has no shared exact scaffold with training; compounds from the same papers
can occur in both. This is a scaffold holdout, **not an independent-study or
prospective validation**. Cross-validation is scaffold-grouped too.

## Required validation

See [RESULTS.md](RESULTS.md) for the generated table and
[metric definitions](docs/METRICS.md) for formulas and limitations.

- Test predictive R², RMSE and MAE in pMIC units.
- Pooled out-of-fold 5-fold Q²; 10-fold R² mean ± sample SD.
- 1,000 held-out-scaffold bootstrap replicates, percentile 95% CI for fixed-model R².
- Golbraikh–Tropsha r0², r0'², k, k', Pearson r² and diagnostic flags.
- 100 full-pipeline Y permutations: actual versus scrambled 5-fold Q², permutation
  p-value and corrected cR²p based on training correlations (not renamed Q²).
- Applicability domain: leverage threshold **and** training-calibrated ECFP4
  nearest-neighbor similarity; percent of test compounds inside both.
- SHAP beeswarm/bar plots, training-only background, held-out explanations and
  explicit additivity verification. The final implementation uses permutation-SHAP
  over the exact predictor; see [explanation validation](docs/SHAP_VALIDATION.md).

OECD's five principles are documented in [OECD_MODEL_CARD.md](docs/OECD_MODEL_CARD.md).
They are a scientific reporting/assessment framework, not a certificate conferred
by passing a numerical checklist. SHAP does not establish an antibacterial mechanism.

## Reproduce and predict

Use an isolated **Python 3.11** environment; do not change the docking environment.
SHAP 0.51 is pinned for Python 3.11 compatibility. Optional DeepChem neural-network
backend warnings are expected; Torch/TensorFlow are not required by this model.

```bash
python -m venv .venv-qsar
# Activate the environment using your platform's normal command.
python -m pip install -r qsar/requirements.txt
python -m unittest discover -s qsar/tests -v
python qsar/scripts/validate.py
python qsar/scripts/predict.py --input qsar/examples/query.csv --output predictions.csv
```

For a full rebuild from the archived snapshot:

```bash
python qsar/scripts/curate.py
python qsar/scripts/train.py --output qsar/results/rerun1
```

New training output directories must not already exist. Raw download is optional:
`fetch_data.py` queries the current ChEMBL release and refuses to overwrite an
existing snapshot. A later release is a new dataset, not a guaranteed identical run.
The original broad snapshot predates final description-level exclusion refinements;
the archived inputs, curation code and fixed split define this release.

The bundled `model.joblib` contains only the model fitted to the 439 training
molecules; the holdout was not folded back into a deployment refit. Load joblib
artifacts only from trusted sources. Prediction output includes `in_AD` and
`seen_in_training`; being inside the AD does not guarantee an accurate MIC.

## Files and attribution

`data/raw/`: compressed unmodified API record payloads, selection audit and provenance.
`data/curated/`: compounds/splits, measurements, rejected records and included assays.
`results/reference/`: model, all predictions, resampling results, explanations and hashes.
`scripts/`: acquisition, curation, modeling, validation and prediction.

ChEMBL data are redistributed under **CC BY-SA 3.0**, separately from the repository's
MIT-licensed original code. See [data attribution](data/LICENSE.md) and
[references](docs/REFERENCES.md). No docking score is used as an activity label.
