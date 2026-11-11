# Execution record

2026-11-11: Python 3.11.16, Windows. Exact installed scientific package versions
are recorded in `results/reference/manifest.json`; the independent QSAR environment
did not modify the docking environment.

- Fixed protocol and data/splits committed as `74bd292` before fitting.
- 217 descriptors computed for 542 compounds; final training selector retained 15.
- 5- and 10-fold scaffold-grouped CV completed with per-fold selection audits.
- 1,000 bootstrap replicates were valid; all 100 Y permutations completed.
- Seven unit tests passed: distinct predictive R²/correlation, GT edge cases,
  cR²p formula, bootstrap, fold-local preprocessing, leverage and chemical-unit rules.
- Archived evidence validator passed: input/output hashes, split integrity,
  recalculated test/CV/GT/bootstrap/permutation metrics, AD coverage and SHAP sums.
- Independent SMILES prediction smoke test passed. The example trimethoprim
  molecule is **outside AD** and absent from training; its numerical prediction
  should not be treated as validated MIC or validation of the docking ligand.
- Figures visually inspected. SHAP maximum reconstruction error over all held-out
  compounds: approximately 6.1e-7 pMIC. See `SHAP_VALIDATION.md` for the corrected
  explanation route; fitting outputs and evaluation metrics were unchanged.

The 10-fold dispersion and cross-publication endpoint limitations remain material.
No test split was changed to improve reported performance. Skill-guided DeepChem
and SHAP procedural sources are credited in `REFERENCES.md`.
