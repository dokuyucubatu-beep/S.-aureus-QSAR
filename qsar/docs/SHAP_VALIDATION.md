# SHAP implementation and additivity audit

The evaluated model is frozen before explanation. Explanations concern predicted
pMIC, with 100 deterministic background compounds drawn exclusively from training.
All 103 held-out molecules are explained. Selected descriptor values are
standardized for the model; plot colors display their corresponding raw values.

On the pinned stack, the initial interventional TreeExplainer reconstruction
differed from the actual Random Forest predictions by up to about 0.0143 pMIC.
The independent strict additivity check caught this despite the library's looser
internal check. Float32/float64 checks did not resolve it; the internal tree
predictor and scikit-learn predictor agreed, so these explanations were not released.

The final implementation uses **PermutationExplainer on the exact RF `predict`
callable**, with an Independent masker holding the same training background,
fixed seed 2026 and ten forward/reverse permutation cycles per row (310 maximum
evaluations for 15 selected descriptors). Attribution values are Monte Carlo
estimates, not exact enumerations of all feature subsets. Their sum plus the
baseline must match the actual model prediction at atol=rtol=1e-5; the observed
maximum error is saved in `shap_audit.json`. No tolerance was relaxed to make an
explanation pass, and no model, split or outcome was changed for this correction.

Correlated-feature masking assumptions still apply after pairwise correlation
filtering. Positive attribution increases predicted pMIC (lower predicted MIC)
relative to the chosen training background. This describes model behavior and
does not prove causal antibacterial action or a biochemical mechanism.
