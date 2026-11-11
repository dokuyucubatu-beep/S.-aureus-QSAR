"""Render the requested metric table directly from the archived JSON report."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'results/reference/metrics.json').read_text())
def f(x): return 'undefined' if x is None else f'{x:.4f}'
g=r['golbraikh_tropsha']; s=r['y_scrambling']; ad=r['AD']; ci=r['bootstrap_test_r2']['ci95']
rows=[
('Molecules : selected descriptors',f"{r['molecules']} : {r['selected_descriptors']} overall; {r['train_molecules']} : {r['selected_descriptors']} fitted ({r['molecules_per_descriptor_train']:.2f}:1)"),
('Descriptor elimination',f"{r['initial_descriptors']} initial -> {r['selected_descriptors']} selected; refitted within each CV fold/permutation"),
('Test R²',f(r['test']['r2'])),('Test RMSE / MAE (pMIC)',f(r['test']['rmse'])+' / '+f(r['test']['mae'])),
('5-fold CV Q² (pooled OOF)',f(r['five_fold_cv_Q2'])),
('10-fold CV R² (mean ± sample SD)',f(r['ten_fold_cv_r2_mean'])+' ± '+f(r['ten_fold_cv_r2_sample_sd'])),
('Bootstrap 95% CI for test R² (n=1000)',f(ci[0])+' to '+f(ci[1])+'; scaffold-cluster bootstrap'),
("Golbraikh–Tropsha r0² / r0'²",f(g['r0_squared'])+' / '+f(g['r0_prime_squared'])),
("Golbraikh–Tropsha k / k'",f(g['k'])+' / '+f(g['k_prime'])),
('Pearson r² (distinct from predictive R²)',f(g['pearson_r_squared'])),
('Y-scrambling cR²p',f(s['cR2p'])),
('Actual Q² / mean scrambled Q²',f(s['actual_q2'])+' / '+f(s['null_q2_mean'])),
('Y-scrambling permutation p (100 repeats)',f(s['permutation_p'])),
('Test compounds within AD',f"{ad['test_inside_percent']:.2f}% ({ad['test_inside_n']}/{r['test_molecules']})")]
text='# Antibacterial QSAR results\n\n'
text+='**S. aureus ATCC 29213; pMIC regression; fixed scaffold-held-out test.**\n\n'
text+='| Requested metric | Result |\n|---|---|\n'+'\n'.join('| '+a+' | '+b+' |' for a,b in rows)+'\n\n'
text+='## Interpretation and limitations\n\n'
text+=f"The RF outperformed the predefined Ridge baseline (test R² {r['baseline_test']['ridge_alpha10']['r2']:.3f}) and training-mean predictor. However, the 10-fold scores range from {min(r['ten_fold_scores']):.3f} to {max(r['ten_fold_scores']):.3f}, with two negative folds: transfer to new scaffolds is uneven. The test confidence interval is conditional on this fitted model.\n\n"
text+=f"The stated directional-OR Golbraikh–Tropsha diagnostic {'passes' if g['prespecified_diagnostic_pass'] else 'fails'}, while the reverse-direction relative discrepancy is {g['relative_difference_prime']:.3f} and does not individually meet 0.1. cR²p {'exceeds' if s['cR2p_above_0_5'] else 'does not exceed'} 0.5. These are conventional diagnostics, not OECD certification.\n\n"
text+=f"The scaffold split has zero scaffold overlap but {r['split']['shared_documents']} shared source documents. Assay protocol metadata are not fully harmonized. The model remains a research baseline; no prospective, independent-study or regulatory validation is claimed.\n\n"
text+='## Figures and evidence\n\n![Held-out predictions](results/reference/test_predictions.png)\n\n![SHAP](results/reference/shap_beeswarm.png)\n\n![Y-scrambling](results/reference/y_scrambling.png)\n\n![Applicability domain](results/reference/applicability_domain.png)\n\n'
text+='[Full metrics](results/reference/metrics.json), [test predictions](results/reference/test_predictions.csv), [CV predictions](results/reference/cv_predictions.csv), [bootstrap replicates](results/reference/bootstrap_r2.csv), [permutations](results/reference/y_scrambling.csv), [selected descriptors](results/reference/selection_audit.json), [SHAP audit](results/reference/shap_audit.json), [OECD evidence map](docs/OECD_MODEL_CARD.md).\n'
(ROOT/'RESULTS.md').write_text(text,encoding='utf-8',newline='\n')
print('Rendered qsar/RESULTS.md from measured outputs')
