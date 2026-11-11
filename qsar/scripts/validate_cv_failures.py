"""Validate diagnostic evidence against frozen CV predictions and curated identities."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def validate(directory=ROOT / 'results/cv_failure_analysis'):
    directory = Path(directory)
    manifest = json.loads((directory / 'manifest.json').read_text())
    for name, sha in manifest['input_sha256'].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == sha, name
    for name, sha in manifest['output_sha256'].items():
        assert hashlib.sha256((directory / name).read_bytes()).hexdigest() == sha, name
    assert hashlib.sha256((ROOT / 'scripts/analyze_cv_failures.py').read_bytes()).hexdigest() == manifest['analysis_source_sha256']
    data = pd.read_csv(ROOT / 'data/curated/compounds.csv').set_index('compound_id')
    cv = pd.read_csv(ROOT / 'results/reference/cv_predictions.csv').set_index('compound_id')
    rows = pd.read_csv(directory / 'compound_diagnostics.csv').set_index('compound_id')
    summary = pd.read_csv(directory / 'fold_summary.csv').set_index('fold')
    scaffolds = pd.read_csv(directory / 'scaffold_diagnostics.csv')
    audit = json.loads((directory / 'fold_audit.json').read_text())
    coords = pd.read_csv(directory / 'negative_fold_pca_coordinates.csv')
    assert rows.index.is_unique and set(rows.index) == set(cv.index)
    assert set(rows.index) == set(data.index[data.split == 'train'])
    assert (rows.scaffold == data.loc[rows.index, 'scaffold']).all()
    np.testing.assert_allclose(rows.pMIC, data.loc[rows.index, 'pMIC'])
    np.testing.assert_allclose(rows.predicted_pMIC, cv.loc[rows.index, 'cv10_pred'], atol=1e-10, rtol=1e-10)
    np.testing.assert_array_equal(rows.fold, cv.loc[rows.index, 'cv10_fold'] + 1)
    np.testing.assert_allclose(rows.residual_pred_minus_obs, rows.predicted_pMIC - rows.pMIC, atol=1e-12)
    np.testing.assert_allclose(rows.absolute_error, abs(rows.residual_pred_minus_obs), atol=1e-12)
    np.testing.assert_allclose(rows.squared_error, rows.residual_pred_minus_obs ** 2, atol=1e-12)
    np.testing.assert_array_equal(rows.in_AD, (rows.leverage <= rows.leverage_threshold) & (rows.nearest_train_tanimoto >= rows.similarity_threshold))
    np.testing.assert_allclose(rows.leverage_ratio, rows.leverage / rows.leverage_threshold)
    np.testing.assert_allclose(rows.descriptor_distance_ratio, rows.descriptor_nn_distance / rows.descriptor_train_LOO_p95)
    for fold, group in rows.groupby('fold'):
        s = summary.loc[fold]
        a = audit[str(fold)]
        train = set(a['train_compound_ids'])
        valid = set(a['validation_compound_ids'])
        assert valid == set(group.index) and train == set(rows.index) - valid
        assert not set(data.loc[list(train), 'scaffold']) & set(group.scaffold)
        assert set(group.closest_structural_train_id) <= train
        assert set(group.closest_descriptor_train_id) <= train
        np.testing.assert_allclose(group.closest_structural_train_pMIC, data.loc[group.closest_structural_train_id, 'pMIC'])
        assert s.validation_n == len(group) and s.train_n == len(train)
        np.testing.assert_allclose(group.leverage_threshold, 3 * (len(a['selected_descriptors']) + 1) / len(train))
        np.testing.assert_allclose(group.similarity_threshold, a['similarity_threshold'])
        np.testing.assert_allclose(group.descriptor_train_LOO_p95, a['descriptor_distance_p95'])
        variance = np.var(group.pMIC)
        expected = {'r2': 1 - group.squared_error.mean() / variance,
                    'rmse': np.sqrt(group.squared_error.mean()), 'mae': group.absolute_error.mean(),
                    'sse': group.squared_error.sum(), 'validation_variance': variance,
                    'within_AD_percent': group.in_AD.mean() * 100,
                    'sse_share_outside_AD': group.loc[~group.in_AD, 'squared_error'].sum() / group.squared_error.sum()}
        for key, value in expected.items():
            np.testing.assert_allclose(s[key], value, atol=1e-12)
        baseline_mse = np.mean((group.pMIC - data.loc[list(train), 'pMIC'].mean()) ** 2)
        np.testing.assert_allclose(s.training_mean_baseline_rmse, np.sqrt(baseline_mse))
        table = scaffolds[scaffolds.fold == fold].set_index('scaffold')
        assert set(table.index) == set(group.scaffold)
        for scaffold, g in group.groupby('scaffold'):
            assert table.loc[scaffold, 'n'] == len(g)
            np.testing.assert_allclose(table.loc[scaffold, 'sse'], g.squared_error.sum())
            np.testing.assert_allclose(table.loc[scaffold, 'fold_sse_percent'], 100 * g.squared_error.sum() / s.sse)
        if s.r2 < 0:
            positions = coords[coords.fold == fold]
            assert set(positions.loc[positions.role == 'fold_train', 'compound_id']) == train
            assert set(positions.loc[positions.role == 'fold_validation', 'compound_id']) == valid
            assert positions.compound_id.is_unique
    assert sorted(summary.index[summary.r2 < 0].tolist()) == manifest['negative_folds_one_based'] == [1, 9]
    print('PASS: diagnostic hashes, frozen predictions, fold-local identities, AD flags, distances, metrics and scaffold error totals')


if __name__ == '__main__':
    validate()
