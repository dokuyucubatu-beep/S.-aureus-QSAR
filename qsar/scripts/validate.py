"""Independently regenerate requested metrics from archived predictions, without fitting."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from modeling import metrics,gt_metrics,corrected_rp

ROOT=Path(__file__).resolve().parents[1]

def validate(directory):
    directory=Path(directory)
    manifest=json.loads((directory/'manifest.json').read_text())
    for rel,sha in manifest['input_sha256'].items():
        if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=sha: raise ValueError(f'Input hash mismatch: {rel}')
    for name,sha in manifest['output_sha256'].items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest()!=sha: raise ValueError(f'Output hash mismatch: {name}')
    data=pd.read_csv(ROOT/'data/curated/compounds.csv')
    train=data[data.split=='train']; test=data[data.split=='test']
    assert not set(train.scaffold)&set(test.scaffold)
    assert not set(train.smiles)&set(test.smiles)
    report=json.loads((directory/'metrics.json').read_text())
    pred=pd.read_csv(directory/'test_predictions.csv')
    assert set(pred.compound_id)==set(test.compound_id)
    truth=test.set_index('compound_id').loc[pred.compound_id,'pMIC'].to_numpy()
    np.testing.assert_allclose(truth,pred.pMIC)
    for key,value in metrics(truth,pred.predicted_pMIC.to_numpy()).items():
        np.testing.assert_allclose(value,report['test'][key],rtol=1e-10)
    gt=gt_metrics(truth,pred.predicted_pMIC.to_numpy())
    for key in ['r0_squared','r0_prime_squared','k','k_prime']:
        np.testing.assert_allclose(gt[key],report['golbraikh_tropsha'][key])
    cv=pd.read_csv(directory/'cv_predictions.csv')
    np.testing.assert_allclose(metrics(cv.observed_pMIC.values,cv.cv5_pred.values)['r2'],report['five_fold_cv_Q2'])
    fold_scores=[metrics(g.observed_pMIC.values,g.cv10_pred.values)['r2'] for _,g in cv.groupby('cv10_fold')]
    np.testing.assert_allclose(np.mean(fold_scores),report['ten_fold_cv_r2_mean'])
    np.testing.assert_allclose(np.std(fold_scores,ddof=1),report['ten_fold_cv_r2_sample_sd'])
    scaffold_map=train.set_index('compound_id').scaffold
    for col in ['cv5_fold','cv10_fold']:
        assigned=cv.assign(scaffold=cv.compound_id.map(scaffold_map))
        assert assigned.groupby('scaffold')[col].nunique().max()==1
    boots=pd.read_csv(directory/'bootstrap_r2.csv')
    assert len(boots)==1000
    np.testing.assert_allclose(boots.test_r2.dropna().quantile([.025,.975]).values,report['bootstrap_test_r2']['ci95'])
    scrambled=pd.read_csv(directory/'y_scrambling.csv')
    assert len(scrambled)==report['y_scrambling']['replicates']
    np.testing.assert_allclose(corrected_rp(report['y_scrambling']['training_pearson_r'],scrambled.training_r),report['y_scrambling']['cR2p'])
    expected_p=(1+(scrambled.q2>=report['five_fold_cv_Q2']).sum())/(len(scrambled)+1)
    np.testing.assert_allclose(expected_p,report['y_scrambling']['permutation_p'])
    np.testing.assert_allclose(pred.in_AD.mean()*100,report['AD']['test_inside_percent'])
    shap=np.load(directory/'shap_values.npz')
    np.testing.assert_allclose(shap['base_values']+shap['values'].sum(1),pred.predicted_pMIC,rtol=1e-5,atol=1e-5)
    shap_audit=json.loads((directory/'shap_audit.json').read_text())
    assert set(shap_audit['background_train_compound_ids'])<=set(train.compound_id)
    assert set(shap_audit['explained_test_compound_ids'])==set(test.compound_id)
    print('PASS: evidence hashes, scaffold separation, metrics, 1000 bootstrap replicates, scrambling, AD and SHAP additivity')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--results',type=Path,default=ROOT/'results/reference')
    validate(p.parse_args().results)
