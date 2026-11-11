"""Fixed-protocol DeepChem descriptor QSAR with scaffold-held-out validation."""
import argparse, json, hashlib, sys, platform, importlib.metadata, warnings
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import deepchem as dc
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from modeling import pipeline, metrics, grouped_cv, gt_metrics, corrected_rp, bootstrap_r2, leverage

ROOT=Path(__file__).resolve().parents[1]

def save_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def fingerprint(smiles):
    generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    return [generator.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]

def nearest_similarity(query,reference,leave_self=False):
    values=[]
    for i,fp in enumerate(query):
        sims=np.asarray(DataStructs.BulkTanimotoSimilarity(fp,reference))
        if leave_self: sims[i]=-1
        values.append(float(sims.max()))
    return np.array(values)

def run(output):
    output=Path(output); output.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((ROOT/'config.json').read_text())
    save_json(output/'protocol.json',cfg)
    data=pd.read_csv(ROOT/'data/curated/compounds.csv')
    train=np.flatnonzero(data.split.values=='train'); test=np.flatnonzero(data.split.values=='test')
    if set(data.iloc[train].scaffold)&set(data.iloc[test].scaffold): raise ValueError('Scaffold leakage')
    feat=dc.feat.RDKitDescriptors(is_normalized=False,use_fragment=True,ipc_avg=True,use_bcut2d=True)
    names=list(feat.reqd_properties)
    X=np.asarray(feat.featurize(data.smiles.tolist()),dtype=float)
    if X.shape!=(len(data),len(names)): raise ValueError('Featurization changed rows/columns')
    y=data.pMIC.to_numpy(); groups=data.scaffold.to_numpy()
    np.savez_compressed(output/'descriptor_matrix.npz',X=X,compound_ids=data.compound_id.to_numpy(dtype=str))
    save_json(output/'descriptor_names.json',names)
    Xtr,ytr,gtr=X[train],y[train],groups[train]
    print(f'Features {X.shape}; train={len(train)}, test={len(test)}',flush=True)
    cv5,fold5,scores5,audit5=grouped_cv(Xtr,ytr,gtr,cfg,5)
    cv10,fold10,scores10,audit10=grouped_cv(Xtr,ytr,gtr,cfg,10)
    q2=metrics(ytr,cv5)['r2']
    pd.DataFrame({'compound_id':data.iloc[train].compound_id.values,'observed_pMIC':ytr,
                  'cv5_pred':cv5,'cv5_fold':fold5,'cv10_pred':cv10,'cv10_fold':fold10}).to_csv(output/'cv_predictions.csv',index=False)
    save_json(output/'cv_selection_audit.json',{'five_fold':audit5,'ten_fold':audit10})
    # The primary RF was fixed in config before any outcomes. No model search.
    core=pipeline(cfg)
    dataset=dc.data.NumpyDataset(X=Xtr,y=ytr[:,None],ids=data.iloc[train].compound_id.values)
    wrapper=dc.models.SklearnModel(core,model_dir=str(output/'deepchem'),use_weights=False)
    wrapper.fit(dataset)
    ptrain=wrapper.predict(dataset).reshape(-1)
    ptest=wrapper.predict(dc.data.NumpyDataset(X[test])).reshape(-1)
    select=core.named_steps['select']; chosen=[names[i] for i in select.selected_indices_]
    save_json(output/'selection_audit.json',dict(select.audit_,selected_descriptors=chosen,
                                                feature_f_scores=select.feature_scores_.tolist()))
    ridge=pipeline(cfg,'ridge').fit(Xtr,ytr)
    baseline={'training_mean':metrics(y[test],np.full(len(test),ytr.mean())),
              'ridge_alpha10':metrics(y[test],ridge.predict(X[test]))}
    Ztr=core[:-1].transform(Xtr); Ztest=core[:-1].transform(X[test])
    h,limit=leverage(Ztr,Ztest); htrain,_=leverage(Ztr,Ztr)
    fps_train=fingerprint(data.iloc[train].smiles); fps_test=fingerprint(data.iloc[test].smiles)
    loo_sim=nearest_similarity(fps_train,fps_train,leave_self=True)
    sim_limit=float(np.quantile(loo_sim,cfg['ad_similarity_quantile']))
    sims=nearest_similarity(fps_test,fps_train)
    in_ad=(h<=limit)&(sims>=sim_limit)
    train_docs=set(';'.join(data.iloc[train].document_ids).split(';'))
    test_docs=set(';'.join(data.iloc[test].document_ids).split(';'))
    predictions=data.iloc[test].copy()
    predictions['predicted_pMIC']=ptest; predictions['residual']=y[test]-ptest
    predictions['leverage']=h; predictions['nearest_train_tanimoto']=sims; predictions['in_AD']=in_ad
    predictions.to_csv(output/'test_predictions.csv',index=False)
    bootstrap,ci=bootstrap_r2(y[test],ptest,groups[test],cfg['bootstrap_replicates'],cfg['seed'])
    pd.DataFrame({'replicate':range(1,len(bootstrap)+1),'test_r2':bootstrap}).to_csv(output/'bootstrap_r2.csv',index=False)
    print(f'5-fold Q2={q2:.4f}; test metrics={metrics(y[test],ptest)}',flush=True)

    def scramble(i):
        yp=np.random.default_rng(cfg['seed']+10000+i).permutation(ytr)
        pp,_,_,_=grouped_cv(Xtr,yp,gtr,cfg,5)
        model=pipeline(cfg).fit(Xtr,yp)
        fit=model.predict(Xtr)
        return {'replicate':i+1,'seed':cfg['seed']+10000+i,'q2':metrics(yp,pp)['r2'],
                'training_r':float(np.corrcoef(yp,fit)[0,1]),'training_r2':metrics(yp,fit)['r2']}
    shuffled=[]
    with ThreadPoolExecutor(max_workers=cfg['scrambling_workers']) as pool:
        futures=[pool.submit(scramble,i) for i in range(cfg['scrambling_replicates'])]
        for future in as_completed(futures):
            shuffled.append(future.result())
            if len(shuffled)%10==0: print(f'Y-scrambling {len(shuffled)}/{cfg["scrambling_replicates"]}',flush=True)
    shuffled=sorted(shuffled,key=lambda r:r['replicate'])
    pd.DataFrame(shuffled).to_csv(output/'y_scrambling.csv',index=False)
    actual_r=float(np.corrcoef(ytr,ptrain)[0,1]); nullq=np.array([r['q2'] for r in shuffled])
    crp=corrected_rp(actual_r,[r['training_r'] for r in shuffled])
    report={'endpoint':cfg['endpoint'],'model':'DeepChem SklearnModel wrapping fixed RandomForest descriptor pipeline',
            'molecules':len(data),'train_molecules':len(train),'test_molecules':len(test),
            'initial_descriptors':len(names),'selected_descriptors':len(chosen),
            'molecules_per_descriptor_all':len(data)/len(chosen),'molecules_per_descriptor_train':len(train)/len(chosen),
            'test':metrics(y[test],ptest),'training':metrics(ytr,ptrain),'baseline_test':baseline,
            'five_fold_cv_Q2':q2,'ten_fold_cv_r2_mean':float(np.mean(scores10)),
            'ten_fold_cv_r2_sample_sd':float(np.std(scores10,ddof=1)),'ten_fold_scores':scores10,
            'bootstrap_test_r2':ci,'golbraikh_tropsha':gt_metrics(y[test],ptest),
            'y_scrambling':{'replicates':len(shuffled),'actual_q2':q2,'null_q2_mean':float(nullq.mean()),
                            'null_q2_sd':float(nullq.std(ddof=1)),'null_q2_max':float(nullq.max()),
                            'permutation_p':float((1+(nullq>=q2).sum())/(len(nullq)+1)),
                            'training_pearson_r':actual_r,'mean_scrambled_training_r':float(np.mean([r['training_r'] for r in shuffled])),
                            'cR2p':crp,'cR2p_above_0_5':bool(crp is not None and crp>0.5)},
            'AD':{'rule':'leverage <= 3(p+1)/n AND nearest training ECFP4 Tanimoto >= training LOO 5th percentile',
                  'leverage_threshold':limit,'similarity_threshold':sim_limit,'test_inside_n':int(in_ad.sum()),
                  'test_inside_percent':100*float(in_ad.mean()),'in_domain_test_metrics':metrics(y[test][in_ad],ptest[in_ad]) if in_ad.any() else None},
            'split':{'type':'fixed scaffold-group holdout; CV scaffold grouped; preprocessing refitted in each fold',
                     'scaffold_overlap':0,'train_scaffolds':len(set(gtr)),'test_scaffolds':len(set(groups[test])),
                     'train_documents':len(train_docs),'test_documents':len(test_docs),'shared_documents':len(train_docs&test_docs)},
            'limitations':['Cross-publication assay details not fully harmonized','Holdout is not an independent laboratory or study test',
                           'Exact MIC-only filtering excludes censored low-activity compounds','No regulatory acceptance or biochemical mechanism claimed']}
    save_json(output/'metrics.json',report)
    joblib.dump({'pipeline':core,'descriptor_names':names,'protocol':cfg,'training_scaled':Ztr,
                 'training_smiles':data.iloc[train].smiles.tolist(),'similarity_threshold':sim_limit},output/'model.joblib',compress=3)

    finalize(output)

def finalize(output):
    """Generate verified explanations/reports from a frozen model without refitting."""
    output=Path(output)
    cfg=json.loads((output/'protocol.json').read_text())
    report=json.loads((output/'metrics.json').read_text())
    data=pd.read_csv(ROOT/'data/curated/compounds.csv')
    train=np.flatnonzero(data.split.values=='train'); test=np.flatnonzero(data.split.values=='test')
    y=data.pMIC.to_numpy(); ytr=y[train]
    bundle=joblib.load(output/'model.joblib'); core=bundle['pipeline']
    X=np.load(output/'descriptor_matrix.npz')['X']
    select=core.named_steps['select']; chosen=[bundle['descriptor_names'][i] for i in select.selected_indices_]
    Ztr=bundle['training_scaled']; Ztest=core[:-1].transform(X[test])
    pred=pd.read_csv(output/'test_predictions.csv'); ptest=pred.predicted_pMIC.to_numpy()
    assert list(pred.compound_id)==data.iloc[test].compound_id.tolist()
    np.testing.assert_allclose(core.predict(X[test]),ptest,atol=1e-12)
    in_ad=pred.in_AD.to_numpy(); h=pred.leverage.to_numpy(); limit=report['AD']['leverage_threshold']
    htrain,_=leverage(Ztr,Ztr)
    cv5=pd.read_csv(output/'cv_predictions.csv').cv5_pred.to_numpy()
    nullq=pd.read_csv(output/'y_scrambling.csv').q2.to_numpy(); q2=report['five_fold_cv_Q2']
    # Explain fixed held-out predictions. Background comes solely from training rows.
    import shap
    rng=np.random.default_rng(cfg['seed'])
    bg_idx=np.sort(rng.choice(len(train),min(cfg['shap_background_size'],len(train)),replace=False))
    masker=shap.maskers.Independent(Ztr[bg_idx],max_samples=len(bg_idx))
    explainer=shap.Explainer(core.named_steps['model'].predict,masker,algorithm='permutation',seed=cfg['seed'])
    explanation=explainer(Ztest,max_evals=10*(2*Ztest.shape[1]+1))
    reconstructed=explanation.base_values+explanation.values.sum(axis=1)
    error=float(np.max(np.abs(reconstructed-ptest)))
    np.testing.assert_allclose(reconstructed,ptest,atol=1e-5,rtol=1e-5)
    raw_test=select.transform(X[test])
    display=shap.Explanation(values=explanation.values,base_values=explanation.base_values,data=raw_test,feature_names=chosen)
    np.savez_compressed(output/'shap_values.npz',values=explanation.values,base_values=explanation.base_values,
                        raw_descriptors=raw_test,background_scaled=Ztr[bg_idx])
    pd.DataFrame({'descriptor':chosen,'mean_abs_SHAP':np.mean(np.abs(explanation.values),axis=0)}).sort_values('mean_abs_SHAP',ascending=False).to_csv(output/'shap_importance.csv',index=False)
    save_json(output/'shap_audit.json',{'algorithm':'PermutationExplainer over exact RF predict callable; 10 forward/reverse permutation cycles per row','output':'predicted pMIC, log10 molar units','background_train_compound_ids':data.iloc[train[bg_idx]].compound_id.tolist(),
              'explained_test_compound_ids':data.iloc[test].compound_id.tolist(),'max_additivity_error':error,
              'interpretation':'Model attribution, not causal mechanism; correlated-feature masking assumptions apply'})
    shap.plots.beeswarm(display,max_display=15,show=False)
    plt.gcf().savefig(output/'shap_beeswarm.png',dpi=170,bbox_inches='tight'); plt.close('all')
    shap.plots.bar(display,max_display=15,show=False)
    plt.gcf().savefig(output/'shap_bar.png',dpi=170,bbox_inches='tight'); plt.close('all')
    fig,ax=plt.subplots(figsize=(6,5),layout='constrained')
    ax.scatter(y[test],ptest,c=np.where(in_ad,'#117f81','#c95a37'),alpha=.75)
    bounds=[min(y[test].min(),ptest.min()),max(y[test].max(),ptest.max())]
    ax.plot(bounds,bounds,'k--'); ax.set(xlabel='Observed pMIC',ylabel='Predicted pMIC',title='Scaffold-held-out test: teal within AD, orange outside')
    fig.savefig(output/'test_predictions.png',dpi=180); plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4),layout='constrained')
    ax.hist(nullq,bins=15,color='#627d98'); ax.axvline(q2,color='#c95a37',label='Actual 5-fold Q2')
    ax.set(xlabel='Scaffold 5-fold Q2',ylabel='Permutation count',title='Y-scrambling: selection refitted for every permutation'); ax.legend()
    fig.savefig(output/'y_scrambling.png',dpi=180); plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4),layout='constrained')
    cv_rmse=metrics(ytr,cv5)['rmse']
    ax.scatter(htrain,(ytr-cv5)/cv_rmse,s=15,alpha=.5,label='Train OOF residual')
    ax.scatter(h,(y[test]-ptest)/cv_rmse,s=20,label='Test residual')
    ax.axvline(limit,color='gray',linestyle='--'); ax.axhline(3,color='gray',linestyle=':'); ax.axhline(-3,color='gray',linestyle=':')
    ax.set(xlabel='Leverage',ylabel='Residual / training OOF RMSE',title='Williams-style diagnostic (AD also uses structural similarity)'); ax.legend()
    fig.savefig(output/'applicability_domain.png',dpi=180); plt.close(fig)
    provenance={'python':platform.python_version(),'versions':{x:importlib.metadata.version(x) for x in ['deepchem','rdkit','numpy','pandas','scikit-learn','shap','scipy','matplotlib','joblib']},
                'input_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'config.json',ROOT/'data/curated/compounds.csv']},
                'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'scripts').glob('*.py')},
                'output_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir() if p.is_file() and p.name!='manifest.json'}}
    save_json(output/'manifest.json',provenance)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output',type=Path,default=ROOT/'results/reference')
    p.add_argument('--finalize-only',action='store_true')
    a=p.parse_args()
    finalize(a.output) if a.finalize_only else run(a.output)
