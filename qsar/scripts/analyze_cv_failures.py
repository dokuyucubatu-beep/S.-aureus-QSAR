"""Post-hoc diagnosis of frozen 10-fold CV: fold-local AD, distances and scaffolds."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, Draw
from modeling import pipeline, metrics, leverage

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'results/reference'

def nearest_fingerprints(query, reference, leave_self=False):
    scores, indices=[],[]
    for i,fp in enumerate(query):
        s=np.asarray(DataStructs.BulkTanimotoSimilarity(fp,reference))
        if leave_self: s[i]=-1
        j=int(np.argmax(s)); scores.append(float(s[j])); indices.append(j)
    return np.array(scores),np.array(indices)

def run(output):
    output=Path(output); output.mkdir(parents=True,exist_ok=False)
    data=pd.read_csv(ROOT/'data/curated/compounds.csv')
    cv=pd.read_csv(REFERENCE/'cv_predictions.csv')
    data=data.set_index('compound_id').loc[cv.compound_id].reset_index()
    assert (data.split=='train').all()  # The external 103-molecule holdout is not used.
    archive=np.load(REFERENCE/'descriptor_matrix.npz')
    index={s:i for i,s in enumerate(archive['compound_ids'])}
    X=archive['X'][[index[s] for s in data.compound_id]]
    names=json.loads((REFERENCE/'descriptor_names.json').read_text())
    cfg=json.loads((REFERENCE/'protocol.json').read_text())
    y=data.pMIC.to_numpy(); fold_ids=cv.cv10_fold.to_numpy()
    generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    mols=[Chem.MolFromSmiles(s) for s in data.smiles]
    fps=[generator.GetFingerprint(m) for m in mols]
    summaries,all_rows,details,positions=[],[],{},[]
    for fold in sorted(set(fold_ids)):
        tr=np.flatnonzero(fold_ids!=fold); va=np.flatnonzero(fold_ids==fold)
        assert not set(data.iloc[tr].scaffold)&set(data.iloc[va].scaffold)
        model=pipeline(cfg).fit(X[tr],y[tr])
        p=model.predict(X[va]); expected=cv.iloc[va].cv10_pred.to_numpy()
        np.testing.assert_allclose(p,expected,atol=1e-10,rtol=1e-10)
        Ztr=model[:-1].transform(X[tr]); Zva=model[:-1].transform(X[va])
        h,hstar=leverage(Ztr,Zva)
        fp_tr=[fps[i] for i in tr]; fp_va=[fps[i] for i in va]
        loo_sim,_=nearest_fingerprints(fp_tr,fp_tr,True)
        sim_limit=float(np.quantile(loo_sim,cfg['ad_similarity_quantile']))
        sim,analog=nearest_fingerprints(fp_va,fp_tr)
        in_ad=(h<=hstar)&(sim>=sim_limit)
        # Additional descriptor-distance diagnostic; it does not redefine the frozen AD.
        nn=NearestNeighbors(n_neighbors=2).fit(Ztr)
        loo_d=nn.kneighbors(Ztr)[0][:,1]/np.sqrt(Ztr.shape[1])
        d_limit=float(np.quantile(loo_d,.95))
        dist,idx=nn.kneighbors(Zva,n_neighbors=1)
        dist=dist[:,0]/np.sqrt(Ztr.shape[1]); idx=idx[:,0]
        residual=p-y[va]; squared=residual**2; m=metrics(y[va],p)
        row=data.iloc[va].copy(); row['fold']=int(fold)+1
        row['predicted_pMIC']=p; row['residual_pred_minus_obs']=residual; row['absolute_error']=abs(residual)
        row['squared_error']=squared; row['leverage']=h; row['leverage_threshold']=hstar; row['leverage_ratio']=h/hstar
        row['nearest_train_tanimoto']=sim; row['similarity_threshold']=sim_limit
        row['in_AD']=in_ad; row['leverage_pass']=h<=hstar; row['similarity_pass']=sim>=sim_limit
        row['descriptor_nn_distance']=dist; row['descriptor_train_LOO_p95']=d_limit
        row['descriptor_distance_ratio']=dist/d_limit
        row['closest_descriptor_train_id']=data.iloc[tr[idx]].compound_id.to_numpy()
        row['closest_structural_train_id']=data.iloc[tr[analog]].compound_id.to_numpy()
        row['closest_structural_train_smiles']=data.iloc[tr[analog]].smiles.to_numpy()
        row['closest_structural_train_pMIC']=y[tr[analog]]
        row['activity_gap_to_structural_analog']=y[va]-y[tr[analog]]
        row['same_source_document_as_analog']=[bool(set(a.split(';'))&set(b.split(';'))) for a,b in zip(row.document_ids,data.iloc[tr[analog]].document_ids)]
        all_rows.append(row)
        baseline=metrics(y[va],np.full(len(va),y[tr].mean()))
        summaries.append(dict(fold=int(fold)+1,train_n=len(tr),validation_n=len(va),validation_scaffolds=data.iloc[va].scaffold.nunique(),
            r2=m['r2'],rmse=m['rmse'],mae=m['mae'],bias=float(residual.mean()),
            validation_pMIC_mean=float(y[va].mean()),validation_pMIC_sd=float(y[va].std(ddof=0)),
            validation_pMIC_min=float(y[va].min()),validation_pMIC_max=float(y[va].max()),
            train_mean_pMIC=float(y[tr].mean()),train_pMIC_sd=float(y[tr].std(ddof=0)),
            mse=float(np.mean(squared)),validation_variance=float(np.var(y[va])),sse=float(squared.sum()),
            training_mean_baseline_r2=baseline['r2'],training_mean_baseline_rmse=baseline['rmse'],
            within_AD_n=int(in_ad.sum()),within_AD_percent=float(in_ad.mean()*100),
            leverage_outside_n=int((h>hstar).sum()),similarity_outside_n=int((sim<sim_limit).sum()),
            descriptor_distance_over_p95_n=int((dist>d_limit).sum()),
            within_AD_mae=float(np.mean(abs(residual[in_ad]))) if in_ad.any() else None,
            outside_AD_mae=float(np.mean(abs(residual[~in_ad]))) if (~in_ad).any() else None,
            sse_share_outside_AD=float(squared[~in_ad].sum()/squared.sum()),
            max_reproduction_difference=float(np.max(abs(p-expected)))))
        selected=[names[i] for i in model.named_steps['select'].selected_indices_]
        audit={'train_compound_ids':data.iloc[tr].compound_id.tolist(),'validation_compound_ids':data.iloc[va].compound_id.tolist(),
               'selected_descriptors':selected,'selector':model.named_steps['select'].audit_,
               'leverage_threshold':hstar,'similarity_threshold':sim_limit,'descriptor_distance_p95':d_limit}
        if m['r2']<0:
            pca=PCA(n_components=2).fit(Ztr)
            ptr,pva=pca.transform(Ztr),pca.transform(Zva)
            audit['pca_explained_variance_ratio']=pca.explained_variance_ratio_.tolist()
            loadings=pd.DataFrame(pca.components_.T,index=selected,columns=['PC1','PC2'])
            loadings.to_csv(output/f'fold_{int(fold)+1}_pca_loadings.csv',index_label='descriptor')
            for ids,xy,role in [(tr,ptr,'fold_train'),(va,pva,'fold_validation')]:
                for i,(a,b) in zip(ids,xy): positions.append({'fold':int(fold)+1,'compound_id':data.iloc[i].compound_id,'role':role,'PC1':a,'PC2':b})
            fig,axs=plt.subplots(1,3,figsize=(16,4.8),layout='constrained')
            ax=axs[0]; ax.scatter(ptr[:,0],ptr[:,1],color='lightgray',s=13,label='Fold training')
            vmax=max(abs(residual).max(),1)
            for mask,marker,label in [(in_ad,'o','Validation: inside AD'),(~in_ad,'^','Validation: outside AD')]:
                sc=ax.scatter(pva[mask,0],pva[mask,1],c=abs(residual[mask]),cmap='magma_r',vmin=0,vmax=vmax,marker=marker,s=45,label=label)
            fig.colorbar(sc,ax=ax,label='Absolute pMIC error')
            ax.set(xlabel=f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)',ylabel=f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)',title='Training-fitted PCA (projection only)'); ax.legend(fontsize=7)
            axs[1].hist(y[tr],bins=15,density=True,alpha=.5,label='Fold training')
            axs[1].hist(y[va],bins=10,density=True,alpha=.65,label='Fold validation')
            axs[1].set(xlabel='Observed pMIC',ylabel='Density',title=f'R2={m["r2"]:.3f}; RMSE={m["rmse"]:.3f}'); axs[1].legend(fontsize=8)
            axs[2].scatter(h/hstar,abs(residual),c=sim,cmap='viridis',vmin=0,vmax=1,s=40)
            axs[2].axvline(1,color='gray',ls='--'); axs[2].set(xlabel='Leverage / fold threshold',ylabel='Absolute pMIC error',title='Fold-local extrapolation diagnostic')
            for j in np.argsort(-squared)[:3]: axs[2].annotate(row.iloc[j].compound_id[:6],(h[j]/hstar,abs(residual[j])),fontsize=7)
            fig.suptitle(f'Frozen 10-fold CV: fold {int(fold)+1}')
            fig.savefig(output/f'fold_{int(fold)+1}_diagnostic.png',dpi=170); plt.close(fig)
        details[str(int(fold)+1)]=audit
        print(f'Fold {int(fold)+1}: R2={m["r2"]:.3f}, n={len(va)}, AD={in_ad.mean()*100:.1f}%',flush=True)
    summary=pd.DataFrame(summaries); predictions=pd.concat(all_rows,ignore_index=True)
    summary.to_csv(output/'fold_summary.csv',index=False); predictions.to_csv(output/'compound_diagnostics.csv',index=False)
    pd.DataFrame(positions).to_csv(output/'negative_fold_pca_coordinates.csv',index=False)
    (output/'fold_audit.json').write_text(json.dumps(details,indent=2)+'\n')
    aggregates=[]
    for (fold,scaffold),g in predictions.groupby(['fold','scaffold'],sort=True):
        fold_sse=float(summary.set_index('fold').loc[fold,'sse'])
        representative=g.sort_values('squared_error',ascending=False).iloc[0]
        aggregates.append({'fold':int(fold),'scaffold':scaffold,'n':len(g),'mae':g.absolute_error.mean(),
            'rmse':np.sqrt(g.squared_error.mean()),'bias':g.residual_pred_minus_obs.mean(),
            'pMIC_mean':g.pMIC.mean(),'pMIC_min':g.pMIC.min(),'pMIC_max':g.pMIC.max(),
            'outside_AD_percent':100*(~g.in_AD).mean(),'mean_nearest_tanimoto':g.nearest_train_tanimoto.mean(),
            'max_leverage_ratio':g.leverage_ratio.max(),'median_descriptor_distance_ratio':g.descriptor_distance_ratio.median(),
            'sse':g.squared_error.sum(),'fold_sse_percent':100*g.squared_error.sum()/fold_sse,
            'representative_id':representative.compound_id,'representative_smiles':representative.smiles,
            'representative_analog_smiles':representative.closest_structural_train_smiles,
            'representative_analog_id':representative.closest_structural_train_id,
            'document_ids':';'.join(sorted(set(';'.join(g.document_ids).split(';'))))})
    scaffold_table=pd.DataFrame(aggregates).sort_values(['fold','sse'],ascending=[True,False])
    scaffold_table.to_csv(output/'scaffold_diagnostics.csv',index=False)
    negative=summary.loc[summary.r2<0,'fold'].tolist()
    for fold in negative:
        top=scaffold_table[scaffold_table.fold==fold].head(4)
        ms=[]; labels=[]
        for rank,(_,r) in enumerate(top.iterrows(),1):
            ms.extend([Chem.MolFromSmiles(r.scaffold) if r.scaffold!='ACYCLIC' else Chem.MolFromSmiles(r.representative_smiles),
                       Chem.MolFromSmiles(r.representative_smiles),Chem.MolFromSmiles(r.representative_analog_smiles)])
            labels.extend([f'S{rank}: n={r.n}; SSE share={r.fold_sse_percent:.1f}%',f'Worst representative {r.representative_id[:8]}',f'Closest training analog {r.representative_analog_id[:8]}'])
        Draw.MolsToGridImage(ms,molsPerRow=3,subImgSize=(350,270),legends=labels).save(str(output/f'fold_{fold}_scaffolds_and_analogs.png'))
    inputs=[ROOT/'data/curated/compounds.csv',REFERENCE/'cv_predictions.csv',REFERENCE/'descriptor_matrix.npz',REFERENCE/'protocol.json']
    manifest={'scope':'Post-hoc diagnostic of archived 10-fold training CV; no test labels, model selection or redefinition of AD',
              'negative_folds_one_based':negative,'input_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
              'output_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir() if p.is_file()},
              'analysis_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output',type=Path,default=ROOT/'results/cv_failure_analysis')
    run(p.parse_args().output)
