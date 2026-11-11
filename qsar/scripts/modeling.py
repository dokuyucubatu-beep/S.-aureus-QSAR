"""Leakage-resistant descriptor selection and explicit validation definitions."""
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.feature_selection import f_regression
from sklearn.metrics import r2_score, mean_absolute_error, root_mean_squared_error
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


class DescriptorSelector(BaseEstimator, TransformerMixin):
    def __init__(self,max_features=15,correlation_cutoff=0.9,max_missing_fraction=0.05,min_ratio=10):
        self.max_features=max_features; self.correlation_cutoff=correlation_cutoff
        self.max_missing_fraction=max_missing_fraction; self.min_ratio=min_ratio

    def fit(self,X,y):
        X=np.asarray(X,dtype=float).copy(); X[~np.isfinite(X)]=np.nan
        self.n_features_in_=X.shape[1]
        self.valid_=np.where(np.mean(np.isnan(X),axis=0)<=self.max_missing_fraction)[0]
        if not len(self.valid_): raise ValueError('No descriptors pass missingness filter')
        self.medians_=np.nanmedian(X[:,self.valid_],axis=0)
        Z=np.where(np.isnan(X[:,self.valid_]),self.medians_,X[:,self.valid_])
        variable=np.where(np.var(Z,axis=0)>1e-12)[0]
        if not len(variable): raise ValueError('No variable descriptors')
        correlation=np.abs(np.corrcoef(Z[:,variable],rowvar=False)) if len(variable)>1 else np.ones((1,1))
        kept=[]
        for j in range(len(variable)):
            if not any(correlation[j,k]>self.correlation_cutoff for k in kept): kept.append(j)
        self.decorrelated_=variable[kept]
        scores=f_regression(Z[:,self.decorrelated_],y,force_finite=True)[0]
        k=min(self.max_features,max(1,len(y)//self.min_ratio),len(scores))
        order=np.argsort(-scores,kind='stable')[:k]
        self.selected_local_=self.decorrelated_[order]
        self.selected_indices_=self.valid_[self.selected_local_]
        self.feature_scores_=scores[order]
        self.audit_={'initial':X.shape[1],'after_missingness':len(self.valid_),
                     'after_variance':len(variable),'after_correlation':len(kept),'selected':k,'fit_n':len(y)}
        return self

    def transform(self,X):
        X=np.asarray(X,dtype=float).copy()
        if X.shape[1]!=self.n_features_in_: raise ValueError('Descriptor count changed')
        X[~np.isfinite(X)]=np.nan
        Z=np.where(np.isnan(X[:,self.valid_]),self.medians_,X[:,self.valid_])
        return Z[:,self.selected_local_]


def pipeline(cfg,kind='rf'):
    model=RandomForestRegressor(**cfg['rf'],random_state=cfg['seed'],n_jobs=1) if kind=='rf' else Ridge(alpha=10)
    return Pipeline([('select',DescriptorSelector(cfg['max_descriptors'],cfg['correlation_cutoff'],cfg['max_missing_fraction'],cfg['min_molecules_per_descriptor'])),
                     ('scale',StandardScaler()),('model',model)])


def metrics(y,p):
    y,p=np.asarray(y),np.asarray(p)
    if y.shape!=p.shape or not np.isfinite(y).all() or not np.isfinite(p).all(): raise ValueError('Invalid predictions')
    return {'r2':float(r2_score(y,p)) if len(y)>1 and np.var(y)>0 else None,
            'rmse':float(root_mean_squared_error(y,p)),'mae':float(mean_absolute_error(y,p))}


def grouped_cv(X,y,groups,cfg,folds=5):
    pred=np.full(len(y),np.nan); fold_id=np.full(len(y),-1); audits=[]; scores=[]
    for i,(tr,va) in enumerate(GroupKFold(folds,shuffle=True,random_state=cfg['seed']+folds).split(X,y,groups)):
        if set(groups[tr])&set(groups[va]): raise AssertionError('Scaffold leakage')
        model=pipeline(cfg).fit(X[tr],y[tr])
        pred[va]=model.predict(X[va]); fold_id[va]=i
        audits.append(dict(fold=i,**model.named_steps['select'].audit_,selected_indices=model.named_steps['select'].selected_indices_.tolist()))
        scores.append(metrics(y[va],pred[va])['r2'])
    return pred,fold_id,scores,audits


def gt_metrics(y,p):
    """y=observed, p=predicted; r0 uses y~k*p; prime uses p~kprime*y."""
    y,p=np.asarray(y,float),np.asarray(p,float)
    if np.var(y)==0 or np.var(p)==0 or np.dot(y,y)==0 or np.dot(p,p)==0:
        return {'defined':False,'reason':'zero_variance_or_norm'}
    r=float(np.corrcoef(y,p)[0,1]); r2=r*r
    k=float(np.dot(y,p)/np.dot(p,p)); kp=float(np.dot(y,p)/np.dot(y,y))
    r0=float(1-np.sum((y-k*p)**2)/np.sum((y-y.mean())**2))
    r0p=float(1-np.sum((p-kp*y)**2)/np.sum((p-p.mean())**2))
    d=(r2-r0)/r2 if r2>0 else None; dp=(r2-r0p)/r2 if r2>0 else None
    passed=r>0 and r2>0.6 and ((0.85<k<1.15 and abs(d)<0.1) or (0.85<kp<1.15 and abs(dp)<0.1)) and abs(r0-r0p)<0.3
    return {'defined':True,'pearson_r':r,'pearson_r_squared':r2,'k':k,'k_prime':kp,'r0_squared':r0,
            'r0_prime_squared':r0p,'relative_difference':d,'relative_difference_prime':dp,
            'absolute_r0_difference':abs(r0-r0p),'prespecified_diagnostic_pass':bool(passed)}


def corrected_rp(r, scrambled_r):
    """Roy cR^2_p = R * sqrt(R^2 - mean(Rr)^2); training correlations, not CV Q2."""
    delta=r*r-float(np.mean(scrambled_r))**2
    return float(r*np.sqrt(delta)) if delta>=0 else None


def bootstrap_r2(y,p,groups,n=1000,seed=2026):
    rng=np.random.default_rng(seed); unique=np.unique(groups); scores=[]
    members={g:np.flatnonzero(groups==g) for g in unique}
    for _ in range(n):
        idx=np.concatenate([members[g] for g in rng.choice(unique,size=len(unique),replace=True)])
        scores.append(metrics(y[idx],p[idx])['r2'])
    finite=[x for x in scores if x is not None and np.isfinite(x)]
    return scores,{'replicates':n,'valid_replicates':len(finite),'unit':'held-out scaffold clusters',
                   'ci95':np.quantile(finite,[0.025,0.975]).tolist() if finite else None,
                   'scope':'Fixed fitted model; conditional on held-out compounds, not training-set uncertainty'}


def leverage(train,query):
    A=np.column_stack([np.ones(len(train)),train]); B=np.column_stack([np.ones(len(query)),query])
    inverse=np.linalg.pinv(A.T@A,rcond=1e-12)
    return np.einsum('ij,jk,ik->i',B,inverse,B),3*A.shape[1]/len(A)
