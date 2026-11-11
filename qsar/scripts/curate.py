"""Reproducible curation: no censoring substitution, no salt-weight assumptions."""
import gzip, json, math, hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupShuffleSplit
from fetch_data import assay_reason

ROOT=Path(__file__).resolve().parents[1]

def normalize_record(r, assays):
    if r['assay_chembl_id'] not in assays: return None,'assay_not_selected'
    if r.get('target_organism') != 'Staphylococcus aureus' or r.get('target_chembl_id') != 'CHEMBL352': return None,'target'
    if r.get('standard_type')!='MIC' or r.get('standard_relation')!='=': return None,'not_exact_MIC'
    if r.get('data_validity_comment'): return None,'validity_flag'
    if r.get('potential_duplicate'): return None,'chembl_potential_duplicate'
    if r.get('standard_flag') != 1: return None,'not_standardized'
    if not r.get('canonical_smiles'): return None,'missing_smiles'
    mol=Chem.MolFromSmiles(r['canonical_smiles'])
    if mol is None: return None,'invalid_smiles'
    if len(Chem.GetMolFrags(mol))!=1: return None,'salt_or_mixture'
    if any(a.GetAtomicNum() not in {1,6,7,8,9,15,16,17,35,53} for a in mol.GetAtoms()): return None,'unsupported_element'
    if not any(a.GetAtomicNum()==6 for a in mol.GetAtoms()): return None,'nonorganic'
    try: value=float(r['standard_value'])
    except (ValueError,TypeError): return None,'invalid_value'
    if not math.isfinite(value) or value<=0: return None,'nonpositive_or_nonfinite'
    mw=Descriptors.MolWt(mol)
    units=r.get('standard_units')
    factor={'nM':1e-9,'uM':1e-6,'mM':1e-3,'M':1}.get(units)
    if units=='ug.mL-1': molar=value*1e-3/mw
    elif factor: molar=value*factor
    else: return None,'unsupported_units'
    smiles=Chem.MolToSmiles(mol,isomericSmiles=True)
    return {'activity_id':r['activity_id'],'molecule_id':r['molecule_chembl_id'],
            'assay_id':r['assay_chembl_id'],'document_id':r['document_chembl_id'],
            'smiles':smiles,'molecular_weight':mw,'pMIC':-math.log10(molar),
            'original_value':value,'original_units':units},''

def main():
    cfg=json.loads((ROOT/'config.json').read_text())
    raw=ROOT/'data/raw'; out=ROOT/'data/curated'; out.mkdir(parents=True,exist_ok=True)
    provenance=json.loads((raw/'provenance.json').read_text())
    for name,sha in provenance['sha256'].items():
        if hashlib.sha256((raw/name).read_bytes()).hexdigest()!=sha: raise ValueError('Raw snapshot hash mismatch')
    all_assays=json.loads(gzip.decompress((raw/'assays.json.gz').read_bytes()))
    assays={a['assay_chembl_id']:a for a in all_assays if not assay_reason(a)}
    records=json.loads(gzip.decompress((raw/'activities.json.gz').read_bytes()))
    kept,rejected=[],[]
    for record in records:
        row,reason=normalize_record(record,assays)
        if reason: rejected.append({'activity_id':record['activity_id'],'reason':reason})
        else: kept.append(row)
    df=pd.DataFrame(kept)
    aggregates=[]
    for smiles, group in df.groupby('smiles',sort=True):
        spread=float(group.pMIC.max()-group.pMIC.min())
        if spread>cfg['replicate_max_range_log10']:
            rejected.extend({'activity_id':int(i),'reason':'replicate_disagreement_over_1_log'} for i in group.activity_id)
            continue
        mol=Chem.MolFromSmiles(smiles)
        # Non-stereochemical scaffold groups keep stereo analogues together.
        scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=mol,includeChirality=False) or 'ACYCLIC'
        aggregates.append({'compound_id':hashlib.sha256(smiles.encode()).hexdigest()[:16],
                           'smiles':smiles,'scaffold':scaffold,'pMIC':float(group.pMIC.median()),
                           'replicate_n':len(group),'replicate_range':spread,
                           'molecule_ids':';'.join(sorted(set(group.molecule_id))),
                           'activity_ids':';'.join(map(str,sorted(group.activity_id))),
                           'assay_ids':';'.join(sorted(set(group.assay_id))),
                           'document_ids':';'.join(sorted(set(group.document_id)))})
    data=pd.DataFrame(aggregates)
    if len(data)<100 or data.scaffold.nunique()<25: raise ValueError('Insufficient compounds/scaffolds for planned validation')
    splitter=GroupShuffleSplit(n_splits=1,test_size=cfg['test_scaffold_fraction'],random_state=cfg['seed'])
    train,test=next(splitter.split(data,groups=data.scaffold))
    data['split']='train'; data.loc[test,'split']='test'
    if len(test)<20 or data.iloc[train].scaffold.nunique()<10: raise ValueError('Insufficient holdout or CV support')
    assert not set(data.iloc[train].scaffold)&set(data.iloc[test].scaffold)
    data.to_csv(out/'compounds.csv',index=False)
    df.to_csv(out/'accepted_measurements_before_replicate_filter.csv',index=False)
    pd.DataFrame(rejected).to_csv(out/'rejected_measurements.csv',index=False)
    pd.DataFrame([{k:a[k] for k in ['assay_chembl_id','document_chembl_id','assay_strain','description']} for a in assays.values()]).to_csv(out/'included_assays.csv',index=False)
    summary={'raw_activities':len(records),'accepted_before_aggregation':len(kept),'curated_compounds':len(data),
             'train_n':len(train),'test_n':len(test),'scaffolds':data.scaffold.nunique(),
             'train_scaffolds':data.iloc[train].scaffold.nunique(),'test_scaffolds':data.iloc[test].scaffold.nunique(),
             'rejection_counts':pd.Series([r['reason'] for r in rejected]).value_counts().to_dict(),
             'endpoint':cfg['endpoint'],'pMIC_range':[float(data.pMIC.min()),float(data.pMIC.max())],
             'assay_selection_not_full_protocol_harmonization':True}
    (out/'curation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__': main()
