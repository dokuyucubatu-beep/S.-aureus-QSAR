"""Predict single-component organic SMILES with the trusted bundled QSAR model."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import joblib
import deepchem as dc
from rdkit import Chem
from rdkit.Chem import Descriptors
from modeling import leverage
from train import fingerprint,nearest_similarity

ROOT=Path(__file__).resolve().parents[1]

def predict(input_file,output_file,model_dir):
    model_dir=Path(model_dir)
    manifest=json.loads((model_dir/'manifest.json').read_text())
    if hashlib.sha256((model_dir/'model.joblib').read_bytes()).hexdigest()!=manifest['output_sha256']['model.joblib']:
        raise ValueError('Model hash mismatch')
    bundle=joblib.load(model_dir/'model.joblib') # Only load this locally generated/trusted artifact.
    df=pd.read_csv(input_file)
    if 'smiles' not in df: raise ValueError('Input CSV needs a smiles column')
    canonical=[]; weights=[]
    for i,s in enumerate(df.smiles):
        mol=Chem.MolFromSmiles(str(s))
        if mol is None or len(Chem.GetMolFrags(mol))!=1: raise ValueError(f'Invalid or multicomponent SMILES at row {i+1}')
        if any(a.GetAtomicNum() not in {1,6,7,8,9,15,16,17,35,53} for a in mol.GetAtoms()) or not any(a.GetAtomicNum()==6 for a in mol.GetAtoms()):
            raise ValueError(f'Unsupported chemistry at row {i+1}')
        canonical.append(Chem.MolToSmiles(mol,isomericSmiles=True)); weights.append(Descriptors.MolWt(mol))
    feat=dc.feat.RDKitDescriptors(descriptors=bundle['descriptor_names'],is_normalized=False,ipc_avg=True)
    if list(feat.reqd_properties)!=bundle['descriptor_names']: raise ValueError('Descriptor implementation/order changed')
    X=np.asarray(feat.featurize(canonical),float)
    model=bundle['pipeline']; pred=model.predict(X); Z=model[:-1].transform(X)
    h,cutoff=leverage(bundle['training_scaled'],Z)
    sims=nearest_similarity(fingerprint(canonical),fingerprint(bundle['training_smiles']))
    df['canonical_smiles']=canonical; df['predicted_pMIC']=pred
    df['predicted_MIC_uM']=1e6*10**(-pred); df['predicted_MIC_ug_mL']=1000*np.asarray(weights)*10**(-pred)
    df['leverage']=h; df['nearest_train_tanimoto']=sims
    df['in_AD']=(h<=cutoff)&(sims>=bundle['similarity_threshold'])
    df['seen_in_training']=[s in set(bundle['training_smiles']) for s in canonical]
    df.to_csv(output_file,index=False)
    print(df.to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True,type=Path); p.add_argument('--output',required=True,type=Path)
    p.add_argument('--model-dir',type=Path,default=ROOT/'results/reference')
    a=p.parse_args(); predict(a.input,a.output,a.model_dir)
