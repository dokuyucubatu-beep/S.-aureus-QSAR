import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from modeling import DescriptorSelector,metrics,gt_metrics,corrected_rp,bootstrap_r2,leverage
from curate import normalize_record

class MetricsTests(unittest.TestCase):
    def test_r2_is_not_squared_correlation(self):
        y=np.arange(1.,11.); p=y+10
        self.assertLess(metrics(y,p)['r2'],0)
        self.assertAlmostEqual(gt_metrics(y,p)['pearson_r_squared'],1)
        self.assertFalse(gt_metrics(y,p)['prespecified_diagnostic_pass'])

    def test_gt_perfect_and_degenerate(self):
        y=np.arange(1.,11.)
        report=gt_metrics(y,y)
        for k in ['r0_squared','r0_prime_squared','k','k_prime']: self.assertAlmostEqual(report[k],1)
        self.assertTrue(report['prespecified_diagnostic_pass'])
        self.assertFalse(gt_metrics(y,np.ones(10))['defined'])

    def test_crp_definition(self):
        self.assertAlmostEqual(corrected_rp(.8,[.1,.3]),.8*np.sqrt(.64-.04))
        self.assertIsNone(corrected_rp(.1,[.9]))

    def test_cluster_bootstrap(self):
        y=np.arange(1.,21.); groups=np.repeat(np.arange(10),2)
        a,ci=bootstrap_r2(y,y,groups,n=1000)
        self.assertEqual(len(a),1000)
        self.assertEqual(ci['valid_replicates'],1000)
        np.testing.assert_allclose(ci['ci95'],[1,1])

    def test_selector_train_only_and_collinearity(self):
        rng=np.random.default_rng(4); x=rng.normal(size=100)
        X=np.column_stack([x,x*2,np.ones(100),rng.normal(size=100),np.full(100,np.nan)])
        s=DescriptorSelector(max_features=2).fit(X,x)
        self.assertNotIn(1,s.selected_indices_)
        self.assertNotIn(2,s.selected_indices_)
        self.assertNotIn(4,s.selected_indices_)
        before=s.medians_.copy()
        s.transform(np.full((3,5),1e9))
        np.testing.assert_array_equal(before,s.medians_)

    def test_leverage_extrapolation(self):
        X=np.arange(20.).reshape(-1,1)
        h,limit=leverage(X,np.array([[1000.]]))
        self.assertGreater(h[0],limit)

    def test_units_censoring_and_salts(self):
        record={'assay_chembl_id':'A','target_organism':'Staphylococcus aureus','target_chembl_id':'CHEMBL352',
                'standard_type':'MIC','standard_relation':'=','standard_flag':1,'canonical_smiles':'CCO',
                'standard_value':'1','standard_units':'uM','activity_id':1,'molecule_chembl_id':'M','document_chembl_id':'D'}
        a,reason=normalize_record(record,{'A':{}})
        self.assertEqual(reason,''); self.assertAlmostEqual(a['pMIC'],6)
        mass=dict(record,standard_value=str(a['molecular_weight']/1000),standard_units='ug.mL-1')
        b,_=normalize_record(mass,{'A':{}}); self.assertAlmostEqual(b['pMIC'],6)
        self.assertEqual(normalize_record(dict(record,standard_relation='>'),{'A':{}})[1],'not_exact_MIC')
        self.assertEqual(normalize_record(dict(record,canonical_smiles='CCO.[Na+]'),{'A':{}})[1],'salt_or_mixture')

if __name__=='__main__': unittest.main()
