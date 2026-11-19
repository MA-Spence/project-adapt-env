#!/usr/bin/env python3
"""Freeze original TrpB replicate/control inputs and the previously fitted RUN-166 panel."""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
BUNDLE=ROOT/'data/processed/trpb-replicate-observation-exp065'
PRIOR=ROOT/'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-166/outputs'

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def prepare():
    provenance=json.loads((BUNDLE/'source_manifest.json').read_text())
    for f in provenance['files']:assert digest(BUNDLE/f['file'])==f['sha256']
    original=pd.read_csv(BUNDLE/'source/four-site_summary_AA_data.csv')
    assert not original.AAs.duplicated().any()
    out=original[['AAs','mu_1-bg/max','mu_2-bg/max','avg_mu-bg/max','avg_InputCount','# Stop','active']].rename(columns={
        'AAs':'genotype','mu_1-bg/max':'replicate_1','mu_2-bg/max':'replicate_2','avg_mu-bg/max':'published_mean',
        'avg_InputCount':'input_count','# Stop':'stop_count','active':'published_active'})
    np.testing.assert_allclose(out[['replicate_1','replicate_2']].mean(axis=1),out.published_mean,rtol=0,atol=1e-14)
    assert np.isfinite(out[['replicate_1','replicate_2','input_count']]).all().all()
    assert (out.input_count>=10).all()
    prior=pd.read_csv(PRIOR/'empirical_cubes.csv');lookup=out.set_index('genotype')
    largest=0.;nodes=set()
    for _,row in prior.iterrows():
        keys=[''.join(n) for n in json.loads(row.node_keys)];nodes.update(keys)
        error=np.max(abs(lookup.loc[keys,'published_mean'].to_numpy()-json.loads(row.raw_values)))
        largest=max(largest,float(error));assert error<1e-7
    out.to_csv(BUNDLE/'replicates.csv',index=False)
    (BUNDLE/'prior').mkdir(exist_ok=True)
    for name in ['empirical_cubes.csv','selected_parameters.json','manifest.json','cube_values.csv']:
        shutil.copyfile(PRIOR/name,BUNDLE/'prior'/name)
    earlier=ROOT/'experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/inputs/run164_trpb_nodes.json'
    previous=json.loads(earlier.read_text())
    nodes.update(''.join(x) for x in previous['node_keys']+previous['engineering_node_keys'])
    (BUNDLE/'excluded_donor_genotypes.json').write_text(json.dumps(sorted(nodes),indent=2)+'\n')
    audit={'source_rows':len(out),'nonstop_rows':int((out.stop_count==0).sum()),'stop_rows':int((out.stop_count>0).sum()),
           'run166_genotypes_verified':768,'maximum_published_vs_imported_fitness_error':largest,
           'noise_donor_exclusions':len(nodes),'shared_baseline':'Both published replicates use avg_InputCount; per-timepoint stop background and maximum normalization precede averaging across 16,20,24,28,36 hours.',
           'source_code_executed':False,'prior_run':'RUN-166','scope':'Diagnostic reanalysis of exposed prior cubes, not fresh empirical confirmation.'}
    files=['replicates.csv','excluded_donor_genotypes.json','source_manifest.json']+['prior/'+x for x in ['empirical_cubes.csv','selected_parameters.json','manifest.json','cube_values.csv']]
    audit['input_sha256']={name:digest(BUNDLE/name) for name in files}
    (BUNDLE/'input_manifest.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='input_sha256'},indent=2))
if __name__=='__main__':prepare()
