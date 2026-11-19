from pathlib import Path
import sys,json,hashlib,importlib.util,importlib
import numpy as np,pandas as pd
from scipy.special import expit
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root/'external/Adapt-Env'))
from adaptenv import ThermodynamicConfig,ThermodynamicLandscape,build_landscape_config_from_preset
spec=importlib.util.spec_from_file_location('adaptenv_frozen175',root/'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-175/source/external/Adapt-Env/adaptenv/__init__.py',submodule_search_locations=[str(root/'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-175/source/external/Adapt-Env/adaptenv')]);frozen=importlib.util.module_from_spec(spec);sys.modules[spec.name]=frozen;spec.loader.exec_module(frozen);old=importlib.import_module('adaptenv_frozen175.thermodynamic')
bundle=root/'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-175/outputs';seeds=pd.read_csv(bundle/'seed_audit.csv');settings=json.loads((bundle/'selected_parameters.json').read_text());count=0;maximum=0.
for protein in ['GB1','TrpB']:
 codes=np.load(bundle/f'{protein}_empirical.npz')['genotypes'][::211]
 for row in seeds[seeds.protein.eq(protein)].itertuples():
  field=np.load(bundle/'latent_fields'/f'{protein}_{row.seed}.npz');a,g=field['a'][::211],field['g'][::211]
  for setting in [s for s in settings if s['protein']==protein]:
   cfg=ThermodynamicConfig(length=4,seed=int(row.seed),reference_rank=.5,representation_rank=8,order_weights=(.35,.45,.20),interaction_ratio=setting['ratio'],reference_occupancy=float(expit(setting['position'])),steepness=setting['steepness'],fitness_min=setting['offset'],fitness_max=setting['offset']+setting['scale'])
   m=ThermodynamicLandscape(cfg);aa,gg=m.latent_components(codes);np.testing.assert_allclose(aa,a,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(gg,g,rtol=1e-12,atol=1e-12)
   expected=setting['offset']+setting['scale']*expit(setting['position']+setting['steepness']*(a+setting['ratio']*g));actual=m.evaluate_batch(codes)
   baseline=old.ThermodynamicLandscape(old.ThermodynamicConfig(**vars(cfg)));np.testing.assert_array_equal(actual,baseline.evaluate_batch(codes));assert m.fingerprint==baseline.fingerprint
   np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-13);maximum=max(maximum,float(abs(actual-expected).max()));count+=len(codes)
   if protein=='TrpB' and setting['ratio']==.5 and setting['observation']=='shared_stop_sensitivity':
    preset=ThermodynamicLandscape(build_landscape_config_from_preset('trpb_run175_r05_shared_stop_intrinsic',seed=int(row.seed)));np.testing.assert_array_equal(preset.evaluate_batch(codes),actual)
print(json.dumps({'seed_fields_replayed':len(seeds),'settings':len(settings),'configuration_seed_pairs':24*6+24*18,'fitness_values_replayed':count,'maximum_absolute_difference':maximum,'stored_latent_values_agree_at_1e12_tolerance':True,'local_frozen_source_fitness_and_fingerprints_match_exactly':True,'run175_preset_replayed_on_all_24_trpb_seeds':True},indent=2))
