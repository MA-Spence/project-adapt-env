import sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from production_landscapes_exp066 import cube_statistics,sample_cubes,config_for,readout
from adaptenv import ThermodynamicConfig,ThermodynamicLandscape


def test_wht_and_slopes_have_known_values_and_unavailable_orders():
    bits=((np.arange(8)[:,None]>>np.arange(3))&1)
    y=(bits@np.array([1.,2.,4.]))[None,:]
    stat,mut,_=cube_statistics(y,pd.DataFrame({'cube':[0]}),[1e-10])
    assert abs(stat.epistatic_fraction.iloc[0])<1e-14
    assert abs(stat.slope.iloc[0])<1e-14
    assert np.isclose(stat.power_1.iloc[0],np.var(y))
    assert np.isnan(stat.power_4.iloc[0])
    parity=(-1.)**bits.sum(axis=1)
    stat,_,_=cube_statistics(parity[None,:],pd.DataFrame({'cube':[0]}),[1e-10])
    assert np.isclose(stat.ho_fraction.iloc[0],1)
    stat,_,_=cube_statistics(np.array([[0,1,1,0.]]),pd.DataFrame({'cube':[0]}),[1e-10])
    assert np.isnan(stat.ho_fraction.iloc[0])
    assert stat.reciprocal_sign_epistasis_fraction.iloc[0]==1


def test_sampling_is_outcome_blind_and_genotype_disjoint():
    codes=np.array(list(np.ndindex(5,5,5)),dtype=np.uint8)
    config=dict(calibration_cubes=3,confirmation_cubes=3,minimum_cubes_per_split=2,max_dimension=3,max_attempts_per_dimension=10000)
    a,_,dim,_=sample_cubes(codes,np.arange(len(codes)),np.zeros(3,dtype=np.uint8),False,config,np.random.default_rng(6))
    b,_,_,_=sample_cubes(codes,-np.arange(len(codes)),np.zeros(3,dtype=np.uint8),False,config,np.random.default_rng(6))
    np.testing.assert_array_equal(a,b)
    assert dim==3
    assert len(set(r.tobytes() for r in a.reshape(-1,3)))==a.shape[0]*8


def test_raw_unit_export_reproduces_public_model():
    base=dict(representation_rank=8,reference_rank=.5,order_weights=(.35,.45,.2))
    model=ThermodynamicLandscape(ThermodynamicConfig(length=7,seed=43,**base))
    nodes=np.random.default_rng(32).integers(20,size=(16,7));a,g=model.latent_components(nodes)
    p=[-.3,np.log(2.),-2.,np.log(.7)];origin=8.;scale=3.
    setting=dict(offset=origin+scale*p[0],scale=scale*np.exp(p[1]),position=p[2],steepness=np.exp(p[3]),interaction_ratio=4)
    public=ThermodynamicLandscape(config_for(setting,7,43,base))
    np.testing.assert_allclose(public.evaluate_batch(nodes),origin+scale*readout('r4_sigmoid',p,a,g),rtol=1e-12)


def test_plateaus_are_not_finite_epistasis_or_dr_evidence():
    stat,mut,_=cube_statistics(np.ones((1,16)),pd.DataFrame({'cube':[0]}),[1e-3])
    assert np.isnan(stat.epistatic_fraction.iloc[0]) and np.isnan(stat.slope.iloc[0])
    assert stat.slope_n.iloc[0]==0
    assert stat.neutral_edge_fraction.iloc[0]==1
    assert stat.strict_peak_fraction.iloc[0]==0


def test_publication_checks_content_and_preserves_previous_bundle(tmp_path):
    import hashlib,json,pytest
    from publish_production_landscapes import publish
    source=tmp_path/'runtime';source.mkdir();(source/'values.csv').write_text('fitness\n1\n')
    manifest=dict(mode='PRODUCTION',execution_errors=[],files={'values.csv':hashlib.sha256((source/'values.csv').read_bytes()).hexdigest()})
    (source/'manifest.json').write_text(json.dumps(manifest));target=tmp_path/'processed'/'RUN-TEST'
    publish(source,target);publish(source,target)
    assert (target/'values.csv').read_text()=='fitness\n1\n'
    (source/'values.csv').write_text('fitness\n2\n')
    with pytest.raises(ValueError,match='checksum'):publish(source,target)
    assert (target/'values.csv').read_text()=='fitness\n1\n'
