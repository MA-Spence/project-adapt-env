"""Numerical controls for measurement corrections, not scientific success gates."""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

spec = importlib.util.spec_from_file_location('replicate_observation_exp065', Path(__file__).parents[1]/'scripts/replicate_observation_exp065.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_exact_spectral_components_and_noiseless_slopes():
    # A known additive + pair + fourth-order landscape in Walsh coordinates.
    c = np.zeros(16); c[0] = 3; c[1] = 2; c[2] = .8; c[3] = .6; c[15] = .4
    f = (c @ m.H)[None,:]
    comp=m.pair_components(f,f,np.zeros_like(f))
    expected=np.array([4.64,.36,0,.16])
    for j in range(3):
        np.testing.assert_allclose(comp['power'][0,j],expected,atol=1e-14)
        for k in range(4):
            x=f[0,m.LOW[k]]; y=f[0,m.HIGH[k]]-x
            np.testing.assert_allclose(comp['numerator'][0,j,k]/comp['denominator'][0,j,k], np.polyfit(x,y,1)[0],atol=1e-13)


def test_heteroscedastic_shared_noise_corrections_recover_signal_in_expectation():
    rng=np.random.default_rng(72904)
    n=60000
    f=(2*m.H[1]+.4*m.H[3]+.2*m.H[7])[None,:]
    q=np.linspace(.1,.8,16)
    shared=rng.normal(size=(n,16))*np.sqrt(q)
    a=f+shared+rng.normal(size=(n,16))*.3
    b=f+shared+rng.normal(size=(n,16))*.5
    out=m.pair_components(a,b,q)
    truth=m.pair_components(f,f,0)
    np.testing.assert_allclose(out['power'][:,2].mean(0),truth['power'][0,0],atol=.008)
    np.testing.assert_allclose(out['numerator'][:,2].mean(0),truth['numerator'][0,0],atol=.012)
    np.testing.assert_allclose(out['denominator'][:,2].mean(0),truth['denominator'][0,0],atol=.012)
    # Differential-error subtraction alone cannot remove shared covariance.
    assert (out['power'][:,1].mean(0)-truth['power'][0,0]).sum() > .3


def test_negative_power_and_unresolved_slopes_are_not_clipped_or_imputed():
    f=np.ones((1,16))
    comp=m.pair_components(f,f,.4)
    frame, mutations=m.statistics(comp,pd.DataFrame({'cube_id':['x']}),[.001],mutations=True)
    assert frame.stop_adjusted_total_power.iloc[0] < 0
    assert np.isnan(frame.stop_adjusted_ho_fraction.iloc[0])
    assert frame.stop_adjusted_slope_n.iloc[0] == 0
    assert mutations.slope.isna().all()
    assert (mutations.loc[mutations.method=='stop_adjusted','denominator']<0).all()


def test_symmetry_and_affine_units():
    rng=np.random.default_rng(119)
    a,b=rng.normal(size=(2,7,16)); q=np.linspace(.01,.03,16)
    ab=m.pair_components(a,b,q); ba=m.pair_components(b,a,q)
    trans=m.pair_components(3*a+9,3*b+9,9*q)
    for key in ('power','numerator','denominator'):
        np.testing.assert_allclose(ab[key],ba[key],atol=1e-13)
        np.testing.assert_allclose(trans[key],9*ab[key],atol=1e-12)


def test_activity_is_two_replicate_definition():
    a=np.ones((3,16)); b=a.copy(); a[0]=0; b[1,:8]=0
    count,label=m.activity(a,b,[.1,.2])
    assert count.tolist()==[0,8,16]
    assert label.tolist()==['no_active','some_active','all_active']


def test_observation_pairing_variances_and_shared_component():
    rng=np.random.default_rng(514)
    # Exact symmetric empirical distributions with known variance.
    noise=dict(fitness_breaks=np.array([0,.1,.5]),depth_breaks=np.array([10,20,30]),
               donors={(f,d):np.array([-2.,2.]) for f in range(4) for d in range(4)},
               depths={f:np.array([5,15,25,35]) for f in range(4)},
               common=[np.array([-3.,3.]) for _ in range(4)],q=np.array([9.]*4))
    uniforms=rng.random((10000,16,6))
    latent=np.full((16,),.2)
    a,b,depth,q=m.observe(latent,noise,uniforms)
    c,d,depth2,q2=m.observe(latent,noise,uniforms,shared=True)
    np.testing.assert_allclose(c-a,d-b,atol=1e-14)
    np.testing.assert_array_equal(depth,depth2)
    np.testing.assert_array_equal(q,q2)
    assert abs(np.var(a-latent)-2)<.025
    assert abs(np.var(a-b)-4)<.05
    assert abs(np.mean((c-.2)*(d-.2))-9)<.07
    assert (a<0).any()  # Never clip negative observed fitness.


def test_noise_training_excludes_exposed_variants_and_holdout(tmp_path):
    n=4000; rng=np.random.default_rng(99)
    mean=np.tile(np.array([-.1,.02,.2,.8]),n//4)
    a=mean+rng.normal(0,.01,n); b=mean+rng.normal(0,.01,n)
    table=pd.DataFrame(dict(genotype=[f'g{i}' for i in range(n)],replicate_1=a,replicate_2=b,
                            published_mean=(a+b)/2,input_count=rng.integers(10,1000,n),
                            stop_count=np.arange(n)%3==0,published_active=mean>.054))
    config=dict(sampling=dict(donor_split_seed=123,training_fraction=.8),
                observation=dict(count_quantiles=[.25,.5,.75],minimum_donors_per_cell=30),
                measurement=dict(activity_cutoffs=[.054,.053]))
    excluded=['g1','g2','g4']
    noise,assign=m.calibrate_noise(table,excluded,config,tmp_path)
    assert set(assign.loc[assign.split=='training','genotype']).isdisjoint(excluded)
    assert (assign.loc[assign.genotype.isin(excluded),'split']=='excluded').all()
    # Arbitrarily change held-out values while preserving the original split.
    changed=table.copy()
    ids=set(assign.loc[assign.split=='holdout','genotype'])
    changed.loc[changed.genotype.isin(ids),['replicate_1','replicate_2','published_mean']]+=100
    again=tmp_path/'again'; again.mkdir()
    noise2,_=m.calibrate_noise(changed,excluded,config,again)
    for name in ('q','fitness_breaks','depth_breaks'):
        np.testing.assert_array_equal(noise[name],noise2[name])
    for key in noise['donors']:
        np.testing.assert_array_equal(noise['donors'][key],noise2['donors'][key])


def test_aggregate_power_ratio_uses_sums_and_retains_empty_strata():
    rng=np.random.default_rng(932)
    f=rng.normal(size=(2,16)); meta=pd.DataFrame({'split':['confirmation']*2,'activity_stratum':['no_active']*2})
    stats,_=m.statistics(m.pair_components(f,f,.05),meta,[.001])
    agg=m.aggregate_table(stats,['split','resolution'])
    assert len(agg)==12
    assert (agg.loc[agg.stratum=='all_active','n_cubes']==0).all()
    for method in m.METHODS:
        row=agg[(agg.stratum=='all')&(agg.method==method)].iloc[0]
        expected=sum(stats[f'{method}_power_{k}'].sum() for k in (3,4))/stats[f'{method}_total_power'].sum()
        assert row.aggregate_ho_fraction == pytest.approx(expected)
