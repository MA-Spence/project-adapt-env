"""Numerical and leakage controls for observation-aware calibration."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from observation_calibration_exp065 import observation_bank,observed,readout,target_from,loss_components,empirical_panel
from replicate_observation_exp065 import observe

def toy_noise():
    return dict(fitness_breaks=np.array([0,.05,.3]),depth_breaks=np.array([20,40,60]),
        depths={f:np.array([10,30,50,70]) for f in range(4)},
        donors={(f,d):np.array([-.01,.01])*(f+d+1) for f in range(4) for d in range(4)},
        common=[np.array([-.02,.02])*(d+1) for d in range(4)],q=np.arange(4)*.0004)

def test_cached_observation_matches_original_at_bin_boundaries():
    noise=toy_noise();u=np.random.default_rng(2).random((3,2,16,6))
    raw=np.resize(np.array([-.1,0,.01,.05,.1,.3,.7]),(2,16))
    for shared in [False,True]:
        a=observe(raw,noise,u,shared);b=observed(raw,observation_bank(noise,u,shared))
        for x,y in zip(a,b):np.testing.assert_allclose(x,y,atol=3e-16)

def test_objective_does_not_see_epistatic_assignment():
    y=np.random.default_rng(3).normal(size=(4,16))
    cfg={'fit':{'quantiles':101},'measurement':{'activity_cutoffs':[.05,.05]}}
    target=target_from(y,cfg);shuffled=np.array([row[np.random.default_rng(i).permutation(16)] for i,row in enumerate(y)])
    other=target_from(shuffled,cfg)
    for k in target:np.testing.assert_allclose(target[k],other[k],atol=1e-14)

def test_sigmoid_export_matches_public_generator():
    from adaptenv import ThermodynamicConfig,ThermodynamicLandscape
    from scipy.special import expit
    m=ThermodynamicLandscape(ThermodynamicConfig(length=4,seed=652301,reference_rank=.5))
    nodes=np.random.default_rng(5).integers(20,size=(16,4));a,g=m.latent_components(nodes)
    params=np.array([-.02,np.log(.7),-3,np.log(1.2)])
    for arm,r in [('r4_sigmoid',4),('r0_sigmoid',0)]:
        model=ThermodynamicLandscape(ThermodynamicConfig(length=4,seed=652301,reference_rank=.5,
            interaction_ratio=r,reference_occupancy=expit(params[2]),steepness=np.exp(params[3]),
            fitness_min=params[0],fitness_max=params[0]+np.exp(params[1])))
        np.testing.assert_allclose(readout(arm,params,a,g),model.evaluate_batch(nodes),atol=1e-14)

def test_panel_disjointness_and_exclusion():
    from itertools import product
    genotypes=[''.join(g) for g in product('ABCD',repeat=4)]
    table=pd.DataFrame(dict(genotype=genotypes,stop_count=0,replicate_1=1.,replicate_2=2.,input_count=10,published_mean=1.5))
    cfg={'sampling':{'empirical_seed':31,'calibration_cubes':1,'confirmation_cubes':1}}
    nodes,used=empirical_panel(table,{'AAAA'},cfg)
    assert len(nodes)==32 and nodes.genotype.nunique()==32
    assert 'AAAA' not in set(nodes.genotype)
    assert set(nodes.genotype).issubset(used)
