"""Controls for the paired sweep: API identity, topology and independent graph checks."""
from pathlib import Path
import sys
import numpy as np
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import interaction_strength_exp065 as sweep
import landscape_graph as graph

def test_public_readout_at_all_sweep_strengths():
    codes=np.random.default_rng(1).integers(0,20,(250,4))
    base=dict(length=4,seed=9812,reference_rank=.5,representation_rank=8,order_weights=(.35,.45,.20))
    m=sweep.ThermodynamicLandscape(sweep.ThermodynamicConfig(**base));a,g=m.latent_components(codes)
    p=np.array([-.1,np.log(1.7),-.4,np.log(.6)])
    for r in [0,.5,1,2,4,8]:
        model=sweep.ThermodynamicLandscape(sweep.ThermodynamicConfig(**base,interaction_ratio=r,reference_occupancy=expit(p[2]),steepness=np.exp(p[3]),fitness_min=p[0],fitness_max=p[0]+np.exp(p[1])))
        np.testing.assert_allclose(sweep.predict(p,a,g,r),model.evaluate_batch(codes),rtol=1e-12,atol=1e-12)

def test_exact_walks_invariant_to_strict_monotone_readout():
    f=np.random.default_rng(74).normal(size=27);valid=np.ones(27,dtype=bool);valid[4]=False
    transformed=-2+3*expit(.5*f)
    a,_,_=graph.graph(f,valid,0,np.quantile(f[valid],.9),f.min()-1,q=3,dim=3)
    b,_,_=graph.graph(transformed,valid,0,np.quantile(transformed[valid],.9),transformed.min()-1,q=3,dim=3)
    np.testing.assert_allclose(a[valid,:16],b[valid,:16],atol=1e-12)

def test_graph_against_independent_small_graph_algorithms():
    assert graph.validate_kernel()['small_graph_cases']==20
