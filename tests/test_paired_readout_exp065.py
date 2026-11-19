import importlib.util
import itertools
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("paired_readout", ROOT / "scripts/paired_readout_exp065.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def test_common_affine_rule_matches_training_scale_and_freezes_for_test():
    raw = np.arange(48, dtype=float).reshape(3, 16) / 11
    empirical = 7 + 0.2 * raw
    p = runner.affine_calibration(raw, empirical)
    transformed = p["offset"] + p["scale"] * raw
    np.testing.assert_allclose(transformed.mean(), empirical.mean())
    np.testing.assert_allclose(transformed.std(axis=1).mean(), empirical.std(axis=1).mean())
    test = p["offset"] + p["scale"] * (2 * raw)
    assert test.std(axis=1).mean() == pytest.approx(2 * empirical.std(axis=1).mean())


def test_affine_rescaling_preserves_resolved_fraction_and_slope_but_scales_power():
    x = runner.BITS.astype(float)
    values = x @ np.array([1., 2., 3., 4.]) + 2 * x[:, 0] * x[:, 1]
    original = runner.describe(values, 4, 1e-12)
    mapped = runner.describe(9 + .13 * values, 4, 1e-12)
    for name in ['walsh_epistatic_fraction', 'walsh_high_order_fraction', 'beneficial_diminishing_returns_slope']:
        assert original[name] == pytest.approx(mapped[name], abs=1e-12)
    for order in range(1,5):
        assert mapped[f'walsh_energy_order_{order}'] == pytest.approx(.13**2 * original[f'walsh_energy_order_{order}'], abs=1e-12)


def test_linear_sigmoid_use_identical_latents_and_preserve_edge_signs():
    rng = np.random.default_rng(78)
    a, g = rng.normal(size=(2,16))
    p = {'ratio':1., 'beta':.1, 'p0':.5, 'scale':.3, 'offset':-.1}
    linear = runner.observe(a,g,{**p, 'readout':'linear'})
    nonlinear = runner.observe(a,g,{**p, 'readout':'sigmoid'})
    for i in range(16):
        for bit in range(4):
            j = i ^ (1 << bit)
            assert np.sign(linear[j]-linear[i]) == np.sign(nonlinear[j]-nonlinear[i])


def test_embedded_cube_holds_inactive_sites_fixed_and_orders_bits():
    anchor = np.arange(12) % 20
    alt = (anchor + 1) % 20
    sites = [2,4,8,11]
    nodes = runner.cube_nodes(anchor,alt,sites)
    np.testing.assert_array_equal(nodes[0],anchor)
    np.testing.assert_array_equal(nodes[-1,sites],alt[sites])
    for bit,site in enumerate(sites):
        assert np.flatnonzero(nodes[1 << bit] != anchor).tolist() == [site]
    assert len({tuple(n) for n in nodes}) == 16


def test_stable_rms_retains_tiny_effects():
    values=np.array([1e-250, -1e-250])
    assert np.sqrt(np.mean(values**2)) == 0
    assert runner.stable_rms(values) == pytest.approx(1e-250, rel=1e-12, abs=0)


def test_fresh_cubes_exclude_prior_genotypes_and_each_other(tmp_path):
    states=list(itertools.product(range(5),repeat=4))
    frame=pd.DataFrame(states,columns=['pos1','pos2','pos3','pos4'])
    frame['fitness']=np.arange(len(frame),dtype=float)
    path=tmp_path/'variants.csv'; frame.to_csv(path,index=False)
    prior={tuple(map(str,s)) for s in states[:16]}
    values,records,_=runner.fresh_empirical(path,prior,count=4,seed=123)
    nodes=[tuple(n) for row in records for n in row['node_keys']]
    assert values.shape == (4,16)
    assert len(set(nodes))==64
    assert not prior.intersection(nodes)


def test_fitting_never_uses_confirmation_panel():
    cfg={'measurement':{'resolution':.001},'fit':{'interaction_ratios':[0.,1.],
         'steepness':[.125], 'reference_occupancy':[.5]}}
    rng=np.random.default_rng(391)
    panels=[{'a':rng.normal(size=16),'g':rng.normal(size=16)} for _ in range(4)]
    empirical=rng.normal(size=(4,16))
    selected,_,_=runner.fit_pairs(cfg,panels,empirical)
    assert len(selected)==4
    for p in selected:
        raw=np.array([runner.raw_readout(row['a'],row['g'],p) for row in panels])
        fitted=p['offset']+p['scale']*raw
        assert fitted.std(axis=1).mean()==pytest.approx(empirical.std(axis=1).mean())
        assert fitted.mean()==pytest.approx(empirical.mean())
