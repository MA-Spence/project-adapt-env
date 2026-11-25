"""Known-answer scientific controls for the paired interaction ablation."""
import sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from interaction_ablation_metrics_exp065 import H, replicate_metrics, summarize_replicates, static_graph
from interaction_ablation_exp065 import fit_arm


def metrics(a,b,q=0.):
    return replicate_metrics(np.atleast_2d(a),np.atleast_2d(b),q,[-100,-100],[0.])


def test_known_pairwise_power_and_parseval():
    y=2*H[1]+3*H[3]+1
    c,f=metrics(y,y)
    assert c.mean_power_1.item()==4
    assert c.cross_power_2.item()==9
    assert c.cross_power_3.item()==0
    assert c.cross_power_4.item()==0
    assert c.filter(regex='^mean_power').sum(axis=1).item()==pytest.approx(y.var())


def test_sigmoid_additive_has_higher_orders_but_no_sign_reversals():
    x=H[1]+2*H[2]+3*H[4]+4*H[8]+1
    y=1/(1+np.exp(-x))
    c,f=metrics(y,y)
    assert c.mean_power_3.item()+c.mean_power_4.item()>0
    assert not f.mean_sign.any()


def test_cross_power_removes_balanced_independent_errors_but_not_shared():
    base=2*H[1]
    a=np.stack([base+u*H[3] for u in [-1,-1,1,1]])
    b=np.stack([base+v*H[3] for v in [-1,1,-1,1]])
    c,_=metrics(a,b)
    assert c.cross_power_2.mean()==0
    assert c.mean_power_2.mean()>0
    shared=3*H[7]
    c,_=metrics(a+shared,b+shared)
    assert c.cross_power_3.mean()==9


def test_opposite_orientations_are_not_reproducible_sign_epistasis():
    c,f=metrics(H[3],-H[3])
    assert f.replicate1_sign.any() and f.replicate2_sign.any()
    assert not f.reproducible_sign.any()
    assert c.cross_power_2.item()==-1


def test_signed_stop_corrections_and_empty_strata_are_preserved():
    c,f=metrics(H[1],H[1],q=100.)
    s=summarize_replicates(c,f)
    assert c.stop_adjusted_power_2.item()<0
    empty=s[s.stratum=='no_active'].iloc[0]
    assert empty.n_cubes==0 and empty.n_faces==0
    assert np.isnan(empty.cross_power_2)


@pytest.mark.parametrize('epsilon',[0.,.1,1.])
def test_graph_fiber_calculation_matches_bruteforce(epsilon):
    q,dim=3,3
    x=np.random.default_rng(92).normal(size=q**dim)
    valid=np.ones(len(x),bool);valid[[0,1,2,5,8]]=False
    codes=np.array(np.unravel_index(np.arange(len(x)),(q,)*dim)).T
    pairs=[(i,j) for i in range(len(x)) for j in range(i+1,len(x))
           if valid[i] and valid[j] and np.sum(codes[i]!=codes[j])==1]
    a,b=np.array([(x[i],x[j]) for i,j in pairs]).T
    u,v=np.r_[a,b],np.r_[b,a]
    traps=[not any(x[j]>x[i]+epsilon for j in range(len(x))
        if valid[j] and np.sum(codes[i]!=codes[j])==1) for i in np.flatnonzero(valid)]
    got=static_graph(x,valid,epsilon,q,dim)
    assert got['edges']==len(pairs)
    assert got['neighbor_correlation']==pytest.approx(np.corrcoef(u,v)[0,1])
    assert got['edge_rms_native']==pytest.approx(np.sqrt(np.mean((a-b)**2)))
    assert got['no_beneficial_fraction']==np.mean(traps)


def test_optimizer_refits_null_and_nests_it_in_full(tmp_path):
    names=['offset','scale','position','steepness']
    inputs=SimpleNamespace(config={'optimizer':{'quick_trials':24,'trials':24,'seed':17,'n_startup_trials':2},
        'ratio_bounds':[0.,8.],'fit':{}},fitting=SimpleNamespace(PARAMETER_NAMES=names,
        parameter_bounds=lambda _: [[0.,8.]]*4),origin=0.,scale=1.,floor=0.,
        historical=[{'ratio':r,'parameters':[1.,2.,3.,4.]} for r in [0.,.5,1.,2.,4.,8.]])
    inputs.loss=lambda r,p: {'total':float((r-1.5)**2+np.sum((np.array(p)-2)**2))}
    null=fit_arm(inputs,'r0',tmp_path,quick=True)
    full=fit_arm(inputs,'full',tmp_path,quick=True,null_setting=null)
    t0=pd.read_csv(tmp_path/'r0_trials.csv.gz');t1=pd.read_csv(tmp_path/'full_trials.csv.gz')
    assert len(t0)==len(t1)==24
    assert (t0.r==0).all() and null['ratio']==0
    assert full['objective']['total']<=null['objective']['total']
    assert (t1.candidate=='paired_refitted_r0').sum()==1
    assert ((t1.r>0)&(t1.r<8)&~t1.r.isin([.5,1,2,4])).any()
    import optuna
    s=optuna.load_study(study_name='TrpB_r0_shared_stop',storage=f'sqlite:///{tmp_path}/r0.sqlite3')
    assert all(set(t.params)==set(names) for t in s.trials)
    with pytest.raises(ValueError,match='overwrite'):
        fit_arm(inputs,'r0',tmp_path,quick=True)
