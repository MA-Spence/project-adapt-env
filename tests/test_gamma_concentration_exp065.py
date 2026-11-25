"""Controls for the new concentration coordinate and matched scientific design."""
import copy
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import gamma_concentration_exp065 as gamma


def test_joint_optimizer_samples_all_six_coordinates_and_logs_alpha(tmp_path):
    names=['p0','p1','p2','p3']
    old=[dict(ratio=r,parameters=[1.,2.,3.,4.]) for r in [0.,.5,1.,2.,4.,8.]]
    def loss(r,p,a):return {'total':float((np.log(a)-np.log(.3))**2+(r-1.7)**2+np.sum((np.array(p)-2.3)**2))}
    null=dict(arm='r0',**old[0],objective=loss(0,old[0]['parameters'],1.))
    full=dict(arm='full',**old[1],objective=loss(.5,old[1]['parameters'],1.))
    x=SimpleNamespace(config=dict(optimizer=dict(trials=32,quick_trials=32,seed=18,n_startup_trials=10),fit={},
        ratio_bounds=[0.,8.],concentration_bounds=[.05,100.]),
        fitting=SimpleNamespace(PARAMETER_NAMES=names,parameter_bounds=lambda _: [[0.,8.]]*4),
        historical=old,controls=[null,full],origin=0.,scale=1.,floor=0.,gamma_loss=loss)
    result=gamma.fit_gamma(x,tmp_path,quick=True)
    trace=pd.read_csv(tmp_path/'gamma_trials.csv.gz')
    assert len(trace)==32 and trace.alpha.between(.05,100).all()
    assert trace.iloc[9:].alpha.nunique()>10
    for n in ['r',*names]:assert trace.iloc[9:][n].nunique()>10
    assert result['objective']['total']<=null['objective']['total']
    import optuna
    study=optuna.load_study(study_name='TrpB_gamma_shared_stop',storage=f'sqlite:///{tmp_path}/gamma.sqlite3')
    assert study.trials[-1].distributions['alpha'].log
    assert result['gamma_shape']==result['gamma_rate']
    with pytest.raises(ValueError,match='overwrite'):gamma.fit_gamma(x,tmp_path,quick=True)


@pytest.fixture(scope='module')
def real_inputs():
    root=Path(__file__).resolve().parents[1]
    config=yaml.safe_load((root/gamma.EXP_REL/'gamma_concentration.yaml').read_text())
    return gamma.Inputs(root/'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-188',config)


def test_real_fields_preserve_profile_support_reference_and_mutation_scale(real_inputs):
    x=real_inputs;seed=int(x.calibration_index.seed.iloc[0]);baseline=x.model(seed,None)
    fields=[]
    for alpha in [.05,1.,100.]:
        model=x.model(seed,alpha)
        np.testing.assert_array_equal(model.reference,baseline.reference)
        np.testing.assert_array_equal(model.field.pair_sites,baseline.field.pair_sites)
        assert model.sigma_a==baseline.sigma_a
        assert model.field.expected_mutation_energy()/model.sigma_a**2==pytest.approx(1.)
        assert np.max(np.abs(model.field.tables.mean(axis=1)))<1e-12
        assert np.max(np.abs(model.field.tables.mean(axis=2)))<1e-12
        fields.append(x.calibration_g(alpha))
    assert not np.allclose(fields[0],fields[-1])
    np.testing.assert_array_equal(fields[1],x.calibration_g(1.))
    control=x.controls[0]
    for alpha in [.05,100.]:
        assert x.gamma_loss(0.,control['parameters'],alpha)==x.loss(0.,control['parameters'])


def test_config_rejects_changed_loss_or_noise(real_inputs):
    x=real_inputs;old=yaml.safe_load((x.control_root/'outputs/config.yaml').read_text())
    for key in ['fit','measurement','model','observation','ratio_bounds','optimizer']:
        changed=copy.deepcopy(x.config);changed[key]={}
        with pytest.raises(ValueError):gamma.validate_config(changed,old)
