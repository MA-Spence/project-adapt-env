"""ANA-064: independent audits and descriptive analysis of frozen RUN-166.

No model fitting or synthetic re-sampling occurs here. Bootstrap uncertainty is
conditional on calibration: empirical cubes and synthetic seed clusters are
resampled, with common weights for paired arms. Numerical assertions verify the
record; they are not scientific acceptance criteria.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import itertools
import json
import math
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ANALYSIS = Path(__file__).resolve().parents[1]
ROOT = ANALYSIS.parents[1]
BUNDLE = ROOT / 'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-166'
FIGURES, TABLES = ANALYSIS / 'figures', ANALYSIS / 'tables'
FEATURES = ['fitness_sd', 'walsh_epistatic_fraction', 'walsh_high_order_fraction',
            'beneficial_diminishing_returns_slope', 'reference_beneficial_fraction',
            'neutral_edge_fraction', 'sign_epistasis_fraction', 'reciprocal_sign_epistasis_fraction',
            'local_maxima_density', 'accessible_path_fraction']
COLORS = {'empirical': '#242424', 'linear': '#0072B2', 'sigmoid': '#D55E00'}
LABELS = {'empirical': 'Measured TrpB', 'linear': 'Linear', 'sigmoid': 'Sigmoid'}
ENSEMBLES = ['distance_0', 'distance_4', 'distance_16', 'uniform']
ECOLORS = dict(zip(ENSEMBLES, ['#0072B2', '#009E73', '#D55E00', '#CC79A7']))
ELABELS = dict(zip(ENSEMBLES, ['Reference anchor', '4 mutations away', '16 mutations away', 'Uniform sequences']))
POLICIES = ['uniform_sites', 'shared_support']
POLABELS = ['Uniform sites', 'One shared order-4 support']
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                    'axes.spines.right': False, 'axes.titlesize': 12, 'axes.labelsize': 10,
                    'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'pdf.fonttype': 42})


def dump_json(path, data):
    def clean(x):
        if isinstance(x, dict): return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple, np.ndarray)): return [clean(v) for v in x]
        if isinstance(x, (float, np.floating)): return float(x) if np.isfinite(x) else None
        if isinstance(x, np.integer): return int(x)
        if isinstance(x, np.bool_): return bool(x)
        return x
    Path(path).write_text(json.dumps(clean(data), indent=2, allow_nan=False) + '\n')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''): h.update(block)
    return h.hexdigest()


def key(frame):
    return (frame.phase + '|' + frame.cube_id + '|' + frame.readout + '|' +
            frame.ratio.fillna(-1).astype(str) + '|' + frame.p0.fillna(-1).astype(str))


def raw_readout(a, g, p):
    z = a + p['ratio'] * g
    return z if p['readout'] == 'linear' else expit(np.log(p['p0'] / (1-p['p0'])) + p['beta'] * z)


def direct_statistics(f, resolution):
    """Vectorized direct contrasts and explicit paths; independent of RUN code."""
    n = np.arange(16)
    h = np.array([[(-1.) ** (int(i & j).bit_count()) for j in n] for i in n])
    coefficients = f @ h / 16
    energy = np.column_stack([(coefficients[:, [i for i in range(1, 16) if i.bit_count() == k]]**2).sum(1) for k in range(1, 5)])
    total = energy.sum(1)
    valid = total > np.maximum(1e-18, resolution**2)
    fractions = np.divide(energy, total[:, None], out=np.full_like(energy, np.nan), where=valid[:, None])
    low = np.array([[i for i in n if not i & (1 << bit)] for bit in range(4)])
    high = low | (1 << np.arange(4))[:, None]
    x, delta = f[:, low], f[:, high] - f[:, low]
    xc, dc = x - x.mean(2, keepdims=True), delta - delta.mean(2, keepdims=True)
    variance = (xc**2).mean(2)
    ok = (variance > np.maximum(1e-18, resolution[:, None]**2)) & (delta.mean(2) > resolution[:, None])
    slopes = np.divide((xc*dc).mean(2), variance, out=np.zeros_like(variance), where=ok)
    slope = np.divide(slopes.sum(1), ok.sum(1), out=np.full(len(f), np.nan), where=ok.sum(1)>0)
    neighbours = n[:, None] ^ (1 << np.arange(4))
    peaks = (f[:, :, None] - f[:, neighbours] > resolution[:, None, None]).all(2).mean(1)
    single, reciprocal = np.zeros(len(f)), np.zeros(len(f))
    for a, b in itertools.combinations(range(4), 2):
        for start in n[(n & ((1 << a) | (1 << b))) == 0]:
            v00, v10, v01, v11 = [f[:, i] for i in [start, start | (1 << a), start | (1 << b), start | (1 << a) | (1 << b)]]
            da, db = v10-v00, v11-v01
            dc1, dd = v01-v00, v11-v10
            flip1 = ((da>resolution)&(db < -resolution)) | ((da < -resolution)&(db>resolution))
            flip2 = ((dc1>resolution)&(dd < -resolution)) | ((dc1 < -resolution)&(dd>resolution))
            single += flip1 | flip2; reciprocal += flip1 & flip2
    oriented = np.where((f[:, -1]<f[:, 0])[:, None], f[:, ::-1], f)
    paths = np.zeros(len(f))
    for order in itertools.permutations(range(4)):
        nodes = [0]
        for bit in order: nodes.append(nodes[-1] | (1 << bit))
        paths += (np.diff(oriented[:, nodes], axis=1) > resolution[:, None]).all(1)
    result = {'fitness_sd': f.std(1), 'total_variance': total,
              'spectral_resolved': valid.astype(float), 'walsh_epistatic_fraction': fractions[:, 1:].sum(1),
              'walsh_high_order_fraction': fractions[:, 2:].sum(1),
              'beneficial_diminishing_returns_slope': slope, 'beneficial_slope_n_mutations': ok.sum(1),
              'reference_beneficial_fraction': ((f[:, 1 << np.arange(4)] - f[:, :1]) > resolution[:, None]).mean(1),
              'neutral_edge_fraction': (abs(delta) <= resolution[:, None, None]).mean((1, 2)),
              'sign_epistasis_fraction': single/24, 'reciprocal_sign_epistasis_fraction': reciprocal/24,
              'local_maxima_density': peaks, 'accessible_path_fraction': paths/24}
    result.update({f'walsh_energy_order_{k}': energy[:, k-1] for k in range(1,5)})
    result.update({f'walsh_fraction_order_{k}': fractions[:, k-1] for k in range(1,5)})
    return result


def embedding(reference, rows):
    ref, x = reference[FEATURES].to_numpy(float), rows[FEATURES].to_numpy(float)
    center = np.array([np.mean(c[np.isfinite(c)]) if np.isfinite(c).any() else 0 for c in ref.T])
    scale = np.array([max(np.std(c[np.isfinite(c)]), .05) if np.isfinite(c).any() else .05 for c in ref.T])
    valid = np.isfinite(x)
    z = np.divide(x-center, scale, out=np.zeros_like(x), where=valid)
    return np.c_[z, valid.astype(float)] / np.sqrt(len(FEATURES))


def distance(x, y):
    return max(0., 2*cdist(x,y).mean()-cdist(x,x).mean()-cdist(y,y).mean())


def load_and_audit():
    FIGURES.mkdir(exist_ok=True); TABLES.mkdir(exist_ok=True)
    if not (BUNDLE/'outputs/manifest.json').exists():
        raise FileNotFoundError('Restore RUN-166.dvc before executing this notebook; see the ANA README.')
    c = {n: pd.read_csv(BUNDLE/'outputs'/f'{n}.csv', float_precision='round_trip') for n in
         ['cube_statistics', 'cube_values', 'latent_cubes', 'empirical_cubes', 'calibration',
          'mutation_effects', 'readout_comparison', 'conditional_summary']}
    c.update({n: json.loads((BUNDLE/'outputs'/f'{n}.json').read_text()) for n in ['manifest', 'selected_parameters', 'summary']})
    m, cfg = c['manifest'], c['manifest']['config']
    assert m['run_id']=='RUN-166' and not m['quick'] and c['summary']['status']=='DESCRIPTIVE_FOLLOWUP_COMPLETE'
    assert m['metric_features']==FEATURES and (BUNDLE/'logs/labproj-exit-code.txt').read_text().strip()=='0'
    remote_files = (BUNDLE/'remote-output-sha256.txt').read_text().splitlines()
    assert len(remote_files)==20
    for line in remote_files:
        expected, filename = line.split(None, 1)
        assert sha256(BUNDLE/filename)==expected, filename
    for filename, expected in m['source_sha256'].items(): assert sha256(BUNDLE/'source'/filename)==expected, filename
    prior_path = BUNDLE/'source'/cfg['data']['prior_nodes']
    prior = json.loads(prior_path.read_text())
    for filename, expected in m['input_sha256'].items():
        if filename.endswith('.json'): local=prior_path
        else: local=ROOT/cfg['data']['root']/'datasets'/cfg['data']['dataset']/'variants.parquet'
        assert sha256(local)==expected, str(local)
    empirical, values, latent, stats = (c[n] for n in ['empirical_cubes','cube_values','latent_cubes','cube_statistics'])
    assert [len(empirical),len(values),len(latent),len(stats),len(c['mutation_effects'])]==[48,3312,1248,9936,105984]
    assert not key(values).duplicated().any() and not latent.cube_id.duplicated().any()
    assert (stats.groupby(key(stats)).size()==3).all()
    assert not set(cfg['fit']['calibration_seeds']) & set(cfg['fit']['confirmation_seeds'])
    assert not set(cfg['fit']['calibration_seeds']+cfg['fit']['confirmation_seeds']) & set(cfg['embedding']['seeds'])
    nodes = [tuple(n) for row in empirical.node_keys.map(json.loads) for n in row]
    excluded = {tuple(n) for n in prior['node_keys']+prior['engineering_node_keys']}
    assert len(nodes)==len(set(nodes))==768 and len(excluded)==896 and not set(nodes)&excluded
    # Confirm the saved measurements against the prepared assay, not just its hash.
    table = pd.read_parquet(ROOT/cfg['data']['root']/'datasets'/cfg['data']['dataset']/'variants.parquet')
    measured = table.groupby(['pos1','pos2','pos3','pos4']).fitness.mean().to_dict()
    raw = np.array(empirical.raw_values.map(json.loads).tolist())
    np.testing.assert_allclose(raw.ravel(), [measured[n] for n in nodes], rtol=0, atol=1e-14)
    eoffset, escale = raw[:24].min(), np.ptp(raw[:24])
    np.testing.assert_allclose(empirical.raw_offset, eoffset); np.testing.assert_allclose(empirical.raw_scale, escale)
    v = np.array(values['values'].map(json.loads).tolist())
    assert np.isfinite(v).all()
    positions = pd.Series(np.arange(len(values)),index=key(values))
    vi = positions.loc[key(stats)].to_numpy()
    independent = direct_statistics(v[vi], stats.resolution.to_numpy())
    errors={}
    for name, reproduced in independent.items():
        observed=stats[name].to_numpy(float)
        np.testing.assert_allclose(reproduced, observed, rtol=2e-7, atol=1e-10, equal_nan=True, err_msg=name)
        errors[name]=float(np.nanmax(abs(reproduced-observed)))
    # All 32 oriented edges per observed cube, including bit/background identity.
    edges=c['mutation_effects']; ei=positions.loc[key(edges)].to_numpy()
    expected=v[ei,edges.target_node.to_numpy()]-v[ei,edges.source_node.to_numpy()]
    assert (edges.target_node.to_numpy()==(edges.source_node.to_numpy() | (1 << edges.mutation.to_numpy()))).all()
    assert (edges.groupby(key(edges)).size()==32).all()
    np.testing.assert_allclose(edges.effect,expected,rtol=0,atol=2e-14)
    np.testing.assert_allclose(edges.background_fitness,v[ei,edges.source_node.to_numpy()],rtol=0,atol=2e-14)
    aa=np.array(latent.additive.map(json.loads).tolist()); gg=np.array(latent.gnk.map(json.loads).tolist())
    li=pd.Series(np.arange(len(latent)),index=latent.cube_id)
    replay_error=0.
    for index, row in values.iterrows():
        if row.phase=='empirical':
            j=int(row.cube_id.split('_')[1]); rebuilt=(raw[j]-eoffset)/escale
        else:
            j=li[row.cube_id]; rebuilt=row.offset+row.scale*raw_readout(aa[j],gg[j],row)
        replay_error=max(replay_error,float(np.max(abs(rebuilt-v[index]))))
        np.testing.assert_allclose(rebuilt,v[index],rtol=3e-13,atol=2e-14)
    # All forty affine calibrations and within-ratio selections use calibration only.
    cal_mask=(latent.phase=='active_four_site')&(latent.split=='calibration')
    emp=(raw[:24]-eoffset)/escale
    cal_reference=stats[(stats.phase=='empirical')&(stats.split=='calibration')&(stats.resolution==.001)]
    for _, p in c['calibration'].iterrows():
        rr=raw_readout(aa[cal_mask],gg[cal_mask],p)
        scale=emp.std(1).mean()/rr.std(1).mean(); offset=emp.mean()-scale*rr.mean()
        np.testing.assert_allclose([p.scale,p.offset],[scale,offset],rtol=1e-12,atol=1e-14)
        cal_stats=pd.DataFrame(direct_statistics(offset+scale*rr,np.full(len(rr),.001)))
        rebuilt_loss=distance(embedding(cal_reference,cal_reference),embedding(cal_reference,cal_stats))
        np.testing.assert_allclose(rebuilt_loss,p.loss,rtol=0,atol=1e-11)
    for p in c['selected_parameters']['within_ratio_pairs']:
        arm=c['calibration'][(c['calibration'].ratio==p['ratio'])&(c['calibration'].readout==p['readout'])]
        winner=arm.loc[arm.loss.idxmin()]
        np.testing.assert_allclose([p.get(k,np.nan) if p.get(k) is not None else np.nan for k in ['beta','p0','offset','scale']],winner[['beta','p0','offset','scale']].to_numpy(float),rtol=1e-12,atol=1e-14,equal_nan=True)
        assert np.isclose(p['calibration_loss'],arm.loss.min())
        arm_values=values[(values.phase=='active_four_site')&(values.ratio==p['ratio'])&(values.readout==p['readout'])]
        for name in ['beta','p0','offset','scale']:
            np.testing.assert_allclose(arm_values[name],np.nan if p[name] is None else p[name],rtol=1e-12,atol=1e-14,equal_nan=True)
    selected=min((p for p in c['selected_parameters']['within_ratio_pairs'] if p['readout']=='sigmoid'),key=lambda p:p['calibration_loss'])
    assert selected==c['selected_parameters']['embedded_pair'][1]
    for frozen in c['selected_parameters']['embedded_pair']:
        transferred=values[(values.phase=='embedded')&(values.readout==frozen['readout'])]
        for name in ['ratio','beta','p0','offset','scale']:
            np.testing.assert_allclose(transferred[name],np.nan if frozen[name] is None else frozen[name],rtol=1e-12,atol=1e-14,equal_nan=True)
    intervened=values[values.phase=='operating_position']
    for name in ['ratio','beta','offset','scale']:
        np.testing.assert_allclose(intervened[name],selected[name],rtol=1e-12,atol=1e-14)
    assert set(intervened.p0)==set(cfg['conditional']['operating_occupancies'])
    # Pairing between policies really holds background and alternate sequence fixed.
    embedded=latent[latent.phase=='embedded']
    pg=['length','reference_rank','seed','library','ensemble','cube']
    assert (embedded.groupby(pg).size()==2).all()
    assert (embedded.groupby(pg).anchor.nunique()==1).all() and (embedded.groupby(pg).alternatives.nunique()==1).all()
    for _, row in latent.iterrows():
        sites=json.loads(row.sites); a=json.loads(row.anchor); alt=json.loads(row.alternatives)
        assert len(sites)==len(set(sites))==4 and len(a)==len(alt)==row.length
        assert all(a[s]!=alt[s] for s in sites)
        if row.phase=='embedded' and row.ensemble!='uniform': assert row.hamming_distance==int(row.ensemble.split('_')[1])
    c['primary']=stats[stats.resolution==.001].copy()
    c['primary']['high_order_power']=c['primary'].walsh_energy_order_3+c['primary'].walsh_energy_order_4
    c['primary']['resolved']=c['primary'].walsh_high_order_fraction.notna().astype(float)
    c['primary']['slope_available']=c['primary'].beneficial_diminishing_returns_slope.notna().astype(float)
    c['ratio']=selected['ratio']
    # Strict monotonicity: inspect actual finite-precision paired active readouts.
    strict=[]
    for ratio in cfg['fit']['interaction_ratios']:
        part=values[(values.phase=='active_four_site')&(values.split=='confirmation')&(values.ratio==ratio)]
        ids=[positions[k] for k in key(part[part.readout=='linear'].sort_values('cube_id'))]
        js=[positions[k] for k in key(part[part.readout=='sigmoid'].sort_values('cube_id'))]
        low=np.concatenate([[n for n in range(16) if not n&(1<<bit)] for bit in range(4)])
        high=low | (1 << np.repeat(np.arange(4),8))
        d1=v[ids][:,high]-v[ids][:,low]; d2=v[js][:,high]-v[js][:,low]
        strict.append({'ratio':ratio,'edges':d1.size,'opposite_sign_edges':int((d1*d2<0).sum()),
                       'sign_or_tie_differences':int((np.sign(d1)!=np.sign(d2)).sum())})
    c['audit']={'run':'RUN-166','remote_output_log_checksums':len(remote_files),'source_hashes':len(m['source_sha256']),
                'input_hashes':len(m['input_sha256']),'unique_empirical_nodes':len(nodes),'excluded_prior_and_engineering_nodes':len(excluded),
                'all_empirical_values_verified_against_prepared_table':True,'latent_cubes':len(latent),
                'replayed_observed_cubes':len(values),'independent_statistic_rows':len(stats),'verified_oriented_edges':len(edges),
                'readout_replay_max_absolute_error':replay_error,'metric_max_absolute_errors':errors,
                'calibration_candidates_affine_and_loss_recomputed':len(c['calibration']),'paired_strict_edges':strict,
                'replay_scope':'Saved latents and frozen readouts; no claim of regenerating GNK coefficients across platforms',
                'analysis_packages':{n:importlib.metadata.version(n) for n in ['numpy','pandas','scipy','matplotlib','nbformat','nbclient','ipykernel']}}
    dump_json(TABLES/'audit.json',c['audit'])
    return c


def weighted_mean(values, weights):
    a=np.asarray(values,float); mask=np.isfinite(a)
    denom=weights@mask.astype(float)
    return np.divide(weights@np.nan_to_num(a,nan=0.),denom,out=np.full(len(weights),np.nan),where=denom>0)


def summarize(c, draws=2000, seed=664001):
    rng=np.random.default_rng(seed)
    we=rng.multinomial(24,np.full(24,1/24),size=draws)/24
    ws=np.repeat(rng.multinomial(8,np.full(8,1/8),size=draws)/64,8,axis=1)
    ref=c['primary'].query("phase=='empirical' and split=='calibration'")
    rows=[]; contrasts=[]; distributions={}
    for resolution in [1e-12,.001,.01]:
        panel=c['cube_statistics'].query("split=='confirmation' and resolution==@resolution")
        emp=panel[panel.phase=='empirical'].sort_values('cube')
        ex=embedding(ref,emp); dee=cdist(ex,ex)
        d_emp=np.sum((we@dee)*we,axis=1)
        for ratio in [0.,.25,1.,4.]:
            for readout in ['linear','sigmoid']:
                syn=panel[(panel.phase=='active_four_site')&(panel.ratio==ratio)&(panel.readout==readout)].sort_values(['seed','cube'])
                sy=embedding(ref,syn)
                assert len(emp)==24 and len(syn)==64 and syn.groupby('seed').size().eq(8).all()
                ed=distance(ex,sy)
                reported=c['readout_comparison'].query('resolution==@resolution and ratio==@ratio and readout==@readout').iloc[0].confirmation_joint_distance
                np.testing.assert_allclose(ed,reported,rtol=0,atol=1e-12)
                dist=np.maximum(0,2*np.sum((we@cdist(ex,sy))*ws,axis=1)-d_emp-np.sum((ws@cdist(sy,sy))*ws,axis=1))
                distributions[resolution,ratio,readout]=dist
                lo,hi=np.quantile(dist,[.025,.975])
                rows.append(dict(resolution=resolution,ratio=ratio,readout=readout,distance=ed,low=lo,high=hi))
            delta=distributions[resolution,ratio,'sigmoid']-distributions[resolution,ratio,'linear']
            pair=[r for r in rows if r['resolution']==resolution and r['ratio']==ratio]
            lo,hi=np.quantile(delta,[.025,.975])
            contrasts.append(dict(resolution=resolution,ratio=ratio,delta_sigmoid_minus_linear=pair[1]['distance']-pair[0]['distance'],low=lo,high=hi))
    c['distances']=pd.DataFrame(rows); c['contrasts']=pd.DataFrame(contrasts)
    metrics=FEATURES+['fitness_mean','high_order_power','spectral_resolved','beneficial_slope_n_mutations']+[f'walsh_{t}_order_{k}' for t in ['energy','fraction'] for k in range(1,5)]
    primary=c['primary'].query("split=='confirmation'")
    summaries=[]
    for readout in ['empirical','linear','sigmoid']:
        ratios=[np.nan] if readout=='empirical' else [0.,.25,1.,4.]
        for ratio in ratios:
            mask=primary.phase.eq('empirical') if readout=='empirical' else (primary.phase.eq('active_four_site')&primary.ratio.eq(ratio)&primary.readout.eq(readout))
            part=primary[mask].sort_values(['seed','cube']); weights=we if readout=='empirical' else ws
            for metric in metrics:
                a=part[metric].to_numpy(float); boot=weighted_mean(a,weights)
                good=boot[np.isfinite(boot)]; lo,hi=np.quantile(good,[.025,.975]) if len(good) else (np.nan,np.nan)
                summaries.append(dict(readout=readout,ratio=ratio,metric=metric,mean=np.nanmean(a) if np.isfinite(a).any() else np.nan,
                                      median=np.nanmedian(a) if np.isfinite(a).any() else np.nan,available=int(np.isfinite(a).sum()),
                                      total=len(a),low=lo,high=hi,bootstrap_available=len(good)))
    c['descriptors']=pd.DataFrame(summaries)
    groupkeys=['reference_rank','readout','site_policy','ensemble','length']
    metrics2=['edge_effect_rms_stable','fitness_sd','high_order_power','walsh_high_order_fraction',
              'walsh_epistatic_fraction','neutral_edge_fraction','zero_edge_fraction','resolved','slope_available','shared_supports']
    emb=c['primary'].query("phase=='embedded'")
    c['embedded_means']=emb.groupby(groupkeys)[metrics2].mean().reset_index()
    c['embedded_seed_means']=emb.groupby(groupkeys+['seed'])[metrics2].mean().reset_index()
    avail=emb.groupby(groupkeys).agg(n_cubes=('cube_id','size'),resolved_cubes=('resolved','sum'),slope_cubes=('slope_available','sum')).reset_index()
    c['embedded_availability']=avail
    c['operating_summary']=primary.query("phase=='operating_position'").groupby('p0').agg(
        mean_fitness=('fitness_mean','mean'),high_order_power=('high_order_power','mean'),
        high_order_fraction=('walsh_high_order_fraction','mean'),beneficial_slope=('beneficial_diminishing_returns_slope','mean'),
        slope_cubes=('slope_available','sum'),resolved_cubes=('resolved','sum'),neutral_edges=('neutral_edge_fraction','mean')).reset_index()
    for name in ['distances','contrasts','descriptors','embedded_means','embedded_seed_means','embedded_availability','operating_summary']:
        c[name].to_csv(TABLES/f'{name}.csv',index=False)
    c['audit']['joint_distances_independently_recomputed']=24
    c['audit']['bootstrap']={'draws':draws,'seed':seed,'units':'24 empirical cubes and 8 synthetic seed clusters of 8 cubes; common resampling weights across readouts and ratios; calibration frozen',
                              'limitations':'Percentile descriptive intervals; no refitting, measurement-error correction, or empirical allele/site-dependence model; embedded panel has only 3 seeds, no bootstrap interval'}
    dump_json(TABLES/'audit.json',c['audit'])
    return c['contrasts']


def save(fig, name):
    for ext in ['png','pdf']: fig.savefig(FIGURES/f'{name}.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)
    return FIGURES/f'{name}.png'


def figure_pairs(c):
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    x=np.arange(4); ratios=[0.,.25,1.,4.]
    for readout,dx in [('linear',-.045),('sigmoid',.045)]:
        p=c['distances'].query('resolution==.001 and readout==@readout').sort_values('ratio')
        ax[0,0].vlines(x+dx,p.low,p.high,color=COLORS[readout],alpha=.55)
        ax[0,0].plot(x+dx,p.distance,'o-',color=COLORS[readout],label=LABELS[readout])
        for a,metric in [(ax[1,0],'fitness_sd'),(ax[1,1],'walsh_high_order_fraction')]:
            p=c['descriptors'].query('readout==@readout and metric==@metric').sort_values('ratio')
            a.vlines(x+dx,p.low,p.high,color=COLORS[readout],alpha=.55);a.plot(x+dx,p['mean'],'o-',color=COLORS[readout])
    p=c['contrasts'].query('resolution==.001').sort_values('ratio')
    ax[0,1].vlines(x,p.low,p.high,color='#6A3D9A',linewidth=2)
    ax[0,1].plot(x,p.delta_sigmoid_minus_linear,'o',color='#6A3D9A')
    ax[0,1].axhline(0,color='grey',ls=':');ax[0,1].text(.43,.93,'Below zero favours sigmoid',transform=ax[0,1].transAxes)
    for a,metric in [(ax[1,0],'fitness_sd'),(ax[1,1],'walsh_high_order_fraction')]:
        p=c['descriptors'].query("readout=='empirical' and metric==@metric").iloc[0]
        a.axhspan(p.low,p.high,color='black',alpha=.07);a.axhline(p['mean'],color='black',ls='--',label='Measured TrpB')
    labels=['A  Joint descriptor distance (lower = closer)','B  Paired difference: sigmoid − linear',
            'C  Within-cube fitness SD','D  WHT fraction at orders ≥3']
    for a,label in zip(ax.flat,labels):
        a.set(title=label,xlabel='Latent interaction/additive ratio r');a.set_xticks(x,labels=['0','0.25','1','4']);a.grid(alpha=.16)
    ax[0,0].legend();ax[1,0].legend()
    fig.suptitle('TrpB confirmation: interaction strength dominates the readout contrast\n24 measured cubes; 64 synthetic cubes / 8 seeds; conditional 95% bootstrap intervals',fontsize=14)
    return save(fig,'01_paired_readouts')


def figure_conditional(c):
    p=c['primary'].query("split=='confirmation'");emp=p.query("phase=='empirical'")
    chosen=p[(p.phase=='active_four_site')&(p.ratio==c['ratio'])]
    fig,ax=plt.subplots(2,3,figsize=(14,9.4),layout='constrained')
    specs=[('high_order_power','Absolute WHT power, orders ≥3'),('walsh_high_order_fraction','WHT fraction, orders ≥3'),
           ('beneficial_diminishing_returns_slope','Beneficial-mutation slope')]
    sets=[[(LABELS[r],chosen[chosen.readout==r],COLORS[r]) for r in ['linear','sigmoid']],
          [(f'p₀ = {p0:g}',p[(p.phase=='operating_position')&(p.p0==p0)],color) for p0,color in zip([.1,.5,.9],['#009E73','#CC79A7','#E69F00'])]]
    for row,groups in enumerate(sets):
        for col,(metric,label) in enumerate(specs):
            a=ax[row,col];counts=[]
            for name,part,color in [*groups,('Measured',emp,'#242424')]:
                a.scatter(part.fitness_mean,part[metric],s=22,alpha=.6,color=color,label=name,marker='x' if name=='Measured' else 'o')
                counts.append(f'{name}: {part[metric].notna().sum()}/{len(part)}')
            for cut in c['selected_parameters']['conditional_bin_boundaries']: a.axvline(cut,color='grey',lw=.6,alpha=.3)
            a.set(xlabel='Mean cube fitness',ylabel=label,title=f'{chr(65+3*row+col)}  '+('Frozen paired readouts' if row==0 else 'Operating-position intervention'))
            a.grid(alpha=.12)
            if col==0:a.set_yscale('log')
            if col==2:a.set_yscale('symlog',linthresh=.3)
            short=[t.replace('Measured','M').replace('Linear','L').replace('Sigmoid','S') for t in counts]
            a.text(.5,-.20,'Available: '+'; '.join(short),ha='center',va='top',transform=a.transAxes,fontsize=8)
        ax[row,0].legend(fontsize=8,loc='upper left')
    fig.suptitle('Conditional statistics at r = 4: curvature moves the pattern without reproducing it\nOne point per cube; vertical guides are empirical calibration quartiles; log power axis; symmetric-log slope axis; M = measured',fontsize=14)
    return save(fig,'02_conditional_statistics')


def figure_fingerprint(c):
    fig,ax=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
    chosen=c['descriptors'][(c['descriptors'].readout=='empirical')|(c['descriptors'].ratio==c['ratio'])]
    x=np.arange(1,5)
    for a,kind in [(ax[0,0],'energy'),(ax[0,1],'fraction')]:
        for i,readout in enumerate(['empirical','linear','sigmoid']):
            part=chosen[(chosen.readout==readout)&chosen.metric.str.startswith(f'walsh_{kind}_order_')].sort_values('metric')
            a.bar(x+(i-1)*.23,part['mean'],width=.22,color=COLORS[readout],label=LABELS[readout],alpha=.85)
            a.vlines(x+(i-1)*.23,part.low,part.high,color=COLORS[readout],lw=1.2)
        a.set_xticks(x);a.set_xlabel('WHT interaction order');a.grid(axis='y',alpha=.15)
    ax[0,0].set(title='A  Mean absolute WHT power',ylabel='Squared normalized fitness');ax[0,0].legend(fontsize=9)
    ax[0,1].set(title='B  Mean per-cube WHT fraction',ylabel='Fraction of nonconstant variance')
    keys=['beneficial_diminishing_returns_slope','sign_epistasis_fraction','reciprocal_sign_epistasis_fraction','local_maxima_density','accessible_path_fraction','neutral_edge_fraction']
    names=['Beneficial slope','Sign epistasis','Reciprocal sign','Strict-peak density','Accessible paths','Neutral edges']
    reference=c['primary'].query("phase=='empirical' and split=='calibration'")
    matrix=[]
    for readout in ['linear','sigmoid']:
        row=[]
        for k in keys:
            em=chosen.query("readout=='empirical' and metric==@k").iloc[0]['mean']
            sm=chosen.query('readout==@readout and metric==@k').iloc[0]['mean']
            scale=max(reference[k].std(ddof=0),.05)
            row.append((sm-em)/scale)
        matrix.append(row)
    limit=max(1.,np.max(abs(np.array(matrix))))
    im=ax[1,0].imshow(matrix,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto')
    ax[1,0].set_yticks([0,1],labels=['Linear','Sigmoid']);ax[1,0].set_xticks(range(6),labels=names,rotation=30,ha='right')
    for i in range(2):
        for j in range(6): ax[1,0].text(j,i,f'{matrix[i][j]:+.2f}',ha='center',va='center',fontsize=9,color='white' if abs(matrix[i][j])>limit*.65 else 'black')
    ax[1,0].set_title('C  Mean differences from measured TrpB');fig.colorbar(im,ax=ax[1,0],label='Difference / frozen calibration scale',shrink=.8)
    for readout in ['linear','sigmoid']:
        part=c['distances'].query('ratio==4 and readout==@readout').sort_values('resolution')
        ax[1,1].plot(range(3),part.distance,'o-',color=COLORS[readout],label=LABELS[readout])
        ax[1,1].vlines(range(3),part.low,part.high,color=COLORS[readout],alpha=.55)
    ax[1,1].set_xticks(range(3),labels=['10⁻¹²','10⁻³ (primary)','10⁻²'])
    ax[1,1].set(title='D  Resolution changes the apparent agreement',xlabel='Descriptive resolution (not measured noise)',ylabel='Joint descriptor distance');ax[1,1].grid(alpha=.15)
    fig.suptitle('At r = 4, aggregate resemblance remains incomplete\n95% conditional bootstrap intervals; absolute and relative WHT summaries answer different questions',fontsize=14)
    return save(fig,'03_aggregate_fingerprint')


def figure_embedded(c):
    p=c['embedded_means'].query('reference_rank==.5')
    fig,ax=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    specs=[('edge_effect_rms_stable','Mean within-cube edge RMS'),('walsh_high_order_fraction','Mean WHT fraction, orders ≥3'),('resolved','Fraction of cubes with resolved WHT')]
    for row,readout in enumerate(['linear','sigmoid']):
        for col,(metric,label) in enumerate(specs):
            a=ax[row,col]
            for ensemble in ENSEMBLES:
                for policy,ls,marker in zip(POLICIES,['-','--'],['o','s']):
                    part=p[(p.readout==readout)&(p.ensemble==ensemble)&(p.site_policy==policy)].sort_values('length')
                    a.plot(part.length,part[metric],color=ECOLORS[ensemble],ls=ls,marker=marker,ms=4,alpha=.85)
            a.set(xlabel='Sequence length',ylabel=label,title=f'{chr(65+3*row+col)}  {LABELS[readout]}')
            a.set_xscale('log',base=2);a.set_xticks([32,128,512],labels=['32','128','512']);a.grid(alpha=.15)
            if col==0:a.set_yscale('log');a.set_ylim(.01,.12)
            if col==2:a.set_ylim(-.05,1.05)
    handles=[Line2D([0],[0],color=color,label=ELABELS[e]) for e,color in ECOLORS.items()]
    handles += [Line2D([0],[0],color='grey',ls=ls,marker=m,label=label) for ls,m,label in zip(['-','--'],['o','s'],POLABELS)]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,fontsize=9)
    fig.suptitle('Embedded four-site libraries: activity near the reference survives increasing length\nReference rank 0.5; frozen r = 4 readouts; 24 cubes / 3 seeds per point; WHT fractions use resolved cubes only',fontsize=14)
    return save(fig,'04_embedded_neighbourhoods')


def figure_saturation(c):
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    p=c['embedded_means'].query("readout=='sigmoid'");avail=c['embedded_availability'].query("readout=='sigmoid'")
    for row,rank in enumerate([.5,1.]):
        for col,policy in enumerate(POLICIES):
            part=p[(p.reference_rank==rank)&(p.site_policy==policy)]
            matrix=part.pivot(index='ensemble',columns='length',values='neutral_edge_fraction').reindex(ENSEMBLES)
            counts=avail[(avail.reference_rank==rank)&(avail.site_policy==policy)].pivot(index='ensemble',columns='length',values='resolved_cubes').reindex(ENSEMBLES)
            zeros=part.pivot(index='ensemble',columns='length',values='zero_edge_fraction').reindex(ENSEMBLES)
            im=ax[row,col].imshow(matrix,cmap='YlOrRd',vmin=0,vmax=1,aspect='auto')
            ax[row,col].set_xticks(range(3),labels=['32','128','512']);ax[row,col].set_yticks(range(4),labels=[ELABELS[e] for e in ENSEMBLES])
            ax[row,col].set(xlabel='Sequence length',title=f'{chr(65+row*2+col)}  Rank {rank:g} · {POLABELS[col]}')
            for i in range(4):
                for j in range(3):
                    val=matrix.iloc[i,j]
                    ax[row,col].text(j,i,f'{100*val:.0f}% neutral\n{int(counts.iloc[i,j])}/24 resolved\n{100*zeros.iloc[i,j]:.0f}% exact zero',ha='center',va='center',fontsize=8,color='white' if val>.65 else '#242424')
    fig.colorbar(im,ax=ax,label='Mean neutral-edge fraction at resolution 0.001',shrink=.8)
    fig.suptitle('Saturation depends on the sampled background and reference convention\nSigmoid readout; rank 1 chooses the most favourable additive allele at each reference site',fontsize=14)
    return save(fig,'05_saturation_and_availability')
